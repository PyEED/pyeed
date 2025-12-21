from __future__ import annotations

import asyncio
import random
import time
from collections.abc import AsyncIterable, AsyncIterator, Iterable
from dataclasses import dataclass

import numpy as np
import torch
from loguru import logger
from rich.progress import Progress, TaskID
from transformers import EsmModel, EsmTokenizer

from .pooling import PoolingFn, l2_normalize, mean_pooling
from .types import EmbeddingRecord, SequenceItem, TokenizedJob, TorchDTypeName
from .utils import _free_device_memory, _login_hf, silence_transformers_init_only

_TORCH_DTYPE_MAP: dict[TorchDTypeName, torch.dtype] = {
    "float16": torch.float16,
    "float32": torch.float32,
}

_NP_DTYPE_MAP: dict[TorchDTypeName, type[np.floating]] = {
    "float16": np.float16,
    "float32": np.float32,
}


class LaunchGate:
    """Global gate to stagger GPU kernel launches and prevent power spikes.

    Enforces minimum time spacing between forward pass starts across all devices
    to desynchronize power consumption peaks. Supports adaptive spacing adjustment
    and optional jitter to prevent phase locking.

    The gate only delays the launch edge of compute; actual forward passes run
    in parallel after launch, preserving throughput while spacing power ramps.

    Args:
        min_spacing_ms: Minimum milliseconds between consecutive launches (0 = no gating).
        jitter_ms: Random jitter [0, jitter_ms] added to spacing to prevent phase lock.

    Example:
        >>> gate = LaunchGate(min_spacing_ms=50.0, jitter_ms=10.0)
        >>> await gate.acquire()  # Waits if needed, then records launch time
        >>> # ... run forward pass ...
        >>> gate.set_spacing(80.0)  # Increase spacing adaptively
    """

    def __init__(
        self,
        min_spacing_ms: float = 0.0,
        jitter_ms: float = 0.0,
    ) -> None:
        if min_spacing_ms < 0:
            msg = "min_spacing_ms must be non-negative"
            raise ValueError(msg)
        if jitter_ms < 0:
            msg = "jitter_ms must be non-negative"
            raise ValueError(msg)

        self._min_spacing_ms = min_spacing_ms
        self._jitter_ms = jitter_ms
        self._last_launch_time: float = 0.0
        self._lock = asyncio.Lock()

    async def acquire(self) -> None:
        """Wait until sufficient time has passed since last launch across all devices."""
        async with self._lock:
            if self._min_spacing_ms > 0:
                now = time.time()
                elapsed_ms = (now - self._last_launch_time) * 1000
                required_spacing = self._min_spacing_ms

                # Add jitter to prevent phase locking
                if self._jitter_ms > 0:
                    required_spacing += random.uniform(0, self._jitter_ms)

                wait_ms = required_spacing - elapsed_ms

                if wait_ms > 0:
                    await asyncio.sleep(wait_ms / 1000)

            # Record launch time
            self._last_launch_time = time.time()

    def set_spacing(self, min_spacing_ms: float) -> None:
        """Update minimum spacing for adaptive adjustment.

        This can be called at runtime to increase/decrease spacing based on
        observed power behavior or system stability.
        """
        if min_spacing_ms < 0:
            msg = "min_spacing_ms must be non-negative"
            raise ValueError(msg)
        self._min_spacing_ms = min_spacing_ms

    @property
    def current_spacing_ms(self) -> float:
        """Current minimum spacing setting."""
        return self._min_spacing_ms


@dataclass(slots=True, frozen=True)
class ESM2DeviceWorker:
    """Processes ESM2 jobs on a single device and emits results."""

    device: str
    in_q: asyncio.Queue[TokenizedJob | None]
    out_q: asyncio.Queue[list[EmbeddingRecord]]

    model: EsmModel
    torch_dtype: torch.dtype
    np_dtype: type[np.floating]
    pooling: PoolingFn | None
    normalize: bool
    launch_gate: LaunchGate | None = None

    async def run(self) -> None:
        if self.device.startswith("cuda:"):
            dev_i = int(self.device.split(":")[1])
            torch.cuda.set_device(dev_i)

        while True:
            job = await self.in_q.get()
            if job is None:
                self.in_q.task_done()
                break
            try:
                batch_results = await self._process_job(job)
                await self.out_q.put(batch_results)
            finally:
                self.in_q.task_done()

    async def _process_job(
        self,
        job: TokenizedJob,
    ) -> list[EmbeddingRecord]:
        # Acquire launch gate before forward pass to stagger kernel launches
        if self.launch_gate is not None:
            await self.launch_gate.acquire()

        logger.debug(f"Lock aquired at {time.time()}")
        return await asyncio.to_thread(self.embed, job)

    def embed(
        self,
        job: TokenizedJob,
    ) -> list[EmbeddingRecord]:
        t0 = time.time()

        with torch.inference_mode():
            input_ids = job.input_ids.to(self.device)
            attention_mask = job.attention_mask.to(self.device)

            use_autocast = self.device.startswith("cuda") and self.torch_dtype == torch.float16

            if use_autocast:
                with torch.autocast(device_type="cuda", dtype=self.torch_dtype):
                    outputs = self.model(input_ids=input_ids, attention_mask=attention_mask)
            else:
                outputs = self.model(input_ids=input_ids, attention_mask=attention_mask)

            hidden = outputs.last_hidden_state

            if self.pooling is None:
                result_t = hidden
            else:
                result_t = self.pooling(hidden, attention_mask)
                if self.normalize:
                    result_t = l2_normalize(result_t)

        arr = result_t.detach().cpu().numpy().astype(self.np_dtype)

        out = [EmbeddingRecord(id=item.id, vector=arr[i]) for i, item in enumerate(job.batch)]

        elapsed = time.time() - t0
        logger.debug(
            f"Embedding batch completed in {elapsed:.3f}s on Device {self.device}",
            extra={
                "batch_size": len(job.batch),
                "elapsed_s": round(elapsed, 3),
                "device": self.device,
            },
        )
        return out


class ESM2Processor[T, R]:
    """Distributes length-sorted chunk batches across independent devices.

    Args:
        model_name: HuggingFace model identifier for ESM2.
        dtype: Data type for model computation and output embeddings.
        pooling: Pooling function to reduce sequence dimension.
        normalize: Whether to L2-normalize output embeddings.
        max_length: Maximum sequence length (tokens).
        huggingface_token: HuggingFace token for private models.
        devices: List of CUDA device indices. Empty uses all available.
        max_padded_tokens: Max padded tokens per batch for memory efficiency.
        shared_queue_size: Size of shared queue for all GPU workers (default 4).
        source_chunk_size: Items to buffer before packing into batches.
        launch_gate_spacing_ms: Minimum milliseconds between GPU kernel launches (0=disabled).
        launch_gate_jitter_ms: Random jitter added to spacing to prevent phase locking.
        progress: Optional Rich Progress for tracking.
        progress_task_id: Task ID for progress updates.
    """

    def __init__(
        self,
        model_name: str = "facebook/esm2_t33_650M_UR50D",
        dtype: TorchDTypeName = "float16",
        pooling: PoolingFn | None = mean_pooling,
        normalize: bool = True,
        max_length: int = 1024,
        huggingface_token: str | None = None,
        devices: list[int] | None = None,
        max_padded_tokens: int = 16_000,
        shared_queue_size: int = 4,
        source_chunk_size: int = 2048,
        launch_gate_spacing_ms: float = 500,
        launch_gate_jitter_ms: float = 100,
        progress: Progress | None = None,
        progress_task_id: TaskID | None = None,
    ) -> None:
        self.model_name = model_name
        self.dtype = dtype
        self.torch_dtype = _TORCH_DTYPE_MAP[dtype]
        self.np_dtype = _NP_DTYPE_MAP[dtype]
        self.pooling = pooling
        self.normalize = normalize
        self.max_length = max_length

        self.hf_token = huggingface_token or _login_hf()
        self.device_args = devices if devices is not None else []

        self.max_padded_tokens = max_padded_tokens
        self.shared_queue_size = shared_queue_size
        self.source_chunk_size = source_chunk_size
        self.launch_gate_spacing_ms = launch_gate_spacing_ms
        self.launch_gate_jitter_ms = launch_gate_jitter_ms

        if progress is None and progress_task_id is not None:
            msg = "progress_task_id provided without a Progress instance"
            raise ValueError(msg)
        self.progress = progress
        self.progress_task_id = progress_task_id

        self._initialized = False
        self._devices: list[str] = []
        self._tokenizer: EsmTokenizer | None = None
        self._models: dict[str, EsmModel] = {}
        self._launch_gate: LaunchGate | None = None

    async def _initialize(self) -> None:
        if self._initialized:
            return

        self._devices = self._resolve_devices()

        for dev in self._devices:
            self._models[dev] = await asyncio.to_thread(self._load_model, dev)

        self._initialized = True

    def _load_tokenizer(self) -> EsmTokenizer:
        with silence_transformers_init_only():
            return EsmTokenizer.from_pretrained(self.model_name, token=self.hf_token)

    def _load_model(self, device: str) -> EsmModel:
        if device.startswith("cuda:"):
            torch.cuda.set_device(int(device.split(":")[1]))

        with silence_transformers_init_only():
            model = EsmModel.from_pretrained(
                self.model_name,
                token=self.hf_token,
                torch_dtype=self.torch_dtype,
            )

        model.config.output_hidden_states = False
        model = model.to(device=device, dtype=self.torch_dtype)
        model.eval()
        return model

    def _resolve_devices(self) -> list[str]:
        if not torch.cuda.is_available():
            return ["cpu"]

        if not self.device_args:
            count = torch.cuda.device_count()
            return [f"cuda:{i}" for i in range(count)] if count > 0 else ["cpu"]

        max_device = torch.cuda.device_count() - 1
        out: list[str] = []
        for idx in self.device_args:
            if 0 <= idx <= max_device:
                out.append(f"cuda:{idx}")
        return out or ["cpu"]

    async def _build_workers(
        self,
    ) -> tuple[
        list[ESM2DeviceWorker],
        asyncio.Queue[TokenizedJob | None],
        asyncio.Queue[list[EmbeddingRecord]],
        EsmTokenizer,
        LaunchGate | None,
    ]:
        await self._initialize()

        # Single shared queue for all GPU workers
        shared_queue: asyncio.Queue[TokenizedJob | None] = asyncio.Queue(
            maxsize=self.shared_queue_size
        )
        out_q: asyncio.Queue[list[EmbeddingRecord]] = asyncio.Queue()
        workers: list[ESM2DeviceWorker] = []

        # Single tokenizer
        tokenizer = await asyncio.to_thread(self._load_tokenizer)

        # Create launch gate if spacing is configured
        launch_gate: LaunchGate | None = None
        if self.launch_gate_spacing_ms > 0 or self.launch_gate_jitter_ms > 0:
            launch_gate = LaunchGate(
                min_spacing_ms=self.launch_gate_spacing_ms,
                jitter_ms=self.launch_gate_jitter_ms,
            )
            logger.info(
                "Launch gate enabled",
                extra={
                    "spacing_ms": self.launch_gate_spacing_ms,
                    "jitter_ms": self.launch_gate_jitter_ms,
                },
            )

        for device in self._devices:
            worker = ESM2DeviceWorker(
                device=device,
                in_q=shared_queue,
                out_q=out_q,
                model=self._models[device],
                torch_dtype=self.torch_dtype,
                np_dtype=self.np_dtype,
                pooling=self.pooling,
                normalize=self.normalize,
                launch_gate=launch_gate,
            )
            workers.append(worker)

        return workers, shared_queue, out_q, tokenizer, launch_gate

    def _pack_by_padded_budget(self, chunk: list[SequenceItem]) -> list[list[SequenceItem]]:
        """Packs a chunk into batches under a padded-token budget."""
        items = sorted(chunk, key=lambda x: len(x.sequence), reverse=True)

        batches: list[list[SequenceItem]] = []
        batch: list[SequenceItem] = []
        max_len = 0

        for item in items:
            item_len = len(item.sequence)
            if item_len > self.max_length - 2:  # -2 for special tokens
                logger.warning(
                    f"Item {item.id} is too long and will be skipped",
                    extra={"id": item.id, "length": item_len},
                )
                continue

            if not batch:
                batch = [item]
                max_len = item_len
                continue

            new_max = max(max_len, item_len)
            new_padded = (len(batch) + 1) * new_max

            if new_padded > self.max_padded_tokens:
                batches.append(batch)
                batch = [item]
                max_len = item_len
                continue

            batch.append(item)
            max_len = new_max

        if batch:
            batches.append(batch)

        return batches

    async def _tokenization_stage(
        self,
        batch_stream: AsyncIterator[list[SequenceItem]],
        output_queue: asyncio.Queue[TokenizedJob | None],
        tokenizer: EsmTokenizer,
        done_event: asyncio.Event,
    ) -> None:
        """Tokenizes batches on CPU and enqueues for GPU workers."""

        async for batch in batch_stream:
            if not batch:
                continue

            # Tokenize on CPU (run in thread to avoid blocking event loop)
            tokenized = await asyncio.to_thread(self._tokenize_batch, batch, tokenizer)

            if tokenized is not None:
                job = TokenizedJob(
                    batch=batch,
                    input_ids=tokenized[1],
                    attention_mask=tokenized[2],
                    num_tokens=int(tokenized[2].sum().item()),
                )
                await output_queue.put(job)

        # Signal that tokenization is complete
        done_event.set()

    def _tokenize_batch(
        self,
        batch: list[SequenceItem],
        tokenizer: EsmTokenizer,
    ) -> tuple[list[str], torch.Tensor, torch.Tensor] | None:
        """Tokenize a batch of sequences on CPU."""
        if not batch:
            return None

        ids, seqs = zip(*[(item.id, item.sequence) for item in batch], strict=True)

        tokens = tokenizer(
            list(seqs),
            return_tensors="pt",
            truncation=False,
            padding=True,
            add_special_tokens=True,
        )

        # Keep tensors on CPU
        input_ids = tokens["input_ids"]
        attention_mask = tokens["attention_mask"]

        return (list(ids), input_ids, attention_mask)

    async def stream_work(
        self,
        chunks: AsyncIterable[list[SequenceItem]],
    ) -> AsyncIterator[EmbeddingRecord]:
        workers, shared_queue, out_q, tokenizer, launch_gate = await self._build_workers()
        self._launch_gate = launch_gate

        job_counter = 0
        processed_jobs = 0
        tokenization_done = asyncio.Event()

        async def batch_packer() -> AsyncIterator[list[SequenceItem]]:
            """Pack chunks into batches and yield them."""
            nonlocal job_counter

            async for chunk in chunks:
                if not chunk:
                    continue

                for batch in self._pack_by_padded_budget(chunk):
                    job_counter += 1
                    yield batch

        async def shutdown_workers() -> None:
            """Send shutdown signal to all workers."""
            await tokenization_done.wait()
            # Send None to signal shutdown (workers share queue, so send one per worker)
            for _ in self._devices:
                await shared_queue.put(None)

        async with asyncio.TaskGroup() as tg:
            # Stage 1 + 2: Batch packing → Tokenization → Shared queue
            tg.create_task(
                self._tokenization_stage(batch_packer(), shared_queue, tokenizer, tokenization_done)
            )

            # Stage 3: GPU workers consume from shared queue
            for w in workers:
                tg.create_task(w.run())

            # Shutdown coordination
            tg.create_task(shutdown_workers())

            # Stage 4: Result consumer
            while True:
                if tokenization_done.is_set() and processed_jobs >= job_counter:
                    break

                batch_results = await out_q.get()
                processed_jobs += 1
                if self.progress is not None and self.progress_task_id is not None:
                    self.progress.advance(self.progress_task_id, len(batch_results))

                for record in batch_results:
                    yield record

    async def _to_chunk_stream(
        self,
        items: Iterable[T] | AsyncIterable[T],
    ) -> AsyncIterator[list[T]]:
        async def aiter() -> AsyncIterator[T]:
            if hasattr(items, "__aiter__"):
                async for x in items:  # type: ignore[union-attr]
                    yield x
            else:
                for x in items:  # type: ignore[union-attr]
                    yield x

        buf: list[T] = []
        async for x in aiter():
            buf.append(x)
            if len(buf) >= self.source_chunk_size:
                yield buf
                buf = []

        if buf:
            yield buf

    async def work(
        self,
        items: Iterable[T] | AsyncIterable[T],
    ) -> list[EmbeddingRecord]:
        results: list[EmbeddingRecord] = []
        async for record in self.stream_work(self._to_chunk_stream(items)):
            results.append(record)
        return results

    def set_launch_gate_spacing(self, spacing_ms: float) -> None:
        """Adaptively adjust launch gate spacing at runtime.

        Use this to dynamically increase or decrease spacing based on observed
        power behavior or system stability. For example, increase spacing if
        correlated power spikes are detected.

        Args:
            spacing_ms: New minimum spacing in milliseconds (0 = no gating).

        Raises:
            RuntimeError: If processor not initialized or launch gate not enabled.

        Example:
            >>> processor.set_launch_gate_spacing(50.0)  # Increase to 50ms
            >>> # Monitor power traces...
            >>> processor.set_launch_gate_spacing(80.0)  # Increase further if needed
        """
        if self._launch_gate is None:
            msg = "Launch gate not enabled. Set launch_gate_spacing_ms > 0 in constructor."
            raise RuntimeError(msg)

        old_spacing = self._launch_gate.current_spacing_ms
        self._launch_gate.set_spacing(spacing_ms)

        logger.info(
            "Launch gate spacing adjusted",
            extra={
                "old_spacing_ms": old_spacing,
                "new_spacing_ms": spacing_ms,
            },
        )

    async def cleanup(self) -> None:
        device_ids: list[int] = []
        for dev in self._devices:
            if dev.startswith("cuda:"):
                device_ids.append(int(dev.split(":")[1]))

        self._models.clear()
        self._tokenizer = None
        self._devices.clear()
        self._initialized = False
        self._launch_gate = None

        if device_ids:
            await asyncio.to_thread(_free_device_memory, device_ids)

    async def __aenter__(self) -> ESM2Processor[T, R]:
        await self._initialize()
        return self

    async def __aexit__(
        self,
        exc_type: type[BaseException] | None,
        exc_val: BaseException | None,
        exc_tb: object,
    ) -> None:
        await self.cleanup()


if __name__ == "__main__":
    import argparse
    import base64
    import json
    from pathlib import Path

    from ..db.neo4j import get_async_driver

    def vec_to_b64_f16(v: np.ndarray) -> str:
        v = np.asarray(v, dtype=np.float16)
        return base64.b64encode(v.tobytes()).decode("ascii")

    async def stream_neo4j(
        driver,
        query: str,
        batch_size: int,
    ) -> AsyncIterator[list[SequenceItem]]:
        async with driver.session() as session:
            res = await session.run(query)
            buf: list[SequenceItem] = []
            async for r in res:
                buf.append(SequenceItem(id=r["id"], sequence=r["sequence"]))
                if len(buf) >= batch_size:
                    yield buf
                    buf = []
            if buf:
                yield buf

    async def main() -> None:
        ap = argparse.ArgumentParser()
        ap.add_argument("--out", required=True, help="Output JSONL path (append)")
        ap.add_argument("--devices", default="0,1,2", help="CUDA device ids, e.g. 0,1,2")
        ap.add_argument("--read-batch", type=int, default=20000, help="Neo4j fetch batch size")
        ap.add_argument("--max-padded-tokens", type=int, default=30000)
        ap.add_argument("--max-length", type=int, default=1024)
        ap.add_argument("--dtype", choices=["float16", "float32"], default="float16")
        ap.add_argument("--where", default="p.embedding_status = 'pending'")
        args = ap.parse_args()

        devices = [int(x) for x in args.devices.split(",") if x.strip()]

        out_path = Path(args.out)
        out_path.parent.mkdir(parents=True, exist_ok=True)

        query = f"""
        MATCH (p:Protein)
        WHERE {args.where}
        RETURN p.id AS id, p.sequence AS sequence
        """

        db = get_async_driver()

        total_written = 0
        with out_path.open("ab") as f:  # append bytes
            async with ESM2Processor(
                devices=devices,
                max_padded_tokens=args.max_padded_tokens,
                max_length=args.max_length,
                dtype=args.dtype,
                progress=None,
                progress_task_id=None,
            ) as processor:
                chunks = stream_neo4j(db, query=query, batch_size=args.read_batch)

                async for record in processor.stream_work(chunks):
                    rec = {
                        "id": record.id,
                        "dtype": "float16",
                        "dim": int(record.vector.shape[-1]),
                        "embedding_b64": vec_to_b64_f16(record.vector),
                    }
                    f.write((json.dumps(rec) + "\n").encode("utf-8"))
                    total_written += 1
                    if total_written % 10000 == 0:
                        f.flush()

        await db.close()

    asyncio.run(main())
