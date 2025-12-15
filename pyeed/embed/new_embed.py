from __future__ import annotations

import asyncio
import time
from collections.abc import AsyncIterable, AsyncIterator, Callable, Iterable
from dataclasses import dataclass
from typing import Literal, NamedTuple

import numpy as np
import torch
from loguru import logger
from rich.progress import Progress, TaskID
from transformers import EsmModel, EsmTokenizer

from .pooling import PoolingFn, l2_normalize, mean_pooling
from .utils import _free_device_memory, _login_hf, silence_transformers_init_only

type DType = Literal["float16", "float32"]

_TORCH_DTYPE_MAP: dict[DType, torch.dtype] = {
    "float16": torch.float16,
    "float32": torch.float32,
}

_NP_DTYPE_MAP: dict[DType, type[np.floating]] = {
    "float16": np.float16,
    "float32": np.float32,
}


class Job[T](NamedTuple):
    """Represents one batched unit of work."""

    job_id: int
    payload: list[T]


@dataclass(slots=True)
class ESM2DeviceWorker:
    """Processes ESM2 jobs on a single device and emits results."""

    device: str
    in_q: asyncio.Queue[Job[tuple[str, str]] | None]
    out_q: asyncio.Queue[tuple[int, list[tuple[str, np.ndarray]]]]

    model: EsmModel
    tokenizer: EsmTokenizer
    torch_dtype: torch.dtype
    np_dtype: type[np.floating]
    pooling: PoolingFn | None
    normalize: bool
    max_length: int

    async def run(self) -> None:
        if self.device.startswith("cuda:"):
            torch.cuda.set_device(int(self.device.split(":")[1]))

        while True:
            job = await self.in_q.get()
            if job is None:
                self.in_q.task_done()
                break
            try:
                batch_results = await self._process_job(job.payload)
                await self.out_q.put((job.job_id, batch_results))
            finally:
                self.in_q.task_done()

    async def _process_job(
        self,
        batch: list[tuple[str, str]],
    ) -> list[tuple[str, np.ndarray]]:
        return await asyncio.to_thread(self.embed, batch)

    def embed(
        self,
        batch: list[tuple[str, str]],
    ) -> list[tuple[str, np.ndarray]]:
        t0 = time.time()
        ids, seqs = zip(*batch, strict=False)

        with torch.inference_mode():
            tokens = self.tokenizer(
                list(seqs),
                return_tensors="pt",
                truncation=False,
                max_length=self.max_length,
                padding=True,
                add_special_tokens=True,
            )

            input_ids = tokens["input_ids"].to(self.device)
            attention_mask = tokens["attention_mask"].to(self.device)

            use_autocast = self.device.startswith("cuda") and self.torch_dtype in (
                torch.float16,
                torch.bfloat16,
            )

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

        out = list(zip(ids, arr, strict=True))

        elapsed = time.time() - t0
        if elapsed > 5.0:
            logger.debug(
                "Embedding batch slow",
                extra={
                    "batch_size": len(batch),
                    "elapsed_s": round(elapsed, 3),
                    "device": self.device,
                },
            )
        return out


class ESM2Processor[T, R]:
    """Distributes length-sorted chunk batches across independent devices."""

    def __init__(
        self,
        model_name: str = "facebook/esm2_t33_650M_UR50D",
        dtype: DType = "float32",
        pooling: PoolingFn | None = mean_pooling,
        normalize: bool = True,
        max_length: int = 1024,
        huggingface_token: str | None = None,
        devices: list[int] | None = None,
        max_padded_tokens: int = 16_000,
        sort_desc: bool = True,
        per_device_queue_size: int = 3,
        source_chunk_size: int = 2048,
        length_key: Callable[[T], int] | None = None,
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
        self.per_device_queue_size = per_device_queue_size
        self.source_chunk_size = source_chunk_size

        self.length_key: Callable[[T], int] = (
            length_key if length_key is not None else (lambda item: len(item[1]))  # type: ignore[index, arg-type]
        )
        if progress is None and progress_task_id is not None:
            msg = "progress_task_id provided without a Progress instance"
            raise ValueError(msg)
        self.progress = progress
        self.progress_task_id = progress_task_id

        self._initialized = False
        self._devices: list[str] = []
        self._tokenizer: EsmTokenizer | None = None
        self._models: dict[str, EsmModel] = {}

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
        list[ESM2DeviceWorker[T, R]],
        list[asyncio.Queue[Job[T] | None]],
        asyncio.Queue[tuple[int, list[tuple[T, R]]]],
    ]:
        await self._initialize()

        out_q: asyncio.Queue[tuple[int, list[tuple[T, R]]]] = asyncio.Queue()
        device_queues: list[asyncio.Queue[Job[T] | None]] = []
        workers: list[ESM2DeviceWorker[T, R]] = []

        for device in self._devices:
            tokenizer = await asyncio.to_thread(self._load_tokenizer)
            dq: asyncio.Queue[Job[T] | None] = asyncio.Queue(self.per_device_queue_size)
            device_queues.append(dq)

            worker = ESM2DeviceWorker(
                device=device,
                in_q=dq,
                out_q=out_q,
                model=self._models[device],
                tokenizer=tokenizer,
                torch_dtype=self.torch_dtype,
                np_dtype=self.np_dtype,
                pooling=self.pooling,
                normalize=self.normalize,
                max_length=self.max_length,
            )
            workers.append(worker)

        return workers, device_queues, out_q

    def _pack_by_padded_budget(self, chunk: list[T]) -> list[list[T]]:
        """Packs a chunk into batches under a padded-token budget."""
        items = sorted(chunk, key=self.length_key, reverse=True)

        batches: list[list[T]] = []
        batch: list[T] = []
        max_len = 0

        for item in items:
            if len(item[1]) > self.max_length - 2:  # -2 for the special tokens
                # log that the item is too long with id
                logger.warning(
                    f"Item {item[0]} is too long and will be skipped",
                    extra={"id": item[0], "length": len(item)},
                )
                continue

            L = self.length_key(item)

            if not batch:
                batch = [item]
                max_len = L
                continue

            new_max = max(max_len, L)
            new_padded = (len(batch) + 1) * new_max

            if new_padded > self.max_padded_tokens:
                batches.append(batch)
                batch = [item]
                max_len = L
                continue

            batch.append(item)
            max_len = new_max

        if batch:
            batches.append(batch)

        return batches

    async def stream_work(
        self,
        chunks: AsyncIterable[list[T]],
    ) -> AsyncIterator[tuple[T, R]]:
        workers, device_queues, out_q = await self._build_workers()

        job_counter = 0
        processed_jobs = 0
        producer_done = asyncio.Event()

        async def producer() -> None:
            nonlocal job_counter
            job_id = 0

            async for chunk in chunks:
                if not chunk:
                    continue

                for batch in self._pack_by_padded_budget(chunk):
                    job = Job[T](job_id=job_id, payload=batch)
                    job_id += 1

                    target_q = min(device_queues, key=lambda q: q.qsize())
                    await target_q.put(job)
                    job_counter += 1

            for dq in device_queues:
                await dq.put(None)

            producer_done.set()

        async with asyncio.TaskGroup() as tg:
            for w in workers:
                tg.create_task(w.run())
            tg.create_task(producer())

            while True:
                if producer_done.is_set() and processed_jobs >= job_counter:
                    break

                _, batch_results = await out_q.get()
                processed_jobs += 1
                if self.progress is not None and self.progress_task_id is not None:
                    self.progress.advance(self.progress_task_id, len(batch_results))

                for inp, res in batch_results:
                    yield inp, res

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
    ) -> list[tuple[T, R]]:
        results: list[tuple[T, R]] = []
        async for inp, res in self.stream_work(self._to_chunk_stream(items)):
            results.append((inp, res))
        return results

    async def cleanup(self) -> None:
        device_ids: list[int] = []
        for dev in self._devices:
            if dev.startswith("cuda:"):
                device_ids.append(int(dev.split(":")[1]))

        self._models.clear()
        self._tokenizer = None
        self._devices.clear()
        self._initialized = False

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
    import time

    from rich import print

    from ..db.neo4j import get_async_driver

    print("Starting")
    db = get_async_driver()
    query = """
    MATCH (p:Protein)
    RETURN p.id, p.sequence
    LIMIT 4048
    """

    async def main() -> None:
        async with db.session() as session:
            print("Running query")
            res = await session.run(query)
            records = await res.data()
            # Output as list of tuple (convert list of dict to list of tuple)
            tuples = [(rec["p.id"], rec["p.sequence"]) for rec in records]
            total = len(tuples)

            # make list of tuple as async iterable
            async def aiter() -> AsyncIterator[tuple[str, str]]:
                for t in tuples:
                    yield t

            # initialize ESM2Processor
            print("Initializing processor")
            with Progress() as progress:
                task_id = progress.add_task("embedding", total=total)
                processor = ESM2Processor(
                    devices=[0, 2],
                    max_padded_tokens=30000,
                    dtype="float16",
                    progress=progress,
                    progress_task_id=task_id,
                )

                # run processor, timing the embedding only (exclude db query)
                print("Running processor")
                t0 = time.time()
                results = await processor.work(aiter())
                t1 = time.time()

            elapsed = t1 - t0
            rate = len(results) / elapsed if elapsed > 0 else float("inf")
            print("Processor completed")
            print(f"ID: {results[0][0]}")
            print(f"Embedding: {results[0][1]}")
            print(f"Embedded {len(results)} sequences in {elapsed:.3f} seconds, {rate:.2f} seq/s")

            print("Query completed")

    asyncio.run(main())


# BATCH_STATUS = 1000
# BATCH_VECTORS = 1000


# async def set_status(session, ids: list[str], status: str) -> None:
#     if not ids:
#         return
#     q = """
#     UNWIND $ids AS id
#     MATCH (p:Protein {id: id})
#     SET p.embedding_status = $status
#     """
#     await session.run(q, ids=ids, status=status)


# async def write_vectors_bulk(items: list[tuple[str, np.ndarray]]) -> None:
#     # Replace with your real vector DB bulk upsert.
#     await asyncio.to_thread(lambda: None)


# async with db.session() as session:
#     processor = ESM2Processor(devices=[0, 2], max_padded_tokens=30000, max_length=1024)

#     # assume you already:
#     # 1) fetched tuples
#     # 2) split too_long vs valid
#     # 3) marked valid in_progress, too_long failed

#     async def valid_iter(valid: list[tuple[str, str]]):
#         for x in valid:
#             yield x

#     completed_ids: list[str] = []
#     vector_buf: list[tuple[str, np.ndarray]] = []

#     async for pid_seq, emb in processor.stream_work(
#         processor._to_chunk_stream(valid_iter(valid))
#     ):
#         # pid_seq is your input item (id, seq)
#         pid = pid_seq[0]

#         vector_buf.append((pid, emb))
#         completed_ids.append(pid)

#         if len(vector_buf) >= BATCH_VECTORS:
#             await write_vectors_bulk(vector_buf)
#             vector_buf.clear()

#         if len(completed_ids) >= BATCH_STATUS:
#             await set_status(session, completed_ids, "completed")
#             completed_ids.clear()

#     # flush remainder
#     if vector_buf:
#         await write_vectors_bulk(vector_buf)

#     if completed_ids:
#         await set_status(session, completed_ids, "completed")
