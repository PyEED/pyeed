import asyncio
from collections.abc import AsyncIterator, Sequence

import numpy as np
import torch
from loguru import logger
from numpy.typing import NDArray
from rich.progress import (
    BarColumn,
    MofNCompleteColumn,
    Progress,
    ProgressColumn,
    SpinnerColumn,
    TaskProgressColumn,
    TextColumn,
    TimeRemainingColumn,
)
from rich.text import Text
from transformers import EsmModel, EsmTokenizer

from ..logging_setup import setup_logging
from .pooling import PoolingFn, PoolingLike, l2_normalize, mean_pooling, normalize_cast_renorm
from .types import (
    NP_DTYPE_MAP,
    TF_DTYPE_MAP,
    ModelDType,
    ReturnDType,
)
from .utils import _free_device_memory, _login_hf, silence_transformers_init_only

setup_logging()

__all__ = [
    "ESM2Embedder",
    "EmbeddingResult",
]

type Vec16 = NDArray[np.float16]
type Vec32 = NDArray[np.float32]
type Vec = Vec16 | Vec32
type EmbeddingResult = dict[str, str | Vec]


class SpeedColumn(ProgressColumn):
    """Custom column that displays sequences per second, handling None values."""

    def render(self, task) -> Text:
        """Render speed with None handling."""
        if task.speed is not None and task.speed > 0:
            return Text(f"{task.speed:.1f} seq/s", style="cyan")
        return Text("0.0 seq/s", style="dim")


class ESM2Embedder:
    def __init__(
        self,
        model_name: str = "facebook/esm2_t33_650M_UR50D",
        pooling_methods: PoolingLike = None,
        normalize: bool = True,
        model_dtype: ModelDType = "float32",
        return_dtype: ReturnDType = "float32",
        max_length: int = 1024,
        verbose: bool = False,
    ):
        """
        Initialize ESM2 embedder.

        Args:
            model_name: ESM-2 model identifier
            pooling_methods: Single or multiple pooling functions
            normalize: L2-normalize embeddings
            model_dtype: Model precision for computation (float32/bfloat16/float16)
            return_dtype: Output array dtype (float32/bfloat16/float16)
            max_length: Maximum sequence length
            verbose: Show progress bar with rich (sequences/sec rate)
        """
        self.model_name = model_name
        self.normalize = normalize
        self.model_dtype = model_dtype
        self.return_dtype = return_dtype
        self.max_length = max_length
        self.verbose = verbose
        self._initialized = False
        self._tf_model_dtype = TF_DTYPE_MAP[model_dtype]
        self._tf_return_dtype = TF_DTYPE_MAP[return_dtype]
        self._np_return_dtype = NP_DTYPE_MAP[return_dtype]

        # Normalize pooling_methods into a list of (name, function) tuples
        self.pooling_configs = self._normalize_pooling_methods(pooling_methods)

        self.models: list[EsmModel | None] = []
        self.tokenizer: EsmTokenizer | None = None  # Single tokenizer on CPU
        self.devices: list[torch.device] = []
        self.token: str | None = None
        self.device_ids: list[int] = []

        logger.info(
            "ESM2 embedder initialized",
            extra={
                "model_name": model_name,
                "normalize": normalize,
                "model_dtype": model_dtype,
                "return_dtype": return_dtype,
                "max_length": max_length,
                "verbose": verbose,
                "pooling_methods": [name for name, _ in self.pooling_configs],
            },
        )

    def _normalize_pooling_methods(
        self, pooling_methods: PoolingLike
    ) -> list[tuple[str, PoolingFn | None]]:
        """
        Convert PoolingLike input into list of (name, function) tuples.

        Returns:
            List of (pooling_name, pooling_function) tuples
        """
        if pooling_methods is None:
            return [("none", None)]

        # Single pooling function
        if callable(pooling_methods):
            name = getattr(pooling_methods, "__name__", "custom")
            return [(name, pooling_methods)]

        # Sequence of pooling functions
        if isinstance(pooling_methods, Sequence):
            result = []
            for idx, method in enumerate(pooling_methods):
                if method is None:
                    result.append(("none", None))
                elif callable(method):
                    name = getattr(method, "__name__", f"custom_{idx}")
                    result.append((name, method))
                else:
                    raise ValueError(f"Invalid pooling method at index {idx}: {method}")
            return result

        raise ValueError(f"Invalid pooling_methods type: {type(pooling_methods)}")

    async def initialize(self) -> None:
        """Load models on all available GPUs and tokenizer on CPU."""
        if self._initialized:
            logger.debug("ESM2 embedder already initialized")
            return

        self.token = _login_hf()

        self.devices = self._detect_devices()[0:2]

        # Pre-allocate lists to avoid race conditions with asyncio.gather
        self.models = [None] * len(self.devices)

        # Load single tokenizer on CPU (thread-safe, no device affinity)
        self.tokenizer = await asyncio.to_thread(self._load_tokenizer)

        # Load models on all devices in parallel
        await asyncio.gather(*[self._load_on_device(i, d) for i, d in enumerate(self.devices)])

        self._initialized = True
        logger.info("ESM2 embedder ready")

    def _detect_devices(self) -> list[torch.device]:
        """Detect available CUDA devices or fallback to CPU."""
        if not torch.cuda.is_available():
            logger.debug("No cuda available, using cpu")
            return [torch.device("cpu")]

        # Initialize device_ids if not set
        if not self.device_ids:
            self.device_ids = list(range(torch.cuda.device_count()))

        devices = [torch.device(f"cuda:{i}") for i in self.device_ids]

        logger.info(
            "Detected devices",
            extra={
                "num_devices": len(self.device_ids),
                "devices": [str(d) for d in devices],
            },
        )

        return devices

    def _load_tokenizer(self) -> EsmTokenizer:
        """Load tokenizer on CPU (called once, thread-safe)."""
        with silence_transformers_init_only():
            tokenizer = EsmTokenizer.from_pretrained(
                self.model_name,
                token=self.token,
            )
        logger.info("ESM2 tokenizer loaded on CPU")
        return tokenizer

    async def _load_on_device(self, device_idx: int, device: torch.device) -> None:
        """Load model on specific device with specified dtype."""

        def _load() -> EsmModel:
            # Set CUDA device context in this thread to ensure tensors use correct device
            if device.type == "cuda":
                torch.cuda.set_device(device.index)

            # Load model with specified dtype
            with silence_transformers_init_only():
                model = EsmModel.from_pretrained(
                    self.model_name,
                    token=self.token,
                    torch_dtype=self._tf_model_dtype,
                )

            # CRITICAL: Disable all-layer hidden states to save memory
            # ESM2-T33 has ~33 layers; materializing all = ~33x memory per forward pass
            model.config.output_hidden_states = False

            model = model.to(device=device, dtype=self._tf_model_dtype).eval()

            return model

        model = await asyncio.to_thread(_load)
        # Use indexed assignment instead of append to avoid race conditions
        self.models[device_idx] = model
        logger.info(
            "ESM2 model loaded on device",
            extra={
                "device": str(device),
                "model_dtype": self.model_dtype,
                "return_dtype": self.return_dtype,
            },
        )

    async def cleanup(self) -> None:
        """Clean up models and free GPU memory on all used devices."""
        logger.info("Cleaning up ESM2 embedder")

        # Store device IDs before clearing
        device_ids = self.device_ids if self.device_ids else None

        self.models.clear()
        self.tokenizer = None
        self.devices.clear()
        self._initialized = False

        # Free memory on all devices that were used
        await asyncio.to_thread(_free_device_memory, device_ids)
        logger.info("ESM2 embedder cleanup complete")

    # ========================================================================
    # Public API
    # ========================================================================

    def _create_length_sorted_batches(
        self,
        sequences: list[str],
        accessions: list[str],
        batch_size: int,
    ) -> list[tuple[list[str], list[str]]]:
        """
        Sort sequences by length and create batches for efficient GPU processing.

        Sorting by length (descending) minimizes padding within each batch,
        improving GPU efficiency. Batches are pre-formed before distribution.

        Args:
            sequences: Protein sequences
            accessions: Accession IDs
            batch_size: Number of sequences per batch

        Returns:
            List of (batch_sequences, batch_accessions) tuples
        """
        # Create pairs and sort by sequence length (descending - longest first)
        seq_acc_pairs = list(zip(sequences, accessions, strict=True))
        sorted_pairs = sorted(seq_acc_pairs, key=lambda x: len(x[0]), reverse=True)

        # Create batches
        batches = []
        for i in range(0, len(sorted_pairs), batch_size):
            batch_pairs = sorted_pairs[i : i + batch_size]
            batch_seqs = [seq for seq, _ in batch_pairs]
            batch_accs = [acc for _, acc in batch_pairs]
            batches.append((batch_seqs, batch_accs))

        return batches

    async def embed_batch(
        self,
        sequences: list[str],
        accessions: list[str] | None = None,
        batch_size: int = 16,
    ) -> EmbeddingResult:
        """
        Generate embeddings for multiple sequences using work queue for load balancing.

        Args:
            sequences: Protein sequences
            accessions: Optional IDs (auto-generated if None)
            batch_size: Sequences per batch per GPU

        Returns:
            Dict with "ids" key (list[str]) and pooling method keys (list[EmbeddingArray])
            Example: {"ids": ["P001", "P002", ...], "mean_pooling": [arr1, arr2, ...]}
        """
        if not self._initialized:
            await self.initialize()

        accessions = accessions or [f"seq_{i}" for i in range(len(sequences))]
        if len(sequences) != len(accessions):
            raise ValueError("sequences and accessions must have same length")

        logger.info(
            "Embedding sequences with work queue",
            extra={
                "num_sequences": len(sequences),
                "num_devices": len(self.devices),
                "batch_size": batch_size,
            },
        )

        # Sort sequences by length and create pre-formed batches
        batches = self._create_length_sorted_batches(sequences, accessions, batch_size)

        logger.debug(
            "Created batches",
            extra={
                "num_batches": len(batches),
                "sequences_per_batch": batch_size,
            },
        )

        # Create work queue and populate with batches
        batch_queue: asyncio.Queue[tuple[list[str], list[str]] | None] = asyncio.Queue()
        for batch in batches:
            await batch_queue.put(batch)

        # Add sentinel None values to signal workers to stop (one per device)
        for _ in range(len(self.devices)):
            await batch_queue.put(None)

        # Launch workers with optional progress bar
        if self.verbose:
            # Create rich progress bar with custom columns
            with Progress(
                SpinnerColumn(),
                TextColumn("[progress.description]{task.description}"),
                BarColumn(),
                TaskProgressColumn(),
                MofNCompleteColumn(),
                TextColumn("•"),
                SpeedColumn(),
                TimeRemainingColumn(),
            ) as progress:
                task_id = progress.add_task(
                    "[green]Embedding sequences...",
                    total=len(sequences),
                )

                # Launch workers with progress tracking
                workers = [
                    self._process_queue_worker(device_idx, batch_queue, progress, task_id)
                    for device_idx in range(len(self.devices))
                ]

                # Wait for all workers to complete
                results = await asyncio.gather(*workers)
        else:
            # Launch workers without progress bar
            workers = [
                self._process_queue_worker(device_idx, batch_queue)
                for device_idx in range(len(self.devices))
            ]

            # Wait for all workers to complete
            results = await asyncio.gather(*workers)

        # Concatenate results from all workers
        final_result: EmbeddingResult = {name: [] for name, _ in self.pooling_configs}
        final_result["ids"] = []

        for worker_result in results:
            for key, values in worker_result.items():
                final_result[key].extend(values)

        logger.info(
            "Generated embeddings",
            extra={
                "num_embeddings": len(final_result["ids"]),
                "num_batches": len(batches),
            },
        )

        return final_result

    async def embed_single(self, sequence: str, accession: str) -> EmbeddingResult:
        """Embed a single sequence. Returns dict with same structure as embed_batch."""
        return await self.embed_batch([sequence], [accession], batch_size=1)

    async def embed_stream(
        self,
        sequences: list[str],
        accessions: list[str],
        batch_size: int = 8,
        chunk_size: int = 1000,
    ) -> AsyncIterator[EmbeddingResult]:
        """
        Stream embeddings in chunks for large datasets.

        Yields batches of embeddings as they complete, enabling
        immediate database writes without waiting for all sequences.

        Args:
            sequences: Protein sequences
            accessions: Accession IDs
            batch_size: Sequences per batch per GPU
            chunk_size: Sequences per yielded chunk

        Yields:
            Dict mapping pooling method to list of embeddings for the chunk
        """
        if not self._initialized:
            await self.initialize()

        if len(sequences) != len(accessions):
            raise ValueError("sequences and accessions must have same length")

        logger.info(
            "Streaming sequences",
            extra={
                "num_sequences": len(sequences),
                "chunk_size": chunk_size,
                "batch_size": batch_size,
            },
        )

        for chunk_start in range(0, len(sequences), chunk_size):
            chunk_end = min(chunk_start + chunk_size, len(sequences))
            chunk_seqs = sequences[chunk_start:chunk_end]
            chunk_accs = accessions[chunk_start:chunk_end]

            # Await the batch computation
            embeddings = await self.embed_batch(chunk_seqs, chunk_accs, batch_size)

            logger.info(
                "Completed chunk",
                extra={
                    "chunk_start": chunk_start,
                    "chunk_end": chunk_end,
                    "total": len(sequences),
                },
            )

            # Yield results immediately
            yield embeddings

    def _distribute_work(
        self,
        sequences: list[str],
        accessions: list[str],
    ) -> list[tuple[list[str], list[str]]]:
        """
        Distribute sequences round-robin across devices.

        Example with 3 devices and 10 sequences:
        - Device 0: sequences [0, 3, 6, 9]
        - Device 1: sequences [1, 4, 7]
        - Device 2: sequences [2, 5, 8]

        Note: Order is reconstructed in _reorder_results using the same pattern.
        """
        num_devices = len(self.devices)
        return [(sequences[i::num_devices], accessions[i::num_devices]) for i in range(num_devices)]

    def _reorder_results(
        self, device_results: list[EmbeddingResult], total: int
    ) -> EmbeddingResult:
        """
        Reorder round-robin results back to original sequence order.

        Args:
            device_results: List of dicts from each device, each containing embeddings
                           in the order they were processed by that device
            total: Total number of sequences

        Returns:
            Dict with reordered embeddings matching input order

        Example:
            3 devices, 10 sequences distributed [0,3,6,9], [1,4,7], [2,5,8]
            device_results = [
                {"mean": [emb0, emb3, emb6, emb9]},  # Device 0
                {"mean": [emb1, emb4, emb7]},         # Device 1
                {"mean": [emb2, emb5, emb8]},         # Device 2
            ]
            Result: {"mean": [emb0, emb1, emb2, emb3, emb4, emb5, emb6, emb7, emb8, emb9]}
        """
        if not device_results:
            return {}

        pooling_names = list(device_results[0].keys())
        result: EmbeddingResult = {name: [] for name in pooling_names}

        for pool_name in pooling_names:
            # Collect from each device (already in correct order per device)
            device_values = [device_result[pool_name] for device_result in device_results]

            # Reconstruct original order using round-robin pattern
            ordered: list = []
            for seq_idx in range(total):
                device_idx = seq_idx % len(self.devices)
                device_seq_idx = seq_idx // len(self.devices)

                # Check if this device has this sequence
                if device_seq_idx < len(device_values[device_idx]):
                    ordered.append(device_values[device_idx][device_seq_idx])

            result[pool_name] = ordered

        return result

    async def _process_queue_worker(
        self,
        device_idx: int,
        batch_queue: asyncio.Queue[tuple[list[str], list[str]] | None],
        progress: Progress | None = None,
        task_id: int | None = None,
    ) -> EmbeddingResult:
        """
        Worker that pulls pre-formed batches from a queue and processes them.

        This implements dynamic load balancing: faster GPUs automatically
        process more batches. Continues until None sentinel is received.

        Args:
            device_idx: Index of GPU device to use
            batch_queue: Queue containing (sequences, accessions) tuples
            progress: Optional Rich Progress instance for visual feedback
            task_id: Optional Progress task ID to update

        Returns:
            Accumulated results with all pooling methods and IDs
        """
        # Initialize accumulator for each pooling method
        accumulated: EmbeddingResult = {name: [] for name, _ in self.pooling_configs}
        accumulated["ids"] = []

        logger.debug(
            "Worker starting",
            extra={
                "device_idx": device_idx,
                "device": str(self.devices[device_idx]),
            },
        )

        batch_count = 0
        while True:
            # Get next batch from queue
            batch_item = await batch_queue.get()

            # None sentinel signals end of work
            if batch_item is None:
                batch_queue.task_done()
                break

            batch_seqs, batch_accs = batch_item

            try:
                batch_result = await asyncio.to_thread(
                    self._compute_embeddings,
                    device_idx,
                    batch_seqs,
                    batch_accs,
                )
            except torch.cuda.OutOfMemoryError as e:
                logger.error(
                    "GPU OOM error on device",
                    extra={
                        "device_idx": device_idx,
                        "device": str(self.devices[device_idx]),
                        "batch_size": len(batch_seqs),
                        "error": str(e),
                    },
                )
                batch_queue.task_done()
                raise RuntimeError(
                    f"GPU out of memory on device {device_idx}. "
                    f"Try reducing batch_size or using dtype=torch.bfloat16"
                ) from e
            except Exception as e:
                logger.error(
                    "Unexpected error during embedding computation",
                    extra={
                        "device_idx": device_idx,
                        "device": str(self.devices[device_idx]),
                        "batch_count": batch_count,
                        "error_type": type(e).__name__,
                        "error": str(e),
                    },
                )
                batch_queue.task_done()
                raise

            # Accumulate results
            for key, values in batch_result.items():
                accumulated[key].extend(values)

            batch_count += 1
            batch_queue.task_done()

            # Update progress bar with number of sequences processed in this batch
            if progress is not None and task_id is not None:
                progress.update(task_id, advance=len(batch_seqs))

            logger.debug(
                "Worker progress",
                extra={
                    "device_idx": device_idx,
                    "batches_completed": batch_count,
                    "total_sequences": len(accumulated["ids"]),
                },
            )

        logger.debug(
            "Worker finished",
            extra={
                "device_idx": device_idx,
                "total_batches": batch_count,
                "total_sequences": len(accumulated["ids"]),
            },
        )

        return accumulated

    async def _process_device_work(
        self,
        device_idx: int,
        sequences: list[str],
        accessions: list[str],
        batch_size: int,
    ) -> EmbeddingResult:
        """Process all sequences assigned to a device."""
        if not sequences:
            # Return empty dict with all pooling methods and IDs
            result = {name: [] for name, _ in self.pooling_configs}
            result["ids"] = []
            return result

        logger.debug(
            "Processing device work",
            extra={
                "device_idx": device_idx,
                "device": str(self.devices[device_idx]),
                "num_sequences": len(sequences),
                "batch_size": batch_size,
            },
        )

        # Initialize accumulator for each pooling method
        accumulated: EmbeddingResult = {name: [] for name, _ in self.pooling_configs}
        accumulated["ids"] = []

        for batch_start in range(0, len(sequences), batch_size):
            batch_end = min(batch_start + batch_size, len(sequences))
            batch_seqs = sequences[batch_start:batch_end]
            batch_accs = accessions[batch_start:batch_end]

            try:
                batch_embs = await asyncio.to_thread(
                    self._compute_embeddings,
                    device_idx,
                    batch_seqs,
                    batch_accs,
                )
            except torch.cuda.OutOfMemoryError as e:
                logger.error(
                    "GPU OOM error on device",
                    extra={
                        "device_idx": device_idx,
                        "device": str(self.devices[device_idx]),
                        "batch_size": len(batch_seqs),
                        "error": str(e),
                    },
                )
                raise RuntimeError(
                    f"GPU out of memory on device {device_idx}. "
                    f"Try reducing batch_size or using dtype=torch.bfloat16"
                ) from e
            except Exception as e:
                logger.error(
                    "Unexpected error during embedding computation",
                    extra={
                        "device_idx": device_idx,
                        "device": str(self.devices[device_idx]),
                        "batch_start": batch_start,
                        "batch_end": batch_end,
                        "error_type": type(e).__name__,
                        "error": str(e),
                    },
                )
                raise

            # Accumulate embeddings for each pooling method
            for pool_name, embeddings_list in batch_embs.items():
                accumulated[pool_name].extend(embeddings_list)

            logger.debug(
                "Device progress",
                extra={
                    "device_idx": device_idx,
                    "completed": len(next(iter(accumulated.values()))),
                    "total": len(sequences),
                },
            )

        return accumulated

    def _compute_embeddings(
        self,
        device_idx: int,
        sequences: Sequence[str],
        accessions: Sequence[str],
    ) -> EmbeddingResult:
        """
        Compute embeddings for a batch (GPU-bound, runs in thread pool).

        This is the core computation that runs synchronously in a thread.
        Applies all configured pooling methods to the hidden states.

        Args:
            device_idx: Index of GPU device to use
            sequences: Batch of protein sequences
            accessions: Accession IDs (currently unused but kept for API compatibility)

        Returns:
            Dict mapping pooling method name to list of embeddings (one per sequence)
            Each embedding is a numpy array of shape [hidden_dim] or [seq_len * hidden_dim]

        Raises:
            RuntimeError: If model/tokenizer not loaded or device mismatch
            torch.cuda.OutOfMemoryError: If GPU runs out of memory (will propagate to caller)
        """
        model = self.models[device_idx]
        tokenizer = self.tokenizer
        device = self.devices[device_idx]

        # Ensure model and tokenizer were properly loaded
        if model is None or tokenizer is None:
            raise RuntimeError(
                f"Model or tokenizer not loaded for device_idx {device_idx}. "
                f"Ensure initialize() was called successfully."
            )

        # CRITICAL: Set CUDA device context BEFORE any tensor operations
        if device.type == "cuda":
            torch.cuda.set_device(device.index)

        # Sanity check: verify model parameters are on the expected device
        param_dev = next(model.parameters()).device
        if param_dev != device:
            raise RuntimeError(
                f"Model on {param_dev}, but assigned device is {device} (idx {device_idx})"
            )

        # Calculate effective max length for this batch to avoid excessive padding
        # +2 accounts for BOS/EOS tokens added by tokenizer
        batch_max_len = max(len(s) for s in sequences) + 2
        # Round up to nearest multiple of 8 for optimal GPU tensor core utilization
        batch_max_len = ((batch_max_len + 7) // 8) * 8
        eff_max_len = min(self.max_length, batch_max_len)

        # Use autocast for fp16/bf16 to reduce memory footprint
        use_autocast = (device.type == "cuda") and (
            self._tf_model_dtype in (torch.float16, torch.bfloat16)
        )

        # inference_mode() is more efficient than no_grad() for inference
        with torch.inference_mode():
            # Batch tokenize - padding to longest sequence in this batch
            inputs = tokenizer(
                list(sequences),
                return_tensors="pt",
                padding=True,  # Pads to longest in batch, not max_length
                truncation=True,
                max_length=eff_max_len,
            )

            # Pin CPU tensors to enable true async H→D transfers
            # Without pinning, non_blocking=True falls back to sync copy
            for k in inputs:
                if inputs[k].device.type == "cpu":
                    inputs[k] = inputs[k].pin_memory()

            # Move tensors to target device (now truly async thanks to pinning)
            inputs = {k: v.to(device, non_blocking=True) for k, v in inputs.items()}

            # Forward pass - model.config.output_hidden_states=False set at load time
            if use_autocast:
                with torch.autocast(device_type=device.type, dtype=self._tf_model_dtype):
                    outputs = model(**inputs)
            else:
                outputs = model(**inputs)

            # Extract last hidden state only (no intermediate layers)
            hidden_states = outputs.last_hidden_state  # [B, L, D]
            attention_mask = inputs.get("attention_mask")  # [B, L] or None

            # Apply all pooling methods
            result: EmbeddingResult = {}

            for pool_name, pool_func in self.pooling_configs:
                if pool_func is None:
                    # No pooling: return raw hidden states [batch, seq_len, hidden_dim]
                    raw_np = hidden_states.cpu().numpy().astype(np.float32)
                    result[pool_name] = raw_np.astype(self.return_dtype)

                else:
                    # Apply pooling function
                    pooled = pool_func(hidden_states, attention_mask)  # [B, D]

                    if self.normalize:
                        if self.return_dtype == "float32":
                            pooled = l2_normalize(pooled)
                        else:
                            pooled = normalize_cast_renorm(pooled, self._tf_return_dtype)

            result[pool_name] = pooled.cpu().numpy().astype(self._np_return_dtype)

            # Add IDs to maintain sequence-embedding correspondence
            result["ids"] = list(accessions)

            # Explicitly delete GPU tensors to free memory immediately
            del inputs, outputs, hidden_states
            if "pooled" in locals():
                del pooled

        # Clear GPU cache after batch processing
        if device.type == "cuda":
            torch.cuda.empty_cache()

        return result

    async def embed_from_generator(
        self,
        seq_generator: AsyncIterator[tuple[list[str], list[str]]],
        batch_size: int = 8,
    ) -> AsyncIterator[EmbeddingResult]:
        """
        Process sequences from an async generator and stream embeddings.

        This method accepts an async generator that yields (sequences, accessions) tuples
        and processes them incrementally. Results are yielded as they complete, enabling
        immediate downstream processing (e.g., database insertion) without waiting for
        all sequences to be embedded.

        Args:
            seq_generator: Async generator yielding (sequences, accessions) tuples
            batch_size: Sequences per batch per GPU

        Yields:
            Dict mapping pooling method to list of embeddings for each chunk

        Example:
            async def load_sequences():
                # Load from file, database, etc.
                for chunk in load_chunks():
                    yield chunk_seqs, chunk_accs

            async with ESM2Embedder(pooling_methods=[mean_pooling]) as embedder:
                async for batch_embeddings in embedder.embed_from_generator(load_sequences()):
                    # batch_embeddings is {"mean_pooling": [array1, array2, ...]}
                    await db.insert_embeddings(batch_embeddings)
        """
        if not self._initialized:
            await self.initialize()

        logger.info(
            "Starting generator-based embedding stream",
            extra={
                "batch_size": batch_size,
                "num_devices": len(self.devices),
            },
        )

        chunk_count = 0
        async for sequences, accessions in seq_generator:
            if len(sequences) != len(accessions):
                raise ValueError(
                    f"Chunk {chunk_count}: sequences and accessions must have same length "
                    f"(got {len(sequences)} vs {len(accessions)})"
                )

            logger.debug(
                "Processing generator chunk",
                extra={
                    "chunk_idx": chunk_count,
                    "num_sequences": len(sequences),
                },
            )

            # Process this chunk through the standard embed_batch pipeline
            embeddings = await self.embed_batch(
                list(sequences),
                list(accessions),
                batch_size=batch_size,
            )

            chunk_count += 1
            # Count embeddings from first pooling method
            num_embeddings = len(next(iter(embeddings.values()))) if embeddings else 0
            logger.debug(
                "Completed generator chunk",
                extra={
                    "chunk_idx": chunk_count - 1,
                    "num_embeddings": num_embeddings,
                },
            )

            # Yield results immediately for downstream processing
            yield embeddings

        logger.info(
            "Completed generator-based embedding",
            extra={
                "total_chunks": chunk_count,
            },
        )

    async def __aenter__(self) -> "ESM2Embedder":
        await self.initialize()
        return self

    async def __aexit__(self, exc_type, exc_val, exc_tb) -> None:  # type: ignore
        await self.cleanup()


if __name__ == "__main__":
    path = "/home/mha/projects/proteingraph/downloads/1000seq.fasta"

    # read fasta and get sequences and accessions (sequences are multiline)
    sequences = []
    accessions = []
    sequence = ""
    accession = ""
    with open(path) as f:
        for idx, line in enumerate(f):
            if line.startswith(">"):
                if idx != 0:
                    sequences.append(sequence)
                    accessions.append(accession)
                accession = line.strip().split()[0][1:]
                sequence = ""
            else:
                sequence += line.strip()
        # Don't forget last sequence
        if sequence:
            sequences.append(sequence)
            accessions.append(accession)

    # Example: Float16 (recommended for Milvus)
    print("\n=== Float16 Embeddings (Milvus-ready) ===")
    embedder = ESM2Embedder(
        model_dtype="float16",
        return_dtype="float16",
        pooling_methods=[mean_pooling],
        verbose=True,  # Show progress bar with sequences/sec
    )
    init = asyncio.run(embedder.initialize())

    result = asyncio.run(
        embedder.embed_batch(
            sequences,
            accessions,
            batch_size=48,
        )
    )

    for pool_name, values in result.items():
        if pool_name == "ids":
            print(f"ids: {len(values)} IDs")
            print(f"  Sample: {values[:3]}")
        else:
            print(f"{pool_name}: {len(values)} embeddings")
            print(f"  Shape: {values[0].shape}")
            print(f"  Dtype: {values[0].dtype}")
            print(f"  Sample: {values[0][:5]}")
