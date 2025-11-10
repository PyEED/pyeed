import asyncio
from collections import defaultdict
from collections.abc import AsyncIterator, Sequence

import torch
from loguru import logger
from rich.progress import Progress, TaskID
from transformers import EsmModel, EsmTokenizer

from .pooling import PoolingFn, PoolingLike, l2_normalize, normalize_cast_renorm
from .types import (
    NP_DTYPE_MAP,
    TF_DTYPE_MAP,
    EmbeddingBatch,
    ModelDType,
    ReturnDType,
)
from .utils import _free_device_memory, _login_hf, silence_transformers_init_only

__all__ = [
    "ESM2Embedder",
]

# Timing statistics
_timing_stats: dict[str, list[float]] = defaultdict(list)


class ESM2Embedder:
    def __init__(
        self,
        model_name: str = "facebook/esm2_t33_650M_UR50D",
        pooling_methods: PoolingLike = None,
        normalize: bool = True,
        model_dtype: ModelDType = "float32",
        return_dtype: ReturnDType = "float32",
        n_gpus: int = -1,
        max_length: int = 1024,
        verbose: bool = False,
        max_batch_queue_size: int | None = 32,
    ):
        """
        Initialize ESM2 embedder.

        Args:
            model_name: ESM-2 model identifier
            pooling_methods: Single or multiple pooling functions
            normalize: L2-normalize embeddings
            model_dtype: Model precision for computation (float32/float16)
            return_dtype: Output array dtype (float32/float16)
            n_gpus: Number of GPUs to use (-1 for all available)
            max_length: Maximum sequence length
            verbose: Show progress bar with rich (sequences/sec rate)
            max_batch_queue_size: Max batches in internal queues (default: n_gpus * 2)
        """
        self.model_name = model_name
        self.normalize = normalize
        self.model_dtype = model_dtype
        self.return_dtype = return_dtype
        self.max_length = max_length
        self.verbose = verbose
        self.n_gpus = n_gpus
        self._initialized = False
        self._tf_model_dtype = TF_DTYPE_MAP[model_dtype]
        self._tf_return_dtype = TF_DTYPE_MAP[return_dtype]
        self._np_return_dtype = NP_DTYPE_MAP[return_dtype]

        # Calculate queue size based on GPUs (will be finalized in initialize)
        self._max_queue_size = max_batch_queue_size

        # Normalize pooling_methods into a list of (name, function) tuples
        self.pooling_configs = self._normalize_pooling_methods(pooling_methods)

        self.models: list[EsmModel | None] = []
        self.tokenizer: EsmTokenizer | None = None  # Single tokenizer on CPU
        self.devices: list[torch.device] = []
        self.token: str | None = None
        self.device_ids: list[int] = []

        # Persistent workers for continuous processing
        self._batch_queue: asyncio.Queue[tuple[list[str], list[str]] | None] | None = None
        self._result_queue: asyncio.Queue[EmbeddingBatch | None] | None = None
        self._worker_tasks: list[asyncio.Task] | None = None

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

        if self.n_gpus == -1:
            self.devices = self._detect_devices()
        else:
            self.devices = [torch.device(f"cuda:{i}") for i in range(self.n_gpus)]

        self.models = [None] * len(self.devices)

        # Load single tokenizer on CPU (thread-safe, no device affinity)
        self.tokenizer = await asyncio.to_thread(self._load_tokenizer)

        # Load models on all devices in parallel
        await asyncio.gather(*[self._load_on_device(i, d) for i, d in enumerate(self.devices)])

        # Finalize queue size based on actual number of devices
        if self._max_queue_size is None:
            self._max_queue_size = len(self.devices) * 2

        # Create BOUNDED persistent queues to enable backpressure
        self._batch_queue = asyncio.Queue(maxsize=self._max_queue_size)
        self._result_queue = asyncio.Queue(maxsize=self._max_queue_size)

        # Launch persistent workers (one per GPU)
        self._worker_tasks = [
            asyncio.create_task(self._persistent_worker(device_idx))
            for device_idx in range(len(self.devices))
        ]

        self._initialized = True
        logger.info(
            "ESM2 embedder ready with persistent workers",
            extra={"max_queue_size": self._max_queue_size, "num_workers": len(self.devices)},
        )

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

        # Shutdown persistent workers
        if self._worker_tasks:
            # Send shutdown sentinels (one per worker)
            for _ in range(len(self.devices)):
                await self._batch_queue.put(None)

            # Wait for workers to finish
            await asyncio.gather(*self._worker_tasks, return_exceptions=True)

            self._worker_tasks = None
            self._batch_queue = None
            self._result_queue = None

        # Store device IDs before clearing
        device_ids = self.device_ids if self.device_ids else None

        self.models.clear()
        self.tokenizer = None
        self.devices.clear()
        self._initialized = False

        # Free memory on all devices that were used
        await asyncio.to_thread(_free_device_memory, device_ids)
        logger.info("ESM2 embedder cleanup complete")

    @property
    def embedding_dim(self) -> int:
        """Get embedding dimension from model config.

        Returns:
            Hidden size of the model (embedding dimension)

        Raises:
            RuntimeError: If model not initialized
        """
        if not self._initialized or not self.models or self.models[0] is None:
            raise RuntimeError("Model not initialized. Call initialize() first.")
        return self.models[0].config.hidden_size

    # ========================================================================
    # Public API
    # ========================================================================

    def create_length_sorted_batches(
        self,
        sequences: list[str],
        accessions: list[str],
        batch_size: int,
    ) -> list[tuple[list[str], list[str]]]:
        """Sort sequences by length and create batches for efficient GPU processing.

        Sorting by length (descending) minimizes padding within each batch,
        improving GPU efficiency. Batches are pre-formed before distribution.

        This is a pure function with no side effects - can be called from
        pipeline to enable overlapping batch creation while embeddings are
        being computed.

        Args:
            sequences: Protein sequences
            accessions: Accession IDs
            batch_size: Number of sequences per batch

        Returns:
            List of (batch_sequences, batch_accessions) tuples

        Example:
            # Create batches externally in pipeline
            batches = embedder.create_length_sorted_batches(sequences, accessions, 32)
            # Then process them
            async for batch in embedder.embed_prepared_batches(batches):
                await db.insert(batch)
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
        accessions: list[str],
        batch_size: int,
        progress: Progress | None = None,
        task_id: TaskID | None = None,
        prepare_task_id: TaskID | None = None,
    ) -> AsyncIterator[EmbeddingBatch]:
        """Create batches and embed sequences (convenience method).

        This is a convenience method that combines batch creation and embedding.
        For pipeline use, prefer calling create_length_sorted_batches() and
        embed_prepared_batches() separately to enable overlapping batch creation
        for chunk N+1 while chunk N is being embedded.

        Args:
            sequences: Protein sequences
            accessions: Optional IDs (auto-generated if None)
            batch_size: Sequences per batch per GPU
            progress: Optional Progress instance
            task_id: Optional TaskID to update existing progress task for embedding
            prepare_task_id: Optional TaskID to update progress for batch preparation

        Yields:
            EmbeddingBatch for each completed batch. id↔embedding↔sequence
            relation is preserved within each batch, but batch order is
            non-deterministic (depends on GPU completion timing).
        """
        if not self._initialized:
            await self.initialize()

        accessions = accessions or [f"seq_{i}" for i in range(len(sequences))]
        if len(sequences) != len(accessions):
            raise ValueError("sequences and accessions must have same length")

        # Create sorted batches in thread pool
        batches = await asyncio.to_thread(
            self.create_length_sorted_batches,
            sequences,
            accessions,
            batch_size,
        )

        # Update prepare progress - batch creation is complete
        if progress is not None and prepare_task_id is not None:
            progress.update(prepare_task_id, advance=len(sequences))

        # Delegate to embed_prepared_batches for processing
        async for batch in self.embed_prepared_batches(batches, progress, task_id):
            yield batch

    async def embed_stream(
        self,
        batches: AsyncIterator[tuple[list[str], list[str]]],
    ) -> AsyncIterator[EmbeddingBatch]:
        """Stream embeddings from batch iterator.

        Pipeline controls pace - embedder processes batches as they arrive.
        Internal GPU workers and queues are completely hidden.

        Uses bounded queues to apply backpressure when GPU workers are saturated,
        preventing upstream reader from racing ahead and consuming excessive memory.

        Args:
            batches: Async iterator of (sequences, accessions) tuples

        Yields:
            EmbeddingBatch for each processed batch
        """
        if not self._initialized:
            await self.initialize()

        async def producer():
            """Enqueue batches as they arrive (blocks when queue full)."""
            async for batch_seqs, batch_accs in batches:
                # This will block when _batch_queue is full, applying backpressure
                await self._batch_queue.put((batch_seqs, batch_accs))
            # Send shutdown sentinels when iterator exhausted
            num_workers = len(self._worker_tasks or [])
            for _ in range(num_workers):
                await self._batch_queue.put(None)

        # Start producer task (runs concurrently with consumer)
        producer_task = asyncio.create_task(producer())

        # Yield results as they complete
        sentinels_received = 0
        num_workers = len(self._worker_tasks or [])

        try:
            while sentinels_received < num_workers:
                batch_result = await self._result_queue.get()

                if batch_result is None:
                    # Sentinel from GPU worker
                    sentinels_received += 1
                    self._result_queue.task_done()
                    continue

                # Real result - yield immediately
                yield batch_result
                self._result_queue.task_done()
        finally:
            # Ensure producer finishes
            await producer_task

    async def embed_prepared_batches(
        self,
        batches: list[tuple[list[str], list[str]]],
        progress: Progress | None = None,
        task_id: TaskID | None = None,
    ) -> AsyncIterator[EmbeddingBatch]:
        """Process pre-created batches (no batch creation step).

        This method processes batches that have already been created and sorted
        externally. Use this in the pipeline to enable overlapping batch creation
        for chunk N+1 while chunk N is being embedded.

        Args:
            batches: Pre-created batches as (sequences, accessions) tuples
            progress: Optional Progress instance
            task_id: Optional TaskID for embedding progress

        Yields:
            EmbeddingBatch for each completed batch. id↔embedding↔sequence
            relation is preserved within each batch.

        Example:
            # In pipeline: create batches externally
            batches = embedder.create_length_sorted_batches(sequences, accessions, batch_size)
            async for batch in embedder.embed_prepared_batches(batches):
                await db.insert(batch)
        """
        if not self._initialized:
            await self.initialize()

        # Enqueue all batches immediately
        for batch_seqs, batch_accs in batches:
            await self._batch_queue.put((batch_seqs, batch_accs))

        # Wait for results - they can arrive in any order
        collected = 0
        while collected < len(batches):
            batch_result = await self._result_queue.get()

            collected += 1

            if progress is not None and task_id is not None:
                progress.update(task_id, advance=len(batch_result.protein_ids))

            yield batch_result

            self._result_queue.task_done()

    async def _persistent_worker(
        self,
        device_idx: int,
    ) -> None:
        """Persistent worker that processes batches until shutdown.

        Runs continuously from initialize() until cleanup(). Pulls jobs from
        batch_queue, processes them, and puts results in result_queue.

        Args:
            device_idx: Index of GPU device to use
        """
        logger.debug(f"Persistent worker {device_idx} starting")

        batch_count = 0
        while True:
            # Get next job from queue
            job = await self._batch_queue.get()
            # None sentinel signals shutdown
            if job is None:
                self._batch_queue.task_done()
                # Signal consumer that this worker is shutting down
                await self._result_queue.put(None)
                logger.debug(f"Persistent worker {device_idx} shutting down")
                break

            batch_seqs, batch_accs = job

            # Process batch on GPU
            batch_result = await asyncio.to_thread(
                self._compute_embeddings,
                device_idx,
                batch_seqs,
                batch_accs,
            )

            # Put result in result queue
            await self._result_queue.put(batch_result)
            batch_count += 1
            self._batch_queue.task_done()

        logger.debug(f"Persistent worker {device_idx} finished ({batch_count} batches)")

    def _compute_embeddings(
        self,
        device_idx: int,
        sequences: Sequence[str],
        protein_ids: Sequence[str],
    ) -> EmbeddingBatch:
        """
        Compute embeddings for a batch (GPU-bound, runs in thread pool).

        This is the core computation that runs synchronously in a thread.
        Applies all configured pooling methods to the hidden states.

        Args:
            device_idx: Index of GPU device to use
            sequences: Batch of protein sequences
            protein_ids: Protein IDs

        Returns:
            EmbeddingBatch with protein_ids, sequences, and embeddings dict
            Each embedding in embeddings dict is a list of numpy arrays (one per sequence)
            Each array has shape [hidden_dim] or [seq_len * hidden_dim]

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
        inputs = tokenizer(
            list(sequences),
            return_tensors="pt",
            padding=True,  # Pads to longest in batch, not max_length
            truncation=True,
            max_length=eff_max_len,
        )

        # Pin CPU tensors to enable true async H→D transfers
        for k in inputs:
            if inputs[k].device.type == "cpu":
                inputs[k] = inputs[k].pin_memory()

        with torch.inference_mode():
            # Batch tokenize - padding to longest sequence in this batch

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

            # Apply all pooling methods on GPU first (no CPU transfers yet)
            # This batches all GPU operations before any CPU transfer
            pooled_tensors: dict[str, torch.Tensor] = {}
            pool_funcs: dict[str, PoolingFn | None] = {}

            for pool_name, pool_func in self.pooling_configs:
                pool_funcs[pool_name] = pool_func

                if pool_func is None:
                    # No pooling: keep raw hidden states on GPU [B, L, D]
                    pooled_tensors[pool_name] = hidden_states
                else:
                    # Apply pooling function (GPU operation)
                    pooled = pool_func(hidden_states, attention_mask)  # [B, D]

                    if self.normalize:
                        if self.return_dtype == "float32":
                            pooled = l2_normalize(pooled)
                        else:
                            pooled = normalize_cast_renorm(pooled, self._tf_return_dtype)

                    pooled_tensors[pool_name] = pooled  # Keep on GPU

            if device.type == "cuda":
                torch.cuda.synchronize(device=device)  # Sync only this device

            # Now transfer all results to CPU (separate phase after all GPU work)
            result = EmbeddingBatch()
            result.protein_ids = list(protein_ids)
            result.sequences = list(sequences)

            for pool_name, pooled_tensor in pooled_tensors.items():
                pool_func = pool_funcs[pool_name]

                # Single GPU→CPU transfer per pooling method
                # (Can't batch different shapes, but all GPU work is done first)
                if pool_func is None:
                    # 3D tensor: [B, L, D] - no pooling
                    raw_np = pooled_tensor.cpu().numpy().astype(self._np_return_dtype)
                    result.embeddings[pool_name] = [raw_np[i] for i in range(raw_np.shape[0])]
                else:
                    # 2D tensor: [B, D] - pooled
                    pooled_np = pooled_tensor.cpu().numpy().astype(self._np_return_dtype)
                    result.embeddings[pool_name] = [pooled_np[i] for i in range(pooled_np.shape[0])]

            # delete tensors
            del inputs, outputs, hidden_states, pooled_tensors

        return result

    async def __aenter__(self) -> "ESM2Embedder":
        await self.initialize()
        return self

    async def __aexit__(self, exc_type, exc_val, exc_tb) -> None:  # type: ignore
        await self.cleanup()
