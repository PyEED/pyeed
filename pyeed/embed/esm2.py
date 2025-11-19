"""ESM2 protein embedder with GPU parallelism and per-GPU thread safety.

Architecture:
    Uses per-GPU locks to serialize model access. Multiple batches can process
    concurrently across different GPUs, but only one batch processes per GPU
    at a time. This prevents rotary embeddings cache corruption.

Thread Safety:
    ESM2's rotary embeddings use per-model-instance cache. Concurrent access
    from multiple threads corrupts this cache. Per-GPU locks ensure sequential
    access to each model instance while maintaining parallelism across GPUs.

Backpressure Flow:
    GPU processing slow → Task limit reached → embed_stream blocks →
    Pipeline queue fills → Reader blocks → No more data ingested

This prevents memory exhaustion when embedding is the bottleneck, while ensuring
thread-safe execution and maintaining GPU utilization.
"""

import asyncio
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


class ESM2Embedder:
    def __init__(
        self,
        model_dtype: ModelDType,
        return_dtype: ReturnDType,
        pooling_methods: PoolingLike,
        model_name: str = "facebook/esm2_t33_650M_UR50D",
        normalize: bool = True,
        max_length: int = 1024,
        n_gpus: int = -1,
        huggingface_token: str | None = None,
    ):
        """
        Initialize ESM2 embedder.

        Args:
            model_name: ESM-2 model identifier
            model_dtype: Model precision for computation (float32/float16)
            return_dtype: Output array dtype (float32/float16)
            pooling_methods: Single or multiple pooling functions
            normalize: L2-normalize embeddings
            max_length: Maximum sequence length
            n_gpus: Number of GPUs to use (-1 for all available)
            huggingface_token: HuggingFace API token (optional, for gated models)
        """
        self.model_name = model_name
        self.normalize = normalize
        self.model_dtype = model_dtype
        self.return_dtype = return_dtype
        self.max_length = max_length
        self.n_gpus = n_gpus
        self._initialized = False
        self._tf_model_dtype = TF_DTYPE_MAP[model_dtype]
        self._tf_return_dtype = TF_DTYPE_MAP[return_dtype]
        self._np_return_dtype = NP_DTYPE_MAP[return_dtype]
        self.huggingface_token = huggingface_token or _login_hf()
        # Normalize pooling_methods into a list of (name, function) tuples
        self.pooling_configs = self._normalize_pooling_methods(pooling_methods)

        self.models: list[EsmModel | None] = []
        self.tokenizer: EsmTokenizer | None = None  # Single tokenizer on CPU
        self.devices: list[torch.device] = []
        self.device_ids: list[int] = []
        self._device_locks: list[asyncio.Lock] | None = None

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

    async def initialize(self) -> None:
        """Load models on all available GPUs and tokenizer on CPU."""
        if self._initialized:
            logger.debug("ESM2 embedder already initialized")
            return

        if self.n_gpus == -1:
            self.devices = self._detect_devices()
        else:
            self.devices = [torch.device(f"cuda:{i}") for i in range(self.n_gpus)]

        self.models = [None] * len(self.devices)

        # Load single tokenizer on CPU (thread-safe, no device affinity)
        self.tokenizer = await asyncio.to_thread(self._load_tokenizer)

        # Load models on all devices in parallel
        await asyncio.gather(*[self._load_on_device(i, d) for i, d in enumerate(self.devices)])

        # Create per-device locks for thread-safe model access
        self._device_locks = [asyncio.Lock() for _ in range(len(self.devices))]
        logger.debug(f"Created {len(self._device_locks)} device locks for thread safety")

        self._initialized = True
        logger.info(
            "ESM2 embedder ready",
            extra={"num_gpus": len(self.devices), "devices": [str(d) for d in self.devices]},
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
                token=self.huggingface_token,
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
                    token=self.huggingface_token,
                    torch_dtype=self._tf_model_dtype,
                )

            # Disable all-layer hidden states to save memory
            model.config.output_hidden_states = False

            model = model.to(device=device, dtype=self._tf_model_dtype).eval()

            return model

        model = await asyncio.to_thread(_load)
        # Use indexed assignment instead of append to avoid race conditions
        self.models[device_idx] = model
        logger.info(
            f"ESM2 model loaded on {device!r}",
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
        """Stream embeddings from batch iterator with concurrent GPU processing.

        Processes batches concurrently across available GPUs using task-based
        parallelism. Results are yielded as they complete (order non-deterministic).

        Backpressure: Limits concurrent tasks to prevent memory exhaustion when
        GPU processing is slower than batch generation.

        Args:
            batches: Async iterator of (sequences, accessions) tuples

        Yields:
            EmbeddingBatch for each processed batch (order non-deterministic)

        Example:
            ```python
            async def batch_gen():
                for batch in batches:
                    yield batch

            async for result in embedder.embed_stream(batch_gen()):
                process(result)
            ```
        """
        num_gpus = len(self.devices)
        device_idx = 0

        # Sliding window: Limit concurrent tasks for backpressure
        # Strategy: Allow 2 batches per GPU in flight
        max_concurrent_batches = num_gpus * 2

        pending_tasks: set[asyncio.Task] = set()
        batch_count = 0
        results_yielded = 0

        logger.info(
            "Starting embedding stream",
            extra={
                "num_gpus": num_gpus,
                "max_concurrent": max_concurrent_batches,
            },
        )

        try:
            async for batch_seqs, batch_accs in batches:
                batch_count += 1

                # Backpressure: Wait if too many tasks in flight
                while len(pending_tasks) >= max_concurrent_batches:
                    logger.debug(
                        f"Concurrent limit reached "
                        f"({len(pending_tasks)}/{max_concurrent_batches}), "
                        f"waiting for completion"
                    )
                    # Yield one completed result
                    done, pending_tasks = await asyncio.wait(
                        pending_tasks, return_when=asyncio.FIRST_COMPLETED
                    )
                    for task in done:
                        result = task.result()
                        results_yielded += 1
                        logger.debug(
                            f"Yielding result {results_yielded}/{batch_count}",
                            extra={"num_proteins": len(result.protein_ids)},
                        )
                        yield result

                # Create task for this batch (round-robin GPU assignment)
                logger.debug(
                    f"Creating task for batch {batch_count}",
                    extra={"device_idx": device_idx, "batch_size": len(batch_seqs)},
                )

                # Define async function that acquires lock before processing
                async def process_with_lock(
                    gpu_idx: int, seqs: list[str], accs: list[str]
                ) -> EmbeddingBatch:
                    """Process batch with per-GPU lock for thread safety."""
                    async with self._device_locks[gpu_idx]:
                        # Only one batch processes on this GPU at a time
                        return await asyncio.to_thread(
                            self._compute_embeddings,
                            gpu_idx,
                            seqs,
                            accs,
                        )

                task = asyncio.create_task(process_with_lock(device_idx, batch_seqs, batch_accs))
                pending_tasks.add(task)
                device_idx = (device_idx + 1) % num_gpus

            # Yield remaining results
            logger.debug(f"All batches enqueued, waiting for {len(pending_tasks)} pending tasks")
            while pending_tasks:
                done, pending_tasks = await asyncio.wait(
                    pending_tasks, return_when=asyncio.FIRST_COMPLETED
                )
                for task in done:
                    result = task.result()
                    results_yielded += 1
                    logger.debug(
                        f"Yielding result {results_yielded}/{batch_count}",
                        extra={"num_proteins": len(result.protein_ids)},
                    )
                    yield result

        except Exception as e:
            logger.error(f"Error in embed_stream: {e}", exc_info=True)
            # Cancel pending tasks on error
            for task in pending_tasks:
                task.cancel()
            raise
        finally:
            logger.info(
                "Embedding stream complete",
                extra={
                    "batches_processed": batch_count,
                    "results_yielded": results_yielded,
                },
            )

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

        logger.debug(f"Computed embeddings for {len(result.protein_ids)} protein IDs")
        return result

    async def __aenter__(self) -> "ESM2Embedder":
        await self.initialize()
        return self

    async def __aexit__(self, exc_type, exc_val, exc_tb) -> None:  # type: ignore
        await self.cleanup()
