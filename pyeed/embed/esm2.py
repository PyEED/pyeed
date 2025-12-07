from __future__ import annotations

import asyncio
import threading
from collections.abc import AsyncIterator, Sequence
from typing import Literal

import numpy as np
import torch
from loguru import logger
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


class ESM2Embedder:
    """ESM2 protein sequence embedder with multi-device support.

    Args:
        model_name: HuggingFace model identifier for ESM2.
        dtype: Data type for model computation and output embeddings.
        pooling: Pooling function to reduce sequence dimension. None for raw hidden states.
        normalize: Whether to L2-normalize output embeddings.
        max_length: Maximum sequence length (tokens). Sequences exceeding this are skipped.
        devices: List of CUDA device indices to use.
            Empty list []: use all available CUDA devices (default).
            [0, 4]: use specific CUDA devices (cuda:0 and cuda:4).
            Falls back to CPU when CUDA is unavailable.
        huggingface_token: HuggingFace token for private models.

    Example:
        >>> embedder = ESM2Embedder(dtype="float16", devices=[0])
        >>> results = await embedder.embed(["P12345"], ["MVLSPADKTN..."])
        >>> for protein_id, embedding in results:
        ...     print(f"{protein_id}: {embedding.shape}")
    """

    def __init__(
        self,
        model_name: str = "facebook/esm2_t33_650M_UR50D",
        dtype: DType = "float32",
        pooling: PoolingFn | None = mean_pooling,
        normalize: bool = True,
        max_length: int = 1024,
        huggingface_token: str | None = None,
        devices: list[int] | None = None,
    ) -> None:
        self._model_name = model_name
        self._dtype = dtype
        self._torch_dtype = _TORCH_DTYPE_MAP[dtype]
        self._np_dtype = _NP_DTYPE_MAP[dtype]
        self._pooling = pooling
        self._normalize = normalize
        self._max_length = max_length
        self._huggingface_token = huggingface_token or _login_hf()
        self._device_args = devices if devices is not None else []

        self._initialized = False
        self._devices: list[str] = []
        self._tokenizer: EsmTokenizer | None = None
        self._tokenizer_lock = threading.Lock()
        self._models: dict[str, EsmModel] = {}
        self._hf_token: str | None = None
        # Serialize embed() calls to avoid concurrent model forwards on same device
        self._embed_lock = asyncio.Lock()

    async def _initialize(self) -> None:
        """Lazy initialization of tokenizer and device list."""
        if self._initialized:
            return

        self._tokenizer = await asyncio.to_thread(self._load_tokenizer)
        self._devices = self._resolve_devices()
        self._initialized = True

        logger.info(
            "ESM2 embedder initialized",
            extra={
                "model": self._model_name,
                "dtype": self._dtype,
                "devices": self._devices,
                "max_length": self._max_length,
            },
        )

    def _load_tokenizer(self) -> EsmTokenizer:
        """Load tokenizer (called once, blocking)."""
        with silence_transformers_init_only():
            return EsmTokenizer.from_pretrained(
                self._model_name,
                token=self._hf_token,
            )

    def _resolve_devices(self) -> list[str]:
        """Resolve device strings based on configuration.

        Returns:
            List of device strings (e.g., ["cuda:0", "cuda:4"] or ["cpu"]).
        """
        if not torch.cuda.is_available():
            return ["cpu"]

        # Empty list means use all available CUDA devices
        if not self._device_args:
            count = torch.cuda.device_count()
            return [f"cuda:{i}" for i in range(count)] if count > 0 else ["cpu"]

        # Use specific device indices
        max_device = torch.cuda.device_count() - 1
        device_strings: list[str] = []
        for idx in self._device_args:
            if idx < 0:
                logger.warning(f"Invalid device index {idx}, skipping")
                continue
            if idx > max_device:
                logger.warning(
                    f"Device index {idx} exceeds available devices (max: {max_device}), skipping"
                )
                continue
            device_strings.append(f"cuda:{idx}")

        return device_strings if device_strings else ["cpu"]

    def _get_model(self, device: str) -> EsmModel:
        """Get or create model for a specific device."""
        if device not in self._models:
            # Set CUDA device context for proper tensor placement
            if device.startswith("cuda:"):
                torch.cuda.set_device(int(device.split(":")[1]))

            with silence_transformers_init_only():
                model = EsmModel.from_pretrained(
                    self._model_name,
                    token=self._hf_token,
                    torch_dtype=self._torch_dtype,
                )

            model.config.output_hidden_states = False
            model = model.to(device=device, dtype=self._torch_dtype)
            model.eval()
            self._models[device] = model

            logger.debug(f"Loaded ESM2 model on {device}")

        return self._models[device]

    def _prepare_sequences(
        self,
        ids: list[str],
        sequences: list[str],
    ) -> tuple[list[str], list[str]]:
        """Validate sequences and filter those exceeding max length.

        Args:
            ids: Protein identifiers.
            sequences: Amino acid sequences.

        Returns:
            Tuple of (valid_ids, valid_sequences) for sequences that passed
            validation. Sequences exceeding max_length are logged and excluded.
        """
        assert self._tokenizer is not None

        valid_ids: list[str] = []
        valid_sequences: list[str] = []

        for protein_id, sequence in zip(ids, sequences, strict=True):
            # Check token length (with lock for thread safety)
            with self._tokenizer_lock:
                tokens = self._tokenizer(sequence, add_special_tokens=True)
            token_len = len(tokens["input_ids"])

            if token_len > self._max_length:
                logger.warning(
                    f"Skipping protein {protein_id}: tokens exceed max length",
                    extra={
                        "id": protein_id,
                        "token_len": token_len,
                        "max_length": self._max_length,
                        "seq_len": len(sequence),
                    },
                )
                continue

            valid_ids.append(protein_id)
            valid_sequences.append(sequence)

        return valid_ids, valid_sequences

    def _chunk_indices(self, n: int, k: int) -> list[tuple[int, int]]:
        """Split range(n) into k nearly equal contiguous chunks."""
        if k <= 1 or n <= 1:
            return [(0, n)]

        base, rem = divmod(n, k)
        chunks: list[tuple[int, int]] = []
        start = 0
        for i in range(k):
            size = base + (1 if i < rem else 0)
            if size > 0:
                chunks.append((start, start + size))
                start += size
        return chunks

    def _embed_on_device(self, sequences: list[str], device: str) -> np.ndarray:
        """Blocking embedding computation on a single device.

        Args:
            sequences: Amino acid sequences to embed.
            device: Target device string (e.g., "cuda:0", "cpu").

        Returns:
            Embeddings array with shape [batch, hidden_dim] if pooling,
            or [batch, seq_len, hidden_dim] if pooling is None.
        """
        model = self._get_model(device)
        assert self._tokenizer is not None

        # Set CUDA device context
        if device.startswith("cuda:"):
            torch.cuda.set_device(int(device.split(":")[1]))

        with torch.inference_mode():
            # Tokenize with lock for thread safety
            with self._tokenizer_lock:
                tokens = self._tokenizer(
                    sequences,
                    return_tensors="pt",
                    truncation=True,
                    max_length=self._max_length,
                    padding=True,
                )

            input_ids = tokens["input_ids"].to(device)
            attention_mask = tokens["attention_mask"].to(device)

            # Use autocast for fp16/bf16
            use_autocast = device.startswith("cuda") and self._torch_dtype in (
                torch.float16,
                torch.bfloat16,
            )

            if use_autocast:
                with torch.autocast(device_type="cuda", dtype=self._torch_dtype):
                    outputs = model(input_ids=input_ids, attention_mask=attention_mask)
            else:
                outputs = model(input_ids=input_ids, attention_mask=attention_mask)

            hidden_states = outputs.last_hidden_state  # [B, L, D]

            # Apply pooling
            if self._pooling is None:
                result = hidden_states
            else:
                result = self._pooling(hidden_states, attention_mask)  # [B, D]

                # Normalize if requested
                if self._normalize:
                    result = l2_normalize(result)

            # Sync before transfer
            if device.startswith("cuda"):
                torch.cuda.synchronize()

        return result.cpu().numpy().astype(self._np_dtype)

    @property
    def embedding_dim(self) -> int:
        """Get embedding dimension from model config.

        Returns:
            Hidden size of the model (embedding dimension).

        Raises:
            RuntimeError: If model not initialized.
        """
        if not self._initialized or not self._models:
            raise RuntimeError("Model not initialized. Call embed() first.")
        model = next(iter(self._models.values()))
        return model.config.hidden_size

    async def embed(
        self,
        ids: Sequence[str],
        sequences: Sequence[str],
    ) -> list[tuple[str, np.ndarray]]:
        """Compute embeddings for protein sequences.

        Args:
            ids: Protein identifiers (must match sequences length).
            sequences: Amino acid sequences to embed.

        Returns:
            List of (id, embedding) tuples. Sequences exceeding max_length
            are logged and excluded. Order is preserved for valid sequences.

        Raises:
            ValueError: If ids and sequences have different lengths or are empty.
        """
        async with self._embed_lock:
            if len(ids) != len(sequences):
                raise ValueError(
                    f"ids and sequences must have same length: {len(ids)} != {len(sequences)}"
                )
            if len(ids) == 0:
                raise ValueError("ids and sequences cannot be empty")

            await self._initialize()

            # Validate sequences
            valid_ids, valid_sequences = await asyncio.to_thread(
                self._prepare_sequences, list(ids), list(sequences)
            )

            if not valid_ids:
                logger.warning("No valid sequences to embed, returning empty results")
                return []

            # Compute embeddings
            if len(self._devices) == 1:
                embeddings = await asyncio.to_thread(
                    self._embed_on_device, valid_sequences, self._devices[0]
                )
            else:
                chunks = self._chunk_indices(len(valid_sequences), len(self._devices))
                tasks: list[asyncio.Task[np.ndarray]] = []

                for (start, end), device in zip(chunks, self._devices, strict=False):
                    chunk_sequences = valid_sequences[start:end]
                    if chunk_sequences:
                        task = asyncio.create_task(
                            asyncio.to_thread(self._embed_on_device, chunk_sequences, device)
                        )
                        tasks.append(task)

                results = await asyncio.gather(*tasks)
                embeddings = np.vstack(results)

            return [(protein_id, embeddings[i]) for i, protein_id in enumerate(valid_ids)]

    async def embed_stream(
        self,
        batches: AsyncIterator[tuple[list[str], list[str]]],
        max_in_flight: int | None = None,
    ) -> AsyncIterator[list[tuple[str, np.ndarray]]]:
        """Process batches concurrently with bounded in-flight tasks.

        Args:
            batches: Async iterator yielding (ids, sequences) tuples.
            max_in_flight: Maximum concurrent batches. Defaults to 2 * num_devices.

        Yields:
            List of (id, embedding) tuples for each batch. Sequences exceeding
            max_length are logged and excluded. Order within each batch is
            preserved; batch completion order may vary.
        """
        await self._initialize()

        concurrency = max_in_flight or max(1, len(self._devices) * 1)
        pending: set[asyncio.Task[list[tuple[str, np.ndarray]]]] = set()

        async for batch_ids, batch_sequences in batches:
            while len(pending) >= concurrency:
                done, pending = await asyncio.wait(pending, return_when=asyncio.FIRST_COMPLETED)
                for task in done:
                    yield task.result()

            task = asyncio.create_task(self.embed(batch_ids, batch_sequences))
            pending.add(task)

        while pending:
            done, pending = await asyncio.wait(pending, return_when=asyncio.FIRST_COMPLETED)
            for task in done:
                yield task.result()

    async def cleanup(self) -> None:
        """Clean up models and free GPU memory."""
        logger.info("Cleaning up ESM2 embedder")

        # Extract device indices before clearing
        device_ids: list[int] = []
        for device in self._devices:
            if device.startswith("cuda:"):
                device_ids.append(int(device.split(":")[1]))

        self._models.clear()
        self._tokenizer = None
        self._devices.clear()
        self._initialized = False

        # Free GPU memory
        if device_ids:
            await asyncio.to_thread(_free_device_memory, device_ids)

        logger.info("ESM2 embedder cleanup complete")

    async def __aenter__(self) -> ESM2Embedder:
        """Async context manager entry."""
        await self._initialize()
        return self

    async def __aexit__(
        self,
        exc_type: type[BaseException] | None,
        exc_val: BaseException | None,
        exc_tb: object,
    ) -> None:
        """Async context manager exit with cleanup."""
        await self.cleanup()


if __name__ == "__main__":
    import numpy as np
    from rich import print as rprint

    async def main() -> None:
        embedder = ESM2Embedder(dtype="float16", device=0)
        await embedder._initialize()

        results = await embedder.embed(
            ["protein_1", "protein_2", "too_long", "invalid_sequence"],
            [
                "MVLSPADKTNVKAAWGKVGAHAGEYGAEALERMFLSFPTTKTYFPHFDLSH",
                "MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQAPILSRVGDGTQDNLSGAEK",
                "M" * 2000,
                "sequence: MGDTDLPCVTGVLLAAGAGKRLGRGPKALLPYRGRTLVEDAAETMLVGGCHEVVIVLGANAQAVCARANLEPYRIVVNHDWSSGMGSSYLAGDAAAHTKNHILVALVDQPGLSVTTVGRLLVSHRPGRISSAAYSSLDSPRVLRRGHPMVIDAGLRPAVASTVSGDAGARVFLRQKPWLVDLIDCSDESTGEDVDTVEQMYRL",
            ],
        )

        rprint(f"Returned {len(results)} embeddings:")
        await embedder.cleanup()

    asyncio.run(main())
