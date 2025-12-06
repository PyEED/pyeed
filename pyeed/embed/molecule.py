from __future__ import annotations

import asyncio
import threading
from collections.abc import AsyncIterator, Sequence
from typing import Final

import numpy as np
import selfies as sf
import torch
from loguru import logger
from transformers import BartModel, PreTrainedTokenizerFast

_MODEL_NAME: Final[str] = "ibm/materials.selfies-ted"


class MoleculeEmbedder:
    """SELFIES/SMILES embedder with multi-device support.

    Args:
        device: Device selection strategy.
            -1: use all available CUDA devices (default).
            0+: use specific CUDA device index.
            Falls back to CPU when CUDA is unavailable.
    """

    def __init__(self, device: int = -1) -> None:
        self._device_arg = device
        self._initialized = False
        self._devices: list[str] = []
        self._tokenizer: PreTrainedTokenizerFast | None = None
        self._tokenizer_lock = threading.Lock()
        self._models: dict[str, BartModel] = {}
        self._max_length: int = 0

    async def _initialize(self) -> None:
        """Lazy initialization of tokenizer and device list."""
        if self._initialized:
            return

        self._tokenizer = await asyncio.to_thread(
            PreTrainedTokenizerFast.from_pretrained, _MODEL_NAME
        )
        self._max_length = self._tokenizer.model_max_length
        self._devices = self._resolve_devices()
        self._initialized = True

    def _resolve_devices(self) -> list[str]:
        """Resolve device strings based on configuration."""
        if not torch.cuda.is_available():
            return ["cpu"]

        if self._device_arg == -1:
            count = torch.cuda.device_count()
            return [f"cuda:{i}" for i in range(count)] if count > 0 else ["cpu"]

        if self._device_arg >= 0:
            idx = min(self._device_arg, max(torch.cuda.device_count() - 1, 0))
            return [f"cuda:{idx}"]

        return ["cpu"]

    def _get_model(self, device: str) -> BartModel:
        """Get or create model for a specific device."""
        if device not in self._models:
            model = BartModel.from_pretrained(_MODEL_NAME)
            model.to(device)
            model.eval()
            self._models[device] = model
        return self._models[device]

    def _prepare_molecules(
        self,
        ids: list[str],
        smiles_list: list[str],
    ) -> tuple[list[str], list[str]]:
        """Convert SMILES to SELFIES and filter invalid molecules.

        Validates each molecule by:
        1. Converting SMILES to SELFIES (skips on EncoderError)
        2. Checking token length against model max (skips if exceeded)

        Args:
            ids: Molecule identifiers.
            smiles_list: SMILES strings to convert.

        Returns:
            Tuple of (valid_ids, valid_selfies) for molecules that passed
            all validation. Invalid molecules are logged and excluded.
        """
        assert self._tokenizer is not None

        valid_ids: list[str] = []
        valid_selfies: list[str] = []

        for mol_id, smiles in zip(ids, smiles_list, strict=True):
            # Convert SMILES to SELFIES
            try:
                selfies_str = sf.encoder(smiles).replace("][", "] [")
            except sf.exceptions.EncoderError as e:
                logger.warning(
                    f"Skipping molecule {mol_id}: SELFIES encoding failed",
                    extra={"id": mol_id, "error": e},
                )
                continue

            # Check token length (with lock for thread safety)
            with self._tokenizer_lock:
                tokens = self._tokenizer(selfies_str, add_special_tokens=True)
            token_len = len(tokens["input_ids"])

            if token_len > self._max_length:
                logger.warning(
                    f"Skipping molecule {mol_id}: tokens exceed max length",
                    extra={"id": mol_id, "token_len": token_len, "max_length": self._max_length},
                )
                continue

            valid_ids.append(mol_id)
            valid_selfies.append(selfies_str)

        return valid_ids, valid_selfies

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

    def _embed_on_device(self, selfies_list: list[str], device: str) -> np.ndarray:
        """Blocking embedding computation on a single device."""
        model = self._get_model(device)
        assert self._tokenizer is not None

        with torch.inference_mode():
            # Tokenize with lock for thread safety
            with self._tokenizer_lock:
                tokens = self._tokenizer(
                    selfies_list,
                    return_tensors="pt",
                    truncation=True,
                    max_length=self._max_length,
                    padding=True,
                )
            input_ids = tokens["input_ids"].to(device)
            attention_mask = tokens["attention_mask"].to(device)

            outputs = model.encoder(input_ids=input_ids, attention_mask=attention_mask)
            hidden = outputs.last_hidden_state

            mask_expanded = attention_mask.unsqueeze(-1).expand(hidden.size()).float()
            summed = torch.sum(hidden * mask_expanded, dim=1)
            mask_sum = torch.clamp(mask_expanded.sum(dim=1), min=1e-9)
            pooled = summed / mask_sum

        return pooled.cpu().numpy().astype(np.float32)

    async def embed(
        self,
        ids: Sequence[str],
        smiles: Sequence[str],
    ) -> list[tuple[str, np.ndarray]]:
        """Compute embeddings for SMILES strings.

        Args:
            ids: Identifiers for each molecule (must match smiles length).
            smiles: SMILES strings to embed.

        Returns:
            List of (id, embedding) tuples. Invalid molecules (encoding errors
            or exceeding token limit) are logged and excluded. Order is
            preserved for valid molecules.

        Raises:
            ValueError: If ids and smiles have different lengths or are empty.
        """
        if len(ids) != len(smiles):
            raise ValueError(f"ids and smiles must have same length: {len(ids)} != {len(smiles)}")
        if len(ids) == 0:
            raise ValueError("ids and smiles cannot be empty")

        await self._initialize()

        # Validate and convert molecules
        valid_ids, valid_selfies = await asyncio.to_thread(
            self._prepare_molecules, list(ids), list(smiles)
        )

        if not valid_ids:
            logger.warning("No valid molecules to embed, returning empty results")
            return []

        # Compute embeddings
        if len(self._devices) == 1:
            embeddings = await asyncio.to_thread(
                self._embed_on_device, valid_selfies, self._devices[0]
            )
        else:
            chunks = self._chunk_indices(len(valid_selfies), len(self._devices))
            tasks: list[asyncio.Task[np.ndarray]] = []

            for (start, end), device in zip(chunks, self._devices, strict=False):
                chunk_selfies = valid_selfies[start:end]
                if chunk_selfies:
                    task = asyncio.create_task(
                        asyncio.to_thread(self._embed_on_device, chunk_selfies, device)
                    )
                    tasks.append(task)

            results = await asyncio.gather(*tasks)
            embeddings = np.vstack(results)

        return [(mol_id, embeddings[i]) for i, mol_id in enumerate(valid_ids)]

    async def embed_stream(
        self,
        batches: AsyncIterator[tuple[list[str], list[str]]],
        max_in_flight: int | None = None,
    ) -> AsyncIterator[list[tuple[str, np.ndarray]]]:
        """Process batches concurrently with bounded in-flight tasks.

        Args:
            batches: Async iterator yielding (ids, smiles) tuples.
            max_in_flight: Maximum concurrent batches. Defaults to 2 * num_devices.

        Yields:
            List of (id, embedding) tuples for each batch. Invalid molecules
            are logged and excluded. Order within each batch is preserved;
            batch completion order may vary.
        """
        await self._initialize()

        concurrency = max_in_flight or max(1, len(self._devices) * 2)
        pending: set[asyncio.Task[list[tuple[str, np.ndarray]]]] = set()

        async for batch_ids, batch_smiles in batches:
            while len(pending) >= concurrency:
                done, pending = await asyncio.wait(pending, return_when=asyncio.FIRST_COMPLETED)
                for task in done:
                    yield task.result()

            task = asyncio.create_task(self.embed(batch_ids, batch_smiles))
            pending.add(task)

        while pending:
            done, pending = await asyncio.wait(pending, return_when=asyncio.FIRST_COMPLETED)
            for task in done:
                yield task.result()

    async def cleanup(self) -> None:
        """Clean up models and free resources."""
        logger.debug("Cleaning up MoleculeEmbedder")
        self._models.clear()
        self._tokenizer = None
        self._devices.clear()
        self._initialized = False
        logger.debug("MoleculeEmbedder cleanup complete")

    async def __aenter__(self) -> MoleculeEmbedder:
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
    from rich import print as rprint

    async def main() -> None:
        embedder = MoleculeEmbedder()
        results = await embedder.embed(
            ["benzene", "ethanol", "too_long", "invalid_smiles"],
            [
                "c1ccccc1",
                "CCO",
                "C" * 5000,
                "not_a_valid_smiles!!!",
            ],
        )
        rprint(f"Returned {len(results)} embeddings:")
        for mol_id, embedding in results:
            rprint(f"  {mol_id}: shape={embedding.shape}, dtype={embedding.dtype}")

    asyncio.run(main())
