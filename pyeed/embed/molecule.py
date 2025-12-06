"""Async SELFIES-based embedding using IBM materials.selfies-ted.

This module exposes `embed_smiles` to compute embeddings for one or many
SMILES strings, optionally using multiple CUDA devices.

Example:
    import asyncio
    from selfies_embedder import embed_smiles

    async def main() -> None:
        emb = await embed_smiles(["c1ccccc1", "CCO"], device=-1)
        print(emb.shape)

    asyncio.run(main())
"""

from __future__ import annotations

import asyncio
from collections.abc import AsyncIterator, Sequence
from typing import Final

import numpy as np
import selfies as sf
import torch
from transformers import AutoModel, AutoTokenizer

_MODEL_NAME: Final[str] = "ibm/materials.selfies-ted"

# Global tokenizer; model instances are created per device.
_TOKENIZER: Final = AutoTokenizer.from_pretrained(_MODEL_NAME)
_DEVICE_MODELS: dict[str, AutoModel] = {}


def _get_devices(device: int) -> list[str]:
    """Return list of device strings to use.

    Args:
        device: CUDA device index. -1 uses all available devices.
            Any non-negative index uses that CUDA device if available.
            Fallback is CPU if CUDA is not available.

    Returns:
        List of torch device strings.
    """
    if not torch.cuda.is_available():
        return ["cpu"]

    if device == -1:
        count = torch.cuda.device_count()
        return [f"cuda:{i}" for i in range(count)] if count > 0 else ["cpu"]

    if device >= 0:
        idx = min(device, max(torch.cuda.device_count() - 1, 0))
        return [f"cuda:{idx}"]

    return ["cpu"]


def _get_model_for_device(device: str) -> AutoModel:
    """Return a model instance placed on the given device."""
    if device not in _DEVICE_MODELS:
        model = AutoModel.from_pretrained(_MODEL_NAME)
        model.to(device)
        model.eval()
        _DEVICE_MODELS[device] = model
    return _DEVICE_MODELS[device]


def _smiles_to_selfies(smiles_list: Sequence[str]) -> list[str]:
    """Convert SMILES to spaced SELFIES tokens."""
    selfies_list: list[str] = []
    for smi in smiles_list:
        s = sf.encoder(smi)
        s = s.replace("][", "] [")
        selfies_list.append(s)
    return selfies_list


def _chunk_indices(n: int, k: int) -> list[tuple[int, int]]:
    """Split range(n) into k nearly equal contiguous chunks."""
    if k <= 1 or n <= 1:
        return [(0, n)]
    base = n // k
    rem = n % k
    chunks: list[tuple[int, int]] = []
    start = 0
    for i in range(k):
        size = base + (1 if i < rem else 0)
        end = start + size
        if size > 0:
            chunks.append((start, end))
        start = end
    return chunks


def _embed_on_device(
    selfies_list: Sequence[str],
    device: str,
    max_length: int,
) -> torch.Tensor:
    """Blocking embedding computation on a single device."""
    model = _get_model_for_device(device)
    with torch.inference_mode():
        tokens = _TOKENIZER(
            list(selfies_list),
            return_tensors="pt",
            max_length=max_length,
            truncation=True,
            padding="max_length",
        )
        input_ids = tokens["input_ids"].to(device)
        attention_mask = tokens["attention_mask"].to(device)

        outputs = model.encoder(input_ids=input_ids, attention_mask=attention_mask)
        hidden = outputs.last_hidden_state  # (batch, seq, dim)

        mask_expanded = attention_mask.unsqueeze(-1).expand(hidden.size()).float()
        summed = torch.sum(hidden * mask_expanded, dim=1)
        mask_sum = torch.clamp(mask_expanded.sum(dim=1), min=1e-9)
        pooled = summed / mask_sum

    return pooled.to("cpu")


def _tensor_rows_to_numpy(tensor: torch.Tensor) -> list[np.ndarray]:
    """Convert a 2D tensor on CPU to a list of numpy arrays."""
    if tensor.numel() == 0:
        return []
    # Ensure contiguous before view, then convert once
    arr = tensor.contiguous().numpy().astype(np.float32, copy=False)
    return [arr[i] for i in range(arr.shape[0])]


async def embed_smiles(
    smiles: str | Sequence[str],
    device: int = 0,
    max_length: int = 128,
) -> torch.Tensor:
    """Compute SELFIES-TED embeddings for SMILES strings asynchronously.

    Args:
        smiles: Single SMILES string or a sequence of SMILES strings.
        device: CUDA device index to use.
            -1: use all available CUDA devices (data-parallel).
             0+: use the given CUDA device index.
             Any value when CUDA is unavailable falls back to CPU.
        max_length: Maximum token length for the tokenizer.

    Returns:
        A tensor of shape (N, D) with embeddings on CPU, where
        N is the number of input SMILES strings.
    """
    smiles_list = [smiles] if isinstance(smiles, str) else list(smiles)

    if not smiles_list:
        raise AssertionError("Smiles list cannot be empty")

    selfies_list = _smiles_to_selfies(smiles_list)
    devices = _get_devices(device)

    if len(devices) == 1:
        return await asyncio.to_thread(
            _embed_on_device,
            selfies_list,
            devices[0],
            max_length,
        )

    chunks = _chunk_indices(len(selfies_list), len(devices))
    tasks: list[asyncio.Task[torch.Tensor]] = []
    for (start, end), dev in zip(chunks, devices, strict=False):
        sub = selfies_list[start:end]
        if not sub:
            continue
        tasks.append(
            asyncio.to_thread(
                _embed_on_device,
                sub,
                dev,
                max_length,
            )
        )

    results = await asyncio.gather(*tasks)
    return torch.cat(results, dim=0)


class MoleculeEmbedder:
    """SELFIES/SMILES embedder with multi-device support and backpressure."""

    Record = dict[str, object]
    BatchResult = list[Record]

    def __init__(
        self,
        device: int = -1,
        max_length: int = 128,
        max_in_flight: int | None = None,
    ):
        """
        Initialize embedder configuration.

        Args:
            device: -1 for all GPUs, 0+ for specific GPU, falls back to CPU.
            max_length: Maximum token length for tokenizer.
            max_in_flight: Optional cap of concurrent batches (defaults to 2×num_devices).
        """
        self.device = device
        self.max_length = max_length
        self.max_in_flight = max_in_flight
        self._initialized = False
        self._num_devices = 1

    async def initialize(self) -> None:
        """Prepare embedder. Idempotent."""
        if self._initialized:
            return

        if self.device == -1:
            self._num_devices = (torch.cuda.device_count() if torch.cuda.is_available() else 1) or 1
        else:
            self._num_devices = 1
        self._initialized = True

    async def embed_batch(
        self,
        smiles: Sequence[str],
        ids: Sequence[str] | None = None,
    ) -> BatchResult:
        """Embed a single batch and return list-of-dict records."""
        if not self._initialized:
            await self.initialize()

        ids = ids or [f"mol_{i}" for i in range(len(smiles))]
        embeddings = await embed_smiles(smiles, device=self.device, max_length=self.max_length)

        emb_rows = _tensor_rows_to_numpy(embeddings)
        return [
            {
                "id": mol_id,
                "smiles": smi,
                "embedding": emb,
                "embedding_name": "selfies_ted",
            }
            for mol_id, smi, emb in zip(ids, smiles, emb_rows, strict=False)
        ]

    async def embed_stream(
        self,
        batches: AsyncIterator[tuple[list[str], list[str]]],
    ) -> AsyncIterator[BatchResult]:
        """Process batches concurrently with bounded in-flight tasks."""
        if not self._initialized:
            await self.initialize()

        max_tasks = self.max_in_flight or max(1, self._num_devices * 2)
        pending: set[asyncio.Task[BatchResult]] = set()

        async def _submit(ids: list[str], smi: list[str]) -> None:
            pending.add(asyncio.create_task(self.embed_batch(smi, ids)))

        async for ids, smi in batches:
            while len(pending) >= max_tasks:
                done, pending = await asyncio.wait(pending, return_when=asyncio.FIRST_COMPLETED)
                for task in done:
                    yield task.result()
            await _submit(ids, smi)

        while pending:
            done, pending = await asyncio.wait(pending, return_when=asyncio.FIRST_COMPLETED)
            for task in done:
                yield task.result()


if __name__ == "__main__":
    import asyncio

    from rich import print

    async def main() -> None:
        emb = await embed_smiles(["c1ccccc1", "CCO"], device=-1)
        print(emb)
        print(emb.shape)
        print(emb.dtype)

    asyncio.run(main())
