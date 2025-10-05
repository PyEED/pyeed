import asyncio
from collections.abc import AsyncIterator, Callable, Sequence

import torch
from transformers import EsmModel, EsmTokenizer

from pyeed.model.embedding import Embedding

from .pooling import mean_pooling, normalize_embeddings
from .utils import free_memory, get_hf_token


class ESM2Embedder:
    """
    ESM-2 embedder with automatic multi-GPU parallelization.

    Usage:
        async with AsyncESM2Embedder() as embedder:
            embeddings = await embedder.embed_batch(sequences, accessions)
    """

    def __init__(
        self,
        model_name: str = "facebook/esm2_t33_650M_UR50D",
        pooling: Callable = mean_pooling,
        layer_index: int = -1,
        normalize: bool = True,
        device_ids: list[int] | None = None,
    ):
        self.model_name = model_name
        self.pooling = pooling
        self.layer_index = layer_index
        self.normalize = normalize
        self.device_ids = device_ids

        self.models: list[EsmModel] = []
        self.tokenizers: list[EsmTokenizer] = []
        self.devices: list[torch.device] = []
        self._initialized = False

        print(f"[ESM2] Model: {model_name}")
        print(f"[ESM2] Pooling: {pooling.__name__}, Layer: {layer_index}, Normalize: {normalize}")

    async def initialize(self) -> None:
        """Load models on all available GPUs."""
        if self._initialized:
            print("[ESM2] Already initialized")
            return

        self.devices = self._detect_devices()
        print(f"[ESM2] Loading model on {len(self.devices)} device(s): {self.devices}")

        await asyncio.gather(*[self._load_on_device(d) for d in self.devices])

        self._initialized = True
        print(f"[ESM2] Ready")

    def _detect_devices(self) -> list[torch.device]:
        """Detect available CUDA devices or fallback to CPU."""
        if not torch.cuda.is_available():
            print("[ESM2] No CUDA available, using CPU")
            return [torch.device("cpu")]

        if self.device_ids is None:
            self.device_ids = list(range(torch.cuda.device_count()))

        return [torch.device(f"cuda:{i}") for i in self.device_ids]

    async def _load_on_device(self, device: torch.device) -> None:
        """Load model and tokenizer on specific device."""
        print(f"[ESM2] Loading on {device}...")

        def _load() -> tuple[EsmModel, EsmTokenizer]:
            token = get_hf_token()
            name = self.model_name if self.model_name.startswith("facebook/") else f"facebook/{self.model_name}"

            model = EsmModel.from_pretrained(name, token=token)
            tokenizer = EsmTokenizer.from_pretrained(name, token=token)
            model = model.to(device).eval()

            return model, tokenizer

        model, tokenizer = await asyncio.to_thread(_load)
        self.models.append(model)
        self.tokenizers.append(tokenizer)
        print(f"[ESM2] Loaded on {device}")

    async def cleanup(self) -> None:
        """Clean up models and free GPU memory."""
        print("[ESM2] Cleaning up...")
        self.models.clear()
        self.tokenizers.clear()
        self.devices.clear()
        self._initialized = False
        await asyncio.to_thread(free_memory)
        print("[ESM2] Cleanup complete")

    # ========================================================================
    # Public API
    # ========================================================================

    async def embed_batch(
        self,
        sequences: Sequence[str],
        accessions: Sequence[str] | None = None,
        batch_size: int = 8,
    ) -> list[Embedding]:
        """
        Generate embeddings for multiple sequences.

        Args:
            sequences: Protein sequences
            accessions: Optional IDs (auto-generated if None)
            batch_size: Sequences per batch per GPU

        Returns:
            List of Embedding objects in same order as input
        """
        if not self._initialized:
            await self.initialize()

        accessions = accessions or [f"seq_{i}" for i in range(len(sequences))]
        if len(sequences) != len(accessions):
            raise ValueError("sequences and accessions must have same length")

        print(f"[ESM2] Embedding {len(sequences)} sequences on {len(self.devices)} device(s), batch_size={batch_size}")

        # Distribute work round-robin across devices
        device_work = self._distribute_work(sequences, accessions)

        # Process all devices in parallel
        results = await asyncio.gather(*[
            self._process_device_work(device_idx, seqs, accs, batch_size)
            for device_idx, (seqs, accs) in enumerate(device_work)
        ])

        # Reorder results to match input order
        embeddings = self._reorder_results(results, len(sequences))
        print(f"[ESM2] Generated {len(embeddings)} embeddings")

        return embeddings

    async def embed_single(self, sequence: str, accession: str) -> Embedding:
        """Embed a single sequence."""
        results = await self.embed_batch([sequence], [accession], batch_size=1)
        return results[0]

    async def embed_stream(
        self,
        sequences: Sequence[str],
        accessions: Sequence[str],
        batch_size: int = 8,
        chunk_size: int = 1000,
    ) -> AsyncIterator[list[Embedding]]:
        """
        Stream embeddings in chunks for large datasets.

        Yields batches of embeddings as they complete, enabling
        immediate database writes without waiting for all 200k.

        Args:
            sequences: Protein sequences
            accessions: Accession IDs
            batch_size: Sequences per batch per GPU
            chunk_size: Sequences per yielded chunk

        Yields:
            Lists of Embedding objects as they complete
        """
        if not self._initialized:
            await self.initialize()

        if len(sequences) != len(accessions):
            raise ValueError("sequences and accessions must have same length")

        print(f"[ESM2] Streaming {len(sequences)} sequences in chunks of {chunk_size}")

        for chunk_start in range(0, len(sequences), chunk_size):
            chunk_end = min(chunk_start + chunk_size, len(sequences))
            chunk_seqs = sequences[chunk_start:chunk_end]
            chunk_accs = accessions[chunk_start:chunk_end]

            # Await the batch computation
            embeddings = await self.embed_batch(chunk_seqs, chunk_accs, batch_size)
            
            print(f"[ESM2] Completed chunk {chunk_start}-{chunk_end} ({chunk_end}/{len(sequences)})")
            
            # Yield results immediately
            yield embeddings

    def _distribute_work(
        self, sequences: Sequence[str], accessions: Sequence[str]
    ) -> list[tuple[Sequence[str], Sequence[str]]]:
        """Distribute sequences round-robin across devices."""
        num_devices = len(self.devices)
        return [
            (sequences[i::num_devices], accessions[i::num_devices])
            for i in range(num_devices)
        ]

    def _reorder_results(
        self, device_results: list[list[Embedding]], total: int
    ) -> list[Embedding]:
        """Reorder round-robin results back to original sequence order."""
        ordered: list[Embedding | None] = [None] * total
        flat = [emb for device_embs in device_results for emb in device_embs]

        idx = 0
        for device_idx in range(len(self.devices)):
            for seq_idx in range(device_idx, total, len(self.devices)):
                if idx < len(flat):
                    ordered[seq_idx] = flat[idx]
                    idx += 1

        return [e for e in ordered if e is not None]

    async def _process_device_work(
        self,
        device_idx: int,
        sequences: Sequence[str],
        accessions: Sequence[str],
        batch_size: int,
    ) -> list[Embedding]:
        """Process all sequences assigned to a device."""
        if not sequences:
            return []

        print(f"[ESM2] Device {device_idx} ({self.devices[device_idx]}): {len(sequences)} sequences")

        embeddings: list[Embedding] = []

        for batch_start in range(0, len(sequences), batch_size):
            batch_end = min(batch_start + batch_size, len(sequences))
            batch_seqs = sequences[batch_start:batch_end]
            batch_accs = accessions[batch_start:batch_end]

            batch_embs = await asyncio.to_thread(
                self._compute_embeddings,
                device_idx,
                batch_seqs,
                batch_accs,
            )

            embeddings.extend(batch_embs)
            print(f"[ESM2] Device {device_idx}: {len(embeddings)}/{len(sequences)} done")

        return embeddings

    def _compute_embeddings(
        self,
        device_idx: int,
        sequences: Sequence[str],
        accessions: Sequence[str],
    ) -> list[Embedding]:
        """
        Compute embeddings for a batch (GPU-bound, runs in thread pool).

        This is the core computation that runs synchronously in a thread.
        """
        model = self.models[device_idx]
        tokenizer = self.tokenizers[device_idx]
        device = self.devices[device_idx]

        embeddings: list[Embedding] = []

        with torch.no_grad():
            for seq, acc in zip(sequences, accessions, strict=True):
                # Tokenize
                inputs = tokenizer(
                    seq,
                    return_tensors="pt",
                    padding=False,
                    truncation=True,
                    max_length=1024,
                ).to(device)

                # Forward pass
                outputs = model(**inputs, output_hidden_states=True)

                # Extract layer
                hidden_states = (
                    outputs.last_hidden_state
                    if self.layer_index == -1
                    else outputs.hidden_states[self.layer_index]
                )

                # Pool
                attention_mask = inputs.get("attention_mask")
                pooled = self.pooling(hidden_states, attention_mask)

                # To numpy and normalize
                vector = pooled.cpu().numpy()[0]
                if self.normalize:
                    vector = normalize_embeddings(vector)

                # Create Embedding object
                embeddings.append(
                    Embedding(
                        model_name=self.model_name,
                        layer_index=self.layer_index,
                        pooling_method=self.pooling.__name__,
                        vector=vector.tolist(),
                        n_dims=len(vector),
                    )
                )

        return embeddings

    async def __aenter__(self) -> "ESM2Embedder":
        await self.initialize()
        return self

    async def __aexit__(self, exc_type, exc_val, exc_tb) -> None:  # type: ignore
        await self.cleanup()

async def embed_proteins_async(
    sequences: Sequence[str],
    accessions: Sequence[str],
    model_name: str = "facebook/esm2_t33_650M_UR50D",
    pooling_method: Callable = mean_pooling,
    batch_size: int = 8,
    normalize: bool = True,
) -> list[Embedding]:
    """
    High-level async function to embed proteins.

    Args:
        sequences: Protein sequences
        accessions: Optional IDs
        model_name: ESM-2 model identifier
        pooling_method: Pooling function (default: mean_pooling)
        batch_size: Batch size per GPU
        normalize: L2-normalize embeddings

    Returns:
        List of Embedding objects

    Example:
        embeddings = await embed_proteins_async(
            sequences=["MKTV", "ARND"],
            accessions=["P12345", "P67890"],
        )
    """
    assert len(sequences) == len(accessions), "sequences and accessions must have the same length"

    async with ESM2Embedder(
        model_name=model_name,
        pooling=pooling_method,
        normalize=normalize,
    ) as embedder:
        return await embedder.embed_batch(sequences, accessions, batch_size)

if __name__ == "__main__":
    from rich import print

    res = asyncio.run(
        embed_proteins_async(
            ["MKTVMKTVMKTVMKTVMKTVMKTVMKTVMKTV", "ARNDARNDARNDARNDARNDARNDARNDARND"],
            ["P12345", "P67890"],
        )
    )
    print(res)
