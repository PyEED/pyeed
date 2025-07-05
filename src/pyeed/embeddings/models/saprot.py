"""SaProt model implementation for protein embeddings."""

from typing import List, Tuple, cast

import numpy as np
import torch
from numpy.typing import NDArray
from transformers import EsmForMaskedLM, EsmTokenizer

from ..base import BaseEmbeddingModel, normalize_embedding
from ..utils import get_hf_token


class SaProtEmbeddingModel(BaseEmbeddingModel):
    """SaProt model implementation."""

    def __init__(self, model_name: str, device: torch.device):
        super().__init__(model_name, device)

    def load_model(self) -> Tuple[EsmForMaskedLM, EsmTokenizer]:
        """Load SaProt model and tokenizer."""
        token = get_hf_token()

        model = EsmForMaskedLM.from_pretrained(
            self.model_name, use_auth_token=token
        )
        tokenizer = EsmTokenizer.from_pretrained(
            self.model_name, use_auth_token=token
        )

        model = model.to(self.device)

        self.model = model
        self.tokenizer = tokenizer

        return model, tokenizer

    def preprocess_sequence(self, sequence: str) -> str:
        """SaProt doesn't need special preprocessing."""
        return sequence

    def get_batch_embeddings(
        self, sequences: List[str], pool_embeddings: bool = True
    ) -> List[NDArray[np.float64]]:
        """Get embeddings for a batch of sequences using SaProt."""
        if self.model is None or self.tokenizer is None:
            self.load_model()

        model = cast(EsmForMaskedLM, self.model)
        tokenizer = cast(EsmTokenizer, self.tokenizer)

        embeddings = []

        for sequence in sequences:
            inputs = tokenizer(
                sequence, padding=True, truncation=True, return_tensors="pt"
            ).to(self.device)

            with torch.no_grad():
                outputs = model.esm(**inputs, output_hidden_states=True)

            hidden_states = outputs.last_hidden_state.cpu().numpy()

            if pool_embeddings:
                embeddings.append(hidden_states.mean(axis=1)[0])
            else:
                embeddings.append(hidden_states)
        return embeddings

    def get_single_embedding_last_hidden_state(
        self, sequence: str
    ) -> NDArray[np.float64]:
        """Get last hidden state embedding for a single sequence."""
        if self.model is None or self.tokenizer is None:
            self.load_model()

        model = cast(EsmForMaskedLM, self.model)
        tokenizer = cast(EsmTokenizer, self.tokenizer)

        inputs = tokenizer(sequence, return_tensors="pt").to(self.device)

        with torch.no_grad():
            outputs = model.esm(**inputs)

        embedding = outputs.last_hidden_state[0, 1:-1, :].detach().cpu().numpy()
        return np.asarray(embedding, dtype=np.float64)

    def get_single_embedding_all_layers(self, sequence: str) -> NDArray[np.float64]:
        """Get embeddings from all layers for a single sequence."""
        if self.model is None or self.tokenizer is None:
            self.load_model()

        model = cast(EsmForMaskedLM, self.model)
        tokenizer = cast(EsmTokenizer, self.tokenizer)

        inputs = tokenizer(sequence, return_tensors="pt").to(self.device)

        with torch.no_grad():
            outputs = model.esm(**inputs, output_hidden_states=True)

        embeddings_list = []
        hidden_states = outputs.hidden_states

        for layer_tensor in hidden_states:
            emb = layer_tensor[0, 1:-1, :].detach().cpu().numpy()
            emb = normalize_embedding(emb)
            embeddings_list.append(emb)

        return np.array(embeddings_list)

    def get_single_embedding_first_layer(self, sequence: str) -> NDArray[np.float64]:
        """Get first layer embedding for a single sequence."""
        if self.model is None or self.tokenizer is None:
            self.load_model()

        model = cast(EsmForMaskedLM, self.model)
        tokenizer = cast(EsmTokenizer, self.tokenizer)

        inputs = tokenizer(sequence, return_tensors="pt").to(self.device)

        with torch.no_grad():
            outputs = model.esm(**inputs, output_hidden_states=True)

        embedding = outputs.hidden_states[0][0, 1:-1, :].detach().cpu().numpy()
        embedding = normalize_embedding(embedding)
        return embedding

    def get_final_embeddings(self, sequence: str) -> NDArray[np.float64]:
        """Get final embeddings for SaProt with robust fallback."""
        try:
            embeddings = self.get_batch_embeddings([sequence], pool_embeddings=True)
            if embeddings and len(embeddings) > 0:
                return np.asarray(embeddings[0], dtype=np.float64)
            raise ValueError("Batch embeddings method returned empty results")
        except Exception as e:
            raise ValueError(f"SaProt embedding extraction failed: {e}")
