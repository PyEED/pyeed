"""
ProtT5 model implementation for protein embeddings.
"""

from typing import List, Tuple, cast

import numpy as np
import torch
from numpy.typing import NDArray
from transformers import T5Model, T5Tokenizer

from ..base import BaseEmbeddingModel, normalize_embedding
from ..utils import get_hf_token, preprocess_sequence_for_prott5


class ProtT5EmbeddingModel(BaseEmbeddingModel):
    """ProtT5 model implementation."""

    def __init__(self, model_name: str, device: torch.device):
        super().__init__(model_name, device)

    def load_model(self) -> Tuple[T5Model, T5Tokenizer]:
        """Load ProtT5 model and tokenizer."""
        token = get_hf_token()

        full_model_name = (
            self.model_name
            if self.model_name.startswith("Rostlab/")
            else f"Rostlab/{self.model_name}"
        )

        model = T5Model.from_pretrained(full_model_name, use_auth_token=token)
        tokenizer = T5Tokenizer.from_pretrained(
            full_model_name, use_auth_token=token, do_lower_case=False
        )

        # Move to device
        model = model.to(self.device)

        self.model = model
        self.tokenizer = tokenizer

        return model, tokenizer

    def preprocess_sequence(self, sequence: str) -> str:
        """ProtT5 needs space-separated sequences with rare AAs mapped to X."""
        return preprocess_sequence_for_prott5(sequence)

    def get_batch_embeddings(
        self, sequences: List[str], pool_embeddings: bool = True
    ) -> List[NDArray[np.float64]]:
        """Get embeddings for a batch of sequences using ProtT5."""
        if self.model is None or self.tokenizer is None:
            self.load_model()

        # Type cast to ensure type checker knows they're not None
        model = cast(T5Model, self.model)
        tokenizer = cast(T5Tokenizer, self.tokenizer)

        # Preprocess sequences for ProtT5
        processed_sequences = [self.preprocess_sequence(seq) for seq in sequences]

        inputs = tokenizer.batch_encode_plus(
            processed_sequences,
            add_special_tokens=True,
            padding="longest",
            return_tensors="pt",
        )

        # Move inputs to device
        input_ids = inputs["input_ids"].to(self.device)
        attention_mask = inputs["attention_mask"].to(self.device)

        with torch.no_grad():
            # For ProtT5, use encoder embeddings for feature extraction
            # Create dummy decoder inputs (just the pad token)
            batch_size = input_ids.shape[0]
            decoder_input_ids = torch.full(
                (batch_size, 1),
                tokenizer.pad_token_id if tokenizer.pad_token_id is not None else 0,
                dtype=torch.long,
                device=self.device,
            )

            outputs = model(
                input_ids=input_ids,
                attention_mask=attention_mask,
                decoder_input_ids=decoder_input_ids,
            )

            # Get encoder last hidden state (encoder embeddings)
            hidden_states = outputs.encoder_last_hidden_state.cpu().numpy()

        if pool_embeddings:
            # Mean pooling across sequence length, excluding padding tokens
            embedding_list = []
            for i, hidden_state in enumerate(hidden_states):
                # Get actual sequence length (excluding padding)
                attention_mask_np = attention_mask[i].cpu().numpy()
                seq_len = attention_mask_np.sum()
                # Pool only over actual sequence tokens
                pooled_embedding = hidden_state[:seq_len].mean(axis=0)
                embedding_list.append(pooled_embedding)
            return embedding_list
        return list(hidden_states)

    def get_single_embedding_last_hidden_state(
        self, sequence: str
    ) -> NDArray[np.float64]:
        """Get last hidden state embedding for a single sequence."""
        if self.model is None or self.tokenizer is None:
            self.load_model()

        # Type cast to ensure type checker knows they're not None
        model = cast(T5Model, self.model)
        tokenizer = cast(T5Tokenizer, self.tokenizer)

        processed_sequence = self.preprocess_sequence(sequence)
        inputs = tokenizer.encode_plus(
            processed_sequence, add_special_tokens=True, return_tensors="pt"
        )

        input_ids = inputs["input_ids"].to(self.device)
        attention_mask = inputs["attention_mask"].to(self.device)

        # Create dummy decoder inputs
        decoder_input_ids = torch.full(
            (1, 1),
            tokenizer.pad_token_id if tokenizer.pad_token_id is not None else 0,
            dtype=torch.long,
            device=self.device,
        )

        with torch.no_grad():
            outputs = model(
                input_ids=input_ids,
                attention_mask=attention_mask,
                decoder_input_ids=decoder_input_ids,
            )

        # Get encoder last hidden state including special tokens
        embedding = outputs.encoder_last_hidden_state[0].detach().cpu().numpy()
        return np.asarray(embedding, dtype=np.float64)

    def get_single_embedding_all_layers(self, sequence: str) -> NDArray[np.float64]:
        """Get embeddings from all layers for a single sequence."""
        if self.model is None or self.tokenizer is None:
            self.load_model()

        # Type cast to ensure type checker knows they're not None
        model = cast(T5Model, self.model)
        tokenizer = cast(T5Tokenizer, self.tokenizer)

        processed_sequence = self.preprocess_sequence(sequence)
        inputs = tokenizer.encode_plus(
            processed_sequence, add_special_tokens=True, return_tensors="pt"
        )

        input_ids = inputs["input_ids"].to(self.device)
        attention_mask = inputs["attention_mask"].to(self.device)

        # Create dummy decoder inputs
        decoder_input_ids = torch.full(
            (1, 1),
            tokenizer.pad_token_id if tokenizer.pad_token_id is not None else 0,
            dtype=torch.long,
            device=self.device,
        )

        with torch.no_grad():
            outputs = model(
                input_ids=input_ids,
                attention_mask=attention_mask,
                decoder_input_ids=decoder_input_ids,
                output_hidden_states=True,
            )

        embeddings_list = []
        # Get all encoder hidden states
        encoder_hidden_states = outputs.encoder_hidden_states
        for layer_tensor in encoder_hidden_states:
            # Remove batch dimension but keep special tokens
            emb = layer_tensor[0].detach().cpu().numpy()
            emb = normalize_embedding(emb)
            embeddings_list.append(emb)

        return np.array(embeddings_list)

    def get_single_embedding_first_layer(self, sequence: str) -> NDArray[np.float64]:
        """Get first layer embedding for a single sequence."""
        if self.model is None or self.tokenizer is None:
            self.load_model()

        # Type cast to ensure type checker knows they're not None
        model = cast(T5Model, self.model)
        tokenizer = cast(T5Tokenizer, self.tokenizer)

        processed_sequence = self.preprocess_sequence(sequence)
        inputs = tokenizer.encode_plus(
            processed_sequence, add_special_tokens=True, return_tensors="pt"
        )

        input_ids = inputs["input_ids"].to(self.device)
        attention_mask = inputs["attention_mask"].to(self.device)

        # Create dummy decoder inputs
        decoder_input_ids = torch.full(
            (1, 1),
            tokenizer.pad_token_id if tokenizer.pad_token_id is not None else 0,
            dtype=torch.long,
            device=self.device,
        )

        with torch.no_grad():
            outputs = model(
                input_ids=input_ids,
                attention_mask=attention_mask,
                decoder_input_ids=decoder_input_ids,
                output_hidden_states=True,
            )

        # Get first encoder hidden state including special tokens
        embedding = outputs.encoder_hidden_states[0][0].detach().cpu().numpy()

        # Normalize the embedding
        embedding = normalize_embedding(embedding)
        return embedding

    def get_final_embeddings(self, sequence: str) -> NDArray[np.float64]:
        """
        Get final embeddings for ProtT5 with robust fallback.
        """
        try:
            embeddings = self.get_batch_embeddings([sequence], pool_embeddings=True)
            if embeddings and len(embeddings) > 0:
                return np.asarray(embeddings[0], dtype=np.float64)
            else:
                raise ValueError("Batch embeddings method returned empty results")
        except Exception as e:
            raise ValueError(f"ProtT5 embedding extraction failed: {e}")
