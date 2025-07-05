"""SaProt model for mutation effect prediction using Foldseek tokens."""

from __future__ import annotations

from typing import Dict

import torch
from transformers import EsmConfig, EsmForMaskedLM, EsmTokenizer

from ..utils import get_hf_token

AA_LIST = [
    "A",
    "C",
    "D",
    "E",
    "F",
    "G",
    "H",
    "I",
    "K",
    "L",
    "M",
    "N",
    "P",
    "Q",
    "R",
    "S",
    "T",
    "V",
    "W",
    "Y",
]
FOLDSEEK_STRUC_VOCAB = "pynwrqhgdlvtmfsaeikc#"


class SaProtFoldseekMutationModel:
    """Lightweight wrapper for mutation effect prediction with SaProt."""

    def __init__(
        self,
        foldseek_path: str | None = None,
        config_path: str | None = None,
        load_pretrained: bool = True,
    ) -> None:
        self.foldseek_path = foldseek_path
        self.config_path = config_path or "westlake-repl/SaProt_650M_AF2"
        self.load_pretrained = load_pretrained

        token = get_hf_token()
        if load_pretrained:
            self.model = EsmForMaskedLM.from_pretrained(self.config_path, use_auth_token=token)
        else:
            cfg = EsmConfig.from_pretrained(self.config_path)
            self.model = EsmForMaskedLM(cfg)

        self.tokenizer = EsmTokenizer.from_pretrained(self.config_path, use_auth_token=token)
        self.device = torch.device("cpu")
        self.model.eval()

    # ------------------------------------------------------------------
    # Basic helpers
    # ------------------------------------------------------------------
    def to(self, device: str | torch.device) -> "SaProtFoldseekMutationModel":
        self.device = torch.device(device)
        self.model.to(self.device)
        return self

    def eval(self) -> None:  # pragma: no cover - wrapper
        self.model.eval()

    # ------------------------------------------------------------------
    # Mutation effect prediction utilities
    # ------------------------------------------------------------------
    def _mask_sequence(self, seq: str, mut_info: str) -> str:
        tokens = self.tokenizer.tokenize(seq)
        for single in mut_info.split(":"):
            pos = int(single[1:-1])
            tokens[pos - 1] = "#" + tokens[pos - 1][-1]
        return " ".join(tokens)

    def predict_mut(self, seq: str, mut_info: str) -> float:
        """Predict effect of one or more mutations."""
        mask_seq = self._mask_sequence(seq, mut_info)
        inputs = self.tokenizer(mask_seq, return_tensors="pt").to(self.device)
        with torch.no_grad():
            probs = self.model(**inputs).logits.softmax(dim=-1)

        score = 0.0
        for single in mut_info.split(":"):
            ori_aa, pos, mut_aa = single[0], int(single[1:-1]), single[-1]
            ori_st = self.tokenizer.get_vocab()[ori_aa + FOLDSEEK_STRUC_VOCAB[0]]
            mut_st = self.tokenizer.get_vocab()[mut_aa + FOLDSEEK_STRUC_VOCAB[0]]
            ori_prob = probs[0, pos, ori_st : ori_st + len(FOLDSEEK_STRUC_VOCAB)].sum()
            mut_prob = probs[0, pos, mut_st : mut_st + len(FOLDSEEK_STRUC_VOCAB)].sum()
            score += torch.log(mut_prob / ori_prob)
        return float(score.item())

    def predict_pos_mut(self, seq: str, pos: int) -> Dict[str, float]:
        """Predict mutation effect for all amino acids at a position."""
        tokens = self.tokenizer.tokenize(seq)
        ori_aa = tokens[pos - 1][0]
        tokens[pos - 1] = "#" + tokens[pos - 1][-1]
        mask_seq = " ".join(tokens)
        inputs = self.tokenizer(mask_seq, return_tensors="pt").to(self.device)
        with torch.no_grad():
            probs = self.model(**inputs).logits.softmax(dim=-1)[0, pos]

        ori_st = self.tokenizer.get_vocab()[ori_aa + FOLDSEEK_STRUC_VOCAB[0]]
        ori_prob = probs[ori_st : ori_st + len(FOLDSEEK_STRUC_VOCAB)].sum()
        scores = {}
        for mut_aa in AA_LIST:
            mut_st = self.tokenizer.get_vocab()[mut_aa + FOLDSEEK_STRUC_VOCAB[0]]
            mut_prob = probs[mut_st : mut_st + len(FOLDSEEK_STRUC_VOCAB)].sum()
            scores[f"{ori_aa}{pos}{mut_aa}"] = float(torch.log(mut_prob / ori_prob).item())
        return scores

    def predict_pos_prob(self, seq: str, pos: int) -> Dict[str, float]:
        """Return probabilities for all amino acids at a position."""
        tokens = self.tokenizer.tokenize(seq)
        tokens[pos - 1] = "#" + tokens[pos - 1][-1]
        mask_seq = " ".join(tokens)
        inputs = self.tokenizer(mask_seq, return_tensors="pt").to(self.device)
        with torch.no_grad():
            probs = self.model(**inputs).logits.softmax(dim=-1)[0, pos]

        scores = {}
        for aa in AA_LIST:
            st = self.tokenizer.get_vocab()[aa + FOLDSEEK_STRUC_VOCAB[0]]
            prob = probs[st : st + len(FOLDSEEK_STRUC_VOCAB)].sum()
            scores[aa] = float(prob.item())
        return scores

