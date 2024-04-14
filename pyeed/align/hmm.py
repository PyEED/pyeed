import os
import pyhmmer
from typing import List, Union
from pydantic import BaseModel, Field

from pyeed.core.abstractsequence import AbstractSequence
from pyeed.core.alignment import Alignment
from pyeed.core.sequence import Sequence


class HMM(BaseModel):
    class Config:
        arbitrary_types_allowed = True

    name: str = Field(
        description="The name of the HMM",
        default="msa",
    )

    alignment: Alignment = Field(
        description="The alignment object to be used for HMM generation",
        default=None,
    )

    model: pyhmmer.plan7.HMM = Field(
        description="The HMM model",
        default=None,
    )

    def __init__(self, **data):
        super().__init__(**data)
        if self.alignment:
            self.model = self.build_model()

    def build_model(self) -> pyhmmer.easel.MSA:
        """
        Builds and returns an HMM model.

        Returns:
            pyhmmer.easel.MSA: The HMM model.
        """

        sequences = self._text_sequences(self.alignment.aligned_sequences)
        msa = pyhmmer.easel.TextMSA(
            name=self.name.encode(), sequences=sequences
        ).digitize(self.alphabet)
        builder = pyhmmer.plan7.Builder(self.alphabet)
        hmm, _, _ = builder.build_msa(msa, self.background)

        return hmm

    def search(
        self,
        sequence: Union[AbstractSequence, List[AbstractSequence]],
        **pipeline_kwargs,
    ):

        query_seqs = self._prepare_sequences(sequence)
        seq_block = pyhmmer.easel.DigitalSequenceBlock(self.alphabet, query_seqs)

        pipeline = pyhmmer.plan7.Pipeline(
            self.alphabet, self.background, **pipeline_kwargs
        )

        hits = pipeline.search_hmm(self.model, seq_block)

        return hits

    def save_model(self, path: str = None, verbose: bool = True):
        """
        Saves the HMM model to a file.

        Args:
            path (str): The path to save the HMM model.
        """
        if not path:
            path = os.getcwd()

        with open(f"{path}/{self.name}.hmm", "wb") as f:
            self.model.write(f)

        if verbose:
            print(f"💾 HMM saved to {path}/{self.name}.hmm")

    @classmethod
    def from_file(cls, path: str):
        """
        Loads an HMM model from a file.

        Args:
            path (str): The path to load the HMM model.

        Returns:
            HMM: The HMM model.
        """
        with pyhmmer.plan7.HMMFile(path) as hmm_file:
            model = hmm_file.read()
            print(type(model))

        return cls(model=model, name=model.name)

    def _prepare_sequences(
        self, sequences: Union[AbstractSequence, List[AbstractSequence]]
    ):
        if not isinstance(sequences, list):
            sequences = [sequences]

        if not all(isinstance(seq, AbstractSequence) for seq in sequences) and not all(
            isinstance(seq, Sequence) for seq in sequences
        ):
            raise ValueError("All sequences must be of type AbstractSequence")

        return [seq.digitize(self.alphabet) for seq in self._text_sequences(sequences)]

    def _text_sequences(
        self,
        sequences: Union[
            Sequence, List[Sequence], AbstractSequence, List[AbstractSequence]
        ],
    ) -> list[pyhmmer.easel.TextSequence]:
        return [
            pyhmmer.easel.TextSequence(
                name=f"{seq.source_id}".encode(), sequence=seq.sequence
            )
            for seq in sequences
        ]

    @property
    def alphabet(self):
        if any("M" in seq.sequence.upper() for seq in self.alignment.input_sequences):
            return pyhmmer.easel.Alphabet.amino()
        else:
            return pyhmmer.easel.Alphabet.dna()

    @property
    def background(self):
        return pyhmmer.plan7.Background(self.alphabet)
