import os
from typing import List

from Bio import AlignIO
from Bio.Align import MultipleSeqAlignment

from pyeed.container import AbstractContainer, ToolImage


class ClustalOmega(AbstractContainer):
    """
    ClustalOmega is a class that manages a container for running the ClustalOmega tool.

    Attributes:
        _container_info (ToolContainer): The information about the ClustalOmega container.

    Methods:
        create_file(data: List[AbstractSequence]) -> str:
            Sets up the input data for the ClustalOmega container.

        setup_command() -> str:
            Sets up the command to run the ClustalOmega container.

        extract_output_data() -> MultiSequenceAlignment:
            Extracts the output data from the ClustalOmega container.

        align(sequences: Union[MultiSequenceAlignment, List[AbstractSequence]]) -> MultiSequenceAlignment:
            Aligns multiple sequences using the ClustalOmega container and returns the alignment result.
    """

    _container_info: ToolImage = ToolImage.CLUSTALO

    def create_file(self, data: List[str]):
        """
        Sets up the input data for the ClustalOmega container.

        Args:
            data (List[str]): List of FASTA formatted sequences to be aligned.
        """
        data = "\n".join(data)
        path = os.path.join(self._tempdir_path, "input.fasta")

        with open(path, "w") as file:
            file.write(data)

    def setup_command(self):
        """
        Sets up the command to run the ClustalOmega container.

        Returns:
            str: The command to run the ClustalOmega container.
        """
        threads = os.cpu_count()
        return f"clustalo -i /data/input.fasta -o /data/output.clu --outfmt=clu --threads={threads}"

    def extract_output_data(self) -> MultipleSeqAlignment:
        """
        Extracts the output data from the ClustalOmega container.

        Returns:
            MultiSequenceAlignment: The alignment result.
        """
        with open(os.path.join(self._tempdir_path, "output.clu"), "r") as file:
            alignment = AlignIO.read(file, "clustal")
        self._delete_temp_dir()

        return alignment

    def align(self, sequences: List[str]):
        """
        Aligns multiple sequences and returns the alignment result.

        Args:
            sequences (List[str]): List of FASTA formatted sequences to be aligned.

        Returns:
            MultiSequenceAlignment: The alignment result.
        """
        self.run_container(
            command=self.setup_command(),
            data=sequences,
        )
        return self.extract_output_data()
