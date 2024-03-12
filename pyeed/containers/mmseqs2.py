import logging
import random
import os
import re
import string
from typing import List
from pathlib import Path
from docker.models.containers import Container
from pydantic import Field
from pyeed.containers import AbstractContainer, ToolImage
from pyeed.core.abstractsequence import AbstractSequence

path_config = Path(__file__).parent.parent.parent / "logging.conf"
logging.config.fileConfig(path_config)
logger = logging.getLogger("pyeed")


class MMSeqs2(AbstractContainer):
    """
    ## About
    MMSeqs2 interfaces the MMSeqs2 (Many-against-Many sequence searching)
    tool for very fast and sensitive searching and clustering of large
    protein sequence sets.

    ## Attributes:
        identity (float): Minimum sequence identity threshold for clustering. Value
            should be between 0 and 1. Default is 0.
        coverage (float): Minimum fraction of aligned (covered) residues to consider a
            match. Value should be between 0 and 1. Default is 0.8.
        cov_mode (int): Defines the coverage mode used in comparisons. Ranges from 0 to 5
            with different meanings explained in the attribute's description. Default is 0.
        alignment_mode (int): Specifies how the alignment is computed. Ranges from 0 to 4,
            explained in the attribute's description. Default is 3.
        alignment_output_mode (int): Specifies the alignment output mode Ranges from 0 to 5,
            explained in the attribute's description. Default is 0.
        cluster_mode (int): Determines the clustering strategy to be used. Ranges from 0 to 3,
            explained in the attribute's description. Default is 0.
        mode (str): Search mode for MMSeqs2. Default is 'easy-cluster'.

    ## Example
        ```py
        clusterer = MMSeqs2(
            coverage=0.8,
            min_seq_id=0.5
        )
        representatives = clusterer.cluster(sequences)
        ```

    """

    _container_info: ToolImage = ToolImage.MMSEQS2
    mode: str = Field(
        description="Search mode",
        default="easy-cluster",
    )
    identity: float = Field(
        description="Minimum sequence identity threshold for clustering",
        default=0,
        ge=0,
        le=1,
    )
    coverage: float = Field(
        description="Minimum fraction of aligned (covered) residues to consider a match",
        default=0.8,
        gt=0,
        le=1,
    )
    cov_mode: int = Field(
        description="""
        0: coverage of query and target
        1: coverage of target
        2: coverage of query
        3: target seq. length has to be at least x% of query length
        4: query seq. length has to be at least x% of target length
        5: short seq. needs to be at least x% of the other seq. length [0]
        """,
        default=0,
        ge=0,
        le=5,
    )
    alignment_mode: int = Field(
        description="""
        How to compute the alignment:
        0: automatic
        1: only score and end_pos
        2: also start_pos and cov
        3: also seq.id
        4: only ungapped alignment
        """,
        default=3,
        ge=0,
        le=4,
    )
    alignment_output_mode: int = Field(
        description="""
        How to compute the alignment:
        0: automatic
        1: only score and end_pos
        2: also start_pos and cov
        3: also seq.id
        4: only ungapped alignment
        5: score only (output) cluster format
        """,
        default=0,
        ge=0,
        le=5,
    )
    cluster_mode: int = Field(
        description="""
        0: Set-Cover (greedy)
        1: Connected component (BLASTclust)
        2,3: Greedy clustering by sequence length (CDHIT)
        """,
        default=0,
        ge=0,
        le=3,
    )

    def create_file(self, data: str):
        """
        Sets up the input data for the MMSeqs2 container.

        Args:
            data (str): The sequence to be searched.
        """
        path = os.path.join(self._tempdir_path, "input.fasta")

        with open(path, "w") as file:
            file.write(data)

    def run_container(self, command: str, data: str) -> Container:
        """Runs a container from a given image and executes a command.
        Data is mounted to the container."""
        try:
            image = self.get_image()
            self.create_file(data=data)
            print("🏃 Clustering sequences with MMSeqs2...")
            print(f"╭── initial sequences: {data.count('>')}")
            print(f"├── min. coverage: {int(self.coverage*100)} %")
            print(f"╰── min. sequence identity: {int(self.identity*100)} %")

            logger.debug(f"🏃 Running {self._container_info.name}")
            self._client.containers.run(
                detach=False,
                image=image,
                command=command,
                name=self._container_info.name
                + "_"
                + "".join(random.choices(string.ascii_lowercase, k=5)),
                auto_remove=True,
                volumes={self._tempdir_path: {"bind": "/app", "mode": "rw"}},
            )
            # output = container.attach(stdout=True, stream=True, logs=True)
            # Set detatch True to access stream

            return self.extract_output_data()

        except Exception as e:
            logger.error(f"Error running {self._container_info} container: {e}")

    def setup_command(self, mode: str):
        """Sets up the command to run the MMSeqs2 container."""

        command = (
            "mmseqs "
            f"{mode} "
            f"/app/input.fasta "
            f"/app/result "
            f"temp "
            f"-c {self.coverage} "
            f"--cov-mode {self.cov_mode} "
            f"--min-seq-id {self.identity} "
            f"--alignment-mode {self.alignment_mode} "
            f"--alignment-output-mode {self.alignment_output_mode} "
            f"--cluster-mode {self.cluster_mode} "
        )

        return command

    @classmethod
    def easy_cluster(
        cls,
        sequences: List[AbstractSequence],
        identity: float = 0,
        coverage: float = 0.8,
        cov_mode: int = 0,
        alignment_mode: int = 3,
        alignment_output_mode: int = 0,
        cluster_mode: int = 0,
    ) -> List[AbstractSequence]:

        mode = "easy-cluster"
        instance = cls(
            identity=identity,
            coverage=coverage,
            cov_mode=cov_mode,
            alignment_mode=alignment_mode,
            alignment_output_mode=alignment_output_mode,
            cluster_mode=cluster_mode,
        )

        command = instance.setup_command(mode)

        multifasta = "\n".join(
            [f">{seq.source_id}\n{seq.sequence}" for seq in sequences]
        )

        result = instance.run_container(command, multifasta)

        ids = re.findall(r">([A-Z0-9\.]+) ", result["representatives"])

        return list(filter(lambda x: x.source_id in ids, sequences))

    def extract_output_data(self):
        results = {}
        with open(f"{self._tempdir_path}/result_rep_seq.fasta") as f:
            results["representatives"] = f.read()

        with open(f"{self._tempdir_path}/result_cluster.tsv") as f:
            results["cluster"] = f.read()

        n_clusters = results["representatives"].count(">")

        self._delete_temp_dir()
        print(
            f"🎉 Clustered intitial sequences in {n_clusters} representative sequences"
        )

        return results
