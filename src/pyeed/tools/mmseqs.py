import os
from typing import List

import httpx
from loguru import logger
from pydantic import BaseModel, Field
from rich.progress import Progress, SpinnerColumn, TextColumn

from pyeed.dbconnect import DatabaseConnector
from pyeed.tools.datamodels.mmseqs import Cluster


class MMSeqs(BaseModel):
    """
    Performs sequence clustering using MMSeqs2.
    The MMSeqs service is run in a Docker container, which is part
    of the pyeed docker service.
    """

    min_seq_id: float = Field(
        default=0.5, ge=0.0, le=1.0, description="Minimum sequence identity (0-1)"
    )
    coverage: float = Field(
        default=0.8, ge=0.0, le=1.0, description="Minimum coverage (0-1)"
    )
    cov_mode: int = Field(
        default=0,
        ge=0,
        le=2,
        description="Coverage mode (0: bidirectional, 1: query, 2: target)",
    )
    threads: int = Field(
        default_factory=lambda: max(1, int((os.cpu_count() or 1) - 1)),
        ge=1,
        description="Number of CPU threads to use",
    )
    cluster_mode: int = Field(
        default=0,
        ge=0,
        le=2,
        description="Cluster mode (0: set-cover, 1: connected-component, 2: greedy)",
    )
    sensitivity: float = Field(
        default=7.5,
        ge=1.0,
        le=9.0,
        description="Sensitivity: 1.0 faster; 7.5 default; 9.0 more sensitive",
    )
    seq_id_mode: int = Field(
        default=0,
        ge=0,
        le=1,
        description="Sequence identity definition (0: alignment length, 1: shorter sequence)",
    )
    rescore_mode: int = Field(
        default=0,
        ge=0,
        le=1,
        description="Rescore overlapping alignments (0: no, 1: yes)",
    )

    def cluster_from_db(
        self,
        accession_ids: list[str],
        dbConnector: DatabaseConnector,
    ) -> List[Cluster]:
        """Cluster sequences from pyeed database.

        Args:
            accession_ids: List of accession IDs
            dbConnector: DatabaseConnector instance

        Returns:
            Clustering results in Cluster objects
        """
        # Get sequences from db
        query = """
        MATCH (p:Protein)
        WHERE p.accession_id IN $ids
        RETURN p.accession_id AS accession_id, p.sequence AS sequence
        """
        sequence_dict = {
            p["accession_id"]: p["sequence"]
            for p in dbConnector.execute_read(query, {"ids": accession_ids})
        }

        return self.cluster_sequence_dict(sequence_dict)

    def cluster_sequence_dict(
        self,
        sequences: dict[str, str],
    ) -> List[Cluster]:
        """Cluster sequences using MMSeqs2.

        Args:
            sequences: Sequence dictionary with sequence id as key and sequence as value
            dbConnector: DatabaseConnector instance

        Returns:
            Clustering results in Cluster objects
        """
        query = self._dict_to_multifasta(sequences)
        response = self._run_mmseqs_service(query)
        sanitized = self._sanitize_response(response)
        return self._parse_clustering_output(sanitized)

    @staticmethod
    def _dict_to_multifasta(sequences: dict[str, str]) -> str:
        """Convert sequence dictionary to multifasta format.

        Args:
            sequences: Sequence dictionary with sequence id as key and sequence as value

        Returns:
            Multifasta formatted sequences
        """
        return "\n".join(f">{k}\n{v}" + "\n" for k, v in sequences.items())

    def _run_mmseqs_service(self, query: str, timeout: int = 6000) -> httpx.Response:
        """Run the MMSeqs service with the provided parameters."""
        try:
            with Progress(
                SpinnerColumn(),
                TextColumn("[progress.description]{task.description}"),
                transient=True,
            ) as progress:
                progress.add_task(
                    description="Running MMSeqs2 clustering...", total=None
                )
            response = httpx.post(
                "http://localhost:8001/cluster",
                json={"query": query, "params": self.model_dump()},
                timeout=timeout,
            )
            self._check_clustering_success(response)
            return response
        except httpx.ConnectError as e:
            logger.error(f"Connection error: {e}")
            raise httpx.ConnectError("PyEED Docker Service not running") from e

    @staticmethod
    def _sanitize_response(response: httpx.Response) -> str:
        """Sanitize the response to remove any unwanted characters."""
        stripped = response.text.strip('"')
        decoded = stripped.encode("utf-8").decode("unicode_escape")
        return decoded

    @staticmethod
    def _check_clustering_success(response: httpx.Response) -> None:
        """Check if the response is successful."""
        if not response.status_code == 200:
            raise ValueError(f"MMSeqs clustering failed: {response.json()}")

    @staticmethod
    def _parse_clustering_output(cluster_string: str) -> List[Cluster]:
        """Parse MMSeqs clustering output into Cluster objects."""

        # Parse TSV output into clusters
        clusters: List[Cluster] = []

        for line in cluster_string.splitlines():
            if not line.strip():
                continue

            parts = line.strip().split("\t")
            if len(parts) == 1:
                # Single ID means sequence is its own representative
                rep_id = member_id = parts[0]
            else:
                # Two IDs mean rep_id clusters member_id
                rep_id, member_id = parts

            # Find or create cluster
            for cluster in clusters:
                if cluster.representative_id == rep_id:
                    cluster.represented_ids.append(member_id)
                    break
            else:
                clusters.append(
                    Cluster(representative_id=rep_id, represented_ids=[member_id])
                )

        return clusters
