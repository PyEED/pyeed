import io
import os
from typing import Literal

import httpx
import pandas as pd
from loguru import logger
from pydantic import BaseModel, Field, field_validator
from rich.progress import Progress, SpinnerColumn, TextColumn


class Blast(BaseModel):
    """Performs BLAST search on a local database.
    The Blast service is run in a Docker container, which is part
    of the pyeed docker service.

    The serive needs to be configured and started before it can be used.
    """

    mode: Literal["blastp", "blastn"]
    db_path: str = Field(..., description="Path to BLAST database")
    db_name: str = Field(..., description="Name of BLAST database")
    evalue: float = Field(default=0.001, ge=0, description="E-value threshold")
    max_target_seqs: int = Field(
        default=100, gt=0, description="Maximum number of target sequences"
    )
    num_threads: int = Field(
        default_factory=lambda: max(1, int((os.cpu_count() or 1) - 1)),
        gt=0,
        description="Number of threads to use (defaults to CPU count - 1)",
    )

    @field_validator("mode")
    def validate_mode(cls, v: str) -> str:
        """Validate BLAST mode"""
        if v not in ["blastp", "blastn"]:
            raise ValueError("Mode must be either 'blastp' or 'blastn'")
        return v

    def search(
        self,
        sequence: str,
        timeout: int = 3600,
    ) -> pd.DataFrame:
        """Search for a sequence in a local BLAST database.

        Args:
            sequence (str): Sequence to search for
            timeout (int): Timeout for BLAST service. Defaults to 3600 seconds.

        Returns:
            pd.DataFrame: BLAST output
        """
        response = self.run_blast_service(timeout=timeout, sequence=sequence)
        if not response.status_code == 200:
            raise ValueError(f"BLAST failed: {response.json()}")
        return self._parse_blast_output(
            response.json()
        )  # Get string from JSON response

    def run_blast_service(
        self,
        sequence: str,
        timeout: int,
    ) -> httpx.Response:
        """Run the BLAST service with the provided parameters.

        Args:
            sequence (str): Sequence to search for
            timeout (int): Timeout for BLAST service

        Returns:
            httpx.Response: BLAST service response

        Raises:
            httpx.ConnectError: If the BLAST service is not running
        """
        params = self.model_dump()
        params["sequence"] = sequence
        logger.debug(f"Initializing BLAST search with params: {params}")

        with Progress(
            SpinnerColumn(),
            TextColumn("[progress.description]{task.description}"),
            transient=True,
        ) as progress:
            progress.add_task(description=f"Running {self.mode}", total=None)
            try:
                return httpx.post(
                    "http://localhost:6001/blast",
                    json=params,
                    timeout=timeout,
                )
            except httpx.ConnectError as e:
                logger.error(f"Connection error: {e}")
                raise httpx.ConnectError("PyEED Docker Service not running") from e

    @staticmethod
    def _parse_blast_output(result_str: dict[str, str]) -> pd.DataFrame:
        """Parse BLAST tab-delimited output into a DataFrame.

        Args:
            result_str (str): BLAST output as a string

        Returns:
            pd.DataFrame: Parsed BLAST output
        """
        df = pd.read_csv(
            io.StringIO(result_str["result"]),
            sep="\t",
            names=[
                "subject_id",
                "identity",
                "alignment_length",
                "mismatches",
                "gap_opens",
                "query_start",
                "query_end",
                "subject_start",
                "subject_end",
                "evalue",
                "bit_score",
            ],
        )
        # Remove query_sequence index and duplicates
        df = df.loc[df.index.repeat(1)]  # Convert index to simple range
        df.index = range(len(df))
        df = df.drop_duplicates(subset=["subject_id", "query_start", "query_end"])
        return df


if __name__ == "__main__":
    from pyeed.tools import Blast

    # Example protein sequence
    sequence = "MSEQVAAVAKLRAKASEAAKEAKAREAAKKLAEAAKKAKAKEAAKRAEAKLAEKAKAAKRAEAKAAKEAKRAAAKRAEAKLAEKAKAAK"

    # Initialize BLAST search
    blast = Blast(
        mode="blastp",  # Use blastp for protein sequences
        db_path="/usr/local/bin/data/test_db/",  # Path in Docker container
        db_name="protein_db",  # Name of your BLAST database
        evalue=0.1,  # E-value threshold
        max_target_seqs=10,  # Maximum number of hits to return
    )

    # Perform search
    results = blast.search(sequence)
    print(results)
