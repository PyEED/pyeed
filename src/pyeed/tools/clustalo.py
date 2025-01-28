import httpx
from loguru import logger
from pydantic import BaseModel, Field

from pyeed.dbconnect import DatabaseConnector
from pyeed.tools.datamodels.mmseqs import MultipleSequenceAlignment, Sequence
from pyeed.tools.services import ServiceURL
from pyeed.tools.utility import dict_to_fasta


class ClustalOmega(BaseModel):
    """
    Performs multiple sequence alignment using Clustal Omega.
    The Clustal Omega service is run in a Docker container.
    """

    service_url: str = Field(default=ServiceURL.CLUSTALO.value)

    def align(self, sequences: dict[str, str]) -> MultipleSequenceAlignment:
        """
        Aligns multiple sequences using Clustal Omega.

        Args:
            sequences: Dictionary of sequences to align. Keys are sequence IDs,
                values are sequence strings.
        Returns:
            MultipleSequenceAlignment: The alignment result
        """
        try:
            data = dict_to_fasta(sequences)
            response = self._run_clustalo_service(data)
            sanitized = self._sanitize_response(response)
            return self._parse_alignment_output(sanitized)
        except Exception as e:
            logger.error(f"Alignment failed: {e}")
            raise

    def align_from_db(
        self, accession_ids: list[str], db: DatabaseConnector
    ) -> MultipleSequenceAlignment:
        """
        Aligns multiple sequences from a database using Clustal Omega.

        Args:
            accession_ids: List of sequence accession IDs to align
            db: Database connector to retrieve sequences

        Returns:
            MultipleSequenceAlignment: The alignment result

        Raises:
            ValueError: If no sequences found for given accession IDs
        """

        query = """
        MATCH (p:Protein)
        WHERE p.accession_id IN $ids
        RETURN p.accession_id AS accession_id, p.sequence AS sequence
        """
        sequence_dict = {
            p["accession_id"]: p["sequence"]
            for p in db.execute_read(query, {"ids": accession_ids})
        }

        return self.align(sequence_dict)

    def _run_clustalo_service(
        self, sequences: str, timeout: int = 3600
    ) -> httpx.Response:
        """Run the Clustal Omega service with the provided parameters."""
        try:
            files = {"file": ("input.fasta", sequences, "text/plain")}
            response = httpx.post(
                self.service_url,
                files=files,
                timeout=timeout,
            )
            self._check_alignment_success(response)
            return response
        except httpx.ConnectError as e:
            logger.debug(f"Connection error: {e}")
            # Try localhost if container name fails
            if "clustalo" in self.service_url:
                self.service_url = self.service_url.replace("clustalo", "localhost")
                return self._run_clustalo_service(sequences, timeout)
            raise httpx.ConnectError("PyEED Docker Service not running") from e

    @staticmethod
    def _sanitize_response(response: httpx.Response) -> str:
        """Sanitize the response to remove any unwanted characters."""
        stripped = response.text.strip('"')
        decoded = stripped.encode().decode("unicode_escape")
        return decoded

    @staticmethod
    def _check_alignment_success(response: httpx.Response) -> None:
        """Check if the response is successful."""
        if not response.status_code == 200:
            raise ValueError(f"Clustal Omega alignment failed: {response.json()}")

    @staticmethod
    def _parse_alignment_output(alignment_string: str) -> MultipleSequenceAlignment:
        """Parse Clustal Omega output into MultipleSequenceAlignment object."""
        sequences: list[Sequence] = []
        lines = alignment_string.splitlines()

        current_id = None
        current_sequence: list[str] = []

        for line in lines:
            line = line.strip()

            if not line:
                continue

            if line.startswith(">"):
                # Save previous sequence if exists
                if current_id is not None:
                    sequences.append(
                        Sequence(id=current_id, sequence="".join(current_sequence))
                    )
                # Start new sequence
                current_id = line[1:]  # Remove '>' prefix
                current_sequence = []
            else:
                # Add line to current sequence
                current_sequence.append(line)

        # Don't forget to add the last sequence
        if current_id is not None:
            sequences.append(
                Sequence(id=current_id, sequence="".join(current_sequence))
            )

        return MultipleSequenceAlignment(sequences=sequences)
