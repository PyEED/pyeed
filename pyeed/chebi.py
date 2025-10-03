# pyeed/chebi.py
from __future__ import annotations

import httpx
from pydantic import BaseModel, ConfigDict, Field

from .model import Molecule

__all__ = ["ChebiClient", "ChebiError"]


class ChebiError(Exception):
    """ChEBI-specific errors during API communication."""

    def __init__(self, message: str, cause: Exception | None = None) -> None:
        super().__init__(message)
        self.cause = cause


# --- API Response Models ---


class ChebiStructure(BaseModel):
    """Chemical structure data from ChEBI."""

    model_config = ConfigDict(frozen=True)

    id: int
    smiles: str | None = None
    standard_inchi: str | None = None
    standard_inchi_key: str | None = None
    wurcs: str | None = None
    is_r_group: bool


class ChebiName(BaseModel):
    """Individual name/synonym entry."""

    model_config = ConfigDict(frozen=True)

    name: str
    status: str
    type: str
    source: str
    ascii_name: str
    adapted: bool
    language_code: str


class ChebiNames(BaseModel):
    """Names and synonyms. All name types are optional."""

    model_config = ConfigDict(frozen=True)

    SYNONYM: list[ChebiName] | None = None
    IUPAC_NAME: list[ChebiName] | None = Field(None, alias="IUPAC NAME")
    INN: list[ChebiName] | None = None


class ChebiChemicalData(BaseModel):
    """Chemical formula and mass data."""

    model_config = ConfigDict(frozen=True)

    formula: str
    charge: int
    mass: str
    monoisotopic_mass: str


class ChebiEntryData(BaseModel):
    """Core data for a ChEBI entry."""

    model_config = ConfigDict(frozen=True)

    id: int
    chebi_accession: str
    name: str
    ascii_name: str
    stars: int
    definition: str
    names: ChebiNames
    chemical_data: ChebiChemicalData
    default_structure: ChebiStructure | None = None
    modified_on: str
    secondary_ids: list[str]
    is_released: bool


class ChebiEntryResult(BaseModel):
    """Individual ChEBI entry result."""

    model_config = ConfigDict(frozen=True)

    standardized_chebi_id: str
    primary_chebi_id: str
    exists: bool
    id_type: str
    data: ChebiEntryData


# --- Client ---


class ChebiClient:
    """Async client for the ChEBI API."""

    def __init__(
        self,
        user_agent: str = "pyeed/1.0",
        timeout_s: float = 20.0,
    ) -> None:
        """
        Initialize ChEBI client.

        Args:
            user_agent: User-Agent header for requests.
            timeout_s: Timeout in seconds for HTTP requests.
        """
        self._base_url = "https://www.ebi.ac.uk/chebi/backend/api/public/compounds/"
        self._headers = {"User-Agent": user_agent}
        self._timeout = httpx.Timeout(timeout_s)

    @staticmethod
    def _normalize_chebi_id(chebi_id: str) -> str:
        """Ensure ChEBI ID has 'CHEBI:' prefix."""
        return chebi_id if chebi_id.startswith("CHEBI:") else f"CHEBI:{chebi_id}"

    async def _fetch_raw(self, chebi_ids: list[str]) -> dict[str, ChebiEntryResult]:
        """
        Fetch raw API response for one or more ChEBI IDs.

        Args:
            chebi_ids: List of ChEBI IDs (with or without 'CHEBI:' prefix).

        Returns:
            Dict mapping ChEBI IDs to entry results.

        Raises:
            ChebiError: On API errors or connection failures.
        """
        normalized = [self._normalize_chebi_id(cid) for cid in chebi_ids]
        params = {"chebi_ids": ",".join(normalized)}

        try:
            async with httpx.AsyncClient() as client:
                resp = await client.get(
                    self._base_url,
                    params=params,
                    headers=self._headers,
                    timeout=self._timeout,
                )
                resp.raise_for_status()
                raw_data = resp.json()
        except httpx.HTTPStatusError as e:
            raise ChebiError(
                f"HTTP {e.response.status_code} fetching ChEBI IDs {normalized}",
                cause=e,
            ) from e
        except httpx.RequestError as e:
            raise ChebiError(f"Request failed for ChEBI IDs {normalized}", cause=e) from e

        if not raw_data:
            raise ChebiError(f"Empty response for ChEBI IDs {normalized}")

        try:
            return {key: ChebiEntryResult.model_validate(val) for key, val in raw_data.items()}
        except Exception as e:
            raise ChebiError(f"Failed to parse ChEBI response: {e}", cause=e) from e

    @staticmethod
    def _entry_to_molecule(entry: ChebiEntryResult) -> Molecule:
        """Convert a ChebiEntryResult to a Molecule."""
        struct = entry.data.default_structure
        return Molecule(
            chebi_id=entry.standardized_chebi_id,
            name=entry.data.ascii_name or None,
            smiles=struct.smiles if struct else None,
            inchi=struct.standard_inchi if struct else None,
        )

    async def get_molecule(self, chebi_id: str) -> Molecule:
        """
        Fetch a single ChEBI entry and return as Molecule.

        Args:
            chebi_id: ChEBI ID (with or without 'CHEBI:' prefix).

        Returns:
            Molecule instance.

        Raises:
            ChebiError: If the ID is not found or the request fails.
        """
        results = await self._fetch_raw([chebi_id])
        if not results:
            raise ChebiError(f"No data found for ChEBI ID {chebi_id}")
        entry = next(iter(results.values()))
        return self._entry_to_molecule(entry)

    async def get_molecules(self, chebi_ids: list[str]) -> list[Molecule]:
        """
        Fetch multiple ChEBI entries in batch and return as Molecules.

        Args:
            chebi_ids: List of ChEBI IDs (with or without 'CHEBI:' prefix).

        Returns:
            List of Molecule instances in the same order as input IDs.

        Raises:
            ChebiError: If any request fails.
        """
        if not chebi_ids:
            return []

        results = await self._fetch_raw(chebi_ids)
        # Maintain input order
        normalized = [self._normalize_chebi_id(cid) for cid in chebi_ids]
        molecules: list[Molecule] = []
        for cid in normalized:
            if cid in results:
                molecules.append(self._entry_to_molecule(results[cid]))
            else:
                # If a specific ID is missing, raise error or skip depending on use case
                raise ChebiError(f"ChEBI ID {cid} not found in batch response")
        return molecules


# --- Usage Example ---
if __name__ == "__main__":
    import asyncio

    from rich import print

    async def demo() -> None:
        client = ChebiClient()

        # Single fetch
        mol = await client.get_molecule("CHEBI:57844")  # l-Methionine
        print(f"Single: {mol}")

        # Batch fetch
        mols = await client.get_molecules(["CHEBI:57844", "CHEBI:16526"])  # water, CO2
        for m in mols:
            print(f"Batch: {m}")

    asyncio.run(demo())
