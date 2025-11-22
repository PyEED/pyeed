"""UniProt taxonomy adapter for fetching organism taxonomy information."""

from __future__ import annotations

from collections.abc import AsyncIterator, Iterable
from typing import Any

import httpx
from loguru import logger
from tenacity import (
    retry,
    retry_if_exception_type,
    stop_after_attempt,
    wait_exponential_jitter,
)

from ..model import Taxon

TAXON_BASE_URL = "https://rest.uniprot.org/taxonomy"
MAX_BATCH_SIZE = 50  # Maximum number of taxon IDs per batch request


class UniProtTaxonomyAdapter:
    """Adapter for fetching taxonomy data from UniProt REST API.

    Uses the batch search endpoint to fetch multiple taxa in a single request,
    significantly reducing API calls and improving performance.
    """

    def __init__(self) -> None:
        self.headers = {"Accept": "application/json"}
        self.timeout = httpx.Timeout(20.0)

    def _build_taxon_query(self, taxon_ids: list[str]) -> str:
        """Build OR query string for batch search.

        Args:
            taxon_ids: List of taxon ID strings

        Returns:
            Query string in format: (id1 OR id2 OR ... OR idN)
        """
        return f"({' OR '.join(taxon_ids)})"

    @retry(
        wait=wait_exponential_jitter(0.5, 3),
        stop=stop_after_attempt(5),
        retry=retry_if_exception_type(httpx.HTTPError),
    )
    async def _fetch_taxa_batch(
        self,
        client: httpx.AsyncClient,
        taxon_ids: list[str],
    ) -> list[dict[str, Any]]:
        """Fetch a batch of taxa using the search endpoint.

        Args:
            client: httpx.AsyncClient
            taxon_ids: List of taxon ID strings (max 500)

        Returns:
            List of taxon record dictionaries from the results array

        Raises:
            httpx.HTTPError: If request fails after retries
            ValueError: If taxon_ids list is empty or exceeds MAX_BATCH_SIZE
        """
        if not taxon_ids:
            raise ValueError("taxon_ids list cannot be empty")
        if len(taxon_ids) > MAX_BATCH_SIZE:
            raise ValueError(
                f"taxon_ids list exceeds maximum batch size of {MAX_BATCH_SIZE}: "
                f"{len(taxon_ids)} provided"
            )

        query = self._build_taxon_query(taxon_ids)
        url = f"{TAXON_BASE_URL}/search"
        params = {"query": query, "size": MAX_BATCH_SIZE}

        r = await client.get(url, params=params, headers=self.headers, timeout=self.timeout)
        r.raise_for_status()

        data = r.json()
        results = data.get("results", [])
        return results

    async def fetch_taxa(
        self,
        client: httpx.AsyncClient,
        taxon_ids: Iterable[int | str],
    ) -> AsyncIterator[dict[str, Any]]:
        """Fetch multiple taxa using batch search endpoint.

        Uses the UniProt taxonomy search endpoint to fetch multiple taxa in a
        single request. If more than 500 taxon IDs are provided, splits them
        into multiple batch requests.

        Args:
            client: httpx.AsyncClient
            taxon_ids: Iterable of NCBI taxonomy IDs

        Yields:
            Taxon record dictionaries from the API results array

        Raises:
            ValueError: If any taxon_id is invalid (non-numeric)
        """
        # Convert to list and validate
        taxon_id_list: list[str] = []
        for tid in taxon_ids:
            taxon_str = str(tid).strip()
            if not taxon_str or not taxon_str.isdigit():
                logger.warning(f"Invalid taxon_id format: {tid!r} (must be numeric), skipping")
                continue
            taxon_id_list.append(taxon_str)

        if not taxon_id_list:
            return

        logger.debug(f"Fetching {len(taxon_id_list)} taxa in {MAX_BATCH_SIZE}-sized batches")
        # Split into batches of MAX_BATCH_SIZE
        for i in range(0, len(taxon_id_list), MAX_BATCH_SIZE):
            batch = taxon_id_list[i : i + MAX_BATCH_SIZE]
            try:
                results = await self._fetch_taxa_batch(client, batch)
                for result in results:
                    yield result
            except httpx.HTTPError as e:
                logger.error(f"HTTP error fetching taxon batch: {e}")
                # Continue with next batch even if one fails
            except Exception as e:
                logger.error(f"Unexpected error fetching taxon batch: {e}", exc_info=True)
                # Continue with next batch even if one fails

    def _extract_taxon_from_dict(self, d: dict[str, Any]) -> Taxon:
        """Extract Taxon object from API response dictionary.

        Args:
            d: Dictionary containing taxon data from API

        Returns:
            Taxon object
        """
        # Extract synonyms, filtering empty strings
        synonyms = d.get("synonyms", [])
        if not isinstance(synonyms, list):
            synonyms = []
        synonyms = [s for s in synonyms if isinstance(s, str) and s.strip()]

        # Handle otherNames field (used for main taxon)
        other_names = d.get("otherNames", [])
        if isinstance(other_names, list):
            other_names = [s for s in other_names if isinstance(s, str) and s.strip()]
            synonyms.extend(other_names)

        return Taxon(
            id=str(d.get("taxonId")),
            scientific_name=d.get("scientificName"),
            common_name=d.get("commonName"),
            rank=d.get("rank"),
            hidden=bool(d.get("hidden")),
            synonyms=synonyms,
        )

    def map(self, t: dict[str, Any] | None) -> list[Taxon]:
        """Map UniProt taxonomy record to list of Taxon objects representing lineage.

        Extracts the main taxon, parent taxon, and all lineage items into a
        deduplicated list of Taxon objects.

        Args:
            t: UniProt taxonomy record dictionary, or None

        Returns:
            List of Taxon objects (main taxon + parent + lineage)
        """
        if t is None:
            return []

        taxa: list[Taxon] = []
        seen_ids: set[str] = set()

        # Extract main taxon
        main_taxon_id = str(t.get("taxonId"))
        if main_taxon_id is not None:
            main_taxon = self._extract_taxon_from_dict(t)
            taxa.append(main_taxon)
            seen_ids.add(main_taxon_id)

        # Extract parent taxon
        parent = t.get("parent")
        if isinstance(parent, dict):
            parent_id = str(parent.get("taxonId"))
            if parent_id is not None and parent_id not in seen_ids:
                parent_taxon = self._extract_taxon_from_dict(parent)
                taxa.append(parent_taxon)
                seen_ids.add(parent_id)

        # Extract lineage items
        lineage = t.get("lineage", [])
        if isinstance(lineage, list):
            for lineage_item in lineage:
                if not isinstance(lineage_item, dict):
                    continue
                lineage_id = str(lineage_item.get("taxonId"))
                if lineage_id is not None and lineage_id not in seen_ids:
                    lineage_taxon = self._extract_taxon_from_dict(lineage_item)
                    taxa.append(lineage_taxon)
                    seen_ids.add(lineage_id)

        return taxa

    def extract_hierarchy_info(self, t: dict[str, Any] | None) -> dict[str, Any]:
        """Extract main_id, parent_id, lineage_ids for IS_A hierarchy construction.

        Args:
            t: UniProt taxonomy record dictionary, or None

        Returns:
            Dict with keys:
                - main_id: The queried taxon ID (as string)
                - parent_id: Direct parent taxon ID (as string, or None)
                - lineage_ids: List of lineage taxon IDs as strings (from root to subfamily)

        Example:
            For Homo sapiens (9606):
            {
                "main_id": "9606",
                "parent_id": "9605",
                "lineage_ids": ["131567", "2759", ..., "207598"]
            }
        """
        if t is None:
            return {"main_id": None, "parent_id": None, "lineage_ids": []}

        main_id = str(t.get("taxonId")) if t.get("taxonId") is not None else None
        parent = t.get("parent", {})
        parent_id = (
            str(parent.get("taxonId"))
            if isinstance(parent, dict) and parent.get("taxonId") is not None
            else None
        )

        lineage = t.get("lineage", [])
        lineage_ids = []
        if isinstance(lineage, list):
            for item in lineage:
                if isinstance(item, dict) and "taxonId" in item:
                    lineage_ids.append(str(item["taxonId"]))  # Convert to string

        return {
            "main_id": main_id,
            "parent_id": parent_id,
            "lineage_ids": lineage_ids,
        }
