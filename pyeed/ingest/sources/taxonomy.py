"""UniProt taxonomy adapter for fetching organism taxonomy information."""

from __future__ import annotations

import asyncio
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


class UniProtTaxonomyAdapter:
    """Adapter for fetching taxonomy data from UniProt REST API.

    Note: UniProt taxonomy API only supports single requests per endpoint.
    For multiple taxon IDs, requests are made concurrently.
    """

    def __init__(self) -> None:
        self.headers = {"Accept": "application/json"}
        self.timeout = httpx.Timeout(20.0)

    @retry(
        wait=wait_exponential_jitter(0.5, 3),
        stop=stop_after_attempt(5),
        retry=retry_if_exception_type(httpx.HTTPError),
    )
    async def fetch_taxon(
        self,
        client: httpx.AsyncClient,
        taxon_id: int | str,
    ) -> dict[str, Any] | None:
        """Fetch a single taxon by NCBI taxonomy ID.

        Args:
            client: httpx.AsyncClient
            taxon_id: NCBI taxonomy ID (e.g., 9606 for human)

        Returns:
            Taxon record dictionary, or None if not found

        Raises:
            httpx.HTTPError: If request fails after retries
            ValueError: If taxon_id is invalid
        """
        # Validate and sanitize taxon_id
        taxon_str = str(taxon_id).strip()
        if not taxon_str or not taxon_str.isdigit():
            raise ValueError(f"Invalid taxon_id format: {taxon_id!r} (must be numeric)")

        # Construct URL safely (taxon_str is now validated as digits only)
        url = f"{TAXON_BASE_URL}/{taxon_str}"
        r = await client.get(url, headers=self.headers, timeout=self.timeout)

        if r.status_code == 404:
            return None

        r.raise_for_status()
        return r.json()

    async def fetch_taxa(
        self,
        client: httpx.AsyncClient,
        taxon_ids: Iterable[int | str],
        max_concurrent: int = 10,
    ) -> AsyncIterator[dict[str, Any]]:
        """Fetch multiple taxa concurrently.

        Since the API only supports single requests, this method:
        - Processes taxon IDs in batches
        - Makes concurrent requests within each batch
        - Yields results as they become available
        - Handles errors gracefully (logs and continues)

        Args:
            client: httpx.AsyncClient
            taxon_ids: Iterable of NCBI taxonomy IDs
            max_concurrent: Maximum number of concurrent requests

        Yields:
            Taxon record dictionaries (None values and errors are skipped)
        """
        semaphore = asyncio.Semaphore(max_concurrent)

        async def fetch_with_semaphore(
            tid: int | str,
        ) -> tuple[int | str, dict[str, Any] | None, Exception | None]:
            """Fetch with error handling, returns (taxon_id, result, error)."""
            async with semaphore:
                try:
                    result = await self.fetch_taxon(client, tid)
                    return (tid, result, None)
                except httpx.HTTPError as e:
                    logger.warning(f"HTTP error fetching taxon {tid}: {e}")
                    return (tid, None, e)
                except ValueError as e:
                    logger.warning(f"Invalid taxon_id {tid}: {e}")
                    return (tid, None, e)
                except Exception as e:
                    logger.error(f"Unexpected error fetching taxon {tid}: {e}", exc_info=True)
                    return (tid, None, e)

        # Create tasks for all taxon IDs
        tasks = [fetch_with_semaphore(tid) for tid in taxon_ids]

        # Process results as they complete
        for coro in asyncio.as_completed(tasks):
            taxon_id, result, error = await coro
            if result is not None:
                yield result
            elif error is None:
                # 404 - taxon not found (silently skipped)
                logger.debug(f"Taxon {taxon_id} not found (404)")

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
        seen_ids: set[int] = set()

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


if __name__ == "__main__":
    import asyncio

    from rich import print

    async def run() -> None:
        async with httpx.AsyncClient() as client:
            taxa = await UniProtTaxonomyAdapter().fetch_taxon(client, 9606)
            print(UniProtTaxonomyAdapter().map(taxa))

    asyncio.run(run())
