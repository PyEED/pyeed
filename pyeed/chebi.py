# adapter_chebi.py
from __future__ import annotations

import asyncio
from collections.abc import AsyncIterator, Iterable
from typing import Any

import httpx
from tenacity import (
    retry,
    retry_if_exception_type,
    stop_after_attempt,
    wait_exponential_jitter,
)

from .model import Molecule

OLS4_CHEBI_TERMS = "https://www.ebi.ac.uk/ols4/api/ontologies/chebi/terms"


class ChebiAdapter:
    """
    Async ChEBI fetcher using the OLS4 REST API.

    Notes:
    - OLS4 doesn't support multi-ID queries here; we achieve "batching" via
      bounded concurrency (multiple GETs in flight).
    - Fields like SMILES/InChI are exposed under term['annotation'] and are often lists.
    - All parsing is defensive; values are optional.
    """

    def __init__(self, *, timeout_s: float = 15.0) -> None:
        self.headers = {"Accept": "application/json"}
        self.timeout = httpx.Timeout(timeout_s)

    @retry(
        wait=wait_exponential_jitter(0.5, 3.0),
        stop=stop_after_attempt(5),
        retry=retry_if_exception_type(httpx.HTTPError),
        reraise=True,
    )
    async def _fetch_one(
        self,
        client: httpx.AsyncClient,
        chebi_id: str,
    ) -> dict[str, Any] | None:
        """
        Fetch a single ChEBI term by obo_id (e.g., 'CHEBI:16072').
        Returns the first matching term dict or None if not found.
        """
        params = {"obo_id": chebi_id}
        r = await client.get(
            OLS4_CHEBI_TERMS, params=params, headers=self.headers, timeout=self.timeout
        )
        r.raise_for_status()
        data = r.json() or {}

        # OLS4 embeds terms under _embedded.terms
        embedded = data.get("_embedded") or {}
        terms = embedded.get("terms") or []
        if not terms:
            return None

        # Typically there is exactly one matching term for an obo_id
        return terms[0]

    async def fetch_ids(
        self,
        client: httpx.AsyncClient,
        chebi_ids: Iterable[str],
        *,
        concurrency: int = 16,
    ) -> AsyncIterator[dict[str, Any]]:
        """
        Fetch many ChEBI terms concurrently with bounded concurrency.

        Yields raw term dicts (as returned by OLS4) for each ID that resolves.
        """

        sem = asyncio.Semaphore(concurrency)

        async def worker(cid: str) -> dict[str, Any] | None:
            async with sem:
                try:
                    return await self._fetch_one(client, cid)
                except httpx.HTTPStatusError as e:
                    err_code = 404
                    if e.response is not None and e.response.status_code == err_code:
                        return None
                    raise

        tasks: list[asyncio.Task[dict[str, Any] | None]] = []
        for cid in chebi_ids:
            tasks.append(asyncio.create_task(worker(cid)))

        for fut in asyncio.as_completed(tasks):
            term = await fut
            if term:
                yield term

    @staticmethod
    def _first_str(value: Any) -> str | None:
        """OLS 'annotation' values are often lists; return the first string if present."""
        if isinstance(value, str):
            return value
        if isinstance(value, list) and value and isinstance(value[0], str):
            return value[0]
        return None

    def map(self, term: dict[str, Any]) -> Molecule:
        """
        Map an OLS ChEBI term into our Molecule model.

        Fields:
        - chebi_id: term['obo_id'] (e.g., 'CHEBI:16072')
        - name:     term['label']
        - smiles/inchi: from term['annotation'] dict if present
        """
        chebi_id = term.get("obo_id") or term.get("oboId")  # be defensive
        if not isinstance(chebi_id, str):
            raise ValueError("ChEBI term missing 'obo_id'")

        name = term.get("label")
        if name is not None and not isinstance(name, str):
            name = None

        ann = term.get("annotation") or {}
        smiles = self._first_str(ann.get("smiles") or ann.get("SMILES"))
        inchi = self._first_str(ann.get("inchi") or ann.get("InChI") or ann.get("INCHI"))

        return Molecule(
            chebi_id=chebi_id,
            name=name,
            smiles=smiles,
            inchi=inchi,
            embedding=[],
        )


async def main() -> None:
    """Minimal async usage example for ChebiAdapter.

    Example:
        >>> import asyncio
        >>> asyncio.run(main())
    """
    adapter = ChebiAdapter()
    async with httpx.AsyncClient() as client:
        async for term in adapter.fetch_ids(client, ["CHEBI:15377", "CHEBI:15379"]):
            molecule = adapter.map(term)
            print("Fetched molecule", f"chebi_id: {molecule.chebi_id}, name: {molecule.name}")
            print(molecule)


if __name__ == "__main__":
    asyncio.run(main())
