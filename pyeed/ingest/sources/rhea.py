"""Rhea database adapter for fetching reaction and molecule information."""

from __future__ import annotations

import asyncio
import csv
import io
import re
from collections import defaultdict
from collections.abc import AsyncIterator, Callable, Iterable
from typing import Any

import httpx
from loguru import logger

from ..model import Molecule, PyeedBase, Reaction

_RHEA_TABLE_URL = "https://www.rhea-db.org/rhea/"
_RHEA_COLS = "rhea-id,equation,chebi-id"
_HEADER_MAP = {
    "Reaction identifier": "rhea-id",
    "Equation": "equation",
    "ChEBI identifier": "chebi-id",
}

# arrows; treat '=' (Rhea table style) as reversible
_REV = re.compile(r"\s*(<=>|=|⇌|⥨|⥫)\s*")
_FWD = re.compile(r"\s*(=>|->|→|⟶)\s*")
_PLUS_OUTSIDE_PARENS = re.compile(r"\+(?![^()]*\))")


class RheaClient:
    def __init__(self, user_agent: str = "pyeed/1.0", timeout_s: float = 20.0) -> None:
        self._hdr = {"User-Agent": user_agent}
        self._timeout = httpx.Timeout(timeout_s)

    async def _fetch_table_row(self, rhea_id: str) -> dict[str, str] | None:
        rid = rhea_id.split(":")[-1]
        params = {"query": f"RHEA:{rid}", "columns": _RHEA_COLS, "format": "tsv", "limit": "1"}
        async with httpx.AsyncClient() as client:
            r = await client.get(
                _RHEA_TABLE_URL, params=params, headers=self._hdr, timeout=self._timeout
            )
            r.raise_for_status()
            text = r.text
        reader = csv.DictReader(io.StringIO(text), delimiter="\t")
        for raw in reader:
            # normalize headers to canonical keys
            return {
                _HEADER_MAP.get(h, h): (raw.get(h, "") or "").strip()
                for h in (reader.fieldnames or [])
            }
        return None

    async def _fetch_json(self, rhea_id: str) -> dict:
        rid = rhea_id.split(":")[-1]
        params = {
            "query": f"RHEA:{rid}",
            "columns": "rhea-id,equation,transport,balanced",
            "format": "json",
            "limit": "1",
        }
        async with httpx.AsyncClient() as client:
            r = await client.get(
                _RHEA_TABLE_URL, params=params, headers=self._hdr, timeout=self._timeout
            )
            r.raise_for_status()
            data = r.json()
        # shape: {"results":[{"rhea-id":"RHEA:16505","equation":"...",
        #                     "balanced":true,"transport":false}], "count":1}
        return (data.get("results") or [{}])[0]

    @staticmethod
    def _split_equation(eq: str) -> tuple[str, str, bool]:
        if m := _REV.search(eq):
            i, j = m.span()
            return eq[:i].strip(), eq[j:].strip(), True
        if m := _FWD.search(eq):
            i, j = m.span()
            return eq[:i].strip(), eq[j:].strip(), False
        return eq.strip(), "", False

    @staticmethod
    def _split_side(side: str) -> list[str]:
        # split on '+' that are not inside parentheses
        parts = _PLUS_OUTSIDE_PARENS.split(side)
        # strip stoich coefficients like "2 H2O"
        out: list[str] = []
        for part in parts:
            stripped = part.strip()
            if not stripped:
                continue
            out.append(re.sub(r"^\s*\d+\s+", "", stripped).strip())
        return out

    def _extract_reaction(
        self,
        row: dict[str, str] | None,
        meta: dict[str, Any] | None,
        rhea_id: str | None = None,
    ) -> list[Reaction]:
        """Extract Reaction object from Rhea data.

        Args:
            row: Table row dictionary with rhea-id, equation, chebi-id
            meta: JSON metadata dictionary with rhea-id, equation, transport, balanced
            rhea_id: Fallback Rhea ID if not in row

        Returns:
            List containing single Reaction object, or empty list if data invalid
        """
        if not row:
            return []

        eq = row.get("equation", "")
        if not eq:
            return []

        # Parse equation
        left, right, rev_from_eq = self._split_equation(eq)
        lhs = self._split_side(left)
        rhs = self._split_side(right)

        # Extract ChEBI IDs
        chebi_ids_str = row.get("chebi-id", "") or ""
        chebi_ids = [x.strip() for x in chebi_ids_str.split(";") if x.strip()]

        # Order assumption: ids follow equation participants left→right
        n_lhs = len(lhs)
        lhs_ids = chebi_ids[:n_lhs]
        rhs_ids = chebi_ids[n_lhs : n_lhs + len(rhs)]

        # Determine reversibility from metadata if available
        reversible = rev_from_eq
        if meta:
            balanced = meta.get("balanced")
            if isinstance(balanced, bool):
                reversible = balanced

        # Get Rhea ID
        rid = row.get("rhea-id") or meta.get("rhea-id") if meta else None
        if not rid and rhea_id:
            rid = f"RHEA:{rhea_id.split(':')[-1]}"
        if not rid:
            return []

        return [
            Reaction(
                id=rid,
                description=eq or None,
                substrate_ids=lhs_ids,
                product_ids=rhs_ids,
                reversible=reversible,
            )
        ]

    def _extract_molecules(
        self,
        row: dict[str, str] | None,
        chebi_enricher: Callable[[str], dict[str, Any] | None] | None = None,
    ) -> list[Molecule]:
        """Extract Molecule objects from Rhea data.

        Args:
            row: Table row dictionary with chebi-id
            chebi_enricher: Optional callback to enrich molecule data

        Returns:
            List of unique Molecule objects (deduplicated by ChEBI ID)
        """
        if not row:
            return []

        # Extract all ChEBI IDs
        chebi_ids_str = row.get("chebi-id", "") or ""
        chebi_ids = [x.strip() for x in chebi_ids_str.split(";") if x.strip()]

        # Deduplicate
        seen_ids: set[str] = set()
        molecules: list[Molecule] = []

        for cid in chebi_ids:
            if not cid or cid in seen_ids:
                continue
            seen_ids.add(cid)

            # Enrich with callback if provided
            info = chebi_enricher(cid) if chebi_enricher else None

            molecules.append(
                Molecule(
                    id=cid,
                    name=(info or {}).get("name") if info else None,
                    smiles=(info or {}).get("smiles") if info else None,
                    inchi=(info or {}).get("inchi") if info else None,
                )
            )

        return molecules

    async def fetch_reactions(
        self,
        rhea_ids: Iterable[int | str],
        max_concurrent: int = 10,
    ) -> AsyncIterator[tuple[dict[str, str] | None, dict[str, Any]]]:
        """Fetch multiple reactions concurrently.

        Makes concurrent requests for both table data and JSON metadata.
        Handles errors gracefully (logs and continues).

        Args:
            rhea_ids: Iterable of Rhea reaction IDs
            max_concurrent: Maximum number of concurrent requests

        Yields:
            Tuples of (table_row, meta_json) for each reaction
            (None values and errors are skipped)
        """
        semaphore = asyncio.Semaphore(max_concurrent)

        async def fetch_with_semaphore(
            rid: int | str,
        ) -> tuple[int | str, dict[str, str] | None, dict[str, Any], Exception | None]:
            """Fetch with error handling, returns (rhea_id, table_row, meta, error)."""
            async with semaphore:
                try:
                    # Fetch both table row and JSON metadata concurrently
                    table_row, meta_json = await asyncio.gather(
                        self._fetch_table_row(str(rid)),
                        self._fetch_json(str(rid)),
                    )
                    return (rid, table_row, meta_json, None)
                except httpx.HTTPError as e:
                    logger.warning(f"HTTP error fetching Rhea {rid}: {e}")
                    return (rid, None, {}, e)
                except ValueError as e:
                    logger.warning(f"Invalid Rhea ID {rid}: {e}")
                    return (rid, None, {}, e)
                except Exception as e:
                    logger.error(f"Unexpected error fetching Rhea {rid}: {e}", exc_info=True)
                    return (rid, None, {}, e)

        # Create tasks for all Rhea IDs
        tasks = [fetch_with_semaphore(rid) for rid in rhea_ids]

        # Process results as they complete
        for coro in asyncio.as_completed(tasks):
            rhea_id, table_row, meta_json, error = await coro
            if table_row is not None:
                yield (table_row, meta_json)
            elif error is None:
                # 404 or empty response - reaction not found
                logger.debug(f"Rhea reaction {rhea_id} not found")

    def map(
        self,
        row: dict[str, str] | None,
        meta: dict[str, Any] | None,
        chebi_enricher: Callable[[str], dict[str, Any] | None] | None = None,
        rhea_id: str | None = None,
    ) -> dict[str, list[PyeedBase]]:
        """Map Rhea data to dictionary of PyeedBase objects.

        Args:
            row: Table row dictionary from _fetch_table_row()
            meta: JSON metadata dictionary from _fetch_json()
            chebi_enricher: Optional callback to enrich molecule data
            rhea_id: Fallback Rhea ID if not in row

        Returns:
            Dictionary mapping class names to lists of PyeedBase objects
        """
        results: dict[str, list[PyeedBase]] = defaultdict(list)

        # Extract reaction
        reactions = self._extract_reaction(row, meta, rhea_id)
        results[Reaction.__name__].extend(reactions)

        # Extract molecules
        molecules = self._extract_molecules(row, chebi_enricher)
        results[Molecule.__name__].extend(molecules)

        return results

    async def get_reaction(
        self,
        rhea_id: str,
        chebi_enricher: Callable[[str], dict | None] | None = None,
    ) -> Reaction | None:
        """Fetch and return a single Reaction object.

        This is a convenience method that maintains backward compatibility.
        For new code, consider using map() directly.

        Args:
            rhea_id: Rhea reaction ID (e.g., "RHEA:32459")
            chebi_enricher: Optional callback to enrich molecule data

        Returns:
            Reaction object, or None if not found
        """
        row = await self._fetch_table_row(rhea_id)
        meta = await self._fetch_json(rhea_id)

        results = self.map(row, meta, chebi_enricher, rhea_id)
        reactions = results.get(Reaction.__name__, [])

        return reactions[0] if reactions else None


## test rhea client
if __name__ == "__main__":
    from rich import print

    rhea_client = RheaClient()
    rx = asyncio.run(rhea_client.get_reaction("RHEA:32459"))
    print(rx)
