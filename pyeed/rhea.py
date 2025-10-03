# pyeed/rhea.py
from __future__ import annotations

import asyncio
import csv
import io
import re
from collections.abc import Callable

import httpx

from .model import Molecule, Reaction

_RHEA_TABLE_URL = "https://www.rhea-db.org/rhea/"
# We request only what we need; header labels are human-readable.
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
    def __init__(
        self, user_agent: str = "pyeed/1.0 (contact: you@example.org)", timeout_s: float = 20.0
    ) -> None:
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
        # shape: {"results":[{"rhea-id":"RHEA:16505","equation":"...","balanced":true,"transport":false}], "count":1}
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
        for p in parts:
            p = p.strip()
            if not p:
                continue
            out.append(re.sub(r"^\s*\d+\s+", "", p).strip())
        return out

    async def get_reaction(
        self,
        rhea_id: str,
        chebi_enricher: Callable[[str], dict | None] | None = None,
    ) -> Reaction | None:
        row = await self._fetch_table_row(rhea_id)
        if not row:
            return None

        eq = row.get("equation", "")
        left, right, rev_from_eq = self._split_equation(eq)
        lhs = self._split_side(left)
        rhs = self._split_side(right)

        chebi_ids = [x.strip() for x in (row.get("chebi-id", "") or "").split(";") if x.strip()]
        # Order assumption used by Rhea examples: ids follow equation participants left→right
        n_lhs = len(lhs)
        lhs_ids = chebi_ids[:n_lhs]
        rhs_ids = chebi_ids[n_lhs : n_lhs + len(rhs)]

        # Optionally refine reversibility/balanced via JSON
        meta = await self._fetch_json(rhea_id)
        reversible = (
            bool(meta.get("balanced")) if meta else rev_from_eq
        )  # or keep rev_from_eq if you prefer

        def make_mol(cid: str) -> Molecule:
            info = chebi_enricher(cid) if chebi_enricher else None
            return Molecule(
                chebi_id=cid,
                name=(info or {}).get("name") if info else None,
                smiles=(info or {}).get("smiles") if info else None,
                inchi=(info or {}).get("inchi") if info else None,
            )

        substrates = [make_mol(cid) for cid in lhs_ids]
        products = [make_mol(cid) for cid in rhs_ids]

        rid = row.get("rhea-id") or f"RHEA:{rhea_id.split(':')[-1]}"
        return Reaction(
            rhea_id=rid,
            description=eq or None,
            substrates=substrates,
            products=products,
            reversible=reversible,
        )


## test rhea client
if __name__ == "__main__":
    rhea_client = RheaClient()
    rx = asyncio.run(rhea_client.get_reaction("RHEA:32459"))
    print(rx)
