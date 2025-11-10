# adapter_uniprot.py
from __future__ import annotations

import re
from collections.abc import AsyncIterator, Iterable
from typing import Any

import httpx
from tenacity import (
    retry,
    retry_if_exception_type,
    stop_after_attempt,
    wait_exponential_jitter,
)

from ..model import (
    Annotation,
    AnnotationType,
    GOAnnotation,
    Organism,
    Protein,
    Reaction,
)

INTERPRO_PATTERN = re.compile(r"^IPR\d{6}$")
UNIPROT_PATTERN = re.compile(
    r"^(?:[OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9](?:[A-Z][A-Z0-9]{2}[0-9]){1,2})$"
)
UNIPROT_SEARCH = "https://rest.uniprot.org/uniprotkb/search"
RETURN_FIELDS = ",".join(
    [
        "accession",
        "id",
        "sequence",
        "organism_id",
        "ft_domain",
        "ft_site",
        "ft_act_site",
        "ft_binding",
        "ec",
        "go",
        "rhea",
    ]
)


def _build_query_for_accessions(accessions: list[str]) -> str:
    """Concatenates the accessions into a query string"""
    return "(" + " OR ".join(f"accession:{acc}" for acc in accessions) + ")"


class UniProtAdapter:
    def __init__(self) -> None:
        self.headers = {"Accept": "application/json"}
        self.timeout = httpx.Timeout(20.0)

    @retry(
        wait=wait_exponential_jitter(0.5, 3),
        stop=stop_after_attempt(5),
        retry=retry_if_exception_type(httpx.HTTPError),
    )
    async def fetch_search_page(
        self,
        client: httpx.AsyncClient,
        query: str,
        fields: str = RETURN_FIELDS,
        size: int = 500,
    ) -> AsyncIterator[dict[str, Any]]:
        params = {
            "query": query,
            "format": "json",
            "fields": fields,
            "size": str(min(size, 500)),  # UniProt caps at 500/page
        }

        url: str | None = UNIPROT_SEARCH
        while url:
            r = await client.get(url, params=params, headers=self.headers, timeout=self.timeout)
            r.raise_for_status()
            data = r.json()

            for rec in data.get("results") or []:
                yield rec

            # Follow RFC5988 Link header: rel="next"
            nxt: str | None = None
            link = r.headers.get("Link")
            if link:
                for part in link.split(","):
                    if 'rel="next"' in part:
                        i = part.find("<") + 1
                        j = part.find(">", i)
                        if i > 0 and j > i:
                            nxt = part[i:j]
                            break
            # After first page, use the full next URL as-is
            url, params = (nxt, {}) if nxt else (None, {})

    async def fetch_accessions(
        self,
        client: httpx.AsyncClient,
        accessions: Iterable[str],
        chunk_size: int = 50,  # keep URL safe; tune as needed
        size_per_page: int = 50,
    ) -> AsyncIterator[dict[str, Any]]:
        """Fetches the proteins by accessions from UniProt
        Args:
            client: httpx.AsyncClient
            accessions: Iterable[str]
            chunk_size: int = 200,  # keep URL safe; tune as needed
            size_per_page: int = 500,

        Returns:
            AsyncIterator[Dict[str, Any]]
        """
        batch: list[str] = []
        for acc in accessions:
            batch.append(acc)
            if len(batch) >= chunk_size:
                async for rec in self.fetch_search_page(
                    client,
                    query=_build_query_for_accessions(batch),
                    size=size_per_page,
                ):
                    yield rec
                batch.clear()
        if batch:
            async for rec in self.fetch_search_page(
                client,
                query=_build_query_for_accessions(batch),
                size=size_per_page,
            ):
                yield rec

    def map(self, p: dict[str, Any]) -> Protein:
        seq_meta = p.get("sequence") or {}
        sequence = seq_meta.get("value")
        if not isinstance(sequence, str):
            raise ValueError(f"Entry {p.get('primaryAccession', '?')} has no sequence")

        desc = p.get("proteinDescription") or {}
        rec = desc.get("recommendedName") or {}

        name = (rec.get("fullName") or {}).get("value")
        ecs = rec.get("ecNumbers") or []
        ec_numbers = [ec.get("value") for ec in ecs if isinstance(ec.get("value"), str)]

        gos: list[GOAnnotation] = []
        for x in p.get("uniProtKBCrossReferences", []):
            if x.get("database") != "GO":
                continue
            gid = x.get("id", "")
            if not gid.startswith("GO:"):
                continue
            term = next(
                (d.get("value") for d in x.get("properties", []) if d.get("key") == "GoTerm"),
                None,
            )
            if term:
                gos.append(GOAnnotation(go_id=gid, term=term, definition=None))

        fmap = {
            "active site": AnnotationType.ACTIVE_SITE,
            "site": AnnotationType.SITE,
            "domain": AnnotationType.DOMAIN,
            "binding site": AnnotationType.BINDING_SITE,
            "metal ion-binding site": AnnotationType.BINDING_SITE,
            "nucleotide phosphate-binding region": AnnotationType.BINDING_SITE,
        }
        anns: list[Annotation] = []
        for feature in p.get("features", []):
            f_type = (feature.get("type") or "").lower().strip()
            annot_type = fmap.get(f_type)
            if not annot_type:
                continue
            loc = feature.get("location") or {}
            start = (loc.get("start") or {}).get("value")
            end = (loc.get("end") or {}).get("value")
            if not isinstance(start, int) or not isinstance(end, int) or start <= 0 or end < start:
                continue

            description = feature.get("description") or None
            ligand = (feature.get("ligand") or {}).get("name")
            custom = {"ligand": ligand} if ligand else {}

            anns.append(
                Annotation(
                    annotation_type=annot_type,
                    positions=list(range(start, end + 1)),
                    description=description,
                    custom=custom,
                )
            )

        rhea_rx = re.compile(r"^RHEA:\d+$")
        rx: list[Reaction] = []
        for c in p.get("comments", []):
            if (c.get("commentType") or "").lower() != "catalytic activity":
                continue
            rct = c.get("reaction") or {}
            for ref in rct.get("reactionCrossReferences", []):
                rid = ref.get("id", "")
                if ref.get("database") == "Rhea" and rhea_rx.match(rid):
                    rx.append(Reaction(rhea_id=rid, description=rct.get("name")))

        return Protein(
            sequence_id=p["primaryAccession"],
            name=name,
            sequence=sequence,
            seq_length=len(sequence),
            mol_weight=seq_meta.get("molWeight"),
            ec_numbers=ec_numbers,
            organisms=[Organism(tax_id=(p.get("organism") or {}).get("taxonId"))],  # type: ignore
            go_terms=gos,
            annotations=anns,
            reactions=rx,
            structure_ids=[],
            embeddings=[],
            custom={},
        )

    @retry(
        wait=wait_exponential_jitter(0.5, 3),
        stop=stop_after_attempt(5),
        retry=retry_if_exception_type(httpx.HTTPError),
    )
    async def fetch_accessions_by_interpro(
        self,
        client: httpx.AsyncClient,
        ipr: str,
        distinct: bool = True,
        limit: int = 1_000_000_000,
    ) -> list[str]:
        """
        One-shot SPARQL: return all UniProt accessions linked to an InterPro ID.
        No pagination; relies on a very large LIMIT and whatever cap the endpoint applies.

        Args:
            client: shared httpx.AsyncClient
            ipr: e.g. "IPR002133"
            distinct: use SELECT DISTINCT
            limit: upper bound for rows (kept huge to approximate 'all at once')

        Returns:
            List of accessions (strings).
        """
        assert ipr.startswith("IPR"), f"Expected InterPro ID like 'IPRxxxxx', got {ipr!r}"

        base_url = "https://sparql.uniprot.org/sparql"
        accept_headers = {"Accept": "application/sparql-results+json", **self.headers}

        q = (
            "PREFIX up:<http://purl.uniprot.org/core/> "
            "PREFIX uniprotkb:<http://purl.uniprot.org/uniprot/> "
            "PREFIX rdfs:<http://www.w3.org/2000/01/rdf-schema#> "
            f"SELECT {'DISTINCT ' if distinct else ''}"
            "(SUBSTR(STR(?p), STRLEN(STR(uniprotkb:)) + 1) AS ?acc) "
            "WHERE { "
            f"?p a up:Protein ; rdfs:seeAlso <http://purl.uniprot.org/interpro/{ipr}> . "
            "} "
            f"LIMIT {int(limit)}"
        )

        r = await client.get(
            base_url, params={"query": q}, headers=accept_headers, timeout=self.timeout
        )
        r.raise_for_status()

        data = r.json()
        bindings = data.get("results", {}).get("bindings", []) or []
        # Extract 'acc' safely
        accs = []
        for b in bindings:
            v = (b.get("acc") or {}).get("value")
            if isinstance(v, str):
                accs.append(v)

        return accs


async def _amain() -> None:
    adapter = UniProtAdapter()
    ipr = "IPR002133"
    async with httpx.AsyncClient() as client:
        accs = await adapter.fetch_all_accessions_sparql(client, ipr)
        return accs


if __name__ == "__main__":
    import asyncio

    length = len(asyncio.run(_amain()))
    print(f"Length: {length}")
