# adapter_uniprot.py
from __future__ import annotations

import re
from typing import Any, AsyncIterator, Dict, Iterable, List, Optional

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


def _build_query_for_accessions(accessions: List[str]) -> str:
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
    ) -> AsyncIterator[Dict[str, Any]]:
        params = {
            "query": query,
            "format": "json",
            "fields": fields,
            "size": str(min(size, 500)),  # UniProt caps at 500/page
        }

        url: Optional[str] = UNIPROT_SEARCH
        while url:
            r = await client.get(
                url, params=params, headers=self.headers, timeout=self.timeout
            )
            r.raise_for_status()
            data = r.json()

            for rec in data.get("results") or []:
                yield rec

            # Follow RFC5988 Link header: rel="next"
            nxt: Optional[str] = None
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
            url, params = (nxt, None) if nxt else (None, None)

    # --- Convenience: chunk IDs to avoid giant URLs; yields records across chunks ---
    async def fetch_accessions(
        self,
        client: httpx.AsyncClient,
        accessions: Iterable[str],
        chunk_size: int = 50,  # keep URL safe; tune as needed
        size_per_page: int = 50,
    ) -> AsyncIterator[Dict[str, Any]]:
        """Fetches the proteins by accessions from UniProt
        Args:
            client: httpx.AsyncClient
            accessions: Iterable[str]
            chunk_size: int = 200,  # keep URL safe; tune as needed
            size_per_page: int = 500,

        Returns:
            AsyncIterator[Dict[str, Any]]
        """
        batch: List[str] = []
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

    # --- Mapping stays strict & defensive ---
    def map(self, p: Dict[str, Any]) -> Protein:
        seq_meta = p.get("sequence") or {}
        sequence = seq_meta.get("value")
        if not isinstance(sequence, str):
            raise ValueError(f"Entry {p.get('primaryAccession','?')} has no sequence")

        desc = p.get("proteinDescription") or {}
        rec = desc.get("recommendedName") or {}

        name = (rec.get("fullName") or {}).get("value")
        ecs = rec.get("ecNumbers") or []
        ec_numbers = [ec.get("value") for ec in ecs if isinstance(ec.get("value"), str)]

        gos: List[GOAnnotation] = []
        for x in p.get("uniProtKBCrossReferences", []):
            if x.get("database") != "GO":
                continue
            gid = x.get("id", "")
            if not gid.startswith("GO:"):
                continue
            term = next(
                (
                    d.get("value")
                    for d in x.get("properties", [])
                    if d.get("key") == "GoTerm"
                ),
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
        anns: List[Annotation] = []
        for feature in p.get("features", []):
            f_type = (feature.get("type") or "").lower().strip()
            annot_type = fmap.get(f_type)
            if not annot_type:
                continue
            loc = feature.get("location") or {}
            start = (loc.get("start") or {}).get("value")
            end = (loc.get("end") or {}).get("value")
            if (
                not isinstance(start, int)
                or not isinstance(end, int)
                or start <= 0
                or end < start
            ):
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
        rx: List[Reaction] = []
        for c in p.get("comments", []):
            if (c.get("commentType") or "").lower() != "catalytic activity":
                continue
            rct = c.get("reaction") or {}
            for ref in rct.get("reactionCrossReferences", []):
                rid = ref.get("id", "")
                if ref.get("database") == "Rhea" and rhea_rx.match(rid):
                    rx.append(Reaction(rhea_id=rid, description=rct.get("name")))

        return Protein(
            accession_id=p["primaryAccession"],
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
            embeddings={},
            custom={},
        )


from tqdm import tqdm


async def main() -> None:
    # Read in ids.tsv and write Accession in list of str
    ids = []
    with open("ids.tsv", "r") as f:
        next(f)  # skip first line
        for line in f:
            if line.strip() and not line.startswith("#"):
                parts = line.strip().split("\t")
                if len(parts) > 2:
                    ids.append(parts[2])

    adapter = UniProtAdapter()

    found = 0
    mapped = 0

    async with httpx.AsyncClient() as client:
        # we know the target count
        proteins = []
        with tqdm(total=len(ids), desc="Mapped", unit="prot") as pbar:
            async for rec in adapter.fetch_accessions(
                client,
                ids,
                chunk_size=30,
            ):
                found += 1
                try:
                    proteins.append(adapter.map(rec))
                    mapped += 1
                finally:
                    # advance for each intended ID that produced a mapped record
                    pbar.update(1)

    print(proteins[30:40])


if __name__ == "__main__":
    import asyncio

    from rich import print

    asyncio.run(main())
