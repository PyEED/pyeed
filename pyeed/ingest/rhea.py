#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path

from rdflib import RDF, RDFS, Graph, Literal, URIRef


def _is_reaction_uri(u: URIRef) -> bool:
    s = str(u)
    return "rdf.rhea-db.org/Reaction_" in s or s.endswith("/Reaction") or s.endswith("#Reaction")


def _chebi_literal(v) -> str | None:
    if isinstance(v, Literal):
        s = str(v)
        if s.startswith("CHEBI:"):
            return s
    return None


def _is_uniprot_uri(v) -> bool:
    if isinstance(v, URIRef):
        s = str(v)
        return "purl.uniprot.org/uniprot/" in s or "uniprot.org/uniprot/" in s
    return False


def _first_label(g: Graph, node: URIRef) -> str | None:
    for p in (RDFS.label,):
        for o in g.objects(node, p):
            if isinstance(o, Literal) and str(o).strip():
                return str(o)
    # fallback: any literal that looks like a name/label
    for o in g.objects(node, None):
        if isinstance(o, Literal):
            s = str(o).strip()
            if s and len(s) <= 200:
                return s
    return None


def extract_first_reactions(g: Graph, n: int = 3) -> list[dict]:
    # 1) collect candidate reaction subjects
    candidates: list[URIRef] = []

    # Prefer explicit rdf:type (more accurate if present)
    for s, _, o in g.triples((None, RDF.type, None)):
        if isinstance(s, URIRef) and isinstance(o, URIRef) and _is_reaction_uri(s):
            candidates.append(s)

    # Fallback: anything with URI pattern Reaction_*
    if not candidates:
        for s in set(g.subjects()):
            if isinstance(s, URIRef) and "rdf.rhea-db.org/Reaction_" in str(s):
                candidates.append(s)

    # stable-ish order by numeric id if possible
    def reaction_key(u: URIRef):
        s = str(u)
        try:
            return int(s.rsplit("_", 1)[1])
        except Exception:
            return s

    candidates = sorted(set(candidates), key=reaction_key)[:n]

    out: list[dict] = []
    for rxn in candidates:
        name = _first_label(g, rxn)

        # 2) find CHEBI accessions within 2 hops of the reaction
        chebi_hits: set[str] = set()
        role_map: dict[str, set[str]] = {}

        for p1, o1 in g.predicate_objects(rxn):
            # direct CHEBI literal
            lit = _chebi_literal(o1)
            if lit:
                chebi_hits.add(lit)
                role_map.setdefault(str(p1), set()).add(lit)

            # 1-hop entity -> CHEBI literal (common: participant/compound nodes carry rh:accession)
            if isinstance(o1, URIRef):
                for p2, o2 in g.predicate_objects(o1):
                    lit2 = _chebi_literal(o2)
                    if lit2:
                        chebi_hits.add(lit2)
                        role_map.setdefault(f"{p1} -> {p2}", set()).add(lit2)

        # 3) find UniProt IDs within 2 hops of the reaction
        uniprot: set[str] = set()
        for p1, o1 in g.predicate_objects(rxn):
            if _is_uniprot_uri(o1):
                uniprot.add(str(o1).rstrip("/").rsplit("/", 1)[-1])
            if isinstance(o1, URIRef):
                for _, o2 in g.predicate_objects(o1):
                    if _is_uniprot_uri(o2):
                        uniprot.add(str(o2).rstrip("/").rsplit("/", 1)[-1])

        out.append(
            {
                "reaction_uri": str(rxn),
                "reaction_name": name,
                "chebi_accessions": sorted(chebi_hits),
                "chebi_paths": {k: sorted(v) for k, v in role_map.items()},
                "uniprot_ids": sorted(uniprot),
            }
        )

    return out


def main() -> None:
    # hardcode your RDF/XML file here
    rdf_path = Path("/home/mha/projects/pyeed/ressources/rhea/rhea.rdf")  # <-- change this
    if not rdf_path.exists():
        raise SystemExit(f"File not found: {rdf_path}")

    g = Graph()
    # auto-detect works often, but Rhea dumps are usually RDF/XML:
    g.parse(rdf_path.as_posix(), format="xml")

    data = extract_first_reactions(g, n=3)
    print(json.dumps(data, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
