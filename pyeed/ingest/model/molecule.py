from __future__ import annotations

from typing import Annotated, Any

from pydantic import Field

from .pyeedbase import BaseNode, LabelProperty


class Molecule(BaseNode):
    id: Annotated[str, LabelProperty(index=True)] = Field(
        ...,
        description="ChEBI ID, e.g. 'CHEBI:10'",
    )
    name: str | None = Field(
        default=None,
        description="Preferred label",
    )

    # basic physicochemical props
    formula: str | None = Field(
        default=None,
        description="Chemical formula",
    )
    charge: int | None = Field(
        default=None,
        description="Charge",
    )
    mass: float | None = Field(
        default=None,
        description="Mass",
    )
    monoisotopic_mass: float | None = Field(
        default=None,
        description="Monoisotopic mass",
    )

    # structure identifiers
    inchi: str | None = Field(
        default=None,
        description="InChI representation",
    )
    inchikey: str | None = Field(
        default=None,
        description="InChI key representation",
    )
    smiles: str | None = Field(
        default=None,
        description="SMILES representation",
    )

    # metadata
    star_rating: int | None = Field(
        default=None,
        description="Star rating",
    )
    deprecated: bool = False

    synonyms: list[str] = Field(
        default_factory=list,
        description="Synonyms",
    )
    cross_refs: list[str] = Field(
        default_factory=list,
        description="Cross references",
    )

    @classmethod
    def from_chebi_node(cls, node: dict[str, Any]) -> Molecule:
        """
        Build a Molecule from a single ChEBI 'node' entry
        from the OBO Graph JSON (chebi.json).
        """
        meta = node.get("meta", {}) or {}
        bpvs = meta.get("basicPropertyValues", []) or []

        # map pred IRI -> val for quick lookup
        pred_map: dict[str, str] = {}
        for bpv in bpvs:
            pred = bpv.get("pred")
            val = bpv.get("val")
            if pred and val:
                pred_map[pred] = val

        # helpers for chemrof IRIs
        chemrof = "https://w3id.org/chemrof/"
        charge_str = pred_map.get(chemrof + "charge")
        mass_str = pred_map.get(chemrof + "mass")
        mono_mass_str = pred_map.get(chemrof + "monoisotopic_mass")
        formula = pred_map.get(chemrof + "generalized_empirical_formula")
        inchi_key = pred_map.get(chemrof + "inchi_key_string")
        inchi = pred_map.get(chemrof + "inchi_string")
        smiles = pred_map.get(chemrof + "smiles_string")

        # star rating from subset like ".../chebi/2_STAR"
        subsets = meta.get("subsets", []) or []
        star_rating: int | None = None
        for s in subsets:
            if "chebi/" in s and s.endswith("_STAR"):
                tail = s.rsplit("/", 1)[-1]  # "2_STAR"
                star_rating = int(tail.split("_", 1)[0])

        # synonyms: take synonym["val"]
        synonyms_raw = meta.get("synonyms", []) or []
        synonyms = [s.get("val") for s in synonyms_raw if s.get("val")]

        # cross refs: list of simple strings
        xrefs_raw = meta.get("xrefs", []) or []
        cross_refs = [x.get("val") for x in xrefs_raw if x.get("val")]

        # deprecated flag
        deprecated = bool(meta.get("deprecated", False))

        # extract CHEBI ID from full IRI "http://purl.obolibrary.org/obo/CHEBI_10"
        full_id = node["id"]
        local = full_id.rsplit("/", 1)[-1]  # "CHEBI_10"
        chebi_id = local.replace("CHEBI_", "CHEBI:")  # "CHEBI:10"

        # label
        name = node.get("lbl")

        return cls(
            id=chebi_id,
            name=name,
            formula=formula,
            charge=int(charge_str) if charge_str is not None else None,
            mass=float(mass_str) if mass_str is not None else None,
            monoisotopic_mass=float(mono_mass_str) if mono_mass_str is not None else None,
            inchi=inchi,
            inchikey=inchi_key,
            smiles=smiles,
            star_rating=star_rating,
            deprecated=deprecated,
            synonyms=synonyms,
            cross_refs=cross_refs,
        )


if __name__ == "__main__":
    # Example: ingest ChEBI from file
    # See pyeed.ingest.chebi_from_file for full ingestion utilities
    import asyncio

    from rich import print as rprint

    from pyeed.db.neo4j import get_async_driver
    from pyeed.ingest.chebi_from_file import ingest_chebi_from_file

    async def main() -> None:
        driver = get_async_driver()
        path = "/home/mha/projects/pyeed/ressources/ChEBI Ontology.json"

        stats = await ingest_chebi_from_file(driver, path, tx_size=5000)
        rprint(stats)

    asyncio.run(main())
