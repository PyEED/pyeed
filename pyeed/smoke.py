# examples/run_smoke.py
from __future__ import annotations

import asyncio
import logging

from pyeed.database import Database
from pyeed.model import (
    MODEL_CLASSES,
    Annotation,
    AnnotationType,
    Embedding,
    Molecule,
    Organism,
    Protein,
    Reaction,
)

logging.basicConfig(level=logging.INFO)


async def main() -> None:
    db = Database()
    db.verify_connection()
    await db.sync_schema(MODEL_CLASSES)

    # ----- 1) save(): one root + subtree -----
    p1 = Protein(
        accession_id="TEST_P00001",
        name="Demo Protein 1",
        sequence="ACDEFGHIKLMNPQRSTVWY",
    )
    # subtree
    p1.organisms.append(Organism(tax_id=9606, name="Homo sapiens"))
    p1.annotations.append(
        Annotation(annotation_type=AnnotationType.ACTIVE_SITE, positions=[10, 11, 12])
    )
    r = Reaction(rhea_id="RHEA:000001", description="Mock reaction")
    r.substrates.append(Molecule(chebi_id="CHEBI:15377", name="H2O"))
    r.products.append(Molecule(chebi_id="CHEBI:15379", name="H+"))
    p1.reactions.append(r)

    emb = Embedding(
        model_name="esm2_t33_650m_ur50s",
        pooling_method="mean",
        n_dims=3,
        vector=[0.1, 0.2, 0.3],
    )
    p1.embeddings.append(emb)

    await db.save(p1)

    print("Saved one Protein with subtree.")

    # sanity check
    res = db.query(
        "MATCH (p:Protein {accession_id:$id}) "
        "OPTIONAL MATCH (p)-[r:CATALYZES]->(:Reaction) "
        "RETURN p.accession_id AS id, count(r) AS rxn_count",
        id="TEST_P00001",
    )
    print("Query result:", res)

    # ----- 2) save_many(): batch roots -----
    p2 = Protein(accession_id="TEST_P00002", sequence="ACDEFGHIKLMNPQR", name="P2")
    p3 = Protein(accession_id="TEST_P00003", sequence="ACDEFGHIKLMNPQR", name="P3")
    await db.save_many([p2, p3])
    print("Saved two more Proteins.")

    # ----- 3) attach(): satellites to existing parent -----
    # add another embedding + an annotation to TEST_P00001
    extra = [
        Embedding(
            model_name="esm2_t33_650m_ur50s",
            pooling_method="cls",
            n_dims=3,
            vector=[0.9, 0.8, 0.7],
        ),
        Annotation(annotation_type=AnnotationType.BINDING_SITE, positions=[2, 3]),
    ]
    await db.attach(Protein, {"TEST_P00001": extra})
    print("Attached satellites to TEST_P00001.")

    # confirm attachment
    res2 = db.query(
        "MATCH (p:Protein {accession_id:$id}) "
        "OPTIONAL MATCH (p)-[:HAS_EMBEDDING]->(e:Embedding) "
        "OPTIONAL MATCH (p)-[:HAS_ANNOTATION]->(a:Annotation) "
        "RETURN p.accession_id AS id, count(DISTINCT e) AS emb_count, count(DISTINCT a) AS ann_count",
        id="TEST_P00001",
    )
    print("After attach:", res2)

    await db.close()


if __name__ == "__main__":
    asyncio.run(main())
