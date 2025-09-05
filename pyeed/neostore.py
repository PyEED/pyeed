from __future__ import annotations

from typing import Any, Dict, List, Tuple

from neo4j import AsyncGraphDatabase

from pyeed.model import BaseNode


class Neo4jStore:
    def __init__(self, uri: str, user: str, password: str, name: str = "neo4j") -> None:
        self.driver = AsyncGraphDatabase.driver(uri, auth=(user, password))

    async def close(self) -> None:
        await self.driver.close()

    async def ensure_constraints(self) -> None:
        """
        Idempotent unique constraints for upsert-safety.
        """
        uniq = [
            ("Protein", "accession_id"),
            ("Organism", "tax_id"),
            ("GO", "go_id"),
            ("Reaction", "rhea_id"),
            ("Annotation", "ann_id"),
            ("Embedding", "id"),
        ]
        async with self.driver.session() as s:
            for label, key in uniq:
                q = (
                    f"CREATE CONSTRAINT `uniq_{label}_{key}` IF NOT EXISTS "
                    f"FOR (n:`{label}`) REQUIRE n.`{key}` IS UNIQUE"
                )
                await s.run(q)

    async def bulk_upsert(
        self,
        nodes: List[Dict[str, Any]],
        edges: List[Dict[str, Any]],
        tx_size: int = 5000,
    ) -> None:
        async with self.driver.session() as session:
            # Batch MERGE for nodes
            buckets: Dict[str, Tuple[str, List[Dict[str, Any]]]] = {}
            for n in nodes:
                lbl = n["label"]
                key_name, key_val = n["key"]
                buckets.setdefault(lbl, (key_name, []))[1].append(
                    {"key": key_val, "props": n["props"]}
                )

            for lbl, (kname, rows) in buckets.items():
                for i in range(0, len(rows), tx_size):
                    chunk = rows[i : i + tx_size]
                    q = f"""
                    UNWIND $rows AS r
                    MERGE (n:`{lbl}` {{ `{kname}`: r.key }})
                    SET n += r.props
                    """
                    await session.run(q, rows=chunk)

            # Batch MERGE for edges
            groups: Dict[Tuple, List[Dict[str, Any]]] = {}
            for e in edges:
                key = (e["type"], *e["src"][:2], *e["dst"][:2])
                groups.setdefault(key, []).append(
                    {"sv": e["src"][2], "dv": e["dst"][2]}
                )

            for (etype, sl, sk, dl, dk), rows in groups.items():
                for i in range(0, len(rows), tx_size):
                    chunk = rows[i : i + tx_size]
                    q = f"""
                    UNWIND $rows AS r
                    MATCH (s:`{sl}` {{ `{sk}`: r.sv }})
                    MATCH (d:`{dl}` {{ `{dk}`: r.dv }})
                    MERGE (s)-[:`{etype}`]->(d)
                    """
                    await session.run(q, rows=chunk)

    async def upsert_node(self, node: BaseNode) -> None:
        nodes, edges = node.graphify()
        await self.bulk_upsert(nodes, edges)


if __name__ == "__main__":
    from .model import Annotation, AnnotationType, Embedding, Protein

    prot = Protein(
        accession_id="P01234",
        sequence="MALWMRLLPLLALLALWGPDPAAA",
        name="Test Protein",
        seq_length=24,
        mol_weight=1000,
        ec_numbers=["1.2.3.4"],
        custom={
            "mentioned_in": [
                "doi:10.1016/j.xinn.2025.100344",
                "doi:10.1016/j.xinn.2025.100345",
            ],
            "already_characterized": False,
            "dictd": 123,
        },
    )

    prot.annotations.append(
        Annotation(
            annotation_type=AnnotationType.ACTIVE_SITE,
            positions=[1, 2, 3],
            custom={"my_custom_evidence": "literature", "validated": True},
        ),
    )

    # add second annotation
    prot.annotations.append(
        Annotation(
            annotation_type=AnnotationType.BINDING_SITE,
            positions=[4, 5, 6],
        ),
    )

    # add embedding
    prot.embeddings.append(
        Embedding(
            model_name="esm2_t33_650M_UR50S",
            pooling_method="mean",
            vector=[0.1, 0.2, 0.3],
            n_dims=3,
        ),
    )

    # add another embedding
    prot.embeddings.append(
        Embedding(
            model_name="esm2_t33_650M_UR50S",
            pooling_method="mean",
            vector=[0.4, 0.5, 0.6],
            n_dims=3,
        ),
    )
    prot.embeddings.append(
        Embedding(
            model_name="esm2-t33-650M-UR50S",
            pooling_method="range",
            vector=[0.4, 0.5, 0.6, 0.7, 0.8],
            n_dims=5,
        ),
    )

    import asyncio

    db = Neo4jStore(
        uri="bolt://127.0.0.1:7687",
        user="neo4j",
        password="12345678",
    )

    async def main():
        await db.upsert_node(prot)
        await db.close()

    asyncio.run(main())
