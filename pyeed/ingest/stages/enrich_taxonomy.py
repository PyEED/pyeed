from __future__ import annotations

import asyncio
from typing import Any

import httpx
from neo4j import AsyncDriver
from rich.progress import Progress, TaskID

from pyeed.db.queries import query_existing_nodes_by_ids
from pyeed.ingest.core.pipeline import IngestItem
from pyeed.ingest.core.protocol import SENTINEL, PipelineContext
from pyeed.ingest.core.utils import group_by_relation
from pyeed.ingest.model import Taxon
from pyeed.ingest.model.protein import Protein
from pyeed.ingest.sources.taxonomy import UniProtTaxonomyAdapter


class TaxonomyEnrichmentStage:
    """Batch proteins, fetch missing taxons, upsert, and relate."""

    REL_KEY = "taxon_ids"

    def __init__(
        self,
        driver: AsyncDriver,
        batch_size: int = 50,
        max_concurrent: int = 50,
    ) -> None:
        self.driver = driver
        self.batch_size = batch_size
        self.max_concurrent = max_concurrent
        self.adapter = UniProtTaxonomyAdapter()
        self._seen: set[str] = set()

    async def run(
        self,
        input_queues: dict[str, asyncio.Queue[Any]],
        output_queues: dict[str, asyncio.Queue[Any]],
        context: PipelineContext,
        progress: Progress | None,
        task_id: TaskID | None,
    ) -> None:
        q = next(iter(input_queues.values()))
        buffer: list[IngestItem[Protein]] = []

        while True:
            item = await q.get()
            if item is SENTINEL:
                if buffer:
                    await self._process_batch(buffer, progress, task_id)
                break

            buffer.append(item)
            if len(self._unique_taxa(buffer)) >= self.batch_size:
                await self._process_batch(buffer, progress, task_id)
                buffer.clear()

    def _unique_taxa(self, items: list[IngestItem[Protein]]) -> set[str]:
        return set(group_by_relation(items, self.REL_KEY).keys())

    async def _process_batch(
        self,
        items: list[IngestItem[Protein]],
        progress: Progress | None,
        task_id: TaskID | None,
    ) -> None:
        taxon_to_proteins = group_by_relation(items, self.REL_KEY)
        if not taxon_to_proteins:
            return

        taxon_ids = list(taxon_to_proteins.keys())
        # 1) filter out already existing taxons
        async with self.driver.session() as session:
            existing = await query_existing_nodes_by_ids(session, "Taxon", "id", taxon_ids)
        missing = [tid for tid in taxon_ids if tid not in existing]

        # 2) fetch only missing taxons
        taxon_nodes: list[Taxon] = []
        if missing:
            async with httpx.AsyncClient() as client:
                async for resp in self.adapter.fetch_taxa(client, missing):
                    for node in self.adapter.map(resp):
                        taxon_nodes.append(node)

        # 3) upsert fetched taxons
        if taxon_nodes:
            await Taxon._bulk_upsert(self.driver, taxon_nodes)

        # 4) link proteins -> taxons (both existing and newly fetched)
        by_id = {t.id: t for t in taxon_nodes}
        pairs = [
            (protein, by_id.get(tid) or Taxon(id=tid))
            for tid, proteins in taxon_to_proteins.items()
            for protein in proteins
            if tid in existing or tid in by_id
        ]
        if pairs:
            await Protein.bulk_relate_to_taxa(
                driver=self.driver,
                rel_type="ORIGINATES_FROM",
                pairs=pairs,
            )

        # 5) progress: count newly seen taxa
        new = [tid for tid in taxon_ids if tid not in self._seen]
        if new and progress and task_id is not None:
            self._seen.update(new)
            progress.advance(task_id, len(new))
