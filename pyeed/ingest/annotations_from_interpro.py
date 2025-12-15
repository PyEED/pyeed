from __future__ import annotations

import asyncio
import xml.etree.ElementTree as ET
from collections.abc import AsyncIterator, Callable, Iterator
from functools import partial
from pathlib import Path
from typing import Any

from rich import print
from rich.progress import BarColumn, Progress, TextColumn, TimeRemainingColumn

from pyeed.db.neo4j import get_async_driver
from pyeed.db.query_utils import execute_write
from pyeed.ingest.model.annotation import Annotation
from pyeed.ingest.model.annotationcategroy import AnnotationCategory
from pyeed.ingest.model.protein import Protein

XML_PATH = Path("/home/mha/projects/pyeed/ressources/match_complete.xml")

# Default: process full file (set >0 to truncate for dev).
DEV_LINE_LIMIT = 0
# How many parsed rows to echo in the ad-hoc main for smoke-testing.
PRINT_SAMPLE = 5

CATEGORY_MAP = {
    "Domain": AnnotationCategory.DOMAIN,
    "Family": AnnotationCategory.FAMILY,
    "Homologous_superfamily": AnnotationCategory.SUPERFAMILY,
    "Active_site": AnnotationCategory.ACTIVE_SITE,
}

AnnotationTuple = tuple[Annotation, dict[str, Any], str]

# InterPro entry types we keep (must have an <ipr> element).
ALLOWED_IPR_TYPES = frozenset(CATEGORY_MAP.keys())


def _iter_interpro_matches(
    xml_path: Path,
    *,
    line_limit: int = DEV_LINE_LIMIT,
    match_dbs: frozenset[str] | None = None,
    progress_cb: Callable[[int], None] | None = None,
) -> Iterator[AnnotationTuple]:
    """Yield integrated InterPro annotations as (annotation, rel_props, protein_id).

    Streaming parser that feeds the XML pull parser line-by-line to avoid
    retaining already-processed chunks in memory. Only the first ``line_limit``
    lines are consumed for fast dev cycles.
    """
    parser = ET.XMLPullParser(events=("start", "end"))
    current_protein_id: str | None = None
    stopped_early = False

    with xml_path.open("rb") as xml_file:
        for line_no, chunk in enumerate(xml_file, start=1):
            parser.feed(chunk)
            if progress_cb is not None:
                progress_cb(xml_file.tell())

            for event, elem in parser.read_events():
                tag = elem.tag

                if event == "start" and tag == "protein":
                    current_protein_id = elem.attrib.get("id")

                elif event == "end" and tag == "match":
                    db_name = elem.attrib.get("dbname")
                    if match_dbs is not None and db_name not in match_dbs:
                        elem.clear()
                        continue

                    ipr_elem = elem.find("ipr")
                    lcn_elem = elem.find("lcn")

                    # Only keep InterPro-integrated entries.
                    if ipr_elem is None or lcn_elem is None:
                        elem.clear()
                        continue

                    ipr_type = ipr_elem.attrib.get("type")
                    if ipr_type not in ALLOWED_IPR_TYPES:
                        elem.clear()
                        continue

                    match_id = elem.attrib.get("id") or ""
                    ipr_id = ipr_elem.attrib.get("id") or ""
                    ipr_name = ipr_elem.attrib.get("name", "")

                    start_raw = lcn_elem.attrib.get("start")
                    end_raw = lcn_elem.attrib.get("end")
                    score_raw = lcn_elem.attrib.get("score")

                    rel_props: dict[str, Any] = {
                        "start": int(start_raw) if start_raw else None,
                        "end": int(end_raw) if end_raw else None,
                        "score": float(score_raw) if score_raw else None,
                        "match_id": match_id,
                    }

                    annotation = Annotation(
                        id=ipr_id,
                        category=CATEGORY_MAP[ipr_type],
                        name=ipr_name or match_id,
                    )

                    if current_protein_id:
                        yield annotation, rel_props, current_protein_id

                    elem.clear()

                elif event == "end" and tag == "protein":
                    elem.clear()

            if line_limit > 0 and line_no >= line_limit:
                stopped_early = True
                break

    # Only close the parser when we've consumed the entire document; closing early
    # after truncation raises ParseError because the XML is incomplete by design.
    if not stopped_early:
        parser.close()


async def iter_interpro_annotations(
    xml_path: Path,
    *,
    line_limit: int = DEV_LINE_LIMIT,
    match_dbs: frozenset[str] | None = None,
    progress_cb: Callable[[int], None] | None = None,
) -> AsyncIterator[AnnotationTuple]:
    """Async wrapper over the streaming parser to keep the call-site async-friendly."""

    loop = asyncio.get_running_loop()
    iterator = _iter_interpro_matches(
        xml_path,
        line_limit=line_limit,
        match_dbs=match_dbs,
        progress_cb=progress_cb,
    )
    sentinel = object()
    next_item = partial(next, iterator, sentinel)

    while True:
        item = await loop.run_in_executor(None, next_item)
        if item is sentinel:
            break
        yield item  # type: ignore[misc]


async def _process_batch(
    driver,
    batch: list[AnnotationTuple],
    *,
    rel_type: str = "HAS_ANNOTATION",
    rel_tx_size: int = 20_000,
) -> tuple[int, int]:
    """Upsert annotations, relate to existing proteins, and report counts."""
    if not batch:
        return 0, 0

    annotations = [ann for ann, _, _ in batch]
    await Annotation._bulk_upsert(driver, annotations, tx_size=rel_tx_size)

    protein_ids = {protein_id for _, _, protein_id in batch}
    if not protein_ids:
        return len(annotations), 0

    proteins = await Protein.get(driver, ids=list(protein_ids))
    protein_by_id = {p.id: p for p in proteins}

    missing = protein_ids - set(protein_by_id.keys())
    if missing:
        sample = sorted(missing)[:5]
        print(
            f"[warn] missing proteins (skipping relationships): {sample} ... total={len(missing)}"
        )

    rows: list[dict[str, Any]] = []
    for ann, rel_props, protein_id in batch:
        protein = protein_by_id.get(protein_id)
        if protein is None:
            continue
        rows.append(
            {
                "protein_id": protein.id,
                "annotation_id": ann.id,
                "props": rel_props,
            }
        )

    if not rows:
        return len(annotations), 0

    query = f"""
    UNWIND $rows AS row
    MATCH (p:`{Protein.__name__}` {{ `{Protein.get_index_field()}`: row.protein_id }})
    MATCH (a:`{Annotation.__name__}` {{ `{Annotation.get_index_field()}`: row.annotation_id }})
    MERGE (p)-[r:`{rel_type}`]->(a)
    SET r += row.props
    """
    await execute_write(driver, query=query, rows=rows)
    return len(annotations), len(rows)


async def ingest_interpro_annotations(
    *,
    xml_path: Path = XML_PATH,
    match_dbs: frozenset[str] | None = None,
    line_limit: int = DEV_LINE_LIMIT,
    batch_size: int = 20_000,
    rel_tx_size: int = 20_000,
    progress_every_pct: float = 5.0,
) -> None:
    """Stream InterPro XML, upsert annotations, and relate to existing proteins."""
    driver = get_async_driver()
    total_bytes = xml_path.stat().st_size
    bytes_seen = 0

    progress_columns = (
        TextColumn("[bold cyan]annotating proteins"),
        BarColumn(),
        TextColumn("{task.percentage:>5.1f}%"),
        TextColumn("• anns={task.fields[ann_count]}"),
        TimeRemainingColumn(),
    )

    with Progress(*progress_columns) as progress:
        task_id = progress.add_task(
            "annotating proteins",
            total=total_bytes if total_bytes > 0 else 1,
            ann_count=0,
        )

        def _progress(bytes_read: int) -> None:
            nonlocal bytes_seen
            bytes_seen = bytes_read
            progress.update(task_id, completed=bytes_seen)

        batch: list[AnnotationTuple] = []
        total_annotations = 0
        total_rels = 0

        async for item in iter_interpro_annotations(
            xml_path,
            line_limit=line_limit,
            match_dbs=match_dbs,
            progress_cb=_progress,
        ):
            batch.append(item)
            if len(batch) >= batch_size:
                ann_count, rel_count = await _process_batch(
                    driver, batch, rel_type="HAS_ANNOTATION", rel_tx_size=rel_tx_size
                )
                total_annotations += ann_count
                total_rels += rel_count
                batch.clear()
                progress.update(task_id, ann_count=total_annotations)

        if batch:
            ann_count, rel_count = await _process_batch(
                driver, batch, rel_type="HAS_ANNOTATION", rel_tx_size=rel_tx_size
            )
            total_annotations += ann_count
            total_rels += rel_count
            progress.update(task_id, ann_count=total_annotations)

        progress.update(task_id, completed=bytes_seen, ann_count=total_annotations)

    print(
        f"[done] annotations upserted={total_annotations}, relationships created={total_rels}, "
        f"bytes_read={bytes_seen}/{total_bytes}"
    )


async def _demo() -> None:
    """Minimal smoke test that prints the first few parsed rows."""
    count = 0
    async for annotation, rel_props, protein_id in iter_interpro_annotations(
        XML_PATH, line_limit=DEV_LINE_LIMIT, match_dbs=None
    ):
        print(annotation)
        print(rel_props)
        print(protein_id)
        print("-" * 100)
        count += 1
        if count >= PRINT_SAMPLE:
            break


def _parse_match_dbs(raw: str | None) -> frozenset[str] | None:
    if raw is None:
        return None
    parts = [p.strip() for p in raw.split(",") if p.strip()]
    return frozenset(parts) if parts else None


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description="Ingest InterPro annotations into Neo4j.")
    parser.add_argument("--xml-path", type=Path, default=XML_PATH)
    parser.add_argument(
        "--line-limit",
        type=int,
        default=DEV_LINE_LIMIT,
        help="Max lines to read (0 means full file).",
    )
    parser.add_argument("--batch-size", type=int, default=20_000)
    parser.add_argument(
        "--rel-tx-size",
        type=int,
        default=20_000,
        help="Relationship batch size for UNWIND MERGE.",
    )
    parser.add_argument(
        "--match-dbs",
        type=str,
        default=None,
        help="Comma-separated dbnames to include (e.g., PFAM,PANTHER). Empty for all.",
    )
    parser.add_argument(
        "--progress-every-pct",
        type=float,
        default=5.0,
        help="Log progress every N percent of bytes read.",
    )
    parser.add_argument(
        "--demo",
        action="store_true",
        help="Run demo mode (prints first few parsed rows without writing).",
    )

    args = parser.parse_args()
    match_dbs = _parse_match_dbs(args.match_dbs)

    if args.demo:
        asyncio.run(_demo())
    else:
        asyncio.run(
            ingest_interpro_annotations(
                xml_path=args.xml_path,
                line_limit=args.line_limit,
                batch_size=args.batch_size,
                rel_tx_size=args.rel_tx_size,
                match_dbs=match_dbs,
                progress_every_pct=args.progress_every_pct,
            )
        )
