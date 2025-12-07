"""Memory-bounded protein ingestion from JSONL files to Neo4j.

Streams protein data from JSONL (JSON Lines) files and bulk-upserts Protein nodes
with Taxon relationships to Neo4j without loading all data into memory.

Example:
    from pyeed.ingest.proteins_from_jsonl import ingest_proteins_from_jsonl
    from pyeed.db.neo4j import get_async_driver

    async def main():
        driver = get_async_driver()
        stats = await ingest_proteins_from_jsonl(
            driver,
            "trembl.jsonl",
        )
        print(stats)  # {'proteins': 50000, 'taxa': 1200, 'relationships': 48500}
"""

from __future__ import annotations

import asyncio
import json
from collections.abc import AsyncIterator
from dataclasses import dataclass
from pathlib import Path
from typing import Any

from loguru import logger
from neo4j import AsyncDriver
from rich.progress import (
    BarColumn,
    MofNCompleteColumn,
    Progress,
    SpinnerColumn,
    TaskID,
    TextColumn,
    TimeElapsedColumn,
    TimeRemainingColumn,
)

from pyeed.ingest.model import Protein, Taxon

type ProteinTaxonMap = dict[str, str]


@dataclass(frozen=True, slots=True)
class LineIndex:
    """Byte-offset index for efficient JSONL streaming.

    Attributes:
        total_lines: Total number of lines in the file.
        batch_offsets: List of byte offsets at the start of each batch.
        batch_size: Number of lines per batch.
        file_size: Total file size in bytes.
    """

    total_lines: int
    batch_offsets: list[int]
    batch_size: int
    file_size: int


async def build_line_index(
    path: str | Path,
    batch_size: int = 5000,
    progress: Progress | None = None,
    task_id: TaskID | None = None,
) -> LineIndex:
    """Build byte-offset index while counting lines.

    Scans the file once to build an index of byte offsets for efficient
    chunked reading. Shows progress as percentage of file size scanned.

    Args:
        path: Path to JSONL file.
        batch_size: Number of lines per batch (affects index granularity).
        progress: Optional Rich Progress instance for progress tracking.
        task_id: Optional TaskID for progress updates.

    Returns:
        LineIndex with total lines, batch offsets, batch size, and file size.

    Raises:
        FileNotFoundError: If JSONL file doesn't exist.
    """
    path = Path(path)
    if not path.exists():
        raise FileNotFoundError(f"JSONL file not found: {path}")

    file_size = path.stat().st_size

    def _scan() -> LineIndex:
        offsets = [0]  # Start of file
        line_count = 0
        bytes_read = 0

        with open(path, "rb") as f:
            while True:
                line = f.readline()
                if not line:
                    break

                line_count += 1
                bytes_read += len(line)

                # Record offset every batch_size lines
                if line_count % batch_size == 0:
                    offsets.append(f.tell())

                # Update progress every 1000 lines
                if progress and task_id and line_count % 1000 == 0:
                    pct = (bytes_read / file_size) * 100 if file_size > 0 else 0
                    mb_read = bytes_read / (1024 * 1024)
                    mb_total = file_size / (1024 * 1024)
                    desc = (
                        f"[cyan]Scanning JSONL... {pct:.1f}% ({mb_read:.0f} MB / {mb_total:.0f} MB)"
                    )
                    progress.update(
                        task_id,
                        completed=bytes_read,
                        description=desc,
                    )

        # Final progress update
        if progress and task_id:
            mb_total = file_size / (1024 * 1024)
            progress.update(
                task_id,
                completed=file_size,
                description=f"[cyan]Scanning complete ({mb_total:.0f} MB)",
            )

        return LineIndex(
            total_lines=line_count,
            batch_offsets=offsets,
            batch_size=batch_size,
            file_size=file_size,
        )

    return await asyncio.to_thread(_scan)


async def stream_lines_from_index(
    path: str | Path,
    index: LineIndex,
) -> AsyncIterator[str]:
    """Stream lines using byte-offset index.

    Reads file in batches without loading entire file into memory.
    Uses pre-built byte-offset index to jump directly to each batch.

    Args:
        path: Path to JSONL file.
        index: Pre-built LineIndex with byte offsets.

    Yields:
        Individual lines (strings) from the file.
    """
    path = Path(path)

    def _read_batch(offset: int, next_offset: int | None) -> list[str]:
        """Read one batch starting at offset."""
        with open(path, "rb") as f:
            f.seek(offset)
            lines = []

            while len(lines) < index.batch_size:
                raw_line = f.readline()
                if not raw_line:
                    break
                # Don't exceed next batch boundary
                if next_offset and f.tell() > next_offset:
                    break
                lines.append(raw_line.decode("utf-8"))

            return lines

    # Process each batch
    for i, offset in enumerate(index.batch_offsets):
        next_offset = index.batch_offsets[i + 1] if i + 1 < len(index.batch_offsets) else None

        batch = await asyncio.to_thread(_read_batch, offset, next_offset)

        for text_line in batch:
            yield text_line


async def stream_proteins_from_jsonl(
    path: str | Path,
    index: LineIndex,
) -> AsyncIterator[tuple[Protein, str | None]]:
    """Stream validated (protein, taxon_id) tuples from JSONL using index.

    Uses pre-built byte-offset index to read file in chunks without loading
    entire file into memory. Validates each entry and yields individual
    (Protein, taxon_id) tuples.

    JSONL format expected:
        {"id": "A0A7C4XJF3", "name": "Protein name", "organism": "1872626", "seq": "MKAY..."}

    Args:
        path: Path to JSONL file.
        index: Pre-built LineIndex with byte offsets.

    Yields:
        Tuples of (Protein object, taxon_id or None).

    Raises:
        FileNotFoundError: If JSONL file doesn't exist.
    """
    path = Path(path)
    if not path.exists():
        raise FileNotFoundError(f"JSONL file not found: {path}")

    skipped_count = 0
    line_number = 0

    # Stream lines using index (memory-bounded)
    async for raw_line in stream_lines_from_index(path, index):
        line_number += 1
        line = raw_line.strip()

        if not line:
            continue

        try:
            data = json.loads(line)
        except json.JSONDecodeError as e:
            logger.warning(
                f"Invalid JSON at line {line_number}, skipping: {e}",
                extra={"line": line_number, "reason": "invalid_json"},
            )
            skipped_count += 1
            continue

        # Extract fields
        protein_id = data.get("id")
        sequence = data.get("seq")
        name = data.get("name")
        organism = data.get("organism")

        # Validate required fields
        if not protein_id:
            logger.warning(
                f"Missing 'id' field at line {line_number}, skipping",
                extra={"line": line_number, "reason": "missing_id"},
            )
            skipped_count += 1
            continue

        if not sequence or not sequence.strip():
            logger.warning(
                f"Empty sequence for protein {protein_id} at line {line_number}, skipping",
                extra={"line": line_number, "protein_id": protein_id, "reason": "empty_sequence"},
            )
            skipped_count += 1
            continue

        # Create Protein object
        try:
            protein = Protein(
                id=protein_id,
                sequence=sequence,
                seq_length=len(sequence),
                name=name,
            )
            yield protein, organism

        except Exception as e:
            logger.warning(
                f"Failed to create Protein for {protein_id} at line {line_number}: {e}",
                extra={"line": line_number, "protein_id": protein_id, "error": str(e)},
            )
            skipped_count += 1
            continue

    if skipped_count > 0:
        logger.warning(
            f"Skipped {skipped_count} invalid entries during JSONL reading",
            extra={"skipped": skipped_count},
        )


async def batch_proteins(
    protein_stream: AsyncIterator[tuple[Protein, str | None]],
    batch_size: int = 5000,
) -> AsyncIterator[tuple[list[Protein], ProteinTaxonMap]]:
    """Accumulate proteins into fixed-size batches.

    Consumes a stream of (Protein, taxon_id) tuples and yields batches
    of proteins with their taxon mapping.

    Args:
        protein_stream: Async iterator of (Protein, taxon_id) tuples.
        batch_size: Maximum number of proteins per batch.

    Yields:
        Tuples of (proteins_list, protein_taxon_map) for each batch.
    """
    buffer_proteins: list[Protein] = []
    buffer_taxon_map: ProteinTaxonMap = {}

    async for protein, taxon_id in protein_stream:
        buffer_proteins.append(protein)
        if taxon_id:
            buffer_taxon_map[protein.id] = taxon_id

        if len(buffer_proteins) >= batch_size:
            yield buffer_proteins, buffer_taxon_map
            buffer_proteins = []
            buffer_taxon_map = {}

    # Yield remainder
    if buffer_proteins:
        yield buffer_proteins, buffer_taxon_map


async def write_protein_batch(
    driver: AsyncDriver,
    proteins: list[Protein],
    protein_taxon_map: ProteinTaxonMap,
    tx_size: int = 5000,
) -> dict[str, int]:
    """Write a single batch: proteins → taxa → relationships.

    Checks for existing proteins and skips them before upserting.

    Args:
        driver: Neo4j async driver.
        proteins: List of Protein objects to upsert.
        protein_taxon_map: Dict mapping protein_id -> taxon_id for this batch.
        tx_size: Transaction batch size for Neo4j operations.

    Returns:
        Stats dict with keys: proteins, taxa, relationships, skipped.
    """
    stats = {"proteins": 0, "taxa": 0, "relationships": 0, "skipped": 0}

    if not proteins:
        return stats

    # Step 1: Check which proteins already exist
    protein_ids = [p.id for p in proteins]
    existing_proteins = await Protein.get(driver, ids=protein_ids)
    existing_ids = {p.id for p in existing_proteins}

    # Filter out existing proteins
    new_proteins = [p for p in proteins if p.id not in existing_ids]
    skipped_count = len(proteins) - len(new_proteins)
    stats["skipped"] = skipped_count

    if skipped_count > 0:
        logger.debug(
            f"Skipped {skipped_count} proteins that already exist in database",
            extra={"skipped": skipped_count, "total": len(proteins)},
        )

    # Step 2: Upsert only new proteins
    if new_proteins:
        await Protein._bulk_upsert(driver, new_proteins, tx_size=tx_size)
        stats["proteins"] = len(new_proteins)

        # Step 3: Ensure taxa exist and create relationships (only for new proteins)
        new_protein_ids = {p.id for p in new_proteins}
        new_taxon_map = {
            protein_id: taxon_id
            for protein_id, taxon_id in protein_taxon_map.items()
            if protein_id in new_protein_ids
        }

        if new_taxon_map:
            unique_taxon_ids = set(new_taxon_map.values())

            # Create minimal Taxon objects
            taxa = [Taxon(id=taxon_id) for taxon_id in unique_taxon_ids]
            await Taxon._bulk_upsert(driver, taxa, tx_size=tx_size)
            stats["taxa"] = len(unique_taxon_ids)

            # Create relationships
            pairs = [
                (Protein.model_construct(id=protein_id), Taxon.model_construct(id=taxon_id))
                for protein_id, taxon_id in new_taxon_map.items()
            ]
            await Protein.bulk_relate_to_taxa(driver, pairs)
            stats["relationships"] = len(new_taxon_map)

    return stats


async def ingest_proteins_from_jsonl(
    driver: AsyncDriver,
    path: str | Path,
    tx_size: int = 5000,
    batch_size: int = 5000,
) -> dict[str, Any]:
    """Ingest proteins from JSONL with memory-bounded streaming.

    Orchestrates the complete ingestion pipeline using byte-offset indexing:
    1. Build byte-offset index (shows % of file scanned)
    2. Stream proteins using index (memory-bounded)
    3. Batch into fixed-size chunks
    4. Write each batch to Neo4j with taxa and relationships
    5. Display progress with accurate totals

    Memory characteristics:
        - Before: O(N) - entire file in RAM
        - After: O(B) - only batch_size lines (~5000) in RAM at once
        - Index overhead: O(N/B) - one int per batch (~8 bytes per 5000 lines)

    Args:
        driver: Neo4j async driver.
        path: Path to JSONL file.
        tx_size: Transaction batch size for Neo4j operations.
        batch_size: Number of proteins per batch (affects memory usage).

    Returns:
        Stats dict with keys:
            - proteins: Total proteins upserted
            - skipped: Total proteins skipped (already exist in database)
            - taxa: Total unique taxa created
            - relationships: Total relationships created

    Example:
        >>> from pyeed.ingest.proteins_from_jsonl import ingest_proteins_from_jsonl
        >>> from pyeed.db.neo4j import get_async_driver
        >>>
        >>> async def main():
        ...     driver = get_async_driver()
        ...     stats = await ingest_proteins_from_jsonl(
        ...         driver,
        ...         "trembl.jsonl",
        ...     )
        ...     print(stats)
        >>>
        >>> asyncio.run(main())
    """
    path = Path(path)
    if not path.exists():
        raise FileNotFoundError(f"JSONL file not found: {path}")

    file_size = path.stat().st_size

    # Setup progress bar
    with Progress(
        SpinnerColumn(),
        TextColumn("[bold blue]{task.description}"),
        BarColumn(),
        MofNCompleteColumn(),
        TimeElapsedColumn(),
        TimeRemainingColumn(),
    ) as progress:
        # Step 1: Build index (shows % of file scanned)
        task_scan = progress.add_task(
            "[cyan]Scanning JSONL file...",
            total=file_size,
        )
        index = await build_line_index(path, batch_size, progress, task_scan)
        progress.remove_task(task_scan)

        logger.info(
            f"Built index: {index.total_lines} lines, "
            f"{len(index.batch_offsets)} batches, "
            f"{index.file_size / (1024 * 1024):.1f} MB"
        )

        if index.total_lines == 0:
            logger.warning("No lines found in JSONL file")
            return {"proteins": 0, "taxa": 0, "relationships": 0, "skipped": 0}

        # Step 2: Main ingestion task
        task_ingest = progress.add_task(
            "[green]Ingesting proteins",
            total=index.total_lines,
        )

        # Aggregate stats
        total_stats: dict[str, int] = {
            "proteins": 0,
            "taxa": 0,
            "relationships": 0,
            "skipped": 0,
        }
        all_taxon_ids: set[str] = set()

        # Step 3: Stream using index, batch, and write
        protein_stream = stream_proteins_from_jsonl(path, index)
        batches = batch_proteins(protein_stream, batch_size)

        async for proteins, taxon_map in batches:
            # Write batch to Neo4j (includes existence check)
            batch_stats = await write_protein_batch(driver, proteins, taxon_map, tx_size)

            # Update aggregate stats
            total_stats["proteins"] += batch_stats["proteins"]
            total_stats["relationships"] += batch_stats["relationships"]
            total_stats["skipped"] += batch_stats["skipped"]
            all_taxon_ids.update(taxon_map.values())

            # Update progress (advance by total proteins processed, including skipped)
            progress.update(task_ingest, advance=len(proteins))

            logger.debug(
                f"Processed batch: {batch_stats['proteins']} new proteins, "
                f"{batch_stats['skipped']} skipped (already exist), "
                f"{batch_stats['taxa']} taxa, {batch_stats['relationships']} relationships",
                extra=batch_stats,
            )

        # Set final taxa count (unique across all batches)
        total_stats["taxa"] = len(all_taxon_ids)

    logger.info(
        f"Protein ingestion complete: {total_stats['proteins']} proteins added, "
        f"{total_stats['skipped']} skipped (already exist), "
        f"{total_stats['taxa']} taxa, {total_stats['relationships']} relationships",
        extra=total_stats,
    )

    return total_stats


if __name__ == "__main__":
    import asyncio
    import sys

    from rich import print as rprint

    from pyeed.db.neo4j import get_async_driver

    MIN_ARGS = 2

    async def main() -> None:
        if len(sys.argv) < MIN_ARGS:
            rprint("[red]Usage: python -m pyeed.ingest.proteins_from_jsonl <jsonl_path>[/red]")
            sys.exit(1)

        jsonl_path = sys.argv[1]
        driver = get_async_driver()

        rprint(f"[cyan]Ingesting proteins from {jsonl_path}...[/cyan]")

        stats = await ingest_proteins_from_jsonl(
            driver,
            jsonl_path,
            tx_size=5000,
            batch_size=5000,
        )

        rprint("[green]Done![/green]")
        rprint(stats)

    asyncio.run(main())
