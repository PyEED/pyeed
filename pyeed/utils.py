"""Refactored FASTA parser with clean async methods.

All methods are async except the high-level wrapper `read_fasta_chunks()`.
"""

from __future__ import annotations

import asyncio
import mmap
import os
import re
from collections.abc import AsyncIterator, Iterator

from .ingest.progress import ProgressReporter


def build_header_index(path: str) -> list[int]:
    """Build an index of byte offsets for all FASTA headers.

    Single-pass scan using memory-mapped file access.

    Args:
        path: Path to FASTA file

    Returns:
        List of byte offsets where each header begins
    """
    offsets: list[int] = []

    with open(path, "rb") as f, mmap.mmap(f.fileno(), 0, access=mmap.ACCESS_READ) as mm:
        pos = 0
        while pos < len(mm):
            newline_pos = mm.find(b"\n", pos)
            if newline_pos == -1:
                break

            if mm[pos] == ord(b">"):
                offsets.append(pos)

            pos = newline_pos + 1

    return offsets


async def _read_sequence_async(
    mm: mmap.mmap,
    start_offset: int,
    end_offset: int,
    header_pattern: re.Pattern[str] | None = None,
) -> tuple[str, str]:
    """Read a single FASTA entry from memory-mapped file.

    Args:
        mm: Memory-mapped file object
        start_offset: Byte position of sequence header
        end_offset: Byte position where sequence ends
        header_pattern: Optional regex to extract protein_id from header

    Returns:
        Tuple of (protein_id, sequence)
    """
    # Validate header
    if mm[start_offset] != ord(b">"):
        header_preview = mm[start_offset : start_offset + 50]
        raise ValueError(
            f"Expected FASTA header at byte offset {start_offset}, found: {header_preview}"
        )

    # Find end of header line
    header_end = mm.find(b"\n", start_offset)
    if header_end == -1 or header_end >= end_offset:
        header_end = end_offset

    # Extract header (without '>' prefix)
    header_bytes = mm[start_offset + 1 : header_end]

    # Extract sequence (skip header line)
    seq_start = header_end + 1
    if seq_start >= end_offset:
        sequence_bytes = b""
    else:
        sequence_bytes = mm[seq_start:end_offset]
        # Fast C-level newline removal
        sequence_bytes = sequence_bytes.replace(b"\n", b"").replace(b"\r", b"")

    # Decode header
    header_str = header_bytes.decode("ascii", errors="ignore")

    # Extract protein_id from header using pattern if provided
    if header_pattern:
        match = header_pattern.search(header_str)
        if match:
            protein_id = match.group(1) if match.lastindex else match.group(0)
        else:
            protein_id = header_str
    else:
        protein_id = header_str

    return (protein_id, sequence_bytes.decode("ascii", errors="ignore"))


async def _read_batch_async(
    mm: mmap.mmap,
    offsets: list[int],
    end_offsets: list[int],
    start_idx: int,
    batch_size: int,
    header_pattern: re.Pattern[str] | None,
) -> tuple[dict[str, str], int]:
    """Read a single batch of sequences starting at start_idx.

    Args:
        mm: Memory-mapped file object
        offsets: List of header byte offsets
        end_offsets: List of end byte offsets
        start_idx: Starting sequence index
        batch_size: Number of sequences per batch
        header_pattern: Optional regex to extract protein_id

    Returns:
        Tuple of (batch_dict, next_idx)
    """
    n = len(offsets)
    if start_idx >= n:
        return {}, start_idx

    batch: dict[str, str] = {}
    i = start_idx

    while i < n and len(batch) < batch_size:
        protein_id, sequence = await _read_sequence_async(
            mm, offsets[i], end_offsets[i], header_pattern=header_pattern
        )
        batch[protein_id] = sequence
        i += 1

    return batch, i


async def _iterate_batches_async(
    path: str,
    offsets: list[int],
    end_offsets: list[int],
    batch_size: int,
    header_pattern: re.Pattern[str] | None,
    progress_report: ProgressReporter | None,
) -> AsyncIterator[dict[str, str]]:
    """Internal implementation: iterate over FASTA in batches.

    This keeps the file open for the entire iteration, avoiding
    repeated open/close overhead.
    """
    idx = 0
    n = len(offsets)

    with open(path, "rb") as f:
        mm = mmap.mmap(f.fileno(), 0, access=mmap.ACCESS_READ)
        try:
            while idx < n:
                # read a batch synchronously inside a thread to avoid blocking loop
                def read_batch_sync(start: int) -> tuple[dict[str, str], int]:
                    batch: dict[str, str] = {}
                    i = start
                    while i < n and len(batch) < batch_size:
                        pid, seq = _read_sequence_from_offset_mmap(
                            mm, offsets[i], end_offsets[i], header_pattern
                        )
                        batch[pid] = seq
                        i += 1
                    return batch, i

                batch, next_idx = await asyncio.to_thread(read_batch_sync, idx)
                if not batch:
                    break

                if progress_report:
                    progress_report(advance=1)
                yield batch
                idx = next_idx
        finally:
            mm.close()


def _read_sequence_from_offset_mmap(
    mm: mmap.mmap,
    start_offset: int,
    end_offset: int,
    header_pattern: re.Pattern[str] | None = None,
) -> tuple[str, str]:
    """Synchronous version for use in thread pool."""
    if mm[start_offset] != ord(b">"):
        header_preview = mm[start_offset : start_offset + 50]
        raise ValueError(
            f"Expected FASTA header at byte offset {start_offset}, found: {header_preview}"
        )

    header_end = mm.find(b"\n", start_offset)
    if header_end == -1 or header_end >= end_offset:
        header_end = end_offset

    header_bytes = mm[start_offset + 1 : header_end]

    seq_start = header_end + 1
    if seq_start >= end_offset:
        sequence_bytes = b""
    else:
        sequence_bytes = mm[seq_start:end_offset]
        sequence_bytes = sequence_bytes.replace(b"\n", b"").replace(b"\r", b"")

    header_str = header_bytes.decode("ascii", errors="ignore")

    if header_pattern:
        match = header_pattern.search(header_str)
        if match:
            protein_id = match.group(1) if match.lastindex else match.group(0)
        else:
            protein_id = header_str
    else:
        protein_id = header_str

    return (protein_id, sequence_bytes.decode("ascii", errors="ignore"))


async def read_fasta_chunks_async(
    path: str,
    chunk_size: int,
    header_pattern: str | re.Pattern[str] | None = None,
    offsets: list[int] | None = None,
    progress_report: ProgressReporter | None = None,
) -> AsyncIterator[dict[str, str]]:
    """Read FASTA file in chunks asynchronously.

    This is the main async entry point.

    Args:
        path: Path to FASTA file
        chunk_size: Number of entries per chunk
        header_pattern: Optional regex to extract protein_id from header
        progress_report: Progress reporter to use

    Yields:
        Dicts mapping protein_id to sequence

    Example:
        >>> async for chunk in read_fasta_chunks_async("proteins.fasta"):
        ...     for protein_id, sequence in chunk.items():
        ...         print(f"{protein_id}: {len(sequence)} aa")
    """
    if chunk_size <= 0:
        raise ValueError(f"chunk_size must be positive, got {chunk_size}")

    compiled_pattern: re.Pattern[str] | None = (
        re.compile(header_pattern) if isinstance(header_pattern, str) else header_pattern
    )

    if offsets is None:
        offsets = await asyncio.to_thread(build_header_index, path)

    file_size = await asyncio.to_thread(os.path.getsize, path)
    end_offsets = [*offsets[1:], file_size]

    async for chunk in _iterate_batches_async(
        path, offsets, end_offsets, chunk_size, compiled_pattern, progress_report
    ):
        yield chunk


def read_fasta_chunks(
    path: str,
    chunk_size: int = 1000,
    header_pattern: str | re.Pattern[str] | None = None,
    progress_report: ProgressReporter | None = None,
) -> Iterator[dict[str, str]]:
    """Synchronous wrapper for async FASTA reading.

    This is the high-level entry point for non-async code.

    Args:
        path: Path to FASTA file
        chunk_size: Number of entries per chunk (default: 1000)
        header_pattern: Optional regex to extract protein_id from header
        show_progress: Display progress bar (default: True)
        progress: Existing Progress instance to share

    Yields:
        Dicts mapping protein_id to sequence

    Example:
        >>> for chunk in read_fasta_chunks("proteins.fasta", chunk_size=1000):
        ...     for protein_id, sequence in chunk.items():
        ...         print(f"{protein_id}: {len(sequence)} aa")
    """

    # Run the async generator in a new event loop
    async def _run():
        chunks = []
        async for chunk in read_fasta_chunks_async(
            path, chunk_size, header_pattern, progress_report
        ):
            chunks.append(chunk)
        return chunks

    chunks = asyncio.run(_run())
    yield from chunks


if __name__ == "__main__":
    import asyncio
    import re
    from collections.abc import Callable

    from rich.progress import BarColumn, Progress, SpinnerColumn, TextColumn, TimeElapsedColumn

    path = "/home/mha/projects/proteingraph/downloads/uniprot_sprot.fasta"

    # -------- reporter type (simple callable) --------
    Reporter = Callable[..., None]

    def null_reporter(**_: object) -> None:  # no-op
        pass

    # -------- pipeline (clean, no monkey-patch) --------
    class AsyncPipelineDemo:
        SENTINEL = object()

        def __init__(
            self,
            path: str,
            chunk_size: int = 1000,
            header_pattern: str | re.Pattern[str] | None = None,
        ):
            self.path = path
            self.chunk_size = chunk_size
            self.header_pattern = header_pattern

            self.read_queue: asyncio.Queue[dict[str, str] | object] = asyncio.Queue(maxsize=2)
            self.write_queue: asyncio.Queue[list[str] | object] = asyncio.Queue(maxsize=2)

        async def _reader_worker(self, report_read: ProgressReporter) -> None:
            async for chunk in read_fasta_chunks_async(  # uses your parser
                self.path, self.chunk_size, self.header_pattern, report_read
            ):
                await self.read_queue.put(chunk)
            # signal end of stream to processor
            await self.read_queue.put(self.SENTINEL)

        async def _processor_worker(self, report_proc: ProgressReporter) -> None:
            while True:
                batch = await self.read_queue.get()
                if batch is self.SENTINEL:
                    # propagate termination to writer and stop
                    await self.write_queue.put(self.SENTINEL)
                    break
                report_proc(advance=1)
                await asyncio.sleep(1.0)  # simulate heavy processing
                await self.write_queue.put(list(batch.keys()))

        async def _writer_worker(
            self, report_write: ProgressReporter, flush_batches: int = 10
        ) -> None:
            buf: list[list[str]] = []
            while True:
                item = await self.write_queue.get()
                if item is self.SENTINEL:
                    if buf:
                        print(f"flushing {len(buf)} batches")
                        await asyncio.sleep(0.3)  # simulate fast upload
                        # Update progress incrementally for remaining items
                        report_write(advance=len(buf))
                        buf.clear()
                    break

                buf.append(item)
                if len(buf) >= flush_batches:
                    await asyncio.sleep(0.3)  # simulate fast upload
                    report_write(advance=len(buf))
                    buf.clear()

        async def run(self) -> None:
            from rich.progress import MofNCompleteColumn

            progress = Progress(
                SpinnerColumn(),
                TextColumn("[bold]{task.description}"),
                BarColumn(),
                MofNCompleteColumn(),
                TimeElapsedColumn(),
                refresh_per_second=3,
            )

            read_task_id = progress.add_task("read", total=None)
            read_reporter = ProgressReporter(progress, read_task_id)

            proc_task_id = progress.add_task("proc", total=None)
            proc_reporter = ProgressReporter(progress, proc_task_id)

            write_task_id = progress.add_task("write", total=None)
            write_reporter = ProgressReporter(progress, write_task_id)

            with progress:
                await asyncio.gather(
                    self._reader_worker(read_reporter),
                    self._processor_worker(proc_reporter),
                    self._writer_worker(write_reporter),
                )

    # -------- Rich Progress wiring via reporters --------

    demo = AsyncPipelineDemo(
        path=path,
        chunk_size=1000,
    )
    asyncio.run(demo.run())
