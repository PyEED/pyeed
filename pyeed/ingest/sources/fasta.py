"""Refactored FASTA parser with clean async methods.

All methods are async except the high-level wrapper `read_fasta_chunks()`.
"""

from __future__ import annotations

import asyncio
import mmap
import os
from collections.abc import AsyncIterator, Callable, Iterator

from ...utils.progress import ProgressReporter


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
    header_extractor: Callable[[str], str] | None = None,
    taxon_extractor: Callable[[str], str] | None = None,
) -> tuple[str, str, str | None]:
    """Read a single FASTA entry from memory-mapped file.

    Args:
        mm: Memory-mapped file object
        start_offset: Byte position of sequence header
        end_offset: Byte position where sequence ends
        header_extractor: Optional function to extract protein_id from header
                         Signature: (header: str) -> str
        taxon_extractor: Optional function to extract taxon_id from header
                         Signature: (header: str) -> str

    Returns:
        Tuple of (protein_id, sequence, taxon_id)
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

    # Extract protein_id and taxon_id
    protein_id = header_extractor(header_str) if header_extractor else header_str
    taxon_id = taxon_extractor(header_str) if taxon_extractor else None

    return (protein_id, sequence_bytes.decode("ascii", errors="ignore"), taxon_id)


async def _read_batch_async(
    mm: mmap.mmap,
    offsets: list[int],
    end_offsets: list[int],
    start_idx: int,
    batch_size: int,
    header_extractor: Callable[[str], str] | None,
    taxon_extractor: Callable[[str], str] | None,
) -> tuple[list[tuple[str, str, str | None]], int]:
    """Read a single batch of sequences starting at start_idx.

    Args:
        mm: Memory-mapped file object
        offsets: List of header byte offsets
        end_offsets: List of end byte offsets
        start_idx: Starting sequence index
        batch_size: Number of sequences per batch
        header_extractor: Optional function to extract protein_id from header
                         Signature: (header: str) -> str
        taxon_extractor: Optional function to extract taxon_id from header
                         Signature: (header: str) -> str

    Returns:
        Tuple of (batch_list, next_idx) where batch_list contains
        (protein_id, sequence, taxon_id) tuples
    """
    n = len(offsets)
    if start_idx >= n:
        return [], start_idx

    batch: list[tuple[str, str, str | None]] = []
    i = start_idx

    while i < n and len(batch) < batch_size:
        protein_id, sequence, taxon_id = await _read_sequence_async(
            mm,
            offsets[i],
            end_offsets[i],
            header_extractor=header_extractor,
            taxon_extractor=taxon_extractor,
        )
        batch.append((protein_id, sequence, taxon_id))
        i += 1

    return batch, i


async def _iterate_batches_async(
    path: str,
    offsets: list[int],
    end_offsets: list[int],
    batch_size: int,
    header_extractor: Callable[[str], str] | None,
    taxon_extractor: Callable[[str], str] | None,
    progress_report: ProgressReporter | None,
) -> AsyncIterator[list[tuple[str, str, str | None]]]:
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
                batch, next_idx = await _read_batch_async(
                    mm,
                    offsets,
                    end_offsets,
                    idx,
                    batch_size,
                    header_extractor,
                    taxon_extractor,
                )
                if not batch:
                    break

                if progress_report:
                    progress_report(advance=1)
                yield batch
                idx = next_idx
        finally:
            mm.close()


async def read_fasta_chunks_async(
    path: str,
    chunk_size: int,
    header_extractor: Callable[[str], str] | None = None,
    taxon_extractor: Callable[[str], str] | None = None,
    offsets: list[int] | None = None,
    progress_report: ProgressReporter | None = None,
) -> AsyncIterator[list[tuple[str, str, str | None]]]:
    """Read FASTA file in chunks asynchronously.

    This is the main async entry point.

    Args:
        path: Path to FASTA file
        chunk_size: Number of entries per chunk
        header_extractor: Optional function to extract protein_id from header
                         Signature: (header: str) -> str
                         If None, uses entire header as protein_id
        taxon_extractor: Optional function to extract taxon_id from header
                         Signature: (header: str) -> str
                         If None, taxon_id will be None
        offsets: Pre-computed header offsets (optional)
        progress_report: Progress reporter to use

    Yields:
        Lists of tuples (protein_id, sequence, taxon_id)

    Example:
        >>> def extract_uniprot_id(header: str) -> str:
        ...     return header.split("|")[1]
        >>> async for chunk in read_fasta_chunks_async(
        ...     "proteins.fasta",
        ...     chunk_size=1000,
        ...     header_extractor=extract_uniprot_id
        ... ):
        ...     for protein_id, sequence, taxon_id in chunk:
        ...         print(f"{protein_id}: {len(sequence)} aa, taxon: {taxon_id}")
    """
    if chunk_size <= 0:
        raise ValueError(f"chunk_size must be positive, got {chunk_size}")

    if offsets is None:
        offsets = await asyncio.to_thread(build_header_index, path)

    file_size = await asyncio.to_thread(os.path.getsize, path)
    end_offsets = [*offsets[1:], file_size]

    async for batch in _iterate_batches_async(
        path, offsets, end_offsets, chunk_size, header_extractor, taxon_extractor, progress_report
    ):
        yield batch


def read_fasta_chunks(
    path: str,
    chunk_size: int = 1000,
    header_extractor: Callable[[str], str] | None = None,
    taxon_extractor: Callable[[str], str] | None = None,
    progress_report: ProgressReporter | None = None,
) -> Iterator[list[tuple[str, str, str | None]]]:
    """Synchronous wrapper for async FASTA reading.

    This is the high-level entry point for non-async code.

    Args:
        path: Path to FASTA file
        chunk_size: Number of entries per chunk (default: 1000)
        header_extractor: Optional function to extract protein_id from header
                         Signature: (header: str) -> str
        taxon_extractor: Optional function to extract taxon_id from header
                         Signature: (header: str) -> str
        progress_report: Progress reporter to use

    Yields:
        Lists of tuples (protein_id, sequence, taxon_id)

    Example:
        >>> def extract_uniprot_id(header: str) -> str:
        ...     return header.split("|")[1]
        >>> for chunk in read_fasta_chunks(
        ...     "proteins.fasta",
        ...     chunk_size=1000,
        ...     header_extractor=extract_uniprot_id
        ... ):
        ...     for protein_id, sequence, taxon_id in chunk:
        ...         print(f"{protein_id}: {len(sequence)} aa, taxon: {taxon_id}")
    """

    # Run the async generator in a new event loop
    async def _run():
        chunks = []
        async for chunk in read_fasta_chunks_async(
            path, chunk_size, header_extractor, taxon_extractor, None, progress_report
        ):
            chunks.append(chunk)
        return chunks

    chunks = asyncio.run(_run())
    yield from chunks
