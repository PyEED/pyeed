"""Efficient FASTA file parser with seek-based random access.

This module provides tools for parsing large FASTA files without loading
the entire file into memory. It uses byte offset indexing for fast random
access to individual sequences.

Optimizations:
- Memory-mapped file access (mmap) for faster reads
- C-level translate() for removing newlines (faster than replace())
"""

from __future__ import annotations

import mmap
import os
import re
from collections.abc import Iterator

from rich.progress import Progress, TaskID

from .ingest.progress import create_progress


def build_header_index(path: str) -> list[int]:
    """Build an index of byte offsets for all FASTA headers in a file.

    Performs a single-pass scan using memory-mapped file access for speed.
    Locates all lines starting with '>' (FASTA header lines) and records
    their byte positions.

    Args:
        path: Path to the FASTA file.

    Returns:
        List of byte offsets where each header begins, in file order.

    Raises:
        FileNotFoundError: If the file doesn't exist.
        IOError: If the file cannot be read.
    """
    offsets: list[int] = []

    with open(path, "rb") as f, mmap.mmap(f.fileno(), 0, access=mmap.ACCESS_READ) as mm:
        pos = 0
        while pos < len(mm):
            # Find next newline
            newline_pos = mm.find(b"\n", pos)
            if newline_pos == -1:
                break

            # Check if line starts with '>'
            if mm[pos] == ord(b">"):
                offsets.append(pos)

            pos = newline_pos + 1

    return offsets


def _read_sequence_from_offset_mmap(
    mm: mmap.mmap,
    start_offset: int,
    end_offset: int,
    header_pattern: re.Pattern[str] | None = None,
) -> dict[str, str]:
    """Read a single FASTA entry (protein_id + sequence) from a memory-mapped file.

    Optimized version using mmap and replace() for removing newlines.

    Args:
        mm: Memory-mapped file object.
        start_offset: Byte position of the sequence header.
        end_offset: Byte position where this sequence ends (start of next
            sequence or EOF).
        header_pattern: Optional regex pattern to extract protein_id from header.
            If provided, the first capture group is used as protein_id.
            If None, the entire header (without '>') is used as protein_id.

    Returns:
        Dict with 'protein_id' and 'sequence' keys. Both are strings.

    Raises:
        ValueError: If start_offset doesn't point to a valid header line.
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

    # Extract header (without '>' prefix, without trailing newline)
    header_bytes = mm[start_offset + 1 : header_end]  # Skip '>' prefix

    # Extract sequence (skip header line)
    seq_start = header_end + 1
    if seq_start >= end_offset:
        # Empty sequence
        sequence_bytes = b""
    else:
        # Read sequence bytes
        sequence_bytes = mm[seq_start:end_offset]
        # Remove newlines using replace (fast, C-level implementation)
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

    return {
        "protein_id": protein_id,
        "sequence": sequence_bytes.decode("ascii", errors="ignore"),
    }


def read_fasta_chunks(
    path: str,
    chunk_size: int = 1000,
    header_pattern: str | re.Pattern[str] | None = None,
    show_progress: bool = True,
    progress: Progress | None = None,
) -> Iterator[list[dict[str, str]]]:
    """Read FASTA file in chunks of entries.

    This is a high-level convenience function that automatically builds
    the header index and yields FASTA entries in chunks. Memory-efficient
    for large FASTA files.

    Optimizations:
    - Uses memory-mapped file access (mmap) for faster reads
    - Uses replace() for removing newlines (C-level, fast)
    - Optional regex pattern to extract protein_id from header
    - Optional progress tracking with Rich progress bars

    Args:
        path: Path to the FASTA file.
        chunk_size: Number of entries per chunk (default: 1000).
        header_pattern: Optional regex pattern (string or compiled Pattern) to extract
            protein_id from header. If provided, the first capture group is used.
            Example: r"^sp\\|([A-Z0-9]+)\\|" to extract UniProt accession from
            "sp|P12345|PROTEIN_NAME".
        show_progress: Display progress bar (default: True).
        progress: Existing Progress instance to share progress context.

    Yields:
        Lists of dicts, each with 'protein_id' and 'sequence' keys. Each chunk contains
        up to chunk_size entries. The final chunk may have fewer entries.

    Raises:
        FileNotFoundError: If the file doesn't exist.
        IOError: If the file cannot be read.
        ValueError: If chunk_size is not positive.

    Example:
        >>> # Read as strings (default)
        >>> for chunk in read_fasta_chunks("proteins.fasta", chunk_size=1000):
        ...     for entry in chunk:
        ...         pid = entry['protein_id']
        ...         seq_len = len(entry['sequence'])
        ...         print(f"Protein ID: {pid}, Sequence length: {seq_len}")
        ...
        >>> # Extract UniProt accession from header
        >>> pattern = re.compile(r"^sp\\|([A-Z0-9]+)\\|")
        >>> for chunk in read_fasta_chunks("proteins.fasta", header_pattern=pattern):
        ...     for entry in chunk:
        ...         print(f"Accession: {entry['protein_id']}")
    """
    if chunk_size <= 0:
        raise ValueError(f"chunk_size must be positive, got {chunk_size}")

    # Compile pattern if string provided
    compiled_pattern: re.Pattern[str] | None = None
    if header_pattern:
        if isinstance(header_pattern, str):
            compiled_pattern = re.compile(header_pattern)
        else:
            compiled_pattern = header_pattern

    # Build index to get total count for progress tracking
    offsets = build_header_index(path)
    total_entries = len(offsets)

    # Pass progress parameters to iter_fasta_batches
    yield from iter_fasta_batches(
        path,
        offsets,
        chunk_size,
        header_pattern=compiled_pattern,
        show_progress=show_progress,
        progress=progress,
        total=total_entries,
    )


def iter_fasta_batches(
    path: str,
    offsets: list[int],
    batch_size: int,
    header_pattern: re.Pattern[str] | None = None,
    show_progress: bool = True,
    progress: Progress | None = None,
    total: int | None = None,
) -> Iterator[list[dict[str, str]]]:
    """Iterate over FASTA entries in batches using memory-mapped access.

    Optimized version using mmap and replace() for removing newlines.

    Args:
        path: Path to the FASTA file.
        offsets: List of byte offsets for each sequence header, typically
            obtained from build_header_index().
        batch_size: Number of entries to yield in each batch.
        header_pattern: Optional regex pattern to extract protein_id from header.
        show_progress: Display progress bar (default: True).
        progress: Existing Progress instance to share progress context.
        total: Total number of entries for progress tracking (default: len(offsets)).

    Yields:
        Lists of dicts, each with 'protein_id' and 'sequence' keys. Each batch contains
        up to batch_size entries. The final batch may have fewer entries.

    Raises:
        FileNotFoundError: If the file doesn't exist.
        IOError: If the file cannot be read.
    """
    file_size = os.path.getsize(path)
    end_offsets = [*offsets[1:], file_size]
    total_entries = total if total is not None else len(offsets)

    # Set up progress tracking
    if show_progress:
        progress_instance = create_progress(progress=progress)
        # Only use 'with' context if we created a new progress instance
        if progress is None:
            with progress_instance:
                task_id = progress_instance.add_task("Read FASTA", total=total_entries)
                yield from _iter_fasta_batches_impl(
                    path,
                    offsets,
                    end_offsets,
                    batch_size,
                    header_pattern,
                    progress_instance,
                    task_id,
                )
        else:
            # Shared progress - don't use context manager
            task_id = progress_instance.add_task("Read FASTA", total=total_entries)
            yield from _iter_fasta_batches_impl(
                path,
                offsets,
                end_offsets,
                batch_size,
                header_pattern,
                progress_instance,
                task_id,
            )
    else:
        # No progress tracking
        yield from _iter_fasta_batches_impl(
            path, offsets, end_offsets, batch_size, header_pattern, None, None
        )


def _iter_fasta_batches_impl(
    path: str,
    offsets: list[int],
    end_offsets: list[int],
    batch_size: int,
    header_pattern: re.Pattern[str] | None,
    progress: Progress | None,
    task_id: TaskID | None,
) -> Iterator[list[dict[str, str]]]:
    """Internal implementation of iter_fasta_batches with progress tracking."""
    with open(path, "rb") as f, mmap.mmap(f.fileno(), 0, access=mmap.ACCESS_READ) as mm:
        batch: list[dict[str, str]] = []

        for start, end in zip(offsets, end_offsets, strict=False):
            entry = _read_sequence_from_offset_mmap(mm, start, end, header_pattern=header_pattern)
            batch.append(entry)

            if len(batch) >= batch_size:
                yield batch
                if task_id is not None and progress is not None:
                    progress.update(task_id, advance=len(batch))
                batch = []

        if batch:
            yield batch
            if task_id is not None and progress is not None:
                progress.update(task_id, advance=len(batch))


def read_fasta_range(
    path: str,
    offsets: list[int],
    start: int,
    end: int,
    header_pattern: re.Pattern[str] | None = None,
    show_progress: bool = True,
    progress: Progress | None = None,
) -> list[dict[str, str]]:
    """Read a range of FASTA entries by their ordinal positions.

    Optimized version using mmap and replace() for removing newlines.

    Args:
        path: Path to the FASTA file.
        offsets: List of byte offsets for each sequence header, typically
            obtained from build_header_index().
        start: Starting entry index (inclusive, 0-based).
        end: Ending entry index (exclusive).
        header_pattern: Optional regex pattern to extract protein_id from header.
        show_progress: Display progress bar (default: True).
        progress: Existing Progress instance to share progress context.

    Returns:
        List of dicts, each with 'protein_id' and 'sequence' keys, for indices [start, end).

    Raises:
        ValueError: If the range [start, end) is invalid.
        FileNotFoundError: If the file doesn't exist.
        IOError: If the file cannot be read.
    """
    if not (0 <= start < end <= len(offsets)):
        raise ValueError(
            f"Invalid range [{start}, {end}): must satisfy 0 <= start < end <= {len(offsets)}"
        )

    file_size = os.path.getsize(path)
    end_offsets = [*offsets, file_size]
    total_entries = end - start

    entries: list[dict[str, str]] = []

    # Set up progress tracking
    if show_progress:
        progress_instance = create_progress(progress=progress)
        # Only use 'with' context if we created a new progress instance
        if progress is None:
            with progress_instance:
                task_id = progress_instance.add_task("Read FASTA Range", total=total_entries)
                with open(path, "rb") as f, mmap.mmap(f.fileno(), 0, access=mmap.ACCESS_READ) as mm:
                    for i in range(start, end):
                        entry = _read_sequence_from_offset_mmap(
                            mm, offsets[i], end_offsets[i + 1], header_pattern=header_pattern
                        )
                        entries.append(entry)
                        progress_instance.update(task_id, advance=1)
        else:
            # Shared progress - don't use context manager
            task_id = progress_instance.add_task("Read FASTA Range", total=total_entries)
            with open(path, "rb") as f, mmap.mmap(f.fileno(), 0, access=mmap.ACCESS_READ) as mm:
                for i in range(start, end):
                    entry = _read_sequence_from_offset_mmap(
                        mm, offsets[i], end_offsets[i + 1], header_pattern=header_pattern
                    )
                    entries.append(entry)
                    progress_instance.update(task_id, advance=1)
    else:
        # No progress tracking
        with open(path, "rb") as f, mmap.mmap(f.fileno(), 0, access=mmap.ACCESS_READ) as mm:
            for i in range(start, end):
                entry = _read_sequence_from_offset_mmap(
                    mm, offsets[i], end_offsets[i + 1], header_pattern=header_pattern
                )
                entries.append(entry)

    return entries


if __name__ == "__main__":
    from rich import print

    path = "/home/mha/projects/proteingraph/downloads/uniprot_sprot.fasta"

    offsets = build_header_index(path)
    print(f"Found {len(offsets)} sequences")

    # Test with strings (default)
    sequences = list(read_fasta_chunks(path, chunk_size=100, header_pattern=r"(?<=\|)[^|]+(?=\|)"))
    print(len(sequences))

    # print first sequence
    print(sequences[0][0])
