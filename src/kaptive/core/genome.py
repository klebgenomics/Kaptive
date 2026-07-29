r"""Module for loading, parsing, and managing genomic sequence assemblies and FASTA I/O.

This module provides high-performance FASTA file reading using `rammappy` Rust bindings
and encapsulates genome assemblies in the [`GenomeAssembly`][kaptive.core.genome.GenomeAssembly]
dataclass, supported by transparent decompression (`.gz`, `.bz2`, `.xz`).
"""

from __future__ import annotations

import threading
from bz2 import open as bzopen
from collections.abc import Callable, Iterable, Iterator
from dataclasses import dataclass, field
from gzip import open as gzopen
from lzma import open as lzopen
from pathlib import Path
from re import compile as re_compile
from typing import IO, Any, ClassVar, Self

from kaptive.core.seq import SeqRecord, Sequences


# Classes --------------------------------------------------------------------------------------------------------------
class FastaReader(Iterator):
    r"""High-performance FASTA file iterator.

    Parses raw binary FASTA streams into [`SeqRecord`][kaptive.core.seq.SeqRecord] instances
    using optimized C/Rust bindings via `rammappy`. Assumes the input handle is opened in
    binary mode (`'rb'`).

    Args:
        handle (IO[bytes]): An open binary stream containing FASTA sequence data.
    """

    def __init__(self, handle: IO[bytes]) -> None:
        r"""Initialize the FASTA reader and parse the underlying stream.

        Args:
            handle (IO[bytes]): Open binary file-like object containing FASTA data.
        """
        self._handle = handle
        import rammappy

        # Read the entire stream and parse using rammappy's high-performance Rust parser
        self._parsed = rammappy.fasta.parse_fasta_bytes(self._handle.read())
        self._generator = (SeqRecord(seq=seq, id=name) for name, seq in self._parsed)

    def __enter__(self) -> Self:
        r"""Enter the runtime context for the FASTA reader.

        Returns:
            FastaReader: The current reader context manager instance.
        """
        return self

    def __exit__(self, exc_type: Any, exc_val: Any, exc_tb: Any) -> None:
        r"""Exit the runtime context and close the underlying binary stream.

        Args:
            exc_type (Any): Exception type if an exception was raised inside context.
            exc_val (Any): Exception instance if an exception was raised.
            exc_tb (Any): Traceback object if an exception was raised.
        """
        self._handle.close()

    def __del__(self) -> None:
        r"""Clean up resources by closing the underlying binary handle if still open."""
        self._handle.close()

    def __iter__(self) -> Iterator[SeqRecord]:
        r"""Return the iterator object itself.

        Returns:
            Iterator[SeqRecord]: Iterator yielding parsed [`SeqRecord`][kaptive.core.seq.SeqRecord] objects.
        """
        return self

    def __next__(self) -> SeqRecord:
        r"""Fetch the next sequence record from the parsed FASTA stream.

        Returns:
            SeqRecord: The next FASTA record as a [`SeqRecord`][kaptive.core.seq.SeqRecord].

        Raises:
            StopIteration: When all FASTA records have been consumed.
        """
        return next(self._generator)


@dataclass(slots=True, frozen=True)
class GenomeAssembly:
    r"""Container for a genome assembly with support for transparent decompression.

    Stores contig sequence data as a contiguous memory block using [`Sequences`][kaptive.core.seq.Sequences]
    and supports thread-safe, lazy-cached index creation for high-performance sequence alignment.

    Attributes:
        id (str): Unique genome assembly identifier.
        contigs (Sequences): Structure-of-Arrays container [`Sequences`][kaptive.core.seq.Sequences]
            holding all contig sequences.
        id_map (dict[str, int]): Mapping from contig ID to sequence index offset within `contigs`.
        rammappy_index (Any): Lazy-cached `rammappy.Index` object for sequence lookups.
    """

    _SEQUENCE_FILE_REGEX = re_compile(
        r"\.(" r"f(asta|a|na|fn|as)|" r")\.?(?P<compression>(gz|bz2|xz))?$"
    )
    _OPENERS: ClassVar[dict[str, Callable]] = {"gz": gzopen, "bz2": bzopen, "xz": lzopen}
    id: str
    contigs: Sequences
    id_map: dict[str, int] = field(init=False, repr=False, hash=False, compare=False)
    rammappy_index: Any = field(default=None, init=False, repr=False, hash=False, compare=False)
    _index_lock: threading.Lock = field(
        default_factory=threading.Lock, init=False, repr=False, hash=False, compare=False
    )

    def __post_init__(self) -> None:
        r"""Initialize derived fields such as `id_map` after dataclass creation."""
        object.__setattr__(self, "id_map", {name: i for i, name in enumerate(self.contigs.ids)})

    @classmethod
    def ensure(cls, genome: Self | str | Path | IO[bytes]) -> Self:
        r"""Ensure the input object is coerced into a [`GenomeAssembly`][kaptive.core.genome.GenomeAssembly].

        Args:
            genome (Self | str | Path | IO[bytes]): Existing assembly object, file path,
                or binary FASTA stream.

        Returns:
            GenomeAssembly: Validated or loaded [`GenomeAssembly`][kaptive.core.genome.GenomeAssembly] instance.
        """
        if isinstance(genome, cls):
            return genome
        elif isinstance(genome, (str, Path)):
            return cls.from_file(genome)
        return cls.from_stream(genome)  # type: ignore

    def __len__(self) -> int:
        r"""Calculate the total number of base pairs across all contigs in the assembly.

        Returns:
            int: Cumulative sequence length in base pairs.
        """
        return len(self.contigs.seqs)

    def __iter__(self) -> Iterator[SeqRecord]:
        r"""Iterate over contigs in the assembly.

        Returns:
            Iterator[SeqRecord]: An iterator over [`SeqRecord`][kaptive.core.seq.SeqRecord] contigs.
        """
        return iter(self.contigs)

    def __str__(self) -> str:
        r"""Return the string representation of the assembly.

        Returns:
            str: Assembly identifier string (`id`).
        """
        return self.id

    def __getitem__(self, item: str) -> bytes:
        r"""Retrieve raw sequence bytes for a contig by its identifier.

        Args:
            item (str): Contig identifier string.

        Returns:
            bytes: Contig sequence bytes.

        Raises:
            KeyError: If `item` is not a recognized contig identifier in `id_map`.
        """
        idx = self.id_map[item]
        offset_val = self.contigs.offsets[idx]
        length_val = self.contigs.lengths[idx]
        return self.contigs.seqs[offset_val : offset_val + length_val].tobytes()

    def get_rammappy_index(self) -> Any:
        r"""Lazily build and return a thread-safe cached `rammappy.Index` for the assembly.

        Returns:
            Any: The compiled `rammappy.Index` instance.
        """
        if self.rammappy_index is None:
            with self._index_lock:
                if self.rammappy_index is None:
                    import rammappy

                    contig_seqs = [(c.id.encode(), c.seq) for c in self.contigs]
                    idx = rammappy.Index.build(contig_seqs)
                    object.__setattr__(self, "rammappy_index", idx)
        return self.rammappy_index

    @classmethod
    def from_file(cls, filepath: str | Path) -> Self:
        r"""Load a genome assembly from a FASTA file path.

        Supports plain `.fasta` / `.fa` / `.fna` files as well as compressed
        `.gz`, `.bz2`, and `.xz` files.

        Args:
            filepath (str | Path): Path to the target FASTA file.

        Returns:
            GenomeAssembly: Loaded [`GenomeAssembly`][kaptive.core.genome.GenomeAssembly] instance.

        Raises:
            NotImplementedError: If the file format or compression extension is unsupported.
        """
        filepath = Path(filepath)
        if not (m := cls._SEQUENCE_FILE_REGEX.search(filepath.name)):
            raise NotImplementedError(f"Unsupported format: {filepath}")

        with cls._OPENERS.get(m.group("compression"), open)(filepath, mode="rb") as handle:
            return cls.from_stream(handle, filepath.name.rstrip(m.group()))

    @classmethod
    def from_stream(cls, handle: IO[bytes], id_: str | None = None) -> Self:
        r"""Load a genome assembly from an open binary stream.

        Args:
            handle (IO[bytes]): Open binary FASTA stream handle.
            id_ (str | None): Custom assembly identifier. If None, derived from
                handle name or defaults to `'unknown'`.

        Returns:
            GenomeAssembly: Loaded [`GenomeAssembly`][kaptive.core.genome.GenomeAssembly] instance.
        """
        with FastaReader(handle) as records:
            return cls.from_records(id_ or getattr(handle, "name", "unknown"), records)

    @classmethod
    def from_records(cls, id_: str, records: Iterable[SeqRecord]) -> Self:
        r"""Construct a [`GenomeAssembly`][kaptive.core.genome.GenomeAssembly] from sequence records.

        Args:
            id_ (str): Genome assembly identifier.
            records (Iterable[SeqRecord]): Iterable of [`SeqRecord`][kaptive.core.seq.SeqRecord] objects.

        Returns:
            GenomeAssembly: Constructed [`GenomeAssembly`][kaptive.core.genome.GenomeAssembly] instance.
        """
        return cls(id_, Sequences.from_records(list(records)))