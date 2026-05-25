"""
Module to handle query (contigs) and target (features) IO.
"""
from dataclasses import dataclass, field
from bz2 import open as bzopen
from collections.abc import Iterator
from gzip import open as gzopen
from lzma import open as lzopen
from pathlib import Path
from re import compile as re_compile
from typing import IO, ClassVar, Callable

from kaptive.io import FastaReader, GfaReader
from kaptive.core.seq import SeqRecord, SeqBatch
from kaptive.core.graph import Edge, AssemblyGraph


# Classes --------------------------------------------------------------------------------------------------------------
@dataclass(slots=True, frozen=True)
class GenomeAssembly:
    """
    Container for a genome assembly, including contigs and their graph topology.

    Handles FASTA and GFA formats, with support for transparent decompression.
    """
    _SEQUENCE_FILE_REGEX = re_compile(
        r'\.('
        r'(?P<fasta>f(asta|a|na|fn|as|aa))|'
        r'(?P<gfa>gfa)|'
        r')\.?(?P<compression>(gz|bz2|xz))?$'
    )
    _OPENERS: ClassVar[dict[str, Callable]] = {'gz': gzopen, 'bz2': bzopen, 'xz': lzopen}
    id: str
    contigs: SeqBatch
    edges: list[Edge]
    _id_map: dict[str, int] = field(init=False, repr=False, hash=False, compare=False)

    def __post_init__(self):
        object.__setattr__(self, '_id_map', {name: i for i, name in enumerate(self.contigs.ids)})

    def __len__(self) -> int:
        """Total number of base pairs in the assembly."""
        return len(self.contigs.seqs)

    def __iter__(self) -> Iterator[SeqRecord]:
        return iter(self.contigs)

    def __str__(self):
        return self.id

    def __getitem__(self, item: str) -> bytes:
        """Access a contig sequence by its ID."""
        idx = self._id_map[item]
        s = self.contigs.offsets[idx]
        l = self.contigs.lengths[idx]
        return self.contigs.seqs[s:s + l].tobytes()

    @classmethod
    def from_file(cls, file: str | Path):
        """
        Load an assembly from a FASTA or GFA file.

        Args:
            file: Path to the file. Supports .gz, .bz2, and .xz compression.
        """
        file = Path(file) # type: Path
        if not (m := cls._SEQUENCE_FILE_REGEX.search(file.name)):
            raise NotImplementedError(f'Unsupported format: {file}')
        reader = FastaReader if m.group('fasta') else GfaReader
        with cls._OPENERS.get(m.group('compression'), open)(file, mode='rb') as handle:
            return cls.from_stream(handle, reader, file.name.rstrip(m.group()))

    @classmethod
    def from_stream(cls, handle: IO[bytes], reader, id_: str | None = None):
        """Load an assembly from an open file stream using the specified reader."""
        contigs, edges = [], []
        for record in reader(handle):
            if isinstance(record, SeqRecord):
                contigs.append(record)
            elif isinstance(record, Edge):
                edges.append(record)
        return cls(id_ or handle.name, SeqBatch.from_records(contigs), edges)

    def as_graph(self) -> AssemblyGraph:
        return AssemblyGraph(self)
