"""
Module to handle query (contigs) and target (features) IO.
"""
from bz2 import open as bzopen
from collections.abc import Callable, Iterable, Iterator
from dataclasses import dataclass, field
from gzip import open as gzopen
from lzma import open as lzopen
from pathlib import Path
from re import compile as re_compile
from typing import IO, ClassVar, Self

from kaptive.core.seq import SeqRecord, Sequences


# Classes --------------------------------------------------------------------------------------------------------------
class FastaReader(Iterator):
    """
    A high-performance FASTA file reader.

    Assumes the provided handle is opened in binary mode ('rb') to return the
    sequence as a byte string, which is highly efficient and ideal for wrapped
    FASTA lines.
    """

    def __init__(self, handle: IO[bytes]):
        self._handle = handle
        self._generator = self._parse_records()

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        self._handle.close()

    def __del__(self):
        self._handle.close()

    def __iter__(self) -> Iterator[SeqRecord]:
        return self

    def __next__(self) -> SeqRecord:
        return next(self._generator)

    def _parse_records(self) -> Iterator[SeqRecord]:
        header: str = ""
        seq_chunks: list[bytes] = []

        # Local variable caching speeds up lookups inside the tight loop
        join_bytes = b"".join

        for line in self._handle:
            # rstrip() is faster than strip() and removes \r\n or \n
            line = line.rstrip()
            if not line:
                continue

            # 62 is the ASCII integer value for '>'
            if line[0] == 62:
                if header:
                    # Yield the previous record
                    yield SeqRecord(seq=join_bytes(seq_chunks), id=header)

                seq_chunks.clear()

                decoded_line = line[1:].decode("utf-8", errors="replace")
                header, _, _ = decoded_line.partition(" ")

            else:
                seq_chunks.append(line)

        # Yield the final record once the loop ends
        if header:
            yield SeqRecord(seq=join_bytes(seq_chunks), id=header)


@dataclass(slots=True, frozen=True)
class GenomeAssembly:
    """Container for a genome assembly in FASTA format, with support for transparent decompression.

    Attributes:
        id (str): The genome assembly identifier.
        contigs (Sequences): The contig sequences as a contiguous block of memory.
        id_map (dict[str, int]): A mapping from contig ID to contig sequence index (offset).
    """
    _SEQUENCE_FILE_REGEX = re_compile(
        r'\.('
        r'f(asta|a|na|fn|as)|'
        r')\.?(?P<compression>(gz|bz2|xz))?$'
    )
    _OPENERS: ClassVar[dict[str, Callable]] = {'gz': gzopen, 'bz2': bzopen, 'xz': lzopen}
    id: str
    contigs: Sequences
    id_map: dict[str, int] = field(init=False, repr=False, hash=False, compare=False)

    def __post_init__(self):
        object.__setattr__(self, 'id_map', {name: i for i, name in enumerate(self.contigs.ids)})

    @classmethod
    def ensure(cls, genome: Self | str | Path | IO[bytes]) -> Self:
        """Ensures the input is a GenomeAssembly, loading it from file if necessary.

        Parameters:
            genome: A GenomeAssembly instance, or a path to a FASTA file or a FASTA binary stream.

        Returns:
            GenomeAssembly: The validated or loaded assembly.
        """
        if isinstance(genome, cls):
            return genome
        elif isinstance(genome, (str, Path)):
            return cls.from_file(genome)
        return cls.from_stream(genome)  # ty:ignore[invalid-argument-type]

    def __len__(self) -> int:
        """Total number of base pairs in the assembly."""
        return len(self.contigs.seqs)

    def __iter__(self) -> Iterator[SeqRecord]:
        return iter(self.contigs)

    def __str__(self):
        return self.id

    def __getitem__(self, item: str) -> bytes:
        """Access a contig sequence by its ID."""
        idx = self.id_map[item]
        s = self.contigs.offsets[idx]
        l = self.contigs.lengths[idx]
        return self.contigs.seqs[s:s + l].tobytes()

    @classmethod
    def from_file(cls, file: str | Path) -> Self:
        """Load a genome from a FASTA file.

        Parameters:
            file (Union[str, Path]): Path to the file. Supports .gz, .bz2, and .xz compression.
        """
        file = Path(file)  # type: Path
        if not (m := cls._SEQUENCE_FILE_REGEX.search(file.name)):
            raise NotImplementedError(f'Unsupported format: {file}')

        with cls._OPENERS.get(m.group('compression'), open)(file, mode='rb') as handle:
            return cls.from_stream(handle, file.name.rstrip(m.group()))

    @classmethod
    def from_stream(cls, handle: IO[bytes], id_: str | None = None) -> Self:
        """Load a genome from a binary stream.

        Parameters:
            handle (IO[bytes]): An opened binary stream.
            id_ (Optional[str]): An optional identifier for the genome, otherwise will extract from handle name.
        """
        with FastaReader(handle) as records:
            return cls.from_records(id_ or getattr(handle, 'name', 'unknown'), records)

    @classmethod
    def from_records(cls, id_: str, records: Iterable[SeqRecord]) -> Self:
        return cls(id_, Sequences.from_records(list(records)))