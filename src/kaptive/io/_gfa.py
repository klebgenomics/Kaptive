from typing import IO, Generator, Any
from collections.abc import Iterator

from kaptive.core.seq import SeqRecord
from kaptive.core.interval import Strand
from kaptive.core.genome import Edge
from kaptive.core.alignment import parse_cigar_string, CigarOp


# Readers --------------------------------------------------------------------------------------------------------------
class GfaReader(Iterator):
    """
    Reader for Graphical Fragment Assembly (GFA) files.

    Parses Segment (S) and Link (L) lines into SeqRecord and Edge objects.
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

    def __iter__(self) -> Iterator[SeqRecord | Edge]:
        return self

    def __next__(self) -> SeqRecord | Edge:
        return next(self._generator)

    @staticmethod
    def _parse_tags(parts: list[bytes]) -> Generator[tuple[str, Any], None, None]:
        for item in parts:
            tag, typ, val = item.split(b':', maxsplit=2)
            typ_str = typ.decode()
            if typ_str == 'f':
                val_parsed = float(val)
            elif typ_str == 'i':
                val_parsed = int(val)
            else:
                val_parsed = val.decode()
            yield tag.decode(), val_parsed

    @classmethod
    def _parse_segment(cls, parts: list[bytes]) -> SeqRecord:
        return SeqRecord(seq=parts[1], id=parts[0].decode())  #, list(cls._parse_tags(parts[2:]))

    @staticmethod
    def _parse_link(parts: list[str]) -> Edge:
        u = parts[0]
        u_strand = Strand(parts[1])
        v = parts[2]
        v_strand = Strand(parts[3])
        cigar_arr = parse_cigar_string(parts[4])
        matches = cigar_arr[(cigar_arr & 0xF) == CigarOp.M]
        overlap = int(matches[0] >> 4) if matches.size > 0 else 0
        return Edge(u, u_strand, v, v_strand, overlap)

    @classmethod
    def _parse_line(cls, line: bytes):
        if line.startswith(b'S\t'):
            return cls._parse_segment(line[2:].rstrip().split(b'\t'))
        elif line.startswith(b'L\t'):
            # Decode Link lines to text for easier string parsing downstream
            return cls._parse_link(line[2:].rstrip().decode().split('\t'))
        else:
            return None

    def _parse_records(self) -> Iterator[SeqRecord | Edge]:
        for line in self._handle:
            if (parsed := self._parse_line(line)) is not None:
                yield parsed
