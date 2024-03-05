"""
Copyright 2023 Tom Stanton (tomdstanton@gmail.com)
https://github.com/klebgenomics/Kaptive

This file is part of Kaptive. Kaptive is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by the Free Software Foundation,
either version 3 of the License, or (at your option) any later version. Kaptive is distributed
in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
details. You should have received a copy of the GNU General Public License along with Kaptive.
If not, see <https://www.gnu.org/licenses/>.
"""
from __future__ import annotations

import re
from typing import Iterable, Generator
from operator import attrgetter, le, ge

from kaptive.intrange import range_overlap

# Constants -----------------------------------------------------------------------------------------------------------
# _CS_OPS_REGEX = re.compile(r'([=:\*\+\-~])([0-9]+)?([acgtn]+)?')  # Regex to match cs tag operations
_CIGAR_OP_REGEX = re.compile(r'(?P<num_bases>[0-9]+)(?P<operation>[MIDNSHPX=])')  # Regex to match cigar operations
_BITWISE_FLAGS = [
    "is_paired", "is_proper_pair", "is_unmapped", "mate_is_unmapped", "is_reverse", "mate_is_reverse", "is_read1",
    "is_read2", "is_secondary", "is_qcfail", "is_duplicate", "is_supplementary"
]


# Classes -------------------------------------------------------------------------------------------------------------
class AlignmentError(Exception):
    pass


class Alignment:
    def __init__(
            self, query_name: str | None = None, query_length: int | None = 0, query_start: int | None = 0,
            query_end: int | None = 0, strand: str | None = None, target_name: str | None = None,
            target_length: int | None = 0, target_start: int | None = 0, target_end: int | None = 0,
            matching_bases: int | None = 0, num_bases: int | None = 0, mapping_quality: int | None = 0,
            flag: dict[str, bool] | None = None, tags: dict | None = None, query_sequence: str | None = None,
            target_sequence: str | None = None, base_quality: list[tuple[str]] | None = None,
            percent_identity: float | None = 0.0, percent_query_coverage: float | None = 0.0,
            percent_target_coverage: float | None = 0.0):
        self.query_name = query_name or ''  # Query sequence name
        self.query_length = query_length  # Query sequence length
        self.query_start = query_start  # Query start coordinate (0-based)
        self.query_end = query_end  # Query end coordinate (0-based)
        self.strand = strand or 'unknown'  # ‘+’ if query/target on the same strand; ‘-’ if opposite
        self.target_name = target_name or ''  # Target sequence name
        self.target_length = target_length  # Target sequence length
        self.target_start = target_start  # Target start coordinate on the original strand (0-based)
        self.target_end = target_end  # Target end coordinate on the original strand (0-based)
        self.matching_bases = matching_bases  # Number of matching bases in the alignment
        self.num_bases = num_bases  # Number bases, including gaps, in the alignment
        self.mapping_quality = mapping_quality  # Mapping quality (0-255 with 255 for missing)
        self.flag = flag or {}  # {flag: bool} bitwise flags from SAM
        self.tags = tags or {}  # {tag: value} pairs
        self.query_sequence = query_sequence or ''  # Query sequence string
        self.target_sequence = target_sequence or ''  # Target sequence string
        self.base_quality = base_quality or []  # Tuple of (base, quality) pairs from SAM
        self.percent_identity = percent_identity
        self.percent_query_coverage = percent_query_coverage
        self.percent_target_coverage = percent_target_coverage

    @classmethod
    def from_paf_line(cls, line: str | bytes, parse_cigar: bool = False, **kwargs):
        if not (line := line.decode().strip() if isinstance(line, bytes) else line.strip()):
            raise AlignmentError("Empty line")
        if len(line := line.split('\t')) < 12:
            raise AlignmentError(f"Line has < 12 columns: {line}")
        self = cls(  # Parse standard fields
            query_name=line[0], query_length=int(line[1]), query_start=int(line[2]), query_end=int(line[3]),
            strand=line[4], target_name=line[5], target_length=int(line[6]), target_start=int(line[7]),
            target_end=int(line[8]), matching_bases=int(line[9]), num_bases=int(line[10]),
            mapping_quality=int(line[11]),
            **kwargs
        )
        self.percent_identity = self.matching_bases / self.num_bases * 100.0
        self.percent_query_coverage = 100.0 * (self.query_end - self.query_start) / self.query_length
        self.percent_target_coverage = 100.0 * (self.target_end - self.target_start) / self.target_length
        for tag in line[12:]:  # Parse tags
            x, y, z = tag.split(":", 2)
            self.tags[x] = int(z) if y == "i" else float(z) if y == "f" else z
        if parse_cigar and (cigar_string := self.tags.get('cg', None)):
            self.tags['cg'] = Cigar(cigar_string)  # Create cigar object
        return self

    @classmethod
    def from_sam_line(cls, line: str | bytes, parse_cigar: bool = False, **kwargs):
        if not (line := line.decode().strip() if isinstance(line, bytes) else line.strip()):
            raise AlignmentError("Empty line")
        if len(line := line.split('\t')) < 12:
            raise AlignmentError(f"Line has < 12 columns: {line}")
        self = Alignment(  # Parse standard fields
            query_name=line[0], target_name=line[2], target_start=int(line[3]) - 1,  # SAM is 1-based
            mapping_quality=int(line[4]), base_quality=[(x, ord(y) - 33) for x, y in zip(line[9], line[10])],
            tags={'cg': line[5], 'rnext': line[6], 'pnext': line[7], 'tlen': int(line[8])}, **kwargs
        )
        self.flag |= {k: bool(int(v)) for k, v in zip(_BITWISE_FLAGS, reversed(format(int(line[1]), '012b')))}
        self.strand = '-' if self.flag['is_reverse'] else '+'  # Set strand
        for tag in line[11:]:  # Parse tags
            x, y, z = tag.split(":", 2)
            self.tags[x] = int(z) if y == "i" else float(z) if y == "f" else z
        if parse_cigar and (cigar_string := self.tags.get('cg', None)):
            self.tags['cg'] = (cigar := Cigar(cigar_string))  # Create cigar object
            self.matching_bases = cigar.M - (self.tags['NM'] - cigar.mm) + self.tags['nn']
            self.num_bases = (cigar.M + cigar.I + cigar.D) - self.tags['nn']
            self.query_start = cigar.clip_right if self.flag['is_reverse'] else cigar.clip_left
            self.query_length = cigar.ql
            qlen = cigar.M + cigar.I + cigar.clip_left + cigar.clip_right  # query length
            # assert qlen == self.query_length, AlignmentError(f"Query length {qlen} != {self.query_length}")
            self.query_end = qlen - cigar.clip_left if self.flag['is_reverse'] else qlen - cigar.clip_right
            self.target_end = self.target_start + cigar.M + cigar.D + cigar.N
        return self

    @classmethod
    def from_blast_line(cls, line: str | bytes, **kwargs):
        if not (line := line.decode().strip() if isinstance(line, bytes) else line.strip()):
            raise AlignmentError("Empty line")
        if len(line := line.split('\t')) < 12:
            raise AlignmentError(f"Line has < 12 columns: {line}")
        return cls(
            query_name=line[0], target_name=line[1], percent_identity=float(line[2]), matching_bases=int(line[3]),
            num_bases=int(line[3]), query_start=int(line[6]) - 1, query_end=int(line[7]), target_start=int(line[8]) - 1,
            target_end=int(line[9]), tags={'evalue': float(line[10]), 'bitscore': float(line[11])}, **kwargs
        )

    def __repr__(self):
        return (f'{self.query_name}:{self.query_start}-{self.query_end} '
                f'{self.target_name}:{self.target_start}-{self.target_end}')

    def __len__(self):
        return self.num_bases

    def __getattr__(self, item):
        if item in self.__dict__:  # First check attributes
            return self.__dict__[item]
        elif item in self.tags:  # Then check tags
            return self.tags[item]
        else:
            raise AttributeError(f"{self.__class__.__name__} object has no attribute {item}")

    def conflicts(self, other: Alignment, overlap_fraction: float = 0.5):
        """
        Returns whether this hit conflicts with the other hit on the target sequence.
        A conflict is defined as the hits overlapping by 50% or more of the shortest hit's length.
        A hit is not considered to conflict with itself.
        """
        if self.query_name == other.query_name:
            return False
        if self.target_name != other.target_name:
            return False
        overlap = range_overlap((self.target_start, self.target_end),(other.target_start, other.target_end), skip_sort=True)
        return overlap / self.num_bases > overlap_fraction


class Cigar:
    def __init__(self, cigar: str | None = None):
        """Simple class to parse and store cigar strings
        :param cigar: Cigar string"""
        self._cigar = cigar or ''
        self._operations = []  # list[(op: str, num_bases: int)] Populated when cigar has been iterated over once.
        [setattr(self, i, 0) for i in ('M', 'N', 'ql', 'tl', 'mm', 'clip_left', 'clip_right',
                                       'I', 'D')]
        self.ext_cigar = False
        if len(self._cigar) > 0:
            self.parse()

    def __len__(self):
        return len(self._operations)

    def __str__(self):
        return self._cigar

    def __iter__(self) -> Generator[tuple[str, int]]:
        if self._operations:
            yield from self._operations
        else:
            for match in _CIGAR_OP_REGEX.finditer(self._cigar):
                self._operations.append(x := (match.group('operation'), int(match.group('num_bases'))))
                yield x

    def parse(self):
        for operation, n_bases in self.__iter__():
            if operation == "M":  # Match or mismatch, consumes reference and query
                [setattr(self, i, getattr(self, i) + n_bases) for i in ('M', 'ql', 'tl')]
                self.ext_cigar = False
            elif operation == "I":  # Insertion to the reference, consumes query
                [setattr(self, i, getattr(self, i) + n_bases) for i in ('I', 'ql')]
            elif operation == "D":  # Deletion from the reference, consumes reference
                [setattr(self, i, getattr(self, i) + n_bases) for i in ('D', 'tl')]
            elif operation == "N":  # Skipped region from the reference, consumes reference
                [setattr(self, i, getattr(self, i) + n_bases) for i in ('N', 'tl')]
            elif operation == "S":  # Soft clipping (clipped sequences present in SEQ), consumes query
                setattr(self, 'clip_left' if self.M == 0 else 'clip_right', n_bases)
                self.ql += n_bases
            elif operation == "H":  # Hard clipping (clipped sequences NOT present in SEQ), consumes nothing
                setattr(self, 'clip_left' if self.M == 0 else 'clip_right', n_bases)
            elif operation == "=":  # Sequence match, consumes reference and query
                [setattr(self, i, getattr(self, i) + n_bases) for i in ('M', 'ql', 'tl')]
                self.ext_cigar = True
            elif operation == "X":  # Sequence mismatch, consumes reference and query
                [setattr(self, i, getattr(self, i) + n_bases) for i in ('M', 'ql', 'tl', 'mm')]
                self.ext_cigar = True


class Pileup:
    """Performs samtools mpileup, see: https://www.htslib.org/doc/samtools-mpileup.html"""

    def __init__(self, positions: dict[int, PileupPosition] | None = None):
        self.positions = positions or {}  # {position: PileupPosition}

    @classmethod
    def from_alignments(cls, alignments: list[Alignment]):
        self = cls()
        [self.add_alignment(alignment) for alignment in alignments]
        return self

    def add_alignment(self, alignment: Alignment):
        ref_pos = alignment.target_start  # index for alignment.target_sequence
        total_bases = 0  # index for alignment.base_quality
        for n_op, (operation, n_bases) in enumerate(alignment.tags['cg']):
            if operation == 'S':  # Skip soft-clipped bases
                total_bases += n_bases
                continue
            if operation in 'M=X':
                for i in range(n_bases):
                    base, quality = alignment.base_quality[total_bases]
                    pos = self.positions.get(
                        ref_pos, PileupPosition(ref_pos, reference_base=alignment.target_sequence[ref_pos])
                    )
                    pos.add_instance('.' if base == pos.reference_base else base, quality, alignment.mapping_quality)
                    self.positions[ref_pos] = pos
                    ref_pos += 1  # Increment reference position
                    total_bases += 1  # Increment total bases as both reference and query have moved forward
            elif operation == 'D':
                for i in range(n_bases):
                    base, quality = alignment.base_quality[total_bases]
                    self.positions.get(
                        ref_pos, PileupPosition(ref_pos, reference_base=alignment.target_sequence[ref_pos])
                    ).add_instance("*", quality, alignment.mapping_quality)
                    ref_pos += 1  # Increment reference position, consume reference bases so don't increment total_bases

    def call(self, **kwargs) -> dict[int, str]:
        return {i: j.call(**kwargs) for i, j in self.positions.items()}

    def get_consensus(self, **kwargs):
        # We need to add "-" if the first position is not 0
        return ''.join(self.call(**kwargs).values())


class PileupPosition:
    """Class representing a single line in a samtools mpileup file"""

    def __init__(self, position: int | None = 0, reference_base: str | None = None, bases: list[str] | None = None,
                 base_qualities: list[int] | None = None, mapping_qualities: list[int] | None = None):
        self.position = position
        self.reference_base = reference_base
        self.bases = bases or []
        self.base_qualities = base_qualities or []
        self.mapping_qualities = mapping_qualities or []

    def __repr__(self):
        return f'{self.position}:{self.reference_base}'

    def __len__(self):
        return len(self.bases)

    def __iter__(self) -> Generator[tuple[str, int, int]]:
        return zip(self.bases, self.base_qualities, self.mapping_qualities)

    def add_instance(self, base: str, base_quality: int, mapping_quality: int):
        self.bases.append(base)
        self.base_qualities.append(base_quality)
        self.mapping_qualities.append(mapping_quality)

    def call(self, min_base_quality: int = 0, min_mapping_quality: int = 0) -> str:
        # Currently use the most common base, but could use a more sophisticated method in the future
        if not self.bases:
            return 'N'
        best = max(
            (i for i in self if all([i[1] >= min_base_quality, i[2] >= min_mapping_quality])),
            key=lambda x: self.bases.count(x[0])
        )[0]
        if best == '*':
            return '-'
        if best == '.':
            return self.reference_base


# Functions ------------------------------------------------------------------------------------------------------------
def get_best_alignments(alignments: Iterable[Alignment], group: str = "query_name", metric: str = "matching_bases",
                        reverse_sort: bool = True, threshold: float = 0.0) -> Generator[Alignment, None, None]:
    """
    Get the best alignments for each group and can additionally filter by calling filter_alignments internally
    """
    current_group = None
    for alignment in sorted(
            filter_alignments(alignments, metric, threshold, reverse_sort),
            key=lambda x: attrgetter(group, metric)(x), reverse=reverse_sort
    ):
        if current_group != getattr(alignment, group):
            current_group = getattr(alignment, group)
            yield alignment


def filter_alignments(alignments: Iterable[Alignment], metric: str = "matching_bases", threshold: float = 0.0,
                      reverse_sort: bool = True) -> Generator[Alignment, None, None]:
    """
    Filter alignments based on a metric and threshold
    """
    operator = ge if reverse_sort else le
    for alignment in alignments:
        if not hasattr(alignment, metric):
            raise AttributeError(f"{alignment} does not have attribute {metric}")
        if not isinstance(getattr(alignment, metric), (int, float)):
            raise TypeError(f"{alignment} attribute {metric} is not a number")
        if magic(getattr(alignment, metric), threshold, operator):
            yield alignment


def magic(left: int | float, right: int | float, operator: callable):
    """
    Takes two numbers and an operator and returns the result of the comparison.
    """
    return operator(left, right)


def cull_conflicting_alignments(alignment_to_keep: Alignment, alignments: Iterable[Alignment]
                                ) -> Generator[Alignment, None, None]:
    """
    Returns a generator of alignments that do not conflict with the alignment_to_keep.
    param alignment_to_keep: Alignment object
    param alignments: Iterable of Alignment objects
    return: Generator of Alignment objects
    """
    for x in alignments:
        if not x.conflicts(alignment_to_keep):
            yield x


def cull_all_conflicting_alignments(alignments: Iterable[Alignment], sort_by: str = 'AS',
                                    sort_large_to_small: bool = True) -> list[Alignment]:
    """
    Returns a list of alignments that do not conflict with each other.
    param alignments: Iterable of Alignment objects
    param sort_by: Attribute of Alignment to sort by
    param sort_large_to_small: Whether to sort by the attribute from large to small or small to large
    return: List of Alignment objects
    """
    kept_alignments = []
    sorted_alignments = sorted(list(alignments), key=lambda x: getattr(x, sort_by), reverse=sort_large_to_small)
    # There must be a more efficient implementation of this
    while sorted_alignments:
        kept_alignments.append(sorted_alignments.pop(0))
        sorted_alignments = list(cull_conflicting_alignments(kept_alignments[-1], sorted_alignments))
    return kept_alignments


def alignments_in_range(alignments: Iterable[Alignment], start: int, end: int, query: bool = False
                        ) -> Generator[Alignment, None, None]:
    """
    Return the alignments that are within the given range
    """
    if query:
        return (i for i in alignments if i.query_start >= start and i.query_end <= end)
    else:
        return (i for i in alignments if i.target_start >= start and i.target_end <= end)


def alignments_overlapping_range(alignments: Iterable[Alignment], start: int, end: int, query: bool = False
                                 ) -> Generator[Alignment, None, None]:
    """
    Return the alignments that are within or overlap the given range
    """
    if query:
        return (i for i in alignments if start <= i.query_start <= end or start <= i.query_end <= end)
    else:
        return (i for i in alignments if start <= i.target_start <= end or start <= i.target_end <= end)
