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

from typing import Iterable, Generator

from kaptive.misc import range_overlap


# Classes -------------------------------------------------------------------------------------------------------------
class AlignmentError(Exception):
    pass


class Alignment:
    """
    Class to store alignment information from PAF, SAM or BLAST tabular (--outfmt 6 / m8) output.
    It is purposely designed to be flexible and can be used with any of the three formats.
    """
    def __init__(
            self, query_name: str | None = None, query_length: int | None = 0, query_start: int | None = 0,
            query_end: int | None = 0, strand: str | None = None, target_name: str | None = None,
            target_length: int | None = 0, target_start: int | None = 0, target_end: int | None = 0,
            matching_bases: int | None = 0, num_bases: int | None = 0, mapping_quality: int | None = 0,
            tags: dict | None = None):
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
        self.tags = tags or {}  # {tag: value} pairs


    @classmethod
    def from_paf_line(cls, line: str | bytes):
        """
        Parse a line from a PAF file and return an Alignment object.
        Optionally parse the cigar string into a Cigar object (stored in tags).
        """
        if not (line := line.decode().strip() if isinstance(line, bytes) else line.strip()):
            raise AlignmentError("Empty line")
        if len(line := line.split('\t')) < 12:
            raise AlignmentError(f"Line has < 12 columns: {line}")
        self = cls(  # Parse standard fields
            query_name=line[0], query_length=int(line[1]), query_start=int(line[2]), query_end=int(line[3]),
            strand=line[4], target_name=line[5], target_length=int(line[6]), target_start=int(line[7]),
            target_end=int(line[8]), matching_bases=int(line[9]), num_bases=int(line[10]),
            mapping_quality=int(line[11]))
        self.tags = {x: int(z) if y == "i" else float(z) if y == "f" else z for tag in line[12:] for x, y, z in tag.split(":", 2)}
        return self

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


# Functions ------------------------------------------------------------------------------------------------------------
def cull_conflicting_alignments(
        a: Alignment, alignments: Iterable[Alignment], overlap_fraction: float = 0.1) -> Generator[Alignment, None, None]:
    """
    Returns a generator of alignments that do not conflict with the alignment_to_keep.
    param a: Alignment object to keep
    param alignments: Iterable of Alignment objects
    return: Generator of Alignment objects that don't conflict with a
    """
    for x in alignments:
        if not x.target_name == a.target_name and range_overlap(
                (x.target_start, x.target_end),(a.target_start, a.target_end), skip_sort=True) / x.num_bases > overlap_fraction:
            yield x


def cull_all_conflicting_alignments(alignments: Iterable[Alignment], sort_by: str = 'matching_bases',
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

