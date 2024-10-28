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
from itertools import groupby
from kaptive.utils import range_overlap


# Classes -------------------------------------------------------------------------------------------------------------
class AlignmentError(Exception):
    pass


class Alignment:
    """
    Similar to `mappy.Alignment` but with additional attributes and methods.
    """

    def __init__(
            self, q: str | None = None, q_len: int | None = 0, q_st: int | None = 0,
            q_en: int | None = 0, strand: str | None = None, ctg: str | None = None,
            ctg_len: int | None = 0, r_st: int | None = 0, r_en: int | None = 0,
            mlen: int | None = 0, blen: int | None = 0, mapq: int | None = 0,
            tags: dict | None = None):
        self.q = q or 'unknown'  # Query sequence name
        self.q_len = q_len  # Query sequence length
        self.q_st = q_st  # Query start coordinate (0-based)
        self.q_en = q_en  # Query end coordinate (0-based)
        self.strand = strand or 'unknown'  # ‘+’ if query/target on the same strand; ‘-’ if opposite
        self.ctg = ctg or 'unknown'  # Target sequence name
        self.ctg_len = ctg_len  # Target sequence length
        self.r_st = r_st  # Target start coordinate on the original strand (0-based)
        self.r_en = r_en  # Target end coordinate on the original strand (0-based)
        self.mlen = mlen  # Number of matching bases in the alignment
        self.blen = blen  # Number bases, including gaps, in the alignment
        self.mapq = mapq  # Mapping quality (0-255 with 255 for missing)
        self.tags = tags or {}  # {tag: value} pairs

    @classmethod
    def from_paf_line(cls, line: str):
        """
        Parse a line in PAF format and return an Alignment object.
        """
        if len(line := line.split('\t')) < 12:
            raise AlignmentError(f"Line has < 12 columns: {line}")
        try:
            return Alignment(  # Parse standard fields
                q=line[0], q_len=int(line[1]), q_st=int(line[2]), q_en=int(line[3]), strand=line[4], ctg=line[5],
                ctg_len=int(line[6]), r_st=int(line[7]), r_en=int(line[8]), mlen=int(line[9]), blen=int(line[10]),
                mapq=int(line[11]),
                tags={(x := t.split(":", 2))[0]: int(x[2]) if x[1] == "i" else float(x[2]) if x[1] == "f" else x[2] for
                      t in line[12:]}
            )
        except Exception as e:
            raise AlignmentError(f"Error parsing PAF line: {line}") from e

    def __repr__(self):
        return f'{self.q}:{self.q_st}-{self.q_en} {self.ctg}:{self.r_st}-{self.r_en} {self.strand}'

    def __len__(self):
        return self.num_bases

    def __getattr__(self, item):
        if item in self.__dict__:  # First check attributes
            return self.__dict__[item]
        elif item in self.tags:  # Then check tags
            return self.tags[item]
        else:
            raise AttributeError(f"{self.__class__.__name__} object has no attribute {item}")

    @property
    def partial(self) -> bool:
        if self.q_len > self.ctg_len:
            return True
        if (self.r_en == self.ctg_len or self.r_st == 0) and (self.q_en - self.q_st) < self.q_len:
            return True
        ref_start_diff = self.r_st
        ref_end_diff = self.ctg_len - self.r_en
        query_start_diff = self.q_st if self.strand == '+' else self.q_len - self.q_en
        query_end_diff = self.q_len - self.q_en if self.strand == '+' else self.q_st
        if query_start_diff > ref_start_diff or query_end_diff > ref_end_diff:
            return True
        return False


# Functions ------------------------------------------------------------------------------------------------------------
def group_alns(alignments: Iterable[Alignment], key: str = 'q') -> Generator[tuple[str, Generator[Alignment]]]:
    """Group alignments by a key"""
    yield from groupby(sorted(alignments, key=lambda x: getattr(x, key)), key=lambda x: getattr(x, key))


def cull(keep: Alignment, alignments: Iterable[Alignment],
         overlap_fraction: float = 0.1) -> Generator[Alignment]:
    """Yield alignments that do not conflict with keep alignment"""
    for a in alignments:
        if (a.ctg != keep.ctg or  # Different contig
                range_overlap((a.r_st, a.r_en), (keep.r_st, keep.r_en), skip_sort=True) / a.blen < overlap_fraction):
            yield a


def cull_all(alignments: list[Alignment]) -> list[Alignment]:
    kept_alignments = []
    sorted_alignments = sorted(list(alignments), key=lambda x: x.mlen, reverse=True)
    while sorted_alignments:
        kept_alignments.append(sorted_alignments.pop(0))
        sorted_alignments = list(cull(kept_alignments[-1], sorted_alignments))
    return kept_alignments


def cull_filtered(pred: callable, alignments: Iterable[Alignment]) -> Generator[Alignment]:
    """
    Cull and flatten alignments that don't overlap with alignments matching the predicate.
    E.g. list(cull_filtered(lambda i: i.q in query_genes, alignments))
    """
    keep, other = [], []
    [(keep if pred(a) else other).append(a) for a in alignments]
    other = cull_all(other)  # Remove conflicting other alignments
    for i in keep:  # Remove other alignments overlapping best match gene alignments
        other = list(cull(i, other))
        yield i
    yield from other


# def extend_aln_ranges(q_len: int, q_st: int, q_en: int, ctg_len: int, r_st: int, r_en: int, tolerance: int = 10
#                       ) -> tuple[int, int, int, int]:
#     """
#     Extend query alignment ranges to the contig boundaries if the distance is within the tolerance and less than
#     or equal to the query length
#     """
#     if start_dist := min(q_st, r_st) <= tolerance:
#         q_st -= start_dist
#         r_st -= start_dist
#     if end_dist := min(ctg_len - q_en, q_len) <= tolerance:
#         q_en += end_dist
#         r_en += end_dist
#     return q_st, q_en, r_st, r_en
