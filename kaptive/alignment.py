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

from Bio.Seq import Seq

from kaptive.misc import decode_line
from kaptive.intrange import range_overlap


# Constants -----------------------------------------------------------------------------------------------------------
# _CS_OPS_REGEX = re.compile(r'([=:\*\+\-~])([0-9]+)?([acgtn]+)?')  # Regex to match cs tag operations
# _MAP_REGIONS_REGEX = re.compile(r'([0-9]+)([MIDNSP=X])')  # Regex to match cigar operations


# Classes -------------------------------------------------------------------------------------------------------------
class AlignmentError(Exception):
    pass


class Alignment:
    def __init__(self, query_name: str | None = None, query_length: int | None = 0,
                 query_start: int | None = 0, query_end: int | None = 0, strand: str | None = None,
                 target_name: str | None = None, target_length: int | None = 0, target_start: int | None = 0,
                 target_end: int | None = 0, matching_bases: int | None = 0, num_bases: int | None = 0,
                 percent_identity: float | None = 0, percent_query_coverage: float | None = 0,
                 percent_target_coverage: float | None = None, cigar: str | None = None,
                 alignment_score: int | None = 0, mapping_quality: int | None = 0, cs_tag: str | None = None
                 # query_sequence: Seq | None = None, target_sequence: Seq | None = None
                 ):
        self.query_name = query_name or ''
        self.query_length = query_length
        self.query_start = query_start
        self.query_end = query_end
        self.strand = strand or ''
        self.target_name = target_name or ''
        self.target_length = target_length
        self.target_start = target_start
        self.target_end = target_end
        self.matching_bases = matching_bases
        self.num_bases = num_bases
        self.percent_identity = percent_identity
        self.percent_query_coverage = percent_query_coverage
        self.percent_target_coverage = percent_target_coverage
        self.cigar = cigar or ''
        self.cs_tag = cs_tag or ''
        self.alignment_score = alignment_score
        self.mapping_quality = mapping_quality
        # self.query_sequence = query_sequence or Seq("")
        # self.target_sequence = target_sequence or Seq("")

    # @classmethod
    # def from_m8_line(cls, line: str | bytes):
    #     "qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore"
    #     if not line:
    #         AlignmentError("Empty line")
    #     if len(line := decode_line(line).split('\t')) < 11:
    #         AlignmentError(f"Line has < 12 columns: {line}")
    #     self = cls(
    #         query_name=line[0], target_name=line[1], percent_identity=float(line[2]), num_bases=int(line[3]),
    #         query_start=int(line[6]), query_end=int(line[7]), target_start=int(line[8]), target_end=int(line[9]),
    #         alignment_score=int(line[11])
    #     )
    #     self.matching_bases = int(self.num_bases * self.percent_identity / 100)
    #     if not self.target_start <= self.target_end:
    #         raise AlignmentError(f"{self}: target start > end")
    #     if not self.query_start <= self.query_end:
    #         raise AlignmentError(f"{self}: query start > end")
    #     return self

    @classmethod
    def from_paf_line(cls, line: str | bytes):
        if not line:
            AlignmentError("Empty line")
        if len(line := decode_line(line).split('\t')) < 12:
            AlignmentError(f"Line has < 12 columns: {line}")
        # Parse standard fields
        self = cls(
            query_name=line[0], query_length=int(line[1]), query_start=int(line[2]), query_end=int(line[3]),
            strand=line[4], target_name=line[5], target_length=int(line[6]), target_start=int(line[7]),
            target_end=int(line[8]), matching_bases=int(line[9]), num_bases=int(line[10]), mapping_quality=int(line[11])
        )
        # Calculate percent identity and coverage
        self.percent_identity = 100.0 * self.matching_bases / self.num_bases
        self.percent_query_coverage = 100.0 * (self.query_end - self.query_start) / self.query_length
        self.percent_target_coverage = 100.0 * (self.target_end - self.target_start) / self.target_length

        for tag in line[12:]:
            if tag.startswith('cg:Z:'):
                self.cigar = tag[5:]
            if tag.startswith('AS:i:'):
                self.alignment_score = int(tag[5:])
            if tag.startswith('cs:Z:'):
                self.cs_tag = tag[5:]
        return self

    def __repr__(self):
        return (f'{self.query_name}:{self.query_start}-{self.query_end} '
                f'{self.target_name}:{self.target_start}-{self.target_end}')

    def __len__(self):
        return self.num_bases

    def add_query_sequence(self, query_sequence: Seq):
        if self.query_length and len(query_sequence) != self.query_length:
            raise AlignmentError(f"Query sequence length ({len(query_sequence)}) does not match query length "
                                 f"({self.query_length})")
        self.query_sequence = query_sequence

    def add_target_sequence(self, target_sequence: Seq):
        if self.target_length and len(target_sequence) != self.target_length:
            raise AlignmentError(f"Target sequence length ({len(target_sequence)}) does not match target length "
                                 f"({self.target_length})")
        self.target_sequence = target_sequence

    def target_overlap(self, other: Alignment, allowed_overlap: int = 0) -> int:
        """
        Tests whether this alignment overlaps with the other alignment in the query sequence. A bit
        of overlap can be allowed using the allowed_overlap parameter.
        """
        if self is other:
            return 0
        if self.target_name != other.target_name:
            return 0
        # With Minimap2, query_start is always < query_end, so we can skip the sort
        return range_overlap((self.target_start - allowed_overlap, self.target_end + allowed_overlap),
                             (other.target_start, other.target_end), skip_sort=True)

    def query_overlap(self, other: Alignment, allowed_overlap: int = 0) -> int:
        """
        Tests whether this alignment overlaps with the other alignment in the target sequence. A bit
        of overlap can be allowed using the allowed_overlap parameter.
        """
        if self is other:
            return 0
        if self.query_name != other.query_name:
            return 0
        # With Minimap2, target_start is always < target_end, so we can skip the sort
        return range_overlap((self.query_start - allowed_overlap, self.query_end + allowed_overlap),
                             (other.query_start, other.query_end), skip_sort=True)

    def conflicts(self, other: Alignment, overlap_fraction: float = 0.5):
        """
        Returns whether this hit conflicts with the other hit on the target sequence.
        A conflict is defined as the hits overlapping by 50% or more of the shortest hit's length.
        A hit is not considered to conflict with itself.
        """
        if self is other:
            return False
        if self.target_name != other.target_name:
            return False
        return self.target_overlap(other) / min(self.num_bases, other.num_bases) > overlap_fraction

    # def construct_target_from_query(self):
    #     """
    #     Constructs the target sequence from the cs tag and the query sequence.
    #     """
    #     if not self.cs_tag:
    #         raise AlignmentError(f"{self}: No cs tag found")
    #     if not self.query_sequence or len(self.query_sequence) == 0:
    #         raise AlignmentError(f"{self}: No query sequence found")
    #
    #     target_sequence = ''
    #     query_sequence = self.query_sequence[self.query_start:self.query_end]  # extract the aligned part of the query
    #     # reverse complement if necessary
    #     query_sequence = query_sequence.reverse_complement() if self.strand == '-' else query_sequence
    #     seq_loc = 0
    #     for operation in _CS_OPS_REGEX.finditer(self.cs_tag):
    #         if operation.group().startswith(":"):  # identical bases
    #             op_len = int(operation.group()[1:])  # number of identical bases
    #             target_sequence += query_sequence[seq_loc:seq_loc + op_len]
    #             seq_loc += op_len
    #         elif operation.group().startswith("*"):  # substitution
    #             substitution = operation.group()[1:].upper()
    #             expected_base = query_sequence[seq_loc]
    #             if expected_base in "MRWSYKVHDBN":  # ambiguous base
    #                 expected_base = "N"  # Some K-loci have ambiguous bases (KL50, KL50-1, KL29, KL26, KL15-D1)
    #             if expected_base != substitution[1]:
    #                 raise AlignmentError(f"{self}: Substitution {substitution[1]} does not match base {seq_loc} "
    #                                      f"in query: {expected_base}")
    #             target_sequence += substitution[0]
    #             seq_loc += 1
    #         elif operation.group().startswith("+"):  # insertion in query (deletion in target)
    #             insertion = operation.group()[1:]
    #             seq_loc += len(insertion)
    #             target_sequence += insertion
    #         elif operation.group().startswith("-"):  # deletion in query (insertion in target)
    #             deletion = operation.group()[1:]
    #             target_sequence += deletion.upper()
    #         else:
    #             raise AlignmentError(f"{self}: Unknown cigar operation: {operation.group()}")
    #     if not target_sequence:
    #         raise AlignmentError(f"{self}: Target sequence could not be constructed from cs tag and query sequence")
    #     self.target_sequence = Seq(target_sequence).reverse_complement() if self.strand == '-' else Seq(target_sequence)


#
# class SamFileError(Exception):
#     pass
#
#
# class SamFile:
#     def __init__(self, path: Path):
#         if not path.exists():
#             raise FileNotFoundError(f'Alignment file {path} not found')
#         if path.stat().st_size == 0:
#             raise SamFileError(f'Alignment file {path} is empty')
#         self.path = path
#         self.open_func = open
#         self.open_mode = 'rt'
#
#     def __repr__(self):
#         return f'{self.path.stem}'
#
#     def __len__(self):
#         return self.path.stat().st_size
#
#     def segments(self):
#         """
#         Opens the SamFile and iterates over the alignments with a generator, yields a SamLine per alignment
#         :return: SamLine
#         """
#         with self.open_func(self.path, self.open_mode) as f:
#             for line in f:
#                 if line.startswith('@'):
#                     continue
#                 yield SamLine(line)

#
# class SamLineError(Exception):
#     pass


# class SamLine:
#     """
#     Parses and stores a line from a SAM file into a class object
#     """
#
#     def __init__(self, line: str | bytes):
#         if not line:
#             raise SamLineError("Empty line passed to SamLine")
#         if line.startswith('@'):
#             raise SamLineError("SAM header line passed to SamLine")
#         line = line.decode().strip().split('\t') if isinstance(line, bytes) else line.strip().split('\t')
#         if len(line) < 11:
#             raise SamLineError("SAM line has less than 11 columns")
#         self.qname = line[0]
#         self.flag = int(line[1])
#         self.rname = line[2]
#         self.pos = int(line[3])
#         self.mapq = int(line[4])
#         self.cigar = line[5]
#         self.rnext = line[6]
#         self.pnext = int(line[7])
#         self.tlen = int(line[8])
#         self.seq = line[9]
#         self.qual = line[10]
#         self.tags = line[11:]
#         self.is_paired = False
#         self.is_proper_pair = False
#         self.is_unmapped = False
#         self.mate_is_unmapped = False
#         self.is_reverse = False
#         self.mate_is_reverse = False
#         self.is_read1 = False
#         self.is_read2 = False
#         self.is_secondary = False
#         self.is_qcfail = False
#         self.is_duplicate = False
#         self.is_supplementary = False
#         # Parse the flag into attributes
#         for attr, byte in zip(
#                 ["is_paired", "is_proper_pair", "is_unmapped", "mate_is_unmapped", "is_reverse", "mate_is_reverse",
#                  "is_read1", "is_read2", "is_secondary", "is_qcfail", "is_duplicate", "is_supplementary"],
#                 bin(self.flag)[2:].zfill(11)[::-1]  # reverse the binary string and pad with zeros to 11 digits
#         ):
#             setattr(self, attr, bool(int(byte)))
#
#     def get_map_regions(self):
#         return _MAP_REGIONS_REGEX.findall(self.cigar)
#
#     def get_fastq(self):
#         return f"@{self.qname}\n" \
#                f"{self.seq}\n" \
#                f"+{self.qname} {self.rname} {self.flag} {self.cigar}\n" \
#                f"{self.qual}\n"

# def extract_soft_clips(self, min_size, max_size) -> tuple:
#     """
#     Extracts soft-clipped sequence from a SAM line and returns FASTQ-formatted records for the left and right
#     """
#     if 'S' not in self.cigar:
#         return None, None
#
#     map_regions = self.get_map_regions()
#     left, right = None, None
#
#     if map_regions[0][-1] == 'S':
#         num_soft_clipped = int(map_regions[0][0])
#         if min_size <= num_soft_clipped <= max_size:
#             if self.is_reverse:
#                 left = f"@{self.qname}\n" \
#                        f"{reverse_complement(self.seq[:num_soft_clipped])}\n" \
#                        f"+{self.qname} {self.rname} {self.flag} {self.cigar}\n" \
#                        f"{self.qual[:num_soft_clipped][::-1]}\n"
#             else:
#                 left = f"@{self.qname}\n" \
#                        f"{self.seq[:num_soft_clipped]}\n" \
#                        f"+{self.qname} {self.rname} {self.flag} {self.cigar}\n" \
#                        f"{self.qual[:num_soft_clipped]}\n"
#
#     if map_regions[-1][-1] == 'S':
#         num_soft_clipped = int(map_regions[-1][0])
#         if min_size <= num_soft_clipped <= max_size:
#             if self.is_reverse:
#                 right = f"@{self.qname}\n" \
#                         f"{reverse_complement(self.seq[-num_soft_clipped:])}\n" \
#                         f"+{self.qname} {self.rname} {self.flag} {self.cigar}\n" \
#                         f"{self.qual[-num_soft_clipped:][::-1]}\n"
#             else:
#                 right = f"@{self.qname}\n" \
#                         f"{self.seq[-num_soft_clipped:]}\n" \
#                         f"+{self.qname} {self.rname} {self.flag} {self.cigar}\n" \
#                         f"{self.qual[-num_soft_clipped:]}\n"
#     return left, right


# Functions --------------------------------------------------------------------
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


# def write_alignment_flanks(sam_file: SamFile, left_fastq: Path, right_fastq: Path, min_clip, max_clip):
#     with open(left_fastq, 'w') as left, open(right_fastq, 'w') as right:
#         for segment in sam_file.segments():
#             if segment.is_unmapped:  # equivalent to -f 4
#                 if not segment.mate_is_unmapped and not segment.mate_is_reverse:  # equivalent to -F 40
#                     left.write(segment.get_fastq())
#                 elif segment.mate_is_reverse:  # equivalent to -f 32 (total of 36)
#                     right.write(segment.get_fastq())
#             else:
#                 clip_left, clip_right = segment.extract_soft_clips(min_clip, max_clip)
#                 if clip_left:
#                     left.write(clip_left)
#                 if clip_right:
#                     right.write(clip_right)
#     return left_fastq, right_fastq


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


def get_consensus_sequence(alignments: Iterable[Alignment], query: bool = True):
    """
    """
    # dict to store the pileup of query bases and mapping quality at each position relative to the target sequence
    pileup = {}
    current_target = None
    for alignment in sorted(alignments, key=attrgetter('target_start')):
        if not alignment.query_sequence:
            raise AlignmentError(f"{alignment}: No query sequence found")
        if not alignment.mapping_quality:
            raise AlignmentError(f"{alignment}: No mapping quality found")
        if current_target and alignment.target_name != current_target:
            raise AlignmentError(f"{alignment}: Target name does not match current target")

        # Now iterate over each base in the query sequence and add it to the pileup
        for i, base in enumerate(alignment.query_sequence):
            target_pos = alignment.target_start + i
            pileup.setdefault(target_pos, {}).setdefault(base, 0)
            pileup[target_pos][base] += alignment.mapping_quality
    return pileup

# def coverage(alignments: Iterable[Alignment]
#              ) -> Generator[tuple[str, int, int, int, int, float, float, float, float], None, None]:
#     """
#     Efficient implementation of Samtools coverage that computes the coverage at each position or region
#     Coverage is defined as the percentage of positions within each bin with at least one base aligned against it.
#     Returns the following data:
#         Reference name / chromosome
#         Start position
#         End position
#         Number reads aligned to the region (after filtering)
#         Number of covered bases with depth >= 1
#         Percentage of covered bases [0..100]
#         Mean depth of coverage
#         Mean baseQ in covered region
#         Mean mapQ of selected reads
#     """
#     for target_name, target_alignments in groupby(alignments, attrgetter('target_name')):
#         target_alignments = list(target_alignments)
#         target_length = target_alignments[0].target_length
#         target_alignments = sorted(target_alignments, key=attrgetter('target_start'))
#         target_coverage = [0] * target_length
#         for alignment in target_alignments:
#             for i in range(alignment.target_start, alignment.target_end):
#                 target_coverage[i] += 1
#         yield (
#             target_name, 0, target_length, len(target_alignments), sum(i > 0 for i in target_coverage),
#             100.0 * sum(i > 0 for i in target_coverage) / target_length,
#             sum(target_coverage) / target_length,
#             sum(i for i in target_coverage if i > 0) / sum(i > 0 for i in target_coverage),
#             sum(i for i in target_coverage if i > 0) / len(target_alignments)
#         )
#
