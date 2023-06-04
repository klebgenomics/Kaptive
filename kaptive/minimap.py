"""
Copyright 2023 Tom Stanton (tomdstanton@gmail.com)
https://github.com/tomdstanton/kaptive-mapper

This file is part of kaptive-mapper. kaptive-mapper is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by the Free Software Foundation,
either version 3 of the License, or (at your option) any later version. kaptive-mapper is distributed
in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
details. You should have received a copy of the GNU General Public License along with kaptive-mapper.
If not, see <https://www.gnu.org/licenses/>.
"""
import subprocess
from re import compile
from typing import Union, List, Dict, Tuple
from pathlib import Path
import time

from .misc import reverse_complement, run_command, find_files_with_suffixes, merge_ranges, range_overlap
from .log import quit_with_error, log, warning
from .reads import ReadFile

MAP_REGIONS_REGEX = compile(r'([0-9]+)([MIDNSP=X])')


# Classes -------------------------------------------------------------------------------------------------------------
class SamLineError(Exception):
    pass


class SamFileError(Exception):
    pass


class PafLineError(Exception):
    pass


class SamFile:
    def __init__(self, path: Path):
        if not path.exists():
            raise FileNotFoundError(f'Alignment file {path} not found')
        if path.stat().st_size == 0:
            raise SamFileError(f'Alignment file {path} is empty')
        self.path = path
        self.open_func = open
        self.open_mode = 'rt'

    def __repr__(self):
        return f'{self.path.stem}'

    def __len__(self):
        return self.path.stat().st_size

    def segments(self):
        """
        Opens the SamFile and iterates over the alignments with a generator, yields a SamLine per alignment
        :return: SamLine
        """
        with self.open_func(self.path, self.open_mode) as f:
            for line in f:
                if line.startswith('@'):
                    continue
                yield SamLine(line)


class SamLine:
    """
    Parses and stores a line from a SAM file into a class object
    """

    def __init__(self, line: Union[str, bytes]):
        if not line:
            raise SamLineError("Empty line passed to SamLine")
        if line.startswith('@'):
            raise SamLineError("SAM header line passed to SamLine")
        line = line.decode().strip().split('\t') if isinstance(line, bytes) else line.strip().split('\t')
        if len(line) < 11:
            raise SamLineError("SAM line has less than 11 columns")
        self.qname = line[0]
        self.flag = int(line[1])
        self.rname = line[2]
        self.pos = int(line[3])
        self.mapq = int(line[4])
        self.cigar = line[5]
        self.rnext = line[6]
        self.pnext = int(line[7])
        self.tlen = int(line[8])
        self.seq = line[9]
        self.qual = line[10]
        self.tags = line[11:]
        self.is_paired = False
        self.is_proper_pair = False
        self.is_unmapped = False
        self.mate_is_unmapped = False
        self.is_reverse = False
        self.mate_is_reverse = False
        self.is_read1 = False
        self.is_read2 = False
        self.is_secondary = False
        self.is_qcfail = False
        self.is_duplicate = False
        self.is_supplementary = False
        # Parse the flag into attributes
        for attr, byte in zip(
                ["is_paired", "is_proper_pair", "is_unmapped", "mate_is_unmapped", "is_reverse", "mate_is_reverse",
                 "is_read1", "is_read2", "is_secondary", "is_qcfail", "is_duplicate", "is_supplementary"],
                bin(self.flag)[2:].zfill(11)[::-1]  # reverse the binary string and pad with zeros to 11 digits
        ):
            setattr(self, attr, bool(int(byte)))

    def get_map_regions(self):
        return MAP_REGIONS_REGEX.findall(self.cigar)

    def get_fastq(self):
        return f"@{self.qname}\n" \
               f"{self.seq}\n" \
               f"+{self.qname} {self.rname} {self.flag} {self.cigar}\n" \
               f"{self.qual}\n"

    def extract_soft_clips(self, min_size, max_size) -> tuple:
        """
        Extracts soft-clipped sequence from a SAM line and returns FASTQ-formatted records for the left and right
        """
        if 'S' not in self.cigar:
            return None, None

        map_regions = self.get_map_regions()
        left, right = None, None

        if map_regions[0][-1] == 'S':
            num_soft_clipped = int(map_regions[0][0])
            if min_size <= num_soft_clipped <= max_size:
                if self.is_reverse:
                    left = f"@{self.qname}\n" \
                           f"{reverse_complement(self.seq[:num_soft_clipped])}\n" \
                           f"+{self.qname} {self.rname} {self.flag} {self.cigar}\n" \
                           f"{self.qual[:num_soft_clipped][::-1]}\n"
                else:
                    left = f"@{self.qname}\n" \
                           f"{self.seq[:num_soft_clipped]}\n" \
                           f"+{self.qname} {self.rname} {self.flag} {self.cigar}\n" \
                           f"{self.qual[:num_soft_clipped]}\n"

        if map_regions[-1][-1] == 'S':
            num_soft_clipped = int(map_regions[-1][0])
            if min_size <= num_soft_clipped <= max_size:
                if self.is_reverse:
                    right = f"@{self.qname}\n" \
                            f"{reverse_complement(self.seq[-num_soft_clipped:])}\n" \
                            f"+{self.qname} {self.rname} {self.flag} {self.cigar}\n" \
                            f"{self.qual[-num_soft_clipped:][::-1]}\n"
                else:
                    right = f"@{self.qname}\n" \
                            f"{self.seq[-num_soft_clipped:]}\n" \
                            f"+{self.qname} {self.rname} {self.flag} {self.cigar}\n" \
                            f"{self.qual[-num_soft_clipped:]}\n"
        return left, right


class PafLine:
    """
    Object to store a line from a PAF output file
    """

    def __init__(self, line: Union[str, bytes]):
        if not line:
            raise PafLineError("Empty line passed to PafLine")

        self.line = line.decode().strip() if isinstance(line, bytes) else line.strip()

        if len(line := self.line.split('\t')) < 11:
            raise PafLineError(f"PAF line has less than 11 columns: {self.line}")

        self.query_name = line[0]
        self.query_length = int(line[1])
        self.query_start = int(line[2])
        self.query_end = int(line[3])
        self.strand = line[4]

        self.target_name = line[5]
        self.target_length = int(line[6])
        self.target_start = int(line[7])
        self.target_end = int(line[8])

        self.matching_bases = int(line[9])
        self.num_bases = int(line[10])
        self.percent_identity = 100.0 * self.matching_bases / self.num_bases

        self.query_cov = 100.0 * (self.query_end - self.query_start) / self.query_length

        self.cigar, self.alignment_score = None, None
        for part in line:
            if part.startswith('cg:Z:'):
                self.cigar = part[5:]
            if part.startswith('AS:i:'):
                self.alignment_score = int(part[5:])

        if self.query_name[-2:] == '/1':
            self.read = self.query_name[:-2]
            self.read_direction = 'forward'
        elif self.query_name[-2:] == '/2':
            self.read = self.query_name[:-2]
            self.read_direction = 'reverse'
        else:
            self.read = self.query_name
            self.read_direction = None

    def __repr__(self):
        return f'{self.query_name}, {self.target_name}, {self.percent_identity:.2f}%'

    def __str__(self):
        return self.line

    def __len__(self):
        return self.target_start - self.target_end


class Minimap2Result:
    def __init__(self, paf_lines: []):
        self.targets = {}
        self.queries = {}
        self.paf_lines = paf_lines
        self.target_ranges = {}

        for paf_line in paf_lines:
            if paf_line.target_name not in self.targets:
                self.targets[paf_line.target_name] = []
            self.targets[paf_line.target_name].append(paf_line)

            # if paf_line.query_name not in self.queries:
            #     self.queries[paf_line.query_name] = []
            # self.queries[paf_line.query_name].append(paf_line)

            if paf_line.read not in self.queries:
                self.queries[paf_line.read] = []
            self.queries[paf_line.read].append(paf_line)

        for target, alignments in self.targets.items():
            self.target_ranges[target] = target_ranges_covered_by_alignments(alignments)[target]

    def __repr__(self):
        return f'{len(self.paf_lines)} PAF lines'

    def __len__(self):
        return len(self.paf_lines)

    def __getitem__(self, item):
        return self.paf_lines[item]

    def __iter__(self):
        return iter(self.paf_lines)

    def __contains__(self, item):
        return item in self.paf_lines

    def get_best_alignment_per_query(self, min_identity: float = 0.0, min_coverage: float = 0.0):
        best_alignments = {}
        for query, alignments in self.queries.items():
            best_alignment = None
            for alignment in alignments:
                if alignment.percent_identity >= min_identity and alignment.query_cov >= min_coverage:
                    if best_alignment is None or alignment.alignment_score > best_alignment.alignment_score:
                        best_alignment = alignment
            if best_alignment is not None:
                best_alignments[query] = best_alignment
        return best_alignments

    def get_best_alignment_per_reference(self, min_identity: float = 0.0, min_coverage: float = 0.0):
        best_alignments = {}
        for reference, alignments in self.targets.items():
            best_alignment = None
            for alignment in alignments:
                if alignment.percent_identity >= min_identity and alignment.query_cov >= min_coverage:
                    if best_alignment is None or alignment.alignment_score > best_alignment.alignment_score:
                        best_alignment = alignment
            if best_alignment is not None:
                best_alignments[reference] = best_alignment
        return best_alignments


def target_ranges_covered_by_alignments(alignments: List[PafLine], tolerance: int = 0)\
        -> Dict[str, List[Tuple[int, int]]]:
    """
    Get the regions of a target covered by the alignments
    :param tolerance: The number of bases to allow between alignments to be considered continuous
    :param alignments: A list of PafLine objects
    :return:
    """
    target_ranges = {}
    for alignment in alignments:
        if alignment.target_name not in target_ranges:
            target_ranges[alignment.target_name] = [(alignment.target_start, alignment.target_end)]
        else:
            target_ranges[alignment.target_name].append((alignment.target_start, alignment.target_end))
    return {target: merge_ranges(ranges, tolerance) for target, ranges in target_ranges.items()}


# Functions --------------------------------------------------------------------
def minimap2(query: Union[Path, List[ReadFile], str], target: Union[Path, str], output: Union[Path, None] = None,
             preset: str = "map-ont", threads: int = 1, extra_args: str = "",
             verbose: bool = False) -> Union[SamFile, List[PafLine], None]:
    """
    Run minimap2 with the given arguments
    :param verbose: Bool for logging messages for command and execution time
    :param query: Path or stdin
    :param target: Path or stdin
    :param output: Path if SAM AlignmentFile object, None if PAF Alignment objects
    :param preset: Alignment preset
    :param threads: Number of threads
    :param extra_args: Extra arguments to pass to minimap2
    :return: AlignmentFile object if output is a Path, else a list of Alignment objects
    """
    pipe = None
    start = time.time()
    command = f"minimap2 -t {threads} -x {preset} "
    if output:
        command += f"-a -o {output} "
    else:
        command += "-c "
    if extra_args:
        command += f"{extra_args} "
    if isinstance(query, str) and isinstance(target, str):
        quit_with_error(f"Query and target cannot both be stdin")
    elif isinstance(query, str):
        command += f"{target} -"
        pipe = query
    elif isinstance(target, str):
        command += f"- {' '.join(str(j) for j in query) if isinstance(query, list) else str(query)}"
        pipe = target
    else:
        command += f"{target} {' '.join(str(j) for j in query) if isinstance(query, list) else str(query)}"

    if verbose:
        log(command)

    out = run_command(command, quiet=True, pipe=pipe) if pipe else run_command(command, quiet=True)

    if verbose:
        log(f"minimap2 completed in {time.time() - start:.2f} seconds")

    if output:
        if output.exists() and output.stat().st_size > 0:
            return SamFile(output)
        else:
            warning(f"minimap2 failed to align {query if isinstance(query, Path) else '-'} to "
                    f"{target if isinstance(target, Path) else '-'}")
            return None
    else:
        if out:
            return [PafLine(i) for i in out.splitlines()]
        else:
            warning(f"minimap2 failed to align {query if isinstance(query, Path) else '-'} to "
                    f"{target if isinstance(target, Path) else '-'}")
            return []


def create_minimap2_index(reference: Union[Path, List[Path], str], prefix: Union[Path, None] = None, threads: int = 1,
                          verbose: bool = False) -> Path:
    """
    Create a minimap2 index
    :param reference:
    :param prefix:
    :param threads:
    :param verbose:
    :return:
    """
    if not prefix:
        if isinstance(reference, str):
            quit_with_error(f"Prefix must be specified if reference is a string")
        else:
            prefix = reference

    if minimap2_index_exists(prefix):
        warning("Index files for reference already exist")
        return prefix
    else:
        start = time.time()
        command = f'minimap2 -2 -t {threads} -d {prefix} '

        if isinstance(reference, str):
            run_command(command + '-', pipe=reference, quiet=True)
        else:
            run_command(
                command + ' '.join(str(i) for i in reference) if isinstance(reference, list) else str(reference))

        if minimap2_index_exists(prefix):
            if verbose:
                log(f"Built minimap2 index in {time.time() - start:.2f} seconds")
            return prefix
        else:
            quit_with_error("Failed to create minimap2 index")


def samtools_depth(sam: Path, threads: int = 1) -> Dict[str, int]:
    sort_cmd = f"samtools sort -l 0 --threads {threads} {sam}".split()
    depth_cmd = "samtools depth -".split()
    sort_proc = subprocess.Popen(sort_cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    depth_proc = subprocess.Popen(depth_cmd, stdin=sort_proc.stdout, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    out, err = depth_proc.communicate()
    if not out:
        if err:
            warning(err.decode())
        return {}
    else:
        depth_per_reference = {}
        for line in out.decode().splitlines():
            fields = line.strip().split('\t')
            reference = fields[0]
            depth = int(fields[2])

            if reference in depth_per_reference:
                depth_per_reference[reference] += depth
            else:
                depth_per_reference[reference] = depth
        return depth_per_reference


def paftools_sam2paf(sam: Union[str, SamFile], long_cs_tag: bool = False) -> List[PafLine]:
    """
    Convert a sam file to a paf file
    :param long_cs_tag: output the cs tag in the long form
    :param sam: Path to sam file or AlignmentFile object
    :return: List[Alignment]
    """
    command = "paftools.js sam2paf"
    if long_cs_tag:
        command += " -L"
    if isinstance(sam, str):
        proc = subprocess.Popen(command.split(), stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        out, err = proc.communicate(sam.encode())
    else:
        proc = subprocess.Popen(f"{command} {sam.path}".split(), stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        out, err = proc.communicate()
    if not out:
        if err:
            warning(err.decode())
        return []
    else:
        return [PafLine(i) for i in out.decode().splitlines()]


# def paftools_splice2bed(alignment: Union[List[PafLine], SamFile]) -> List[BedRecord]:
#     """
#     Run `paftools.js splice2bed` on a list of PafLines or a SamFile to get a list of BedRecords
#     :param alignment: Either a list of PafLines or a SamFile
#     :return: list of BedRecords
#     """
#     command = 'paftools.js splice2bed'
#     out = run_command(f'{command} {alignment.path}') if isinstance(alignment, SamFile) else \
#         run_command(f'{command} {alignment.path} -', pipe="\n".join(str(i) for i in alignment))
#     if not out:
#         return []
#     else:
#         return [BedRecord(i) for i in out.splitlines()]


def minimap2_index_exists(reference: Union[Path, str]) -> bool:
    """
    Check if minimap2 index exists for the given reference
    :param reference: Path to reference
    :return: bool
    """
    suffixes = [".mmi"]
    return len(find_files_with_suffixes(reference, suffixes)) == len(suffixes)


def get_overlapping_alignments(paf_lines: List[PafLine], start: int, end: int, query: bool = False) -> List[PafLine]:
    """
    Get all alignments that overlap the given coordinates
    :param query: return query coordinates instead of target coordinates
    :param paf_lines: list of PafLines
    :param start: start coordinate (zero-based)
    :param end: end coordinate (zero-based)
    :return: list of PafLines
    """
    # if query:
    #     return [i for i in paf_lines if i.query_start < end and i.query_end > start]
    # else:
    #     return [i for i in paf_lines if i.target_start < end and i.target_end > start]
    if query:
        return [i for i in paf_lines if range_overlap((i.query_start, i.query_end), (start, end))]
    else:
        return [i for i in paf_lines if range_overlap((i.target_start, i.target_end), (start, end))]


def weighted_identity(paf_lines: List[PafLine], start: int, end: int):
    overlapping_bases = 0
    weighted_identity_sum = 0
    for paf_line in paf_lines:
        if paf_line.target_start <= end and paf_line.target_end >= start:
            overlap_start = max(start, paf_line.target_start)
            overlap_end = min(end, paf_line.target_end)
            overlap_length = overlap_end - overlap_start
            if paf_line.percent_identity == 0:
                continue
            weighted_identity_sum += paf_line.percent_identity * overlap_length
            overlapping_bases += overlap_length
    if overlapping_bases > 0:
        return weighted_identity_sum / overlapping_bases
    else:
        return 0

# def paf_to_bed(paf_line: str) -> str:
#     fields = paf_line.strip().split('\t')
#     chrom = fields[0]
#     start = int(fields[2])
#     end = int(fields[3])
#     name = fields[5]
#     strand = fields[4]
#     block_sizes = fields[10].split(',')
#     block_sizes[-1] = str(int(fields[1]) - int(fields[9]))
#     block_sizes = ','.join(block_sizes)
#     block_starts = fields[11].split(',')
#     block_starts[-1] = '0'
#     for i in range(1, len(block_starts)):
#         block_starts[i] = str(int(block_starts[i]) - start)
#     block_starts = ','.join(block_starts)
#     bed_fields = [chrom, str(start), str(end), name, '0', strand, str(start), str(end), '0', '1', block_sizes,
#                   block_starts]
#     bed_line = '\t'.join(bed_fields)
#     return bed_line


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

