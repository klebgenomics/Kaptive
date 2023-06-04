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
from typing import Union, List
from pathlib import Path

from .log import warning, quit_with_error, log
from .misc import run_command
from .regions import Region


class BedRecordError(Exception):
    """
    Exception for BedRecord
    """
    pass


class BedRecord:
    """
    Class to hold a BED record
    """

    def __init__(self, line: Union[str, bytes]):
        if not line:
            raise BedRecordError("Empty line passed to BedRecord")
        self.line = line.decode().strip() if isinstance(line, bytes) else line.strip()
        if len(line := self.line.split('\t')) < 4:
            raise BedRecordError(f"Invalid BED line: {self.line}")
        self.chrom = line[0]
        self.start = int(line[1])
        self.end = int(line[2])
        self.name = line[3]
        self.score = line[4] if len(line) > 4 else None
        self.strand = line[5] if len(line) > 5 else None
        self.thick_start = int(line[6]) if len(line) > 6 else None
        self.thick_end = int(line[7]) if len(line) > 7 else None
        self.item_rgb = line[8] if len(line) > 8 else None
        self.block_count = float(line[9]) if len(line) > 9 else None
        self.block_sizes = [x for x in line[10].split(',')] if len(line) > 10 else None
        self.block_starts = [x for x in line[11].split(',')] if len(line) > 11 else None
        self.region = Region(self.chrom, self.start, self.end, self.strand)

    def __repr__(self):
        return f"{self.chrom}:{self.start}-{self.end} {self.name}"

    def __str__(self):
        return self.line

    def __len__(self):
        return self.end - self.start

    def __eq__(self, other):
        return self.chrom == other.chrom and self.start == other.start and self.end == other.end


def bedtools_coverage(a: Union[Path, List[BedRecord], str], b: Union[Path, List[BedRecord], str],
                      extra_args: str = '', verbose: bool = False):
    """
    Run bedtools coverage on a sam file wtih a bed file and return a list of BedRecords
    :param extra_args:
    :param verbose: Print command
    :param a: Path to bed file
    :param b: Path to sam file
    :return: List of BedRecords
    """
    pipe = None
    command = f'bedtools coverage '
    command += f'{extra_args} ' if extra_args else ''

    if not isinstance(a, Path) and not isinstance(b, Path):
        quit_with_error(f"A and B cannot both be stdin")
    elif isinstance(a, Path):
        command += f"-a {a} -b -"
        # String representation of BedRecords is the same as the original line, so we can use str(record)
        pipe = '\n'.join(str(j) for j in b) if isinstance(b, list) else b
    elif isinstance(b, Path):
        command += f"-a - -b {b}"
        pipe = '\n'.join(str(j) for j in a) if isinstance(a, list) else a
    else:
        command += f"-a {a} -b {b}"
    if verbose:
        log(command)
    out = run_command(command, pipe=pipe) if pipe else run_command(command)
    if not out:
        warning(f"bedtools coverage returned no results")
        return []
    else:
        # return [BedRecord(line) for line in out.splitlines()]
        return out.splitlines()


def bedtools_genomecov(i: Union[Path, List[BedRecord], str], g: Union[Path, List[BedRecord], str],
                       extra_args: str = '', verbose: bool = False, bam: bool = False) -> List[BedRecord]:
    """
    Run bedtools coverage on a sam file wtih a bed file and return a list of BedRecords
    :param g:
    :param bam:
    :param extra_args:
    :param i:
    :param verbose: Print command
    :return: List of BedRecords
    """
    pipe = None
    command = f'bedtools genomecov '
    command += f'{extra_args} ' if extra_args else ''
    icmd = 'ibam' if bam else 'i'

    if not isinstance(i, Path) and not isinstance(g, Path):
        quit_with_error(f"i and g cannot both be stdin")
    elif not isinstance(g, Path):
        command += f"-{icmd} {i} -g -"
        pipe = '\n'.join(str(j) for j in g) if isinstance(g, list) else g
    elif not isinstance(i, Path):
        command += f"-{icmd} - -g {g}"
        pipe = '\n'.join(str(j) for j in i) if isinstance(i, list) else i
    else:
        command += f"-{icmd} {i} -g {g}"
    if verbose:
        log(command)
    out = run_command(command, pipe=pipe) if pipe else run_command(command)
    if not out:
        warning(f"bedtools genomecov returned no results")
        return []
    else:
        return [BedRecord(line) for line in out.splitlines()]


def bedtools_intersect(a: Union[Path, List[BedRecord], str], b: Union[Path, List[BedRecord], str],
                       extra_args: str = '', verbose: bool = False) -> List[BedRecord]:
    """
    Runs `bedtools intersect` on two bed files/records and returns a list of BedRecords
    :return: List of BedRecords
    """
    pipe = None
    command = f'bedtools intersect '
    command += f'{extra_args} ' if extra_args else ''
    if not isinstance(a, Path) and not isinstance(b, Path):
        quit_with_error(f"A and B cannot both be stdin")
    elif not isinstance(a, Path):
        command += f"-a - -b {b}"
        pipe = '\n'.join(str(j) for j in a) if isinstance(a, list) else a
    elif not isinstance(b, Path):
        command += f"-a {a} -b -"
        pipe = '\n'.join(str(j) for j in b) if isinstance(b, list) else b
    else:
        command += f"-a {a} -b {b}"
    if verbose:
        log(command)
    out = run_command(command, pipe=pipe) if pipe else run_command(command)
    if not out:
        warning(f"bedtools intersect returned no results")
        return []
    else:
        return [line for line in out.splitlines()]


def bedtools_merge(i: Union[Path, List[BedRecord], str], extra_args: str = '', verbose: bool = False
                   ) -> List[BedRecord]:
    """
    Runs `bedtools merge` on two bed files/records and returns a list of BedRecords
    :return: List of BedRecords
    """
    command = f'bedtools merge '
    command += f'{extra_args} ' if extra_args else ''
    command += f"-i {i}" if isinstance(i, Path) else f"-i -"
    if verbose:
        log(command)
    out = run_command(command, pipe='\n'.join(str(j) for j in i)) if not isinstance(i, Path) else run_command(command)
    if not out:
        warning(f"bedtools intersect returned no results")
        return []
    else:
        return [line for line in out.splitlines()]


def sort_bed_records(records: List[BedRecord]) -> List[BedRecord]:
    """
    Sorts a list of BedRecords by chrom, start, end, strand
    :param records: List of BedRecords
    :return: List of BedRecords
    """
    return sorted(records, key=lambda x: (x.chrom, x.start, x.end, x.strand))

