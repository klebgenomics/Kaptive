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

import argparse
import os
import sys
from pathlib import Path
from typing import List, Tuple
import gzip
import subprocess
from textwrap import wrap

from Bio.Seq import Seq

from .log import log, bold_cyan, quit_with_error, warning


def check_file(parser: 'argparse.ArgumentParser', file: Path) -> Path:
    if not file.is_file():
        parser.error(f'{file} does not exist')
    elif file.stat().st_size == 0:
        parser.error(f'{file} is empty')
    else:
        return file.absolute()


def check_dir(parser: 'argparse.ArgumentParser', dirpath) -> Path:
    try:
        dirpath.mkdir(parents=True, exist_ok=True)
        return dirpath
    except Exception as e:
        parser.error(str(e))


def check_programs(progs: List[str], verbose: bool = False):
    """Check if programs are installed and executable"""
    bins = {  # Adapted from: https://unix.stackexchange.com/a/261971/375975
        f: f'{p}/{f}' for p in filter(
            os.path.isdir, os.environ["PATH"].split(os.path.pathsep)
        ) for f in os.listdir(p) if os.access(f'{p}/{f}', os.X_OK)
    }
    for program in progs:
        if program in bins.keys():
            if verbose:
                log(f'{program}: in path {bins[program]}')
        else:
            quit_with_error(f'{program} not found')


def check_python_version():
    if sys.version_info.major < 3 or sys.version_info.minor < 9:
        quit_with_error('kaptive-mapper requires Python 3.9 or later')


def get_gc_content(seq: str, length: int) -> float:
    seq = seq.upper()
    return (seq.count('G') + seq.count('C')) / length


def fasta_wrap(fasta_string: str, width: int = 80) -> str:
    return '\n'.join(wrap(fasta_string, width))


def run_command(command: str, pipe: str = '', cmd_split=' ', quiet: bool = False):
    """Run a command and return the output"""
    if pipe:
        with subprocess.Popen(command.split(cmd_split), stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                              stdin=subprocess.PIPE) as child:
            out, err = child.communicate(input=pipe.encode())
    else:
        with subprocess.Popen(command.split(cmd_split), stdout=subprocess.PIPE, stderr=subprocess.PIPE) as child:
            out, err = child.communicate()

    if err and not quiet:
        warning(err.decode())

    return '' if not out else out.decode()


def write_text_file(path, text: str, trailing_newline: bool = False, compressed=False):
    _open = gzip.open if compressed else open
    with _open(path, 'wt') as f:
        f.write(text)
        if trailing_newline:
            f.write('\n')
    if path.is_file() and path.stat().st_size > 0:
        log(f'Wrote {path}')
        return path
    elif path.is_file() and path.stat().st_size == 0:
        warning(f'Wrote empty file {path}')
        return None
    else:
        quit_with_error(f'Could not write {path}')
        return None


def get_ascii_art():
    return bold_cyan(fr'''
              
''')


REV_COMP_DICT = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G', 'a': 't', 't': 'a', 'g': 'c', 'c': 'g',
                 'R': 'Y', 'Y': 'R', 'S': 'S', 'W': 'W', 'K': 'M', 'M': 'K', 'B': 'V', 'V': 'B',
                 'D': 'H', 'H': 'D', 'N': 'N', 'r': 'y', 'y': 'r', 's': 's', 'w': 'w', 'k': 'm',
                 'm': 'k', 'b': 'v', 'v': 'b', 'd': 'h', 'h': 'd', 'n': 'n', '.': '.', '-': '-',
                 '?': '?'}


def complement_base(base):
    try:
        return REV_COMP_DICT[base]
    except KeyError:
        return 'N'


def reverse_complement(seq):
    return ''.join([complement_base(x) for x in seq][::-1])


def format_str(string, format_number):
    return f"\033[{format_number}m{string}\033[0m"


def get_executable_binaries():
    """Adapted from: https://unix.stackexchange.com/a/261971/375975"""
    paths = os.environ["PATH"].split(os.path.pathsep)
    executables = {}
    for path in filter(os.path.isdir, paths):
        for file_ in os.listdir(path):
            if os.access(os.path.join(path, file_), os.X_OK):
                executables[file_] = f'{path}/{file_}'
    return executables


def seq_to_dict(seq: str) -> dict:
    return {r[0]: ''.join(r[1:]) for record in seq.split('>')[1:] if (r := record.splitlines())}


def dict_to_seq(seq: dict) -> str:
    return '\n'.join([f'>{k}\n{v}' for k, v in seq.items()])


def load_fasta(fasta_file: Path) -> dict:
    """Load a FASTA file into a dictionary"""
    _open = gzip.open if fasta_file.suffix == '.gz' else open
    try:
        with _open(fasta_file, 'rt') as fasta:
            return seq_to_dict(fasta.read())
    except Exception as e:
        log(format_str(f'Error loading {fasta_file.name}', '91'))
        log(str(e))
        return None


def quick_translate(seq: str, to_stop=False, gap='-'):
    """ Quick DNA translation using bacterial translation-table
    'to_stop=True' tests for truncation as a consequence of SNP """
    if not isinstance(seq, Seq):
        if gap in seq:  # Remove blast deletion symbol '-' to avoid error and test effect of deletion
            seq = seq.replace(gap, '')
        seq = Seq(seq)
    return str(seq.translate(table=11, to_stop=to_stop, stop_symbol='*'))


def find_files_with_suffixes(prefix: Path, suffixes: List[str], min_size: int = 1) -> List[Path]:
    """
    Find files with given suffixes for the given prefix, useful for finding blast database or bwa index files
    :param suffixes: List of existing Paths with suffixes
    :param prefix: Prefix as Path
    :return: List[Path]
    """
    return [p for i in suffixes if
            (p := prefix.with_suffix(prefix.suffix + i)).exists() and p.stat().st_size >= min_size]


def merge_ranges(ranges: List[Tuple[int, int]], tolerance=0) -> List[Tuple[int, int]]:
    """
    Merge overlapping ranges
    :param ranges: List of tuples of start and end positions
    :param tolerance: The number of bases to allow between alignments to be considered continuous
    :return: List of merged ranges
    """
    if not ranges:
        return []

    # Sort the ranges by their start values
    sorted_ranges = sorted(ranges, key=lambda x: x[0])

    merged_ranges = []
    current_range = sorted_ranges[0]

    for start, end in sorted_ranges[1:]:
        adjusted_start = start - tolerance
        adjusted_end = end + tolerance

        if adjusted_start <= current_range[1] + tolerance:
            # Ranges overlap, merge them
            current_range = (current_range[0], max(current_range[1], adjusted_end))
        else:
            # No overlap, add the current range to the merged list and start a new range
            merged_ranges.append(current_range)
            current_range = (start, end)

    # Add the last range to the merged list
    merged_ranges.append(current_range)

    return merged_ranges


def range_overlap(range1: Tuple[int, int], range2: Tuple[int, int]) -> int:
    """
    Returns the overlap between two ranges
    :param range1: Tuple of start and end positions
    :param range2: Tuple of start and end positions
    :return: Integer of overlap
    """
    start1, end1 = range1
    start2, end2 = range2
    overlap_start = max(start1, start2)
    overlap_end = min(end1, end2)
    overlap = max(0, overlap_end - overlap_start)
    return overlap


def get_compression_type(filename: Path) -> str:
    """
    Attempts to guess the compression (if any) on a file using the first few bytes.
    http://stackoverflow.com/questions/13044562
    """
    magic_dict = {'gz': (b'\x1f', b'\x8b', b'\x08'),
                  'bz2': (b'\x42', b'\x5a', b'\x68'),
                  'zip': (b'\x50', b'\x4b', b'\x03', b'\x04')}
    max_len = max(len(x) for x in magic_dict)

    with open(filename, 'rb') as unknown_file:
        file_start = unknown_file.read(max_len)

    compression_type = 'plain'
    for filetype, magic_bytes in magic_dict.items():
        if file_start.startswith(magic_bytes):
            compression_type = filetype
    if compression_type == 'bz2':
        quit_with_error('cannot use bzip2 format - use gzip instead')
    if compression_type == 'zip':
        quit_with_error('cannot use zip format - use gzip instead')
    return compression_type


def check_out(parser: 'argparse.ArgumentParser', out_path) -> Path:
    if out_path.endswith('/'):
        try:
            Path(out_path).mkdir(parents=True, exist_ok=True)
            return Path(os.path.join(out_path, 'kaptive_results')).absolute()
        except Exception as e:
            parser.error(str(e))
    else:
        return Path(out_path).absolute()


def float_to_str(float_in):
    """
    This function converts a float to a string in a special manner: if the float is an integer,
    the resulting string has no decimal point. Otherwise, one decimal point is used.
    """
    if float_in == int(float_in):
        return str(int(float_in))
    else:
        return '%.1f' % float_in


def good_start_and_end(start, end, length, allowed_margin):
    """
    Checks whether the given start and end coordinates are within the accepted margin of error.
    """
    good_start = start <= allowed_margin
    good_end = end >= length - allowed_margin
    start_before_end = start < end
    return good_start and good_end and start_before_end
