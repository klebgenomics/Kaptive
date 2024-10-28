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

import os
import sys
from zlib import decompress as gz_decompress
from gzip import open as gz_open
from bz2 import (decompress as bz2_decompress, open as bz2_open)
from lzma import (decompress as xz_decompress, open as xz_open)
from typing import Generator, TextIO, Any, BinaryIO
from operator import itemgetter

from kaptive.log import log, quit_with_error, bold_cyan, warning

# Constants -----------------------------------------------------------------------------------------------------------
_MAX_CPUS = 32
_MAGIC_BYTES = {b'\x1f\x8b': 'gz', b'\x42\x5a': 'bz2', b'\xfd7zXZ\x00': 'xz'}
_OPEN = {'gz': gz_open, 'bz2': bz2_open, 'xz': xz_open}
_DECOMPRESS = {'gz': gz_decompress, 'bz2': bz2_decompress, 'xz': xz_decompress}
_MIN_N_BYTES = max(len(i) for i in _MAGIC_BYTES)  # Minimum number of bytes to read in a file to guess the compression)
_LOGO = r"""  _  __    _    ____ _____ _____     _______ 
 | |/ /   / \  |  _ \_   _|_ _\ \   / / ____|
 | ' /   / _ \ | |_) || |  | | \ \ / /|  _|  
 | . \  / ___ \|  __/ | |  | |  \ V / | |___ 
 |_|\_\/_/   \_\_|    |_| |___|  \_/  |_____|                                   
"""


# Functions -----------------------------------------------------------------------------------------------------------
def check_programs(progs: list[str], verbose: bool = False):
    """Check if programs are installed and executable"""
    bins = {  # Adapted from: https://unix.stackexchange.com/a/261971/375975
        binary: x for path in filter(
            os.path.isdir, os.environ["PATH"].split(os.path.pathsep)
        ) for binary in os.listdir(path) if os.access((x := os.path.join(path, binary)), os.X_OK)
    }
    for program in progs:
        if program in bins:
            log(f'{program}: {bins[program]}', verbose=verbose)
        else:
            quit_with_error(f'{program} not found')


def check_file(file: str | os.PathLike, panic: bool = False) -> os.PathLike | None:
    """Checks a file exists and is non-empty and returns the absolute path"""
    func = quit_with_error if panic else warning
    if not os.path.exists(file):
        return func(f'{file} does not exist')
    if not os.path.isfile(file):
        return func(f'{file} is not a file')
    elif not os.path.getsize(file):
        return func(f'{file} is empty')
    else:
        return os.path.abspath(file)


def check_cpus(cpus: Any = None, max_cpus: int = _MAX_CPUS, verbose: bool = False) -> int:
    avail_cpus = os.cpu_count() or max_cpus
    if isinstance(cpus, str):
        cpus = int(cpus) if cpus.isdigit() else avail_cpus
    elif isinstance(cpus, float):
        cpus = int(cpus)
    else:
        cpus = avail_cpus
    cpus = min(cpus, avail_cpus, max_cpus)
    log(f'Using {cpus=}', verbose)
    return cpus


def check_out(path: str | os.PathLike, mode: str = "at", exist_ok: bool = True) -> os.PathLike | TextIO:
    """
    Check if the user wants to create/append a file or directory.
    If it looks like/is already a file (has an extension), return the file object.
    If it looks like/is already a directory, return the directory path.
    """
    if path == '-':  # If the path is '-', return stdout
        return sys.stdout
    if os.path.splitext(path)[1]:  # If the path has an extension, it's probably a file
        try:
            return open(path, mode)  # Open the file
        except Exception as e:
            quit_with_error(f'Could not open {path}: {e}')
    if not os.path.exists(path):  # Assume directory
        try:
            os.makedirs(path, exist_ok=exist_ok)  # Create the directory if it doesn't exist
        except Exception as e:
            quit_with_error(f'Could not create {path}: {e}')
    return path


def opener(file: str | os.PathLike, verbose: bool = False, *args, **kwargs) -> TextIO | BinaryIO:
    """
    Opens a file with the appropriate open function based on the magic bytes at the beginning of the data
    :param file: File to open
    :param verbose: Print log messages to stderr
    :return: File handle
    """
    try:
        file = check_file(file)
    except FileNotFoundError as e:
        raise e
    basename = os.path.basename(file)
    with open(file, 'rb') as f:  # Open the file to read bytes
        first_bytes = f.read(_MIN_N_BYTES)  # Get the bytes necessary to guess the compression type
    for magic, compression in _MAGIC_BYTES.items():
        if first_bytes.startswith(magic):
            log(f"Assuming {basename} is compressed with {compression}", verbose=verbose)
            try:
                return _OPEN[compression](file, *args, **kwargs)
            except Exception as e:
                return warning(f"Error opening {basename} with {compression}; {first_bytes=}\n{e}")
    log(f"Assuming {basename} is uncompressed", verbose=verbose)
    return open(file, *args, **kwargs)


def get_logo(message: str, width: int = 43) -> str:  # 43 is the width of the logo
    return bold_cyan(f'{_LOGO}\n{message.center(width)}')


def merge_ranges(ranges: list[tuple[int | float, int | float]], tolerance: int | float = 0, skip_sort: bool = False
                 ) -> Generator[tuple[int | float, int | float], None, None]:
    """
    Merge overlapping ranges
    :param ranges: List of tuples of start and end positions
    :param tolerance: Integer or float of tolerance for merging ranges
    :param skip_sort: Skip sorting the ranges before merging
    :return: List of merged ranges
    """
    if not ranges:
        return None
    if len(ranges) == 1:
        yield ranges[0]
        return None
    current_range = (ranges := ranges if skip_sort else sorted(ranges, key=itemgetter(0)))[0]  # Start with the first range
    for start, end in ranges[1:]:  # Iterate through the ranges
        if start - tolerance <= current_range[1]:  # Overlap, merge the ranges
            current_range = (current_range[0], max(current_range[1], end))
        else:  # No overlap, add the current range to the merged list and start a new range
            yield current_range  # Yield the current range
            current_range = (start, end)   # Start a new range
    yield current_range  # Yield the last range


def range_overlap(range1: tuple[int, int], range2: tuple[int, int], skip_sort: bool = False) -> int:
    """
    Returns the overlap between two ranges
    :param range1: Tuple of start and end positions
    :param range2: Tuple of start and end positions
    :param skip_sort: Skip sorting each range before calculating the overlap
    :return: Integer of overlap
    """
    start1, end1 = range1 if skip_sort else sorted(range1)
    start2, end2 = range2 if skip_sort else sorted(range2)
    overlap_start = max(start1, start2)
    overlap_end = min(end1, end2)
    return max(0, overlap_end - overlap_start)


def str2val(value: str) -> str | bool | int | float:
    """
    Convert a string to a bool, integer or float, falling back to the original string if conversion fails
    """
    if value == 'True':
        return True
    elif value == 'False':
        return False
    try:
        return int(value)
    except ValueError:
        try:
            return float(value)
        except ValueError:
            return value

