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
from pathlib import Path
from gzip import open as gzopen
from bz2 import open as bzopen
from typing import Generator, TextIO

from kaptive.log import log, quit_with_error, bold_cyan, warning

# Constants -----------------------------------------------------------------------------------------------------------
_COMPRESSION_MAGIC = {b'\x1f\x8b': 'gz', b'\x42\x5a': 'bz2', b'\x50\x4b': 'zip', b'\x37\x7a': '7z', b'\x78\x01': 'xz'}
_READ_N_BYTES = max(len(x) for x in _COMPRESSION_MAGIC)
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
        f: Path(f'{p}/{f}') for p in filter(
            os.path.isdir, os.environ["PATH"].split(os.path.pathsep)
        ) for f in os.listdir(p) if os.access(f'{p}/{f}', os.X_OK)
    }
    for program in progs:
        if program in bins:
            log(f'{program}: {bins[program]}', verbose=verbose)
        else:
            quit_with_error(f'{program} not found')


def check_file(path: str | Path) -> Path | None:
    path = Path(path) if isinstance(path, str) else path
    if not path.exists():
        return warning(f'{path} does not exist')
    if not path.is_file():
        return warning(f'{path} is not a file')
    elif path.stat().st_size == 0:
        return warning(f'{path} is empty')
    else:
        return path.absolute()


def check_cpus(cpus: int | None = 0, verbose: bool = False) -> int:
    cpus = os.cpu_count() if not cpus else min(cpus, os.cpu_count())
    log(f'Using {cpus} CPUs', verbose)
    return cpus


def check_out(path: str, mode: str = "at", parents: bool = True, exist_ok: bool = True) -> Path | TextIO:
    """
    Check if the user wants to create/append a file or directory.
    If it looks like/is already a file (has an extension), return the file object.
    If it looks like/is already a directory, return the directory path.
    """
    if path == '-':  # If the path is '-', return stdout
        return sys.stdout
    if (path := Path(path)).suffix:  # If the path has an extension, it's probably a file
        try:
            return path.open(mode)  # Open the file
        except Exception as e:
            quit_with_error(f'Could not open {path}: {e}')
    if not path.exists():  # Assume directory
        try:
            path.mkdir(parents=parents, exist_ok=exist_ok)  # Create the directory if it doesn't exist
        except Exception as e:
            quit_with_error(f'Could not create {path}: {e}')
    return path


def check_python_version(major: int = 3, minor: int = 8):
    if sys.version_info.major < major or sys.version_info.minor < minor:
        quit_with_error(f'Python version {major}.{minor} or greater required')


def check_biopython_version(major: int = 1, minor: int = 79):
    try:
        from Bio import __version__ as biopython_version
    except ImportError:
        quit_with_error('BioPython is required')
    if ((major_version := int(biopython_version.split('.')[0])) < major or
            (minor_version := int(biopython_version.split('.')[1])) < minor):
        quit_with_error(f'Biopython version {major}.{minor} or greater required, got {major_version}.{minor_version}')


def opener(file: Path | str, check: bool = True, verbose: bool = False, *args, **kwargs):
    """Opens a file with the appropriate open function based on the compression format of the file"""
    with open(check_file(file) if check else file, 'rb') as f:
        file_start = f.read(_READ_N_BYTES)
    compression = next((comp for magic, comp in _COMPRESSION_MAGIC.items() if file_start.startswith(magic)), 'no')
    log(f'Opening {file} with {compression} compression, file start: {file_start}', verbose)
    if compression == 'no':
        return open(file, *args, **kwargs)  # Use the built-in open function
    elif compression == 'gz':
        return gzopen(file, *args, **kwargs)  # Use the gzip open function
    elif compression == 'bz2':
        return bzopen(file, *args, **kwargs)  # Use the bzip2 open function
    else:
        quit_with_error(f'Unsupported compression format: {compression}')


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
        return
    if len(ranges) == 1:
        yield ranges[0]
        return
    current_range = (ranges := ranges if skip_sort else sorted(ranges, key=lambda x: x[0]))[0]  # Start with the first range
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

