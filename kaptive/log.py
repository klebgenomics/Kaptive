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

import os
import datetime
import textwrap
import sys


END_FORMATTING = '\033[0m'
ITALICS = '\033[3m'
BOLD = '\033[1m'
UNDERLINE = '\033[4m'
RED = '\033[31m'
GREEN = '\033[32m'
MAGENTA = '\033[35m'
YELLOW = '\033[93m'
DIM = '\033[2m'
BLUE = '\033[34m'
CYAN = '\033[36m'
WHITE = '\033[37m'


def bold(text: str):
    return BOLD + text + END_FORMATTING


def italic(text: str):
    return ITALICS + text + END_FORMATTING


def bold_yellow(text: str):
    return YELLOW + BOLD + text + END_FORMATTING


def bold_yellow_underline(text: str):
    return YELLOW + BOLD + UNDERLINE + text + END_FORMATTING


def dim(text: str):
    return DIM + text + END_FORMATTING


def red(text: str):
    return RED + text + END_FORMATTING


def bold_red(text: str):
    return RED + BOLD + text + END_FORMATTING


def cyan(text: str):
    return CYAN + text + END_FORMATTING


def bold_cyan(text: str):
    return CYAN + BOLD + text + END_FORMATTING


def bold_magenta_underline(text: str):
    return MAGENTA + BOLD + UNDERLINE + text + END_FORMATTING


def log(message: str = '', verbose: bool = True, end: str = '\n', sep: str = ' ', flush: bool = True, file=sys.stderr):
    if verbose:
        print(f"{get_timestamp()}: {message}", file=file, flush=flush, end=end, sep=sep)


def warning(text: str):
    width, _ = get_terminal_size_stderr()
    for line in textwrap.wrap(text, width=width-1):
        log(bold_yellow(f"WARNING: {line}"))


def error(text: str):
    width, _ = get_terminal_size_stderr()
    for line in textwrap.wrap(text, width=width-1):
        log(bold_red(f"ERROR: {line}"))


def quit_with_error(text: str):
    error(text)
    sys.exit(1)


def get_timestamp():
    return '{:%Y-%m-%d %H:%M:%S}'.format(datetime.datetime.now())


def get_terminal_size_stderr(fallback=(80, 24)):
    """
    Unlike shutil.get_terminal_size, which looks at stdout, this looks at stderr.
    """
    try:
        size = os.get_terminal_size(sys.__stderr__.fileno())
    except (AttributeError, ValueError, OSError):
        size = os.terminal_size(fallback)
    return size
