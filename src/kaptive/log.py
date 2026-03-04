"""
This module contains functions for logging messages to stderr.

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
from datetime import datetime
import sys
from inspect import stack


def bold(text: str):
    return f"\033[1m{text}\033[0m"


def bold_yellow(text: str):
    return f"\033[1;33m{text}\033[0m"


def bold_red(text: str):
    return f"\033[1;31m{text}\033[0m"


def bold_cyan(text: str):
    return f"\033[1;36m{text}\033[0m"


def log(message: str = '', verbose: bool = True, rjust: int = 20, stack_depth: int = 1):
    """
    Simple function for logging messages to stderr. Only runs if verbose == True.
    Stack depth can be increased if the parent function name needs to be exposed.
    """
    if verbose:  # Only build log if verbosity is requested; simple way of controlling log
        sys.stderr.write(f"{datetime.now():%Y-%m-%d %H:%M:%S} {stack()[stack_depth].function:>{rjust}}] {message}\n")


def warning(message: str):
    for line in message.splitlines():
        log(bold_yellow(f"WARNING] {line}"), verbose=True, stack_depth=2)


def quit_with_error(message: str):
    for line in message.splitlines():
        log(bold_red(f"  ERROR] {line}"), verbose=True, stack_depth=2)
    sys.exit(1)

