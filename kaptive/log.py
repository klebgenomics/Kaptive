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
from inspect import currentframe


def bold(text: str):
    return f"\033[1m{text}\033[0m"


def bold_yellow(text: str):
    return f"\033[1;33m{text}\033[0m"


def bold_red(text: str):
    return f"\033[1;31m{text}\033[0m"


def bold_cyan(text: str):
    return f"\033[1;36m{text}\033[0m"


def log(message: str = '', verbose: bool = True, rjust: int = 15):
    if verbose:
        prefix = f"[{datetime.now().strftime('%Y-%m-%d %H:%M:%S')} {currentframe().f_back.f_code.co_name:>{rjust}}]"
        sys.stderr.write(f"{prefix} {message}\n")
        # sys.stderr.flush()


def warning(message: str):
    # The currentframe() in log() will tell user the message is a warning, so no need to prefix message
    for line in message.splitlines():
        log(bold_yellow(line))


def error(message: str):
    # The currentframe() in log() will tell user the message is an error, so no need to prefix message
    for line in message.splitlines():
        log(bold_red(line))


def quit_with_error(text: str):
    error(text)
    sys.exit(1)

