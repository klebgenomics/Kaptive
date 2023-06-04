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
from re import compile, IGNORECASE
from typing import List
from pathlib import Path
from gzip import open as gzopen

from .log import warning, quit_with_error

# File path regexes ---------------------------------------------------------------------------------------------------
FASTQ_REGEX = compile('|'.join([
    '_R[12]\.(f(?:ast)?q(?:\.gz)?)$',
    '_R[12]_[0-9]+?\.(f(?:ast)?q(?:\.gz)?)$',
    '_R[12].[0-9]+?\.(f(?:ast)?q(?:\.gz)?)$',
    '_[12]\.(f(?:ast)?q(?:\.gz)?)$',
    '_[12]_[0-9]+?\.(f(?:ast)?q(?:\.gz)?)$',
    '_[12].[0-9]+?\.(f(?:ast)?q(?:\.gz)?)$',
    '\.(f(?:ast)?q(?:\.gz)?)$']), IGNORECASE)


# Classes -------------------------------------------------------------------------------------------------------------
class ReadGroup:
    """
    Sample class to hold sample name and reads
    Acts like a ReadGroup but will also be used to hold results
    """

    def __init__(self, reads: List['ReadFile'], sample_name):
        self.reads = reads
        self.name = sample_name
        [setattr(read, 'sample', self) for read in self.reads]
        self.path = self.reads[0].path.parent
        self.prefix = self.path / self.name

    def __repr__(self):
        return self.name


class ReadFile:
    def __init__(self, filepath: Path, extension: str):
        self.path = filepath
        self.extension = extension
        self.sample = None
        self.is_gzipped = self.extension.endswith('.gz')
        self.open_mode = 'rt' if self.is_gzipped else 'r'
        self.open_function = open if not self.is_gzipped else gzopen

    def __repr__(self):
        return str(self.path)


def load_read_files(file_list: List[Path]):
    reads = {}
    for file in file_list:
        extension_match = FASTQ_REGEX.search(str(file))
        if extension_match:
            sample_name = file.name.replace(extension_match[0], '')
            if sample_name in reads.keys():
                reads[sample_name].append(ReadFile(file, extension_match[0]))
            else:
                reads[sample_name] = [ReadFile(file, extension_match[0])]
        else:
            warning(f'{file.name} does not match the fastq extension regex')
    samples = []
    if reads:
        for sample, files in reads.items():
            if len(files) != 2:
                warning(f'Files for {sample} files not paired: {" ".join(str(i) for i in files)}')
            else:
                samples.append(ReadGroup(files, sample))
    if not samples:
        quit_with_error('No files to analyse')
    return samples

