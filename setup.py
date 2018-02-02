#!/usr/bin/env python3
"""
Copyright 2018 Ryan Wick (rrwick@gmail.com)
https://github.com/katholt/Kaptive

This file is part of Kaptive. Kaptive is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by the Free Software Foundation,
either version 3 of the License, or (at your option) any later version. Kaptive is distributed in
the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
details. You should have received a copy of the GNU General Public License along with Kaptive. If
not, see <http://www.gnu.org/licenses/>.
"""

from setuptools import setup

with open('README.md', 'rb') as readme:
    long_description = readme.read()
if not isinstance(long_description, str):
    long_description = long_description.decode()

# Get the version from kaptive.py.
__version__ = '0.0.0'
kaptive_script_lines = open('kaptive.py').readlines()
version_line = [x for x in kaptive_script_lines if x.startswith('__version__')][0]
exec(version_line)

setup(name='Kaptive',
      version=__version__,
      description='K and O locus typing for Klebsiella assemblies',
      long_description=long_description,
      url='http://github.com/katholt/Kaptive',
      author='Ryan Wick',
      author_email='rrwick@gmail.com',
      license='GPLv3',
      scripts=['kaptive.py'],
      install_requires=['biopython'],
      zip_safe=False)
