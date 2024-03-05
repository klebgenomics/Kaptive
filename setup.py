#!/usr/bin/env python3
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

from setuptools import setup, find_packages
from distutils.util import convert_path

with open('README.md', 'rb') as readme:
    long_description = readme.read()
if not isinstance(long_description, str):
    long_description = long_description.decode()

# Get the version from kaptive.__version__.py.
__version__ = '0.0.0'
ver_path = convert_path('kaptive/version.py')
with open(ver_path) as ver_file:
    exec(ver_file.read())

setup(
    name='kaptive',
    version=__version__,
    description='In silico serotyping',
    long_description=long_description,
    url='http://github.com/klebgenomics/Kaptive',
    author='Tom Stanton',
    author_email='tomdstanton@gmail.com',
    license='GPLv3',
    install_requires=['biopython'],
    zip_safe=False,
    packages=find_packages(),
    package_data={'kaptive': ['reference_database/*', 'extras/*']},
    package_dir={'kaptive': 'kaptive'},
    entry_points={
        'console_scripts': [
            'kaptive = kaptive.__main__:main',
        ]
    }
)
