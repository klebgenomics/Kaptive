#!/usr/bin/env python3
"""
Copyright 2018 Ryan Wick (rrwick@gmail.com)
https://github.com/katholt/Kaptive

Kaptive - SLURM job generator

This tool generates SLURM jobs so the Kaptive program can be run in parallel. All instances
of the script will append to the same output table.

This file is part of Kaptive. Kaptive is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by the Free Software Foundation,
either version 3 of the License, or (at your option) any later version. Kaptive is distributed in
the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
details. You should have received a copy of the GNU General Public License along with Kaptive. If
not, see <http://www.gnu.org/licenses/>.
"""

from __future__ import print_function
from __future__ import division
import sys
import argparse
import os
import time


def find_script():
    """Returns the location of kaptive.py."""
    script_dir = os.path.dirname(os.path.realpath(__file__))
    while not os.path.isfile(os.path.join(script_dir, 'kaptive.py')):
        script_dir = os.path.dirname(script_dir)
    if not os.path.isfile(os.path.join(script_dir, 'kaptive.py')):
        script_dir = os.getcwd()
        while not os.path.isfile(os.path.join(script_dir, 'kaptive.py')):
            script_dir = os.path.dirname(script_dir)
    if not os.path.isfile(os.path.join(script_dir, 'kaptive.py')):
        quit_with_error('Could not find kaptive.py')
    else:
        return os.path.join(script_dir, 'kaptive.py')


script_path = find_script()
sys.path.append(os.path.dirname(script_path))
import kaptive


def main():
    """Script execution starts here."""
    args = get_arguments()
    kaptive.fix_paths(args)
    kaptive.check_files_exist(args.assembly + [args.k_refs])

    main_parser = kaptive.get_argument_parser()
    using_threads = (args.threads != main_parser.get_default('threads'))
    using_start_end_margin = (args.start_end_margin != main_parser.get_default('start_end_margin'))
    using_min_gene_cov = (args.min_gene_cov != main_parser.get_default('min_gene_cov'))
    using_min_gene_id = (args.min_gene_id != main_parser.get_default('min_gene_id'))
    using_low_gene_id = (args.low_gene_id != main_parser.get_default('low_gene_id'))
    using_min_assembly_piece = (args.min_assembly_piece != main_parser.get_default('min_assembly_piece'))
    using_gap_fill_size = (args.gap_fill_size != main_parser.get_default('gap_fill_size'))

    for assembly_filename in args.assembly:
        assembly_name = simple_assembly_name(assembly_filename)
        assembly_filename = os.path.abspath(assembly_filename)

        # Build the kaptive.py command.
        kaptive_command = script_path
        kaptive_command += ' --assembly ' + os.path.abspath(assembly_filename)
        kaptive_command += ' --k_refs ' + os.path.abspath(args.k_refs)
        if args.allelic_typing:
            kaptive_command += ' --allelic_typing ' + os.path.abspath(args.allelic_typing)
        kaptive_command += ' --out ' + os.path.abspath(args.out)
        if args.verbose:
            kaptive_command += ' --verbose'
        if using_threads:
            kaptive_command += ' --threads ' + str(args.threads)
        if args.no_seq_out:
            kaptive_command += ' --no_seq_out'
        if args.no_table:
            kaptive_command += ' --no_table'
        if args.no_json:
            kaptive_command += ' --no_json'
        if using_start_end_margin:
            kaptive_command += ' --start_end_margin ' + str(args.start_end_margin)
        if using_min_gene_cov:
            kaptive_command += ' --min_gene_cov ' + str(args.min_gene_cov)
        if using_min_gene_id:
            kaptive_command += ' --min_gene_id ' + str(args.min_gene_id)
        if using_low_gene_id:
            kaptive_command += ' --low_gene_id ' + str(args.low_gene_id)
        if using_min_assembly_piece:
            kaptive_command += ' --min_assembly_piece ' + str(args.min_assembly_piece)
        if using_gap_fill_size:
            kaptive_command += ' --gap_fill_size ' + str(args.gap_fill_size)

        # Build the SLURM script.
        slurm = ['#!/bin/bash',
                 '#SBATCH -p sysgen',
                 '#SBATCH --job-name=kaptive_' + assembly_name,
                 '#SBATCH --ntasks=1',
                 '#SBATCH --mem-per-cpu=' + args.memory,
                 '#SBATCH --time=' + args.walltime,
                 'module load BLAST+/2.2.30-vlsci_intel-2015.08.25-Python-2.7.10',
                 'module load Python/3.5.2-vlsci_intel-2015.08.25',
                 'module load Biopython/1.67-iccifort-2015.2.164-GCC-4.9.2-Python-3.5.2',
                 kaptive_command]

        slurm_string = '\n'.join(slurm)
        print(slurm_string)
        print()
        os.system('echo "' + slurm_string + '" | sbatch')
        time.sleep(0.1)
        print()


def get_arguments():
    """Specifies the command line arguments required by the script."""
    main_parser = kaptive.get_argument_parser()

    parser = argparse.ArgumentParser(description='Kaptive - SLURM job generator',
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--walltime', type=str, required=False,
                        help='wall time', default='0-0:10')
    parser.add_argument('--memory', type=str, required=False,
                        help='memory (Megabytes)', default='2048')
    kaptive.add_arguments_to_parser(parser)
    return parser.parse_args()


def quit_with_error(message): # type: (str) -> None
    """Displays the given message and ends the program's execution."""
    print('Error:', message, file=sys.stderr)
    sys.exit(1)


def simple_assembly_name(assembly_filename):
    """Returns just the assembly name, without a path or file extension."""
    assembly_filename = os.path.basename(assembly_filename)
    assembly_filename = rchop(assembly_filename, '.fasta')
    assembly_filename = rchop(assembly_filename, '.FASTA')
    assembly_filename = rchop(assembly_filename, '.fa')
    return assembly_filename


def rchop(s, e):
    if s.endswith(e):
        return s[:-len(e)]
    else:
        return s


if __name__ == '__main__':
    main()
