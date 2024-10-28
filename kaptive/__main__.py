#!/usr/bin/env python
"""
This is the main entry point for the kaptive package. It is called when the package is run as a script via entry_points.

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

import sys
import re
import argparse
from io import TextIOWrapper
from os import path

from Bio import __version__ as biopython_version

from kaptive.version import __version__
from kaptive.log import bold, quit_with_error, log
from kaptive.utils import get_logo, check_out, check_cpus

# Constants -----------------------------------------------------------------------------------------------------------
_URL = 'https://kaptive.readthedocs.io/en/latest/'


# Functions -----------------------------------------------------------------------------------------------------------
def parse_args(a) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=get_logo('In silico serotyping'), usage="%(prog)s <command>", add_help=False,
        prog="kaptive", formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=f'For more help, visit: {bold(_URL)}')

    subparsers = parser.add_subparsers(title=bold('Command'), dest='subparser_name', metavar="")
    assembly_subparser(subparsers)
    extract_subparser(subparsers)
    convert_subparser(subparsers)
    opts = parser.add_argument_group(bold('Other options'), '')
    other_opts(opts)

    if len(a) == 0:  # No arguments, print help message
        parser.print_help(sys.stderr)
        quit_with_error(f'Please specify a command; choose from {{assembly,extract,convert}}')
    if any(x in a for x in {'-v', '--version'}):  # Version message
        print(__version__)
        sys.exit(0)
    if subparser := subparsers.choices.get(a[0], None):  # Check if the first argument is a subparser
        if len(a) == 1:  # Subparser help message
            subparser.print_help(sys.stderr)
            quit_with_error(f'Insufficient arguments for kaptive {a[0]}')
        if any(x in a[1:] for x in {'-h', '--help'}):  # Subparser help message
            subparser.print_help(sys.stderr)
            sys.exit(0)
    elif any(x in a for x in {'-h', '--help'}):  # Help message
        parser.print_help(sys.stderr)
        sys.exit(0)
    else:  # Unknown command
        parser.print_help(sys.stderr)
        quit_with_error(f'Unknown command "{a[0]}"; choose from {{assembly,extract,convert}}')
    return parser.parse_args(a)


def assembly_subparser(subparsers):
    assembly_parser = subparsers.add_parser(
        'assembly', description=get_logo('In silico serotyping of assemblies'),
        epilog=f'For more help, visit: {bold(_URL)}', add_help=False, formatter_class=argparse.RawTextHelpFormatter,
        help='In silico serotyping of assemblies', usage="kaptive assembly <db> <fasta> [<fasta> ...] [options]")
    opts = assembly_parser.add_argument_group(bold('Inputs'), "")
    opts.add_argument('db', metavar='db path/keyword', help='Kaptive database path or keyword')
    opts.add_argument('input', nargs='+', metavar='fasta', help='Assemblies in fasta(.gz|.xz|.bz2) format')
    opts = assembly_parser.add_argument_group(bold('Output options'), "\nNote, text outputs accept '-' for stdout")
    # Note these are different to the convert output options as TSV is the main output and fna is the main fasta output
    opts.add_argument('-o', '--out', metavar='', default=sys.stdout, type=argparse.FileType('at'),
                      help='Output file to write/append tabular results to (default: stdout)')
    opts.add_argument('-f', '--fasta', metavar='', nargs='?', default=None, const='.', type=check_out,
                      help='Turn on fasta output\n'
                           'Accepts a single file or a directory (default: cwd)')
    opts.add_argument('-j', '--json', metavar='', nargs='?', default=None, const='kaptive_results.json',
                      type=argparse.FileType('at'),
                      help='Turn on JSON lines output\n'
                           'Optionally choose file (can be existing) (default: %(const)s)')
    opts.add_argument('-s', '--scores', metavar='', nargs='?', default=None, const=sys.stdout,
                      type=argparse.FileType('at'),
                      help='Dump locus score matrix to tsv (typing will not be performed!)\n'
                           'Optionally choose file (can be existing) (default: stdout)')
    other_fmt_opts(opts)
    opts = assembly_parser.add_argument_group(bold('Scoring options'), "")
    opts.add_argument('--min-cov', type=float, required=False, default=50.0, metavar='',
                      help='Minimum gene %%coverage (blen/q_len*100) to be used for scoring (default: %(default)s)')
    opts.add_argument("--score-metric", metavar='', default=0, type=int, choices=range(4),
                      help="Metric for scoring each locus (default: %(default)s)\n"
                           "  0: AS (alignment score of genes found)\n"
                           "  1: mlen (matching bases of genes found)\n"
                           "  2: blen (aligned bases of genes found)\n"
                           "  3: q_len (query length of genes found)")
    opts.add_argument("--weight-metric", metavar='', default=3, type=int, choices=range(6),
                      help="Weighting for the 1st stage of the scoring algorithm (default: %(default)s)\n"
                           "  0: No weighting\n"
                           "  1: Number of genes found\n"
                           "  2: Number of genes expected\n"
                           "  3: Proportion of genes found\n"
                           "  4: blen (aligned bases of genes found)\n"
                           "  5: q_len (query length of genes found)")
    opts.add_argument('--n-best', type=int, default=2, metavar='',
                      help='Number of best loci from the 1st round of scoring to be\n'
                           'fully aligned to the assembly (default: %(default)s)')

    opts = assembly_parser.add_argument_group(bold('Confidence options'), "")
    opts.add_argument("--gene-threshold", type=float, metavar='',
                      help="Species-level locus gene identity threshold (default: database specific)")
    opts.add_argument("--max-other-genes", type=int, metavar='', default=1,
                      help="Typeable if <= other genes (default: %(default)s)")
    opts.add_argument("--percent-expected", type=float, metavar='', default=50,
                      help="Typeable if >= %% expected genes (default: %(default)s)")
    opts.add_argument("--below-threshold", type=bool, default=False, metavar='',
                      help="Typeable if any genes are below threshold (default: %(default)s)")
    opts = assembly_parser.add_argument_group(bold('Database options'), "")
    db_opts(opts)
    opts.add_argument('--filter', type=re.compile, metavar='',
                      help='Python regular-expression to select loci to include in the database')
    opts = assembly_parser.add_argument_group(bold('Other options'), "")
    other_opts(opts)
    opts.add_argument('-t', '--threads', type=check_cpus, default=check_cpus(), metavar='',
                      help="Number of alignment threads or 0 for all available (default: 0)")


def convert_subparser(subparsers):
    convert_parser = subparsers.add_parser(
        'convert', description=get_logo('Convert Kaptive results into different formats'),
        epilog=f'For more help, visit: {bold(_URL)}', add_help=False, formatter_class=argparse.RawTextHelpFormatter,
        help='Convert Kaptive results into different formats', usage="kaptive convert <db> <json> [formats] [options]")
    opts = convert_parser.add_argument_group(bold('Inputs'), "")
    opts.add_argument('db', metavar='db path/keyword', help='Kaptive database path or keyword')
    opts.add_argument('input', help='Kaptive JSON lines file or - for stdin', type=argparse.FileType('rt'),
                      metavar='json')
    opts = convert_parser.add_argument_group(bold('Formats'), "\nNote, text outputs accept '-' for stdout")
    opts.add_argument('-t', '--tsv', metavar='', nargs='?', default=None, const='-', type=check_out,
                      help='Convert to tabular format in file (default: stdout)')
    opts.add_argument('-j', '--json', metavar='', nargs='?', default=None, const='-', type=check_out,
                      help='Convert to JSON lines format in file (default: stdout)')
    fmt_opts(opts)
    other_fmt_opts(opts)
    opts = convert_parser.add_argument_group(bold('Filter options'),
                                             "\nNote, filters take precedence in descending order")
    opts.add_argument('-r', '--regex', metavar='', type=re.compile,
                      help='Python regular-expression to select JSON lines (default: All)')
    opts.add_argument('-l', '--loci', metavar='', nargs='+',
                      help='Space-separated list to filter locus names (default: All)')
    opts.add_argument('-s', '--samples', metavar='', nargs='+',
                      help='Space-separated list to filter sample names (default: All)')
    opts = convert_parser.add_argument_group(bold('Database options'), "")
    db_opts(opts)
    # Note, we don't allow users to filter the database here in case the results contain a locus that has been filtered
    # out of the database
    opts = convert_parser.add_argument_group(bold('Other options'), "")
    other_opts(opts)


def extract_subparser(subparsers):
    extract_parser = subparsers.add_parser(
        'extract', description=get_logo('Extract entries from a Kaptive database'),
        epilog=f'For more help, visit: {bold(_URL)}', add_help=False, formatter_class=argparse.RawTextHelpFormatter,
        help='Extract entries from a Kaptive database', usage="kaptive extract <db> [formats] [options]")
    opts = extract_parser.add_argument_group(bold('Inputs'), "\nNote, combine with --filter to select loci")
    opts.add_argument('db', metavar='db path/keyword', help='Kaptive database path or keyword')
    opts = extract_parser.add_argument_group(bold('Formats'), "\nNote, text outputs accept '-' for stdout")
    fmt_opts(opts)
    opts = extract_parser.add_argument_group(bold('Database options'), "")
    db_opts(opts)
    opts.add_argument('--filter', type=re.compile, metavar='',
                      help='Python regular-expression to select loci to include in the database')
    opts = extract_parser.add_argument_group(bold('Other options'), "")
    other_opts(opts)


def fmt_opts(opts: argparse.ArgumentParser):
    """Format opts shared by convert and extract"""
    opts.add_argument('--fna', metavar='', nargs='?', default=None, const='.', type=check_out,
                      help='Convert to locus nucleotide sequences in fasta format\n'
                           'Accepts a single file or a directory (default: cwd)')
    opts.add_argument('--ffn', metavar='', nargs='?', default=None, const='.', type=check_out,
                      help='Convert to locus gene nucleotide sequences in fasta format\n'
                           'Accepts a single file or a directory (default: cwd)')
    opts.add_argument('--faa', metavar='', nargs='?', default=None, const='.', type=check_out,
                      help='Convert to locus gene protein sequences in fasta format\n'
                           'Accepts a single file or a directory (default: cwd)')


def other_fmt_opts(opts: argparse.ArgumentParser):
    """Format opts shared by convert and assembly"""
    opts.add_argument('-p', '--plot', metavar='', nargs='?', default=None, const='.', type=check_out,
                      help='Plot results to "./{assembly}_kaptive_results.{fmt}"\n'
                           'Optionally choose a directory (default: cwd)')
    opts.add_argument('--plot-fmt', default='png', metavar='png/svg', choices={'png', 'svg'},
                      help='Format for locus plots (default: %(default)s)')
    opts.add_argument('--no-header', action='store_true', help='Suppress header line')


def db_opts(opts: argparse.ArgumentParser):
    opts.add_argument('--locus-regex', type=re.compile, metavar='',
                      help=f'Python regular-expression to match locus names in db source note')
    opts.add_argument('--type-regex', type=re.compile, metavar='',
                      help=f'Python regular-expression to match locus types in db source note')


def other_opts(opts: argparse.ArgumentParser):
    opts.add_argument('-V', '--verbose', action='store_true', help='Print debug messages to stderr')
    opts.add_argument('-v', '--version', help='Show version number and exit', metavar='')
    opts.add_argument('-h', '--help', help='Show this help message and exit', metavar='')


# Main -----------------------------------------------------------------------------------------------------------------
def main():
    if sys.version_info.major < 3 or sys.version_info.minor < 9:
        quit_with_error(f'Python version 3.9 or greater required')

    if int(biopython_version.split('.')[0]) < 1 or int(biopython_version.split('.')[1]) < 83:
        quit_with_error('Biopython version 1.83 or greater required')

    args = parse_args(sys.argv[1:])  # Parse the arguments

    # Assembly mode ----------------------------------------------------------------------------------------------------
    if args.subparser_name == 'assembly':
        from kaptive.assembly import typing_pipeline, write_headers
        from kaptive.database import load_database

        args.db = load_database(
            args.db, args.gene_threshold, locus_filter=args.filter, load_locus_seqs=True, verbose=args.verbose,
            extract_translations=False, locus_regex=args.locus_regex, type_regex=args.type_regex)

        write_headers(args.scores or args.out, args.no_header, args.scores)

        for assembly in args.input:
            if result := typing_pipeline(assembly, args.db, args.threads, args.score_metric, args.weight_metric,
                                         args.min_cov, args.n_best, args.max_other_genes, args.percent_expected,
                                         args.below_threshold, args.scores, args.verbose):
                result.write(args.out, args.json, args.fasta, None, None, args.plot, args.plot_fmt)

    # Extract mode -----------------------------------------------------------------------------------------------------
    elif args.subparser_name == 'extract':
        from kaptive.database import parse_database
        for locus in parse_database(args.db, args.filter, args.fna, args.faa, args.verbose,
                                    locus_regex=args.locus_regex, type_regex=args.type_regex):
            locus.write(args.fna, args.ffn, args.faa)

    # Convert mode -----------------------------------------------------------------------------------------------------
    elif args.subparser_name == 'convert':
        from kaptive.database import load_database
        from kaptive.assembly import parse_result, write_headers

        args.db = load_database(  # Load database in memory, we don't need to load the full sequences (False)
            args.db, verbose=args.verbose, load_locus_seqs=False, extract_translations=False,
            locus_regex=args.locus_regex, type_regex=args.type_regex)

        write_headers(args.tsv, args.no_header)

        for line in args.input:
            if result := parse_result(line, args.db, args.regex, args.samples, args.loci):
                result.write(args.tsv, args.json, args.fna, args.ffn, args.faa, args.plot, args.plot_fmt)

    # Cleanup ----------------------------------------------------------------------------------------------------------
    for attr in vars(args):  # Close all open files in the args namespace if they aren't sys.stdout or sys.stdin
        if (x := getattr(args, attr, None)) and isinstance(x, TextIOWrapper) and x not in {sys.stdout, sys.stdin}:
            x.close()  # Close the file

    log("Done!", verbose=args.verbose)
