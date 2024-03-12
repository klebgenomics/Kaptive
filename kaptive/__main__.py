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

import io
import sys
import re
import argparse
from pathlib import Path

from kaptive.version import __version__
from kaptive.log import bold, quit_with_error
from kaptive.misc import check_python_version, check_programs, get_logo, check_cpus, check_dir, check_file
from kaptive.database import Database, get_database

# Constants -----------------------------------------------------------------------------------------------------------
_ASSEMBLY_HEADERS = [
    'Assembly', 'Best match locus', 'Best match type', 'Confidence', 'Problems', 'Identity', 'Coverage',
    'Length discrepancy', 'Expected genes in locus', 'Expected genes in locus, details', 'Missing expected genes',
    'Other genes in locus', 'Other genes in locus, details', 'Expected genes outside locus',
    'Expected genes outside locus, details', 'Other genes outside locus', 'Other genes outside locus, details',
    'Truncated genes, details'
]
_ASSEMBLY_EXTRA_HEADERS = [
    'Extra genes', 'Contigs', 'Pieces', 'Pieces, details', 'Score', 'Zscore', 'All scores', 'Args'
]


# Functions -----------------------------------------------------------------------------------------------------------
def parse_args(a):
    parser = argparse.ArgumentParser(
        description=get_logo('In silico serotyping'), usage="%(prog)s <command>", add_help=False,
        prog="kaptive", formatter_class=argparse.RawDescriptionHelpFormatter, epilog=f"%(prog)s version {__version__}")

    subparsers = parser.add_subparsers(title=bold('Command'), dest='subparser_name', metavar="")
    assembly_subparser(subparsers)
    # reads_subparser(subparsers)
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
        epilog=f'kaptive assembly v{__version__}', add_help=False, formatter_class=argparse.RawTextHelpFormatter,
        help='In silico serotyping of assemblies', usage="kaptive assembly <db> <input> [<input> ...] [options]")
    opts = assembly_parser.add_argument_group(bold('Inputs'), "")
    opts.add_argument('db', type=get_database, help='Kaptive database path or keyword')
    opts.add_argument('input', nargs='+', type=check_file, help='Assemblies in fasta(.gz) format')

    opts = assembly_parser.add_argument_group(bold('Output options'), "")
    output_opts(opts)

    opts = assembly_parser.add_argument_group(bold('Scoring options'), "")
    opts.add_argument("--score-metric", type=str, default='AS', metavar='',
                      help="Alignment metric to use for scoring (default: %(default)s)")
    opts.add_argument("--weight-metric", type=str, metavar='', default='prop_genes_found',
                      help="Weighting for scoring metric (default: %(default)s)\n"
                           " - none: No weighting\n"
                           " - locus_length: length of the locus\n"
                           " - genes_expected: # of genes expected in the locus\n"
                           " - genes_found: # of genes found in the locus\n"
                           " - prop_genes_found: genes_found / genes_expected")
    opts.add_argument("--min-zscore", type=float, metavar='', default=3.0,
                      help="Minimum zscore for confidence (default: %(default)s)")
    opts.add_argument('--min-cov', type=float, required=False, default=50.0, metavar='',
                      help='Minimum gene %%coverage to be used for scoring (default: %(default)s)')
    opts.add_argument("--gene-threshold", type=float, metavar='',
                      help="Species-level locus gene identity threshold (default: database specific)")

    opts = assembly_parser.add_argument_group(bold('Database options'), "")
    db_opts(opts)
    # opts.add_argument('--is-seqs', type=check_file, metavar='',
    #                   help='Fasta file of IS element sequences to include in the database (default: None)')

    opts = assembly_parser.add_argument_group(bold('Other options'), "")
    other_opts(opts)
    # opts.add_argument('-@', '--mp', const=8, nargs='?', type=int, metavar='#',
    #                   help="Process multiple samples in parallel, optionally pass max workers (default: %(const)s)")
    opts.add_argument('-t', '--threads', type=check_cpus, default=1, metavar='',
                      help="Number of threads for alignment (default: %(default)s)")


# def reads_subparser(subparsers):
#     assembly_parser = subparsers.add_parser(
#         'reads', description=bold_cyan(LOGO + f"\n{'In silico serotyping of reads' : ^43}"),
#         epilog=f'kaptive reads v{__version__}', add_help=False, formatter_class=argparse.RawTextHelpFormatter,
#         help='In silico serotyping of reads', usage="kaptive reads <db> <reads> [<reads> ...] [options]")
#     opts = assembly_parser.add_argument_group(bold('Inputs'), "")
#     opts.add_argument('db', type=get_database, help='Kaptive database path or keyword')
#     opts.add_argument('reads', nargs='+', type=ReadFile.from_path, help='Reads in fastq(.gz) format')
#     opts = assembly_parser.add_argument_group(bold('Output options'), "")
#     output_opts(opts)
#     opts = assembly_parser.add_argument_group(bold('Alignment options'), "")
#     alignmnent_opts(opts)
#     opts = assembly_parser.add_argument_group(bold('Database options'), "")
#     db_opts(opts)
#     opts.add_argument('--is-seqs', type=check_file, metavar='',
#                       help='Fasta file of IS element sequences to include in the database (default: None)')
#     opts = assembly_parser.add_argument_group(bold('Other options'), "")
#     other_opts(opts)


def extract_subparser(subparsers):
    extract_parser = subparsers.add_parser(
        'extract', description=get_logo('Extract entries from a Kaptive database'),
        epilog=f'kaptive extract v{__version__}', add_help=False, formatter_class=argparse.RawTextHelpFormatter,
        help='Extract entries from a Kaptive database', usage="kaptive extract <db> <format> [options]")
    opts = extract_parser.add_argument_group(bold('Inputs'), "\n - Note: combine with --filter to select loci")
    opts.add_argument('db', help='Kaptive database path or keyword', type=get_database)
    opts.add_argument('format', choices=['loci', 'genes', 'proteins', 'gbk', 'gff', 'ids'], metavar='format',
                      help='Format to extract database\n - loci: Loci (fasta nucleotide)\n'
                           ' - genes: Genes (fasta nucleotide)\n - proteins: Proteins (fasta amino acid)\n'
                           ' - gbk: Genbank format\n - gff: GFF in NCBI format\n - ids: List of Locus IDs')
    opts = extract_parser.add_argument_group(bold('Output options'), "")
    opts.add_argument('-o', '--out', metavar='', default=sys.stdout, type=Path, help='Output file (default: stdout)')
    opts = extract_parser.add_argument_group(bold('Database options'), "")
    db_opts(opts)
    opts = extract_parser.add_argument_group(bold('Other options'), "")
    other_opts(opts)


def convert_subparser(subparsers):
    convert_parser = subparsers.add_parser(
        'convert', description=get_logo('Convert Kaptive results into different formats'),
        epilog=f'kaptive convert v{__version__}', add_help=False, formatter_class=argparse.RawTextHelpFormatter,
        help='Convert Kaptive results into different formats',
        usage="kaptive convert <result> [<result> ...] [options]")
    opts = convert_parser.add_argument_group(
        bold('Inputs'),
        "\n - Note: If you used --is-seqs during the run, make sure to provide the same fasta file here")
    opts.add_argument('db', help='Kaptive database in genbank format', type=get_database)
    opts.add_argument('json', help='Kaptive result files', type=check_file, nargs='+')

    opts = convert_parser.add_argument_group(bold('Filter options'),
                                             "\n - Note: filters take precedence in descending order")
    opts.add_argument('-r', '--regex', metavar='', type=re.compile,
                      help='Regex to filter the string interpretation of the result (default: All)')
    opts.add_argument('-l', '--loci', metavar='', nargs='+',
                      help='Space-separated list to filter locus names (default: All)')
    opts.add_argument('-s', '--samples', metavar='', nargs='+',
                      help='Space-separated list to filter sample names (default: All)')

    opts = convert_parser.add_argument_group(bold('Output options'), "")
    # opts.add_argument('-o', '--out', metavar='', default=sys.stdout, type=argparse.FileType('wt'),
    #                   help='Output file (default: stdout)')
    opts.add_argument('-f', '--format', metavar='', default='json',
                      choices=['json', 'tsv', 'locus', 'genes', 'proteins'],
                      help='Output format (default: %(default)s)\n - json: JSON format\n - tsv: Tab-separated values\n'
                           ' - locus: Locus nucleotide sequence in fasta format\n'
                           ' - proteins: Proteins in fasta format\n - genes: Genes in fasta format')

    opts = convert_parser.add_argument_group(bold('Database options'), "")
    db_opts(opts)
    opts = convert_parser.add_argument_group(bold('Other options'), "")
    other_opts(opts)


def db_opts(opts: argparse.ArgumentParser):
    opts.add_argument('--locus-regex', type=re.compile, metavar='',
                      help='Pattern to match locus names in db source note, (default: %(default)s)')
    opts.add_argument('--type-regex', type=re.compile, metavar='',
                      help='Pattern to match locus types in db source note, (default: %(default)s)')
    opts.add_argument('--filter', type=re.compile, metavar='',
                      help='Pattern to select loci to include in the database (default: All)')


def output_opts(opts: argparse.ArgumentParser):
    opts.add_argument('-o', '--out', metavar='', default=sys.stdout, type=argparse.FileType('at'),
                      help='Output file (default: stdout)')
    opts.add_argument('--fasta', metavar='path', nargs='?', default=None, const='.', type=check_dir,
                      help='Output locus sequence to "{input}_kaptive_results.fna"\n'
                           'Optionally pass output directory (default: current directory)')
    opts.add_argument('--json', metavar='path', nargs='?', default=None, const='kaptive_results.json',
                      type=argparse.FileType('at'),
                      help='Output results to JSON lines\n'
                           'Optionally pass file name (default: %(const)s)')
    opts.add_argument('--draw', metavar='path', nargs='?', default=None, const='.', type=check_dir,
                      help='Output locus figures to "{input}_kaptive_results.{fmt}"\n'
                           'Optionally pass output directory (default: current directory)')
    opts.add_argument('--draw-fmt', default='png', metavar='png,svg',
                      help='Format for locus figures (default: %(default)s)')
    # opts.add_argument('--bokeh', action='store_true', help='Plot locus figures using bokeh')
    opts.add_argument('--no-header', action='store_true', help='Do not print header line')
    opts.add_argument('--debug', action='store_true', help='Append debug columns to table output')


def other_opts(opts: argparse.ArgumentParser):
    opts.add_argument('-V', '--verbose', action='store_true', help='Print debug messages to stderr')
    opts.add_argument('-v', '--version', help='Show version number and exit', metavar='')
    opts.add_argument('-h', '--help', help='Show this help message and exit', metavar='')


def write_headers(out: io.TextIOWrapper, no_header: bool = False, subparser_name: str = '', debug: bool = False):
    """
    Write headers to output file if not already written
    """
    if out.name != '<stdout>' and out.tell() != 0:  # If file is path and not already written to
        no_header = True  # Headers already written, useful for running on HPC
    if not no_header:
        if subparser_name == 'assembly':
            out.write('\t'.join(_ASSEMBLY_HEADERS + _ASSEMBLY_EXTRA_HEADERS if debug else _ASSEMBLY_HEADERS) + '\n')


# Main -----------------------------------------------------------------------------------------------------------------
def main():
    check_python_version()
    args = parse_args(sys.argv[1:])

    if args.subparser_name == 'assembly':
        check_programs(['minimap2'], verbose=args.verbose)  # Check for minimap2

        # Load database in memory, we don't need to load the full sequences (False)
        db = Database.from_genbank(
            path=args.db, gene_threshold=args.gene_threshold, locus_filter=args.filter, verbose=args.verbose,
            load_seq=False, locus_regex=args.locus_regex, type_regex=args.type_regex)

        # Write headers to output file if not already written
        write_headers(args.out, args.no_header, args.subparser_name, args.debug)

        from kaptive.assembly import typing_pipeline  # Import here to avoid circular import
        if args.draw:
            from matplotlib import pyplot as plt  # Import here to avoid unnecessary imports

        for assembly in args.input:  # type: 'Path'
            if (result := typing_pipeline(
                    assembly, db, args.threads, args.min_cov, args.min_zscore, args.score_metric, args.weight_metric, args.verbose)):
                args.out.write(result.as_table(args.debug))  # type: 'TextIOWrapper'
                if args.fasta:  # Create and write locus pieces to fasta file
                    (args.fasta / f'{result.assembly.name}_kaptive_results.fna').write_text(result.as_fasta())
                if args.json:  # type: 'TextIOWrapper'
                    args.json.write(result.as_json())
                if args.draw:  # type: 'Path'
                    ax = result.as_GraphicRecord().plot()[0]  # Create figure and axis
                    ax.set_title(
                        f"{result.assembly.name} {result.best_match} ({result.phenotype}) - {result.confidence}")  # Add title
                    ax.figure.savefig(args.draw / f'{result.assembly.name}_kaptive_results.{args.draw_fmt}')  # Save figure
                    plt.close()  # Close figure to prevent memory leaks

    elif args.subparser_name == 'extract':
        from kaptive.database import extract
        extract(args)

    elif args.subparser_name == 'convert':
        from kaptive.typing import parse_results
        args.db = Database.from_genbank(  # Load database in memory, we don't need to load the full sequences (False)
            args.db, args.is_seqs,  args.filter, False, locus_regex=args.locus_regex, type_regex=args.type_regex)
        for result_file in args.json:
            for result in parse_results(result_file, args.db, args.regex, args.samples, args.loci):
                if args.format == 'json':
                    sys.stdout.write(result.as_json())
                elif args.format == 'tsv':
                    sys.stdout.write(result.as_table())
                elif args.format == 'locus':
                    sys.stdout.write(result.as_fasta())
                elif args.format == 'genes':
                    sys.stdout.write(result.as_gene_fasta())
                elif args.format == 'proteins':
                    sys.stdout.write(result.as_protein_fasta())

