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
from json import dumps, loads

from Bio import SeqIO

from kaptive.version import __version__
from kaptive.log import bold, quit_with_error
from kaptive.misc import (check_python_version, check_biopython_version, check_programs, get_logo, check_cpus,
                          check_dir, check_file)
from kaptive.database import Database, get_database
from kaptive.assembly import typing_pipeline, TypingResult

# Constants -----------------------------------------------------------------------------------------------------------
_ASSEMBLY_HEADERS = [
    'Assembly', 'Best match locus', 'Best match type', 'Confidence', 'Problems', 'Identity', 'Coverage',
    'Length discrepancy', 'Expected genes in locus', 'Expected genes in locus, details', 'Missing expected genes',
    'Other genes in locus', 'Other genes in locus, details', 'Expected genes outside locus',
    'Expected genes outside locus, details', 'Other genes outside locus', 'Other genes outside locus, details',
    'Truncated genes, details'
]
_ASSEMBLY_EXTRA_HEADERS = ['Extra genes', 'Pieces, details', 'Score', 'Zscore', 'All scores', 'Scoring args',
                           'Confidence args']
_URL = 'https://kaptive.readthedocs.io/en/latest/'


# Functions -----------------------------------------------------------------------------------------------------------
def parse_args(a):
    parser = argparse.ArgumentParser(
        description=get_logo('In silico serotyping'), usage="%(prog)s <command>", add_help=False,
        prog="kaptive", formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=f'For more help, visit: {bold(_URL)}')

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
        epilog=f'For more help, visit: {bold(_URL)}', add_help=False, formatter_class=argparse.RawTextHelpFormatter,
        help='In silico serotyping of assemblies', usage="kaptive assembly <db> <fasta> [<fasta> ...] [options]")
    opts = assembly_parser.add_argument_group(bold('Inputs'), "")
    opts.add_argument('db', type=get_database, metavar='db path/keyword', help='Kaptive database path or keyword')
    opts.add_argument('input', nargs='+', metavar='fasta', type=check_file, help='Assemblies in fasta format')

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
    opts.add_argument('--min-cov', type=float, required=False, default=50.0, metavar='',
                      help='Minimum gene %%coverage to be used for scoring (default: %(default)s)')

    opts = assembly_parser.add_argument_group(bold('Confidence options'), "")
    opts.add_argument("--gene-threshold", type=float, metavar='',
                      help="Species-level locus gene identity threshold (default: database specific)")
    opts.add_argument("--max-other-genes", type=int, metavar='', default=1,
                      help="Typeable if <= other genes (default: %(default)s)")
    opts.add_argument("--percent-expected-genes", type=float, metavar='', default=50,
                      help="Typeable if >= %% expected genes (default: %(default)s)")
    opts.add_argument("--allow-below-threshold", action='store_true', help="Typeable if any genes are below threshold")

    opts = assembly_parser.add_argument_group(bold('Database options'), "")
    db_opts(opts)
    # opts.add_argument('--is-seqs', type=check_file, metavar='',
    #                   help='Fasta file of IS element sequences to include in the database (default: None)')

    opts = assembly_parser.add_argument_group(bold('Other options'), "")
    other_opts(opts)
    # opts.add_argument('-@', '--mp', const=8, nargs='?', type=int, metavar='#',
    #                   help="Process multiple samples in parallel, optionally pass max workers (default: %(const)s)")
    opts.add_argument('-t', '--threads', type=check_cpus, default=check_cpus(None), metavar='',
                      help="Number of threads for alignment (default: maximum available CPUs)")


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
        epilog=f'For more help, visit: {bold(_URL)}', add_help=False, formatter_class=argparse.RawTextHelpFormatter,
        help='Extract entries from a Kaptive database', usage="kaptive extract <db> <format> [options]")
    opts = extract_parser.add_argument_group(bold('Inputs'), "\n - Note: combine with --filter to select loci")
    opts.add_argument('db', type=get_database, metavar='db path/keyword', help='Kaptive database path or keyword')
    opts.add_argument('format', choices=['loci', 'genes', 'proteins', 'genbank'], metavar='format',
                      help='Format to extract database\n - loci: Loci (fasta nucleotide)\n'
                           ' - genes: Genes (fasta nucleotide)\n - proteins: Proteins (fasta amino acid)\n'
                           ' - genbank: Genbank format')
    opts = extract_parser.add_argument_group(bold('Output options'), "")
    opts.add_argument('-o', '--out', metavar='', default=sys.stdout, type=argparse.FileType('at'),
                      help='Output file to write/append loci to (default: stdout)')
    opts.add_argument('-d', '--outdir', metavar='', type=check_dir,
                      help='Output directory for converted results\n'
                           ' - Note: This forces the output to be written to files (instead of stdout)\n'
                           '         and one file will be written per locus.'
                      )
    opts = extract_parser.add_argument_group(bold('Database options'), "")
    db_opts(opts)
    opts = extract_parser.add_argument_group(bold('Other options'), "")
    other_opts(opts)


def convert_subparser(subparsers):
    convert_parser = subparsers.add_parser(
        'convert', description=get_logo('Convert Kaptive results into different formats'),
        epilog=f'For more help, visit: {bold(_URL)}', add_help=False, formatter_class=argparse.RawTextHelpFormatter,
        help='Convert Kaptive results into different formats', usage="kaptive convert <db> <json> <format> [options]")
    opts = convert_parser.add_argument_group(
        bold('Inputs'), ""
        # "\n - Note: If you used --is-seqs during the run, make sure to provide the same fasta file here"
    )
    opts.add_argument('db', type=get_database, metavar='db path/keyword', help='Kaptive database path or keyword')
    opts.add_argument('input', help='Kaptive JSON lines file or - for stdin', type=argparse.FileType('rt'),
                      metavar='json')
    opts.add_argument('format', metavar='format',
                      choices=['json', 'tsv', 'loci', 'genes', 'proteins', 'png', 'svg'],
                      help='Output format\n'
                           ' - json: JSON lines format (same as input but optionally filtered)\n'
                           ' - tsv: Tab-separated values (results table)\n'
                           ' - loci: Locus nucleotide sequence in fasta format\n'
                           ' - proteins: Locus proteins in fasta format\n - genes: Locus genes in fasta format\n'
                           ' - png: Locus plot in PNG format\n - svg: Locus plot in SVG format'
                      )
    opts = convert_parser.add_argument_group(bold('Filter options'),
                                             "\n - Note: filters take precedence in descending order")
    opts.add_argument('-r', '--regex', metavar='', type=re.compile,
                      help='Python regular-expression to select JSON lines (default: All)')
    opts.add_argument('-l', '--loci', metavar='', nargs='+',
                      help='Space-separated list to filter locus names (default: All)')
    opts.add_argument('-s', '--samples', metavar='', nargs='+',
                      help='Space-separated list to filter sample names (default: All)')

    opts = convert_parser.add_argument_group(bold('Output options'), "")
    opts.add_argument('-o', '--out', metavar='', default=sys.stdout, type=argparse.FileType('at'),
                      help='Output file to write/append results to (default: stdout)\n'
                           ' - Note: Only for text formats, figures will be written to files')
    opts.add_argument('-d', '--outdir', metavar='', type=check_dir,
                      help='Output directory for converted results\n'
                           ' - Note: This forces the output to be written to files\n'
                           '         If used with locus, proteins or genes, one file will be written per sample'
                      )
    opts = convert_parser.add_argument_group(bold('Database options'), "")
    db_opts(opts)
    opts = convert_parser.add_argument_group(bold('Other options'), "")
    other_opts(opts)


def db_opts(opts: argparse.ArgumentParser):
    opts.add_argument('--locus-regex', type=re.compile, metavar='',
                      help=f'Python regular-expression to match locus names in db source note')
    opts.add_argument('--type-regex', type=re.compile, metavar='',
                      help=f'Python regular-expression to match locus types in db source note')
    opts.add_argument('--filter', type=re.compile, metavar='',
                      help='Python regular-expression to select loci to include in the database')


def output_opts(opts: argparse.ArgumentParser):
    opts.add_argument('-o', '--out', metavar='file', default=sys.stdout, type=argparse.FileType('at'),
                      help='Output file to write/append tabular results to (default: stdout)')
    opts.add_argument('-f', '--fasta', metavar='dir', nargs='?', default=None, const='.', type=check_dir,
                      help='Turn on fasta output, defaulting "{input}_kaptive_results.fna"\n'
                           ' - Optionally choose the output directory (default: cwd)')
    opts.add_argument('-j', '--json', metavar='file', nargs='?', default=None, const='kaptive_results.json',
                      type=argparse.FileType('at'),
                      help='Turn on JSON lines output\n'
                           ' - Optionally choose file (can be existing) (default: %(const)s)')
    opts.add_argument('-p', '--plot', metavar='dir', nargs='?', default=None, const='.', type=check_dir,
                      help='Turn on plot output, defaulting to "{input}_kaptive_results.{fmt}"\n'
                           ' - Optionally choose the output directory (default: cwd)')
    opts.add_argument('--plot-fmt', default='png', metavar='png,svg', choices=['png', 'svg'],
                      help='Format for locus plots (default: %(default)s)')
    # opts.add_argument('--bokeh', action='store_true', help='Plot locus figures using bokeh')
    opts.add_argument('--no-header', action='store_true', help='Suppress header line')
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


def plot_result(result: TypingResult, outdir: Path, fmt: str):
    """
    Convienence function to plot a TypingResult and save to file
    """
    ax = result.as_graphic_record().plot(figure_width=18)[0]  # type: 'matplotlib.axes.Axes'
    ax.set_title(  # Add title to figure
        f"{result.sample_name} {result.best_match} ({result.phenotype}) - {result.confidence}")
    ax.figure.savefig(outdir / f'{result.sample_name}_kaptive_results.{fmt}', bbox_inches='tight')
    ax.figure.clear()  # Close figure to prevent memory leak


# Main -----------------------------------------------------------------------------------------------------------------
def main():
    check_python_version(3, 9)
    check_biopython_version(1, 83)
    args = parse_args(sys.argv[1:])

    # Assembly mode ----------------------------------------------------------------------------------------------------
    if args.subparser_name == 'assembly':
        check_programs(['minimap2'], verbose=args.verbose)  # Check for minimap2
        args.db = Database.from_genbank(
            path=args.db, gene_threshold=args.gene_threshold, locus_filter=args.filter, verbose=args.verbose,
            load_seq=True, locus_regex=args.locus_regex, type_regex=args.type_regex)

        write_headers(args.out, args.no_header, args.subparser_name, args.debug)
        for sample in args.input:
            if (result := typing_pipeline(
                    sample, args.db, args.threads, args.min_cov, args.score_metric, args.weight_metric,
                    args.max_other_genes, args.percent_expected_genes, args.allow_below_threshold, args.debug,
                    args.verbose

            )):
                args.out.write(result.as_table(args.debug))  # type: 'TextIOWrapper'
                if args.fasta:  # Create and write locus pieces to fasta file
                    (args.fasta / f'{result.sample_name}_kaptive_results.fna').write_text(result.as_fasta())
                if args.json:  # type: 'TextIOWrapper'
                    args.json.write(dumps(result.as_dict()) + '\n')
                if args.plot:  # type: 'Path'
                    plot_result(result, args.plot, args.plot_fmt)

    # Extract mode -----------------------------------------------------------------------------------------------------
    elif args.subparser_name == 'extract':
        from kaptive.database import parse_database, name_from_record
        if args.format == "genbank":
            for record in SeqIO.parse(args.db, 'genbank'):
                locus_name, type_name = name_from_record(record, args.locus_regex, args.type_regex)
                record.name = locus_name  # Set the name of the record to the locus name
                if args.filter and not args.filter.search(locus_name):
                    continue
                if args.outdir:
                    (args.outdir / f'{locus_name.replace("/", "_")}.gbk').write_text(record.format('genbank'))
                else:
                    args.out.write(record.format('genbank'))  # Do we need to add a newline?
        else:  # load seqs == True if format == 'loci' else False
            for locus in parse_database(
                    args.db, args.filter, args.format == 'loci', locus_regex=args.locus_regex, type_regex=args.type_regex):
                if args.format == 'loci':
                    format_func, ext = locus.as_fasta, 'fna'
                elif args.format == 'genes':
                    format_func, ext = locus.as_gene_fasta, 'ffn'
                elif args.format == 'proteins':
                    format_func, ext = locus.as_protein_fasta, 'faa'
                else:
                    quit_with_error(f"Invalid format: {args.format}")
                if args.outdir:  # Perform a replacement "/" with "_" here as some O-loci names have '/'
                    (args.outdir / f'{locus.name.replace("/", "_")}.{ext}').write_text(format_func())
                else:
                    args.out.write(format_func())

    # Convert mode -----------------------------------------------------------------------------------------------------
    elif args.subparser_name == 'convert':
        args.db = Database.from_genbank(  # Load database in memory, we don't need to load the full sequences (False)
            args.db, args.filter, False, verbose=args.verbose, locus_regex=args.locus_regex, type_regex=args.type_regex)
        if args.format in ['png', 'svg'] and not args.outdir:  # Check if output directory is required and not set
            args.outdir = Path('.')  # Set output directory to current directory
        for line in args.input:
            if args.regex and not args.regex.search(line):
                return
            try:
                if (d := loads(line)):
                    if args.samples and d['sample_name'] not in args.samples:
                        return
                    if args.loci and d['best_match'] not in args.loci:
                        return
                    result = TypingResult.from_dict(d, args.db)
                    if args.format == 'loci':  # Create and write locus pieces to fasta file
                        if args.outdir:
                            (args.outdir / f'{result.sample_name}_kaptive_results.fna').write_text(result.as_fasta())
                        else:
                            args.out.write(result.as_fasta())
                    elif args.format == 'genes':
                        if args.outdir:
                            (args.outdir / f'{result.sample_name}_kaptive_results.ffn').write_text(
                                result.as_gene_fasta())
                        else:
                            args.out.write(result.as_gene_fasta())
                    elif args.format == 'proteins':
                        if args.outdir:
                            (args.outdir / f'{result.sample_name}_kaptive_results.faa').write_text(
                                result.as_protein_fasta())
                        else:
                            args.out.write(result.as_protein_fasta())
                    elif args.format == 'json':
                        if args.outdir:
                            (args.outdir / 'kaptive_results.json').write_text(dumps(result.as_dict()))
                        else:
                            args.out.write(dumps(result.as_dict()) + '\n')
                    elif args.format in ['png', 'svg']:
                        plot_result(result, args.outdir, args.format)
            except Exception as e:
                quit_with_error(f"Invalid JSON: {e}\n{line}")

