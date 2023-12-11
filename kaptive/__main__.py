#!/usr/bin/env python
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
import sys
import re
import argparse
from concurrent.futures import ProcessPoolExecutor as Pool
from itertools import chain  # repeat

# Silence Biopython warnings
import warnings
from Bio import BiopythonWarning

warnings.simplefilter('ignore', BiopythonWarning)

from kaptive.version import __version__
from kaptive.log import bold, bold_cyan, log, quit_with_error
from kaptive.misc import check_python_version, check_programs, LOGO, check_cpus, check_dir, check_file
from kaptive.database import Database, LOCUS_REGEX, TYPE_REGEX, extract, get_database, DB_KEYWORDS
from kaptive.assembly import type_assembly, Assembly
from kaptive.result import parse_results

# Constants -----------------------------------------------------------------------------------------------------------
ASSEMBLY_HEADERS = [
    'Assembly', 'Best match locus', 'Best match type', 'Problems', 'Coverage', 'Identity', 'Length discrepancy',
    'Expected genes in locus', 'Expected genes in locus, details', 'Missing expected genes',
    'Other genes in locus', 'Other genes in locus, details', 'Expected genes outside locus',
    'Expected genes outside locus, details', 'Other genes outside locus', 'Other genes outside locus, details'
]
ASSEMBLY_EXTRA_HEADERS = [
    'Truncated proteins', 'Extra genes', 'Contigs', 'Pieces', 'Pieces, details', 'Score', 'Zscore',
    'Score metric', 'Weighting', 'Minimap2 preset', 'Minimap2 args', 'Gene threshold', 'Gene distance',
    'Minimum coverage', 'Not used for scoring', 'All scores'
]


# Functions -----------------------------------------------------------------------------------------------------------
def parse_args(a):
    parser = argparse.ArgumentParser(
        description=bold_cyan(LOGO + f"\n{'In silico serotyping' : ^43}"), usage="%(prog)s <command>", add_help=False,
        prog="kaptive", formatter_class=argparse.RawDescriptionHelpFormatter, epilog=f"%(prog)s version {__version__}")

    subparsers = parser.add_subparsers(title=bold('Command'), dest='subparser_name', metavar="")
    assembly_subparser(subparsers)
    extract_subparser(subparsers)
    convert_subparser(subparsers)
    opts = parser.add_argument_group(bold('Other options'), '')
    other_opts(opts)

    if len(a) == 0:
        parser.print_help(sys.stderr)
        quit_with_error(f'Please specify a command; choose from {{assembly,extract,convert}}')
    elif '-v' in a or '--version' in a:
        print(__version__)
        sys.exit(0)
    elif a[0] in ['-h', '--help']:  # General help message
        parser.print_help(sys.stderr)
        sys.exit(0)
    else:  # Any number of arguments, first should be a subparser
        subparser = next((v for k, v in subparsers.choices.items() if k == a[0]), None)
        if subparser is not None:
            if len(a) == 1:  # Subparser help message
                subparser.print_help(sys.stderr)
                quit_with_error(f'Insufficient arguments for kaptive {a[0]}')
            elif '-h' in a or '--help' in a:  # Subparser help message if -h or --help is in the arguments
                subparser.print_help(sys.stderr)
                sys.exit(0)
        else:
            parser.print_help(sys.stderr)
            quit_with_error(f'Unknown command "{a[0]}"; choose from {{assembly,extract,convert}}')

    return parser.parse_args(a)


def assembly_subparser(subparsers):
    assembly_parser = subparsers.add_parser(
        'assembly', description=bold_cyan(LOGO + f"\n{'In silico serotyping of assemblies' : ^43}"),
        epilog=f'kaptive assembly v{__version__}', add_help=False, formatter_class=argparse.RawTextHelpFormatter,
        help='In silico serotyping of assemblies', usage="kaptive assembly <db> <assembly> [<assembly> ...] [options]")
    opts = assembly_parser.add_argument_group(bold('Inputs'), "")
    opts.add_argument('db', type=get_database, help='Kaptive database path or keyword')
    opts.add_argument('assembly', nargs='+', type=Assembly.from_path, help='Assemblies in fasta(.gz) format')
    opts = assembly_parser.add_argument_group(bold('Output options'), "")
    opts.add_argument('-o', '--out', metavar='', default=sys.stdout, type=argparse.FileType('w'),
                      help='Output file (default: stdout)')
    opts.add_argument('--fasta', metavar='path', nargs='?', default=None, const='.', type=check_dir,
                      help='Output results to fasta, optionally pass output directory (default: current directory)')
    opts.add_argument('--json', metavar='prefix', nargs='?', default=None, const='kaptive_results.json',
                      type=argparse.FileType('at', encoding='utf-8'),
                      help='Output results to json, optionally pass file name (default: %(const)s)')
    opts = assembly_parser.add_argument_group(bold('Alignment options'), "")
    opts.add_argument('-t', '--threads', type=check_cpus, default=1, metavar='',
                      help="Number of threads for minimap2 (default: %(default)s)\n"
                           " - Note, with --mp it is recommended to use 1 thread per sample")
    opts.add_argument('--args', metavar='', default='',
                      help='Additional arguments for minimap2 (default: %(default)s)')
    opts.add_argument('--preset', help='Preset for minimap2 (default: None)', metavar='',
                      choices=['map-pb', 'map-ont', 'map-hifi', 'ava-pb', 'ava-ont', 'asm5', 'asm10', 'asm20',
                               'splice', 'splice:hq', 'sr'])
    opts = assembly_parser.add_argument_group(bold('Scoring options'), "")
    opts.add_argument("--score", type=str, default='alignment_score', metavar='',
                      help="Alignment metric to use for scoring (default: %(default)s)")
    opts.add_argument("--weight", type=str, metavar='',
                      help="Weighting for scoring metric (default: None)\n"
                           " - locus_length: length of the locus\n"
                           " - genes_expected: # of genes expected in the locus\n"
                           " - genes_found: # of genes found in the locus\n"
                           " - ~ Any attribute of class Alignment ~")
    opts = assembly_parser.add_argument_group(bold('Locus reconstruction options'), "")
    opts.add_argument("--gene-distance", type=int, default=10000, metavar='',
                      help="Distance in bp between genes on the same contig to be considered part of the locus "
                           "(default: %(default)s)")
    opts.add_argument("--gene-threshold", type=float, default=95.0, metavar='',
                      help="Mark genes as low identity if they are below this threshold (default: %(default)s)")
    opts.add_argument('--min-cov', type=float, required=False, default=50.0, metavar='',
                      help='Minimum %%coverage for gene alignment to be used for scoring (default: %(default)s)')

    opts = assembly_parser.add_argument_group(bold('Database options'), "")
    db_opts(opts)
    opts = assembly_parser.add_argument_group(bold('Extra options'), "")
    opts.add_argument('-@', '--mp', const=8, nargs='?', type=int, metavar='#',
                      help="Process multiple samples in parallel, optionally pass max workers (default: %(const)s)")
    opts.add_argument('--no-header', action='store_true', help='Do not print header line')
    opts.add_argument('--debug', action='store_true', help='Append debug columns to table output')
    opts = assembly_parser.add_argument_group(bold('Other options'), "")
    other_opts(opts)


def extract_subparser(subparsers):
    extract_parser = subparsers.add_parser(
        'extract', description=bold_cyan(LOGO + f"\n{'Extract entries from a Kaptive database' : ^43}"),
        epilog=f'kaptive extract v{__version__}', add_help=False, formatter_class=argparse.RawTextHelpFormatter,
        help='Extract entries from a Kaptive database', usage="kaptive extract <db> <format> [options]")
    opts = extract_parser.add_argument_group(bold('Inputs'), "\n - Note: combine with --filter to select loci")
    opts.add_argument('db', help='Kaptive database in genbank format', type=get_database)
    opts.add_argument('format', choices=['loci', 'genes', 'proteins', 'gbk', 'gff', 'ids'], metavar='',
                      help='Format to extract database\n - loci: Loci (fasta nucleotide)\n'
                           ' - genes: Genes (fasta nucleotide)\n - proteins: Proteins (fasta amino acid)\n'
                           ' - gbk: Genbank format\n - gff: GFF in NCBI format\n - ids: List of Locus IDs')
    opts = extract_parser.add_argument_group(bold('Output options'), "")
    opts.add_argument('-o', '--out', metavar='', default=sys.stdout, type=argparse.FileType('w'),
                      help='Output file (default: stdout)')
    opts = extract_parser.add_argument_group(bold('Database options'), "")
    db_opts(opts)
    opts = extract_parser.add_argument_group(bold('Other options'), "")
    other_opts(opts)


def convert_subparser(subparsers):
    convert_parser = subparsers.add_parser(
        'convert', description=bold_cyan(LOGO + f"\n{'Convert Kaptive results into different formats' : ^43}"),
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
    opts.add_argument('-o', '--out', metavar='', default=sys.stdout, type=argparse.FileType('w'),
                      help='Output file (default: stdout)')
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
    opts.add_argument('--locus-regex', type=re.compile, default=LOCUS_REGEX, metavar='',
                      help='Pattern to match locus names in db source note, (default: %(default)s)')
    opts.add_argument('--type-regex', type=re.compile, default=TYPE_REGEX, metavar='',
                      help='Pattern to match locus types in db source note, (default: %(default)s)')
    opts.add_argument('--filter', type=re.compile, metavar='',
                      help='Pattern to select loci to include in the database (default: All)')
    opts.add_argument('--is-seqs', type=check_file, metavar='',
                      help='Fasta file of IS element sequences to include in the database (default: None)')


def other_opts(opts: argparse.ArgumentParser):
    opts.add_argument('-V', '--verbose', action='store_true', help='Print debug messages to stderr')
    opts.add_argument('-v', '--version', help='Show version number and exit', metavar='')
    opts.add_argument('-h', '--help', help='Show this help message and exit', metavar='')


# Main -----------------------------------------------------------------------------------------------------------------
def main():
    check_python_version()
    args = parse_args(sys.argv[1:])

    if args.subparser_name == 'assembly':
        check_programs(['minimap2'], verbose=args.verbose)
        # Load database in memory, we don't need to load the full sequences
        args.db = Database.from_genbank(args.db, args.is_seqs, args.locus_regex, args.type_regex, args.filter,
                                        load_seq=False)
        log(f"Loaded:\t{args.db}", args.verbose)
        if not args.no_header:
            args.out.write(
                '\t'.join(ASSEMBLY_HEADERS + ASSEMBLY_EXTRA_HEADERS if args.debug else ASSEMBLY_HEADERS) + '\n'
            )

        if not args.mp:  # No parallel processing
            [type_assembly(assembly, args) for assembly in args.assembly]
        else:  # Parallel processing
            log(f"Processing {len(args.assembly)} assemblies in parallel", args.verbose)
            with Pool(args.mp) as pool:
                [pool.submit(type_assembly, assembly, args) for assembly in args.assembly]

    elif args.subparser_name == 'extract':
        extract(args)

    elif args.subparser_name == 'convert':
        # Load database in memory, we don't need to load the full sequences
        args.db = Database.from_genbank(args.db, args.is_seqs, args.locus_regex, args.type_regex, args.filter,
                                        load_seq=False)
        log(f"Loaded:\t{args.db}", args.verbose)
        for result_file in args.json:
            for result in parse_results(result_file, args.db, args.regex, args.samples, args.loci):
                if args.format == 'json':
                    args.out.write(result.as_json())
                elif args.format == 'tsv':
                    args.out.write('\t'.join(result.as_list()) + '\n')
                elif args.format == 'locus':
                    args.out.write(result.as_fasta())
                elif args.format == 'genes':
                    args.out.write(result.as_gene_fasta())
                elif args.format == 'proteins':
                    args.out.write(result.as_protein_fasta())
