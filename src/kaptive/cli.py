# import sys
# import re
# import argparse
# from io import TextIOWrapper
# from importlib.metadata import version
# from shutil import which

# from kaptive.utils import Resources, write_to_file_or_directory

# # Functions ------------------------------------------------------------------------------------------------------------
# def _bold(text: str) -> str:
#     return f'\033[1m{text}\033[0m'

# # Constants -----------------------------------------------------------------------------------------------------------
# __version__ = version(Resources.package)
# _EPILOG = f'For more help, visit: {_bold("https://klebgenomics.github.io/Kaptive/")}'


# # Subparsers --




# def extract_subparser(subparsers):
#     extract_parser = subparsers.add_parser(
#         'extract', description=_bold(f'Extract entries from a {Resources.package} database'),
#         epilog=_EPILOG, add_help=False, formatter_class=argparse.RawTextHelpFormatter,
#         help=f'Extract entries from a {Resources.package} database', 
#         usage=f"{Resources.package} extract <database> [formats] [options]"
#     )
#     opts = extract_parser.add_argument_group(_bold('Inputs'), "\nNote, combine with --filter to select loci")
#     opts.add_argument('database', metavar='database path/keyword', help='Database path or keyword')
#     opts = extract_parser.add_argument_group(_bold('Formats'), "\nNote, text outputs accept '-' for stdout")
#     fmt_opts(opts)
#     opts = extract_parser.add_argument_group(_bold('Database options'), "")

#     opts.add_argument('--filter', type=re.compile, metavar='',
#                       help='Python regular-expression to select loci to include in the database')
#     opts = extract_parser.add_argument_group(_bold('Other options'), "")
#     other_opts(opts)


# def fmt_opts(opts: argparse.ArgumentParser):
#     """Format opts shared by convert and extract"""
#     opts.add_argument('--fna', metavar='', nargs='?', default=None, const='.', type=write_to_file_or_directory,
#                       help='Convert to locus nucleotide sequences in fasta format\n'
#                            'Accepts a single file or a directory (default: cwd)')
#     opts.add_argument('--ffn', metavar='', nargs='?', default=None, const='.', type=write_to_file_or_directory,
#                       help='Convert to locus gene nucleotide sequences in fasta format\n'
#                            'Accepts a single file or a directory (default: cwd)')
#     opts.add_argument('--faa', metavar='', nargs='?', default=None, const='.', type=write_to_file_or_directory,
#                       help='Convert to locus gene protein sequences in fasta format\n'
#                            'Accepts a single file or a directory (default: cwd)')


# def other_fmt_opts(opts: argparse.ArgumentParser):
#     """Format opts shared by convert and assembly"""
#     opts.add_argument('--no-header', action='store_true', help='Suppress header line')


# def other_opts(opts: argparse.ArgumentParser):
#     opts.add_argument('-v', '--version', help='Show version number and exit', action='version')
#     opts.add_argument('-h', '--help', help='Show this help message and exit', metavar='')


# # Main -----------------------------------------------------------------------------------------------------------------
# def main():
#     parser = argparse.ArgumentParser(
#         description=_bold('In silico serotyping'), usage="%(prog)s <command>", add_help=False,
#         prog=Resources.package, formatter_class=argparse.RawDescriptionHelpFormatter, epilog=_EPILOG
#     )

#     subparsers = parser.add_subparsers(title=_bold('Command'), dest='subparser_name', metavar="")
#     assembly_subparser(subparsers)
#     extract_subparser(subparsers)
#     convert_subparser(subparsers)
#     opts = parser.add_argument_group(_bold('Other options'), '')
#     other_opts(opts)
    
#     args = sys.argv[1:]
    
#     if len(args) == 0:  # No arguments, print help message
#         parser.print_help(sys.stderr)
#         raise RuntimeError(f'Please specify a command; choose from {{assembly,extract,convert}}')
#     if any(x in args for x in {'-v', '--version'}):  # Version message
#         print(__version__)
#         sys.exit(0)
#     if subparser := subparsers.choices.get(args[0], None):  # Check if the first argument is a subparser
#         if len(args) == 1:  # Subparser help message
#             subparser.print_help(sys.stderr)
#             raise RuntimeError(f'Insufficient arguments for {Resources.package} {a[0]}')
#         if any(x in args[1:] for x in {'-h', '--help'}):  # Subparser help message
#             subparser.print_help(sys.stderr)
#             sys.exit(0)
#     elif any(x in args for x in {'-h', '--help'}):  # Help message
#         parser.print_help(sys.stderr)
#         sys.exit(0)
#     else:  # Unknown command
#         parser.print_help(sys.stderr)
#         raise RuntimeError(f'Unknown command "{args[0]}"; choose from {{assembly,extract,convert}}')
    
#     args = parser.parse_args(args)

#     # Assembly mode ----------------------------------------------------------------------------------------------------
#     if args.subparser_name == 'assembly':
        
#         if not which('minimap2'):
#             raise RuntimeError(f'Minimap2 not installed; please install minimap2')
        
#         from kaptive.genomeassembly import typing_pipeline, write_headers
#         from kaptive.database import load_database

#         args.db = load_database(args.db, args.gene_threshold, locus_filter=args.filter)

#         write_headers(args.scores or args.out, args.no_header, args.scores)

#         for assembly in args.input:
#             if result := typing_pipeline(assembly, args.db, args.threads, args.score_metric, args.weight_metric,
#                                          args.min_cov, args.n_best, args.max_other_genes, args.percent_expected,
#                                          args.below_threshold, args.scores):
#                 result.write(args.out, args.json, args.fasta, None, None)

#     # Extract mode -----------------------------------------------------------------------------------------------------
#     elif args.subparser_name == 'extract':
#         from kaptive.database import parse_database
#         for locus in parse_database(args.db, args.filter, args.fna, args.faa):
#             locus.write(args.fna, args.ffn, args.faa)

#     # Convert mode -----------------------------------------------------------------------------------------------------
#     elif args.subparser_name == 'convert':
#         from kaptive.database import load_database
#         from kaptive.genomeassembly import parse_result, write_headers

#         args.db = load_database(args.db)

#         write_headers(args.tsv, args.no_header)

#         for line in args.input:
#             if result := parse_result(line, args.db, args.regex, args.samples, args.loci):
#                 result.write(args.tsv, args.json, args.fna, args.ffn, args.faa)

#     # Cleanup ----------------------------------------------------------------------------------------------------------
#     for attr in vars(args):  # Close all open files in the args namespace if they aren't sys.stdout or sys.stdin
#         if (x := getattr(args, attr, None)) and isinstance(x, TextIOWrapper) and x not in {sys.stdout, sys.stdin}:
#             x.close()  # Close the file
