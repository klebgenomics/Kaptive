# import argparse

# # Functions ------------------------------------------------------------------------------------------------------------
# def assembly_subparser(subparsers):
#     assembly_parser = subparsers.add_parser(
#         'assembly', description=_bold('In silico serotyping of assemblies'),
#         epilog=_EPILOG, add_help=False, formatter_class=argparse.RawTextHelpFormatter,
#         help='In silico serotyping of assemblies',
#         usage=f"{Resources.package} assembly <database> <fasta> [<fasta> ...] [options]"
#     )
#     opts = assembly_parser.add_argument_group(_bold('Inputs'), "")
#     opts.add_argument('database', metavar='database path/keyword', help='Database path or keyword')
#     opts.add_argument('input', nargs='+', metavar='fasta', help='Assemblies in fasta(.gz|.xz|.bz2) format')
#     opts = assembly_parser.add_argument_group(_bold('Output options'), "\nNote, text outputs accept '-' for stdout")
#     # Note these are different to the convert output options as TSV is the main output and fna is the main fasta output
#     opts.add_argument('-o', '--out', metavar='', default=sys.stdout, type=argparse.FileType('at'),
#                       help='Output file to write/append tabular results to (default: stdout)')
#     opts.add_argument('-f', '--fasta', metavar='', nargs='?', default=None, const='.', type=write_to_file_or_directory,
#                       help='Turn on fasta output\n'
#                            'Accepts a single file or a directory (default: cwd)')
#     opts.add_argument('-j', '--json', metavar='', nargs='?', default=None, const=f'{Resources.package}_results.json',
#                       type=argparse.FileType('at'),
#                       help='Turn on JSON lines output\n'
#                            'Optionally choose file (can be existing) (default: %(const)s)')
#     opts.add_argument('-s', '--scores', metavar='', nargs='?', default=None, const=sys.stdout,
#                       type=argparse.FileType('at'),
#                       help='Dump locus score matrix to tsv (serotyping will not be performed!)\n'
#                            'Optionally choose file (can be existing) (default: stdout)')
#     other_fmt_opts(opts)
#     opts = assembly_parser.add_argument_group(_bold('Scoring options'), "")
#     opts.add_argument('--min-cov', type=float, required=False, default=50.0, metavar='',
#                       help='Minimum gene %%coverage (blen/q_len*100) to be used for scoring (default: %(default)s)')
#     opts.add_argument("--score-metric", metavar='', default=0, type=int, choices=range(4),
#                       help="Metric for scoring each locus (default: %(default)s)\n"
#                            "  0: AS (alignment score of genes found)\n"
#                            "  1: mlen (matching bases of genes found)\n"
#                            "  2: blen (aligned bases of genes found)\n"
#                            "  3: q_len (query length of genes found)")
#     opts.add_argument("--weight-metric", metavar='', default=3, type=int, choices=range(6),
#                       help="Weighting for the 1st stage of the scoring algorithm (default: %(default)s)\n"
#                            "  0: No weighting\n"
#                            "  1: Number of genes found\n"
#                            "  2: Number of genes expected\n"
#                            "  3: Proportion of genes found\n"
#                            "  4: blen (aligned bases of genes found)\n"
#                            "  5: q_len (query length of genes found)")
#     opts.add_argument('--n-best', type=int, default=2, metavar='',
#                       help='Number of best loci from the 1st round of scoring to be\n'
#                            'fully aligned to the assembly (default: %(default)s)')

#     opts = assembly_parser.add_argument_group(_bold('Confidence options'), "")
#     opts.add_argument("--gene-threshold", type=float, metavar='',
#                       help="Species-level locus gene identity threshold (default: database specific)")
#     opts.add_argument("--max-other-genes", type=int, metavar='', default=1,
#                       help="Typeable if <= other genes (default: %(default)s)")
#     opts.add_argument("--percent-expected", type=float, metavar='', default=50,
#                       help="Typeable if >= %% expected genes (default: %(default)s)")
#     opts.add_argument("--below-threshold", type=bool, default=False, metavar='',
#                       help="Typeable if any genes are below threshold (default: %(default)s)")

#     opts = assembly_parser.add_argument_group(_bold('Database options'), "")
#     opts.add_argument('--filter', type=re.compile, metavar='',
#                       help='Python regular-expression to select loci to include in the database')
#     opts = assembly_parser.add_argument_group(_bold('Other options'), "")
#     other_opts(opts)
#     opts.add_argument('-t', '--threads', type=int, default=Resources.available_cpus, metavar='',
#                       help="Number of alignment threads or 0 for all available (default: 0)")


# def convert_subparser(subparsers):
#     convert_parser = subparsers.add_parser(
#         'convert', description=_bold(f'Convert {Resources.package} results into different formats'),
#         epilog=_EPILOG, add_help=False, formatter_class=argparse.RawTextHelpFormatter,
#         help=f'Convert {Resources.package} results into different formats',
#         usage=f"{Resources.package} convert <database> <json> [formats] [options]")
#     opts = convert_parser.add_argument_group(_bold('Inputs'), "")
#     opts.add_argument('database', metavar='database path/keyword', help='Database path or keyword')
#     opts.add_argument('inputf', help=f'{Resources.package} JSON lines file or - for stdin', type=argparse.FileType('rt'),
#                       metavar='json')
#     opts = convert_parser.add_argument_group(_bold('Formats'), "\nNote, text outputs accept '-' for stdout")
#     opts.add_argument('-t', '--tsv', metavar='', nargs='?', default=None, const='-', type=write_to_file_or_directory,
#                       help='Convert to tabular format in file (default: stdout)')
#     opts.add_argument('-j', '--json', metavar='', nargs='?', default=None, const='-', type=write_to_file_or_directory,
#                       help='Convert to JSON lines format in file (default: stdout)')
#     fmt_opts(opts)
#     other_fmt_opts(opts)
#     opts = convert_parser.add_argument_group(_bold('Filter options'),
#                                              "\nNote, filters take precedence in descending order")
#     opts.add_argument('-l', '--loci', metavar='', nargs='+',
#                       help='Space-separated list to filter locus names (default: All)')
#     opts.add_argument('-s', '--samples', metavar='', nargs='+',
#                       help='Space-separated list to filter sample names (default: All)')
#     opts = convert_parser.add_argument_group(_bold('Database options'), "")

#     opts = convert_parser.add_argument_group(_bold('Other options'), "")
#     other_opts(opts)