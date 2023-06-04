#!/usr/bin/env python
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
import os
import sys
from pathlib import Path
from argparse import SUPPRESS

from kaptive.version import __version__
from kaptive.log import bold
from kaptive.misc import get_ascii_art, check_file, check_python_version, check_programs, check_dir
from kaptive.help_formatter import MyHelpFormatter, MyParser
from kaptive.map import run_map
from kaptive.database import Database


# from kaptive.align import run_align


def parse_args(args):
    description = 'R|' + get_ascii_art() + '\n' + \
                  bold('kaptive-mapper: IS mapping software for Illumina sequencing data') + '\n'
    parser = MyParser(description=description, formatter_class=MyHelpFormatter, add_help=False)

    subparsers = parser.add_subparsers(title='Commands', dest='subparser_name')
    map_subparser(subparsers)

    longest_choice_name = max(len(c) for c in subparsers.choices)

    subparsers.help = 'R|'
    for choice, choice_parser in subparsers.choices.items():
        padding = ' ' * (longest_choice_name - len(choice))
        subparsers.help += choice + ': ' + padding
        d = choice_parser.description
        subparsers.help += d[0].lower() + d[1:]  # don't capitalise the first letter
        subparsers.help += '\n'

    help_args = parser.add_argument_group('Help')
    help_args.add_argument('-h', '--help', action='help', default=SUPPRESS,
                           help='Show this help message and exit')
    help_args.add_argument('--version', action='version', version='kaptive-mapper v' + __version__,
                           help="Show program's version number and exit")

    # If no arguments were used, print the base-level help which lists possible commands.
    if len(args) == 0:
        parser.print_help(file=sys.stderr)
        sys.exit(1)

    return parser.parse_args(args)


def map_subparser(subparsers):
    group = subparsers.add_parser('map', description='in silico serotyping of reads',
                                  formatter_class=MyHelpFormatter, add_help=False)

    required_args = group.add_argument_group('Required arguments')
    required_args.add_argument('reference', metavar='[genbank]',
                               help='Kaptive reference database in genbank format',
                               type=lambda p: check_file(subparsers, Path(p)))
    required_args.add_argument('reads', nargs='+', metavar='[R_{1,2}.fastq(.gz) ...]',
                               help='Paired end reads for analysing, can be gzipped',
                               type=lambda p: check_file(subparsers, Path(p)))

    output_args = group.add_argument_group('Output arguments', 'Arguments for output files')
    # output_args.add_argument('-o', '--out', type=argparse.FileType('w'), metavar='[stdout]', default=sys.stdout,
    #                          help='Output file for tabular results')
    output_args.add_argument('--keep_sam', action='store_true', default=False, help='Keep SAM file')
    output_args.add_argument('--vcf', type=lambda p: check_dir(subparsers, Path(p).absolute()), metavar='[cwd]',
                             default=None, nargs='?', const=Path.cwd(),
                             help='Output VCF files from the SNP calling to the specified directory\n'
                                  'Default is current directory')

    mapping_args = group.add_argument_group('Mapping arguments', 'Arguments for mapping reads to the reference')
    mapping_args.add_argument('--preset', help='Mapping preset for minimap2', default='sr',
                              choices=['map-pb', 'map-ont', 'map-hifi', 'ava-pb', 'ava-ont', 'asm5', 'asm10', 'asm20',
                                       'splice', 'splice:hq', 'sr'])

    ise_args = group.add_argument_group('IS Element arguments', 'Arguments for mapping IS elements')
    ise_args.add_argument('-i', '--ise', type=lambda p: check_file(subparsers, Path(p)),
                          help='File containing IS element sequences in fasta format;\n'
                               'This will activate IS element mapping')
    ise_args.add_argument('--ise_range_merge_tolerance', metavar='[10]', default=10, type=int,
                          help='Number of bases to extend ISE alignment ranges to chain and overlap with reads that '
                               'mapped to the locus; results in fewer fragmented ISE alignment ranges, especially '
                               'when using the "sr" preset')

    snp_args = group.add_argument_group('SNP calling arguments', 'Arguments for calling SNPs in the locus')
    snp_args.add_argument('--min_depth', metavar='[10]', default=10, type=int)
    snp_args.add_argument('--min_qual', metavar='[40]', default=40, type=int)
    snp_args.add_argument('--min_mq', metavar='[60]', default=60, type=int)
    snp_args.add_argument('--csq_phase', choices=['a', 'm', 'r', 'R', 's'], default='m',
                          help="how to handle unphased heterozygous genotypes:\n"
                               "  a: take GTs as is, create haplotypes regardless of phase (0/1 → 0|1)\n"
                               "  m: merge all GTs into a single haplotype (0/1 → 1, 1/2 → 1)\n"
                               "  r: require phased GTs, throw an error on unphased heterozygous GTs\n"
                               "  R: create non-reference haplotypes if possible (0/1 → 1|1, 1/2 → 1|2)\n"
                               "  s: skip unphased heterozygous GTs\n"),
    snp_args.add_argument('--all_csqs', action='store_true', default=False,
                          help='Report all consequences per gene, rather than the first "harmful" consequence')
    snp_args.add_argument('--synonymous', action='store_true', default=False, help='Report synonymous SNPs')

    ref_args = group.add_argument_group('Reference database arguments')
    ref_settings(ref_args)
    other_args = group.add_argument_group('Other arguments')
    other_settings(other_args)


def align_subparser(subparsers):
    group = subparsers.add_parser('align', description='in silico serotyping of assemblies',
                                  formatter_class=MyHelpFormatter, add_help=False)

    required_args = group.add_argument_group('Required arguments')
    required_args.add_argument('reference', metavar='[genbank]',
                               help='Kaptive reference database in genbank format',
                               type=lambda p: check_file(subparsers, Path(p)))
    required_args.add_argument('assemblies', nargs='+', metavar='[assembly.fasta ...]',
                               help='Assembly fasta for analysing, can be gzipped',
                               type=lambda p: check_file(subparsers, Path(p)))

    output_parser = subparsers.add_argument_group('Output options')
    output_parser.add_argument('--no_seq_out', action='store_true',
                               help='Suppress output files of sequences matching locus')
    output_parser.add_argument('--no_table', action='store_true', help='Suppress output of tab-delimited table')
    output_parser.add_argument('--no_json', action='store_true', help='Suppress output of JSON file')
    output_parser.add_argument('--scores', action='store_true', help='Add locus scores to JSON file (debugging)')

    score_parser = subparsers.add_argument_group('Scoring options')
    score_parser.add_argument('--start_end_margin', type=int, required=False, default=10,
                              help='Missing bases at the ends of locus allowed in a perfect match.')
    score_parser.add_argument('--min_gene_cov', type=float, required=False, default=90.0,
                              help='minimum required %% coverage for genes')
    score_parser.add_argument('--min_gene_id', type=float, required=False, default=80.0,
                              help='minimum required %% identity for genes')
    score_parser.add_argument('--low_gene_id', type=float, required=False, default=95.0,
                              help='genes with a %% identity below this value will be flagged as low identity')
    score_parser.add_argument('--min_assembly_piece', type=int, required=False, default=100,
                              help='minimum locus matching assembly piece to return')
    score_parser.add_argument('--gap_fill_size', type=int, required=False, default=1500,
                              help=f'when separate parts of the assembly are found within this distance, '
                                   f'they will be merged')

    ref_args = group.add_argument_group('Reference database arguments')
    ref_settings(ref_args)
    other_args = group.add_argument_group('Other arguments')
    other_settings(other_args)


def ref_settings(ref_args):
    ref_args.add_argument('--locus_label', type=str, default='',
                          help='In the Genbank file, the source feature must have a note  identifying the locus name, '
                               'starting with this label followed by a colon (e.g. /note="K locus: KL1")')
    ref_args.add_argument('--type_label', type=str, default='',
                          help='In the Genbank file, the source feature must have a note identifying the type name, '
                               'starting with this label followed by a colon (e.g. /note="K type: K1")')


def other_settings(other_args):
    other_args.add_argument('-t', '--threads', type=lambda t: min([os.cpu_count(), int(t)]), default=os.cpu_count(),
                            help="Number of threads; defaults to all")
    other_args.add_argument('-v', '--version', action='version', version=f'%(prog)s {__version__}')
    other_args.add_argument('-h', '--help', action='help', help='Show this help message and exit')
    other_args.add_argument('--verbose', action='store_true', default=False, help='Print debug messages to stderr')
    other_args.add_argument('--header', action='store_true', default=False, help='Print header with output')


def main():
    check_python_version()
    args = parse_args(sys.argv[1:])

    # args = parse_args(
    #     f'map '
    #     f'{Path.cwd() / "reference_database/Klebsiella_k_locus_primary_reference.gbk"} '
    #     f'{" ".join(str(i) for i in Path("/Users/tom/Bioinformatics/kaptive-mapper/test_data/").glob("Klebs-K3_*"))} '
    #     f'-i {Path.cwd() / "is.fna"} --preset sr --verbose'.split()
    # )

    db = Database(args.reference, args.locus_label, args.type_label)

    if args.subparser_name == 'map':
        check_programs(['minimap2', 'samtools', 'bcftools', 'paftools.js'], args.verbose)
        run_map(args, db)

    # elif args.subparser_name == 'align':
    #     check_programs(['blastn', 'makeblastdb'], args.verbose)
    #     run_align(args, db)
    #
    # elif args.subparser_name == 'type':
    #     check_programs(['minimap2'], args.verbose)
    #     run_type(args, db)


if __name__ == '__main__':
    main()
