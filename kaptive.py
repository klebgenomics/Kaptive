#!/usr/bin/env python3
"""
Copyright 2018 Ryan Wick (rrwick@gmail.com)
https://github.com/katholt/Kaptive

Kaptive is a tool which reports information about the K and O types for Klebsiella genome
assemblies. It will help a user to decide whether their Klebsiella sample has a known or novel
locus type, and if novel, how similar it is to a known type.

This script needs the following input files to run:
* A Genbank file with known locus types
* One or more assemblies in FASTA format

Example command:
kaptive.py -a path/to/assemblies/*.fasta -k k_loci_refs.gbk -o output_directory

For each input assembly file, Kaptive will identify the closest known locus type and report
information about the corresponding locus genes.

It generates the following output files:
* A FASTA file for each input assembly with the nucleotide sequences matching the closest locus type
* A table summarising the results for all input assemblies

Character codes indicate problems with the locus match:
* `?` indicates that the match was not in a single piece, possible due to a poor match or
      discontiguous assembly
* `-` indicates that genes expected in the locus were not found
* `+` indicates that extra genes were found in the locus
* `*` indicates that one or more expected genes was found but with low identity

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
import argparse
import sys
import os
import multiprocessing
import subprocess
import json
import fcntl
import gzip
import copy
import random
from collections import OrderedDict
from Bio import SeqIO

__version__ = '0.7.3'


def main():
    """Script execution starts here."""
    args = get_argument_parser().parse_args()

    check_for_blast()
    check_files_exist(args.assembly + [args.k_refs] + [args.allelic_typing])
    check_assembly_format(args.assembly)
    fix_paths(args)

    output_table = not args.no_table
    output_json = not args.no_json

    temp_dir = make_temp_dir(args)
    k_ref_seqs, gene_seqs, k_ref_genes = parse_genbank(args.k_refs, temp_dir, args.locus_label)
    all_gene_dict = {}
    for gene_list in k_ref_genes.values():
        for gene in gene_list:
            all_gene_dict[gene.full_name] = gene

    k_refs = load_k_locus_references(k_ref_seqs, k_ref_genes)
    type_gene_names = get_type_gene_names(args.allelic_typing)

    if output_table:
        create_table_file(args.out, type_gene_names)
    json_list = []

    for fasta_file in args.assembly:
        assembly = Assembly(fasta_file)
        best_k = get_best_k_type_match(assembly, k_ref_seqs, k_refs, args.threads)
        if best_k is None:
            type_gene_results = {}
            best_k = KLocus('None', '', [])
        else:
            find_assembly_pieces(assembly, best_k, args)
            assembly_pieces_fasta = save_assembly_pieces_to_file(best_k, assembly, args.out)
            type_gene_results = type_gene_search(assembly_pieces_fasta, type_gene_names, args)
            if args.no_seq_out and assembly_pieces_fasta is not None:
                os.remove(assembly_pieces_fasta)
            protein_blast(assembly, best_k, gene_seqs, args)
            check_name_for_o1_o2(best_k)

        output(args.out, assembly, best_k, args, type_gene_names, type_gene_results,
               json_list, output_table, output_json, all_gene_dict)

    if output_json:
        write_json_file(args.out, json_list)

    clean_up(k_ref_seqs, gene_seqs, temp_dir)
    sys.exit(0)


def get_argument_parser():
    """Specifies the command line arguments required by the script."""
    parser = argparse.ArgumentParser(description='Kaptive',
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    add_arguments_to_parser(parser)
    return parser


def add_arguments_to_parser(parser):
    parser.add_argument('--version', action='version', version='Kaptive v' + __version__,
                        help="Show Kaptive's version number and exit")
    parser.add_argument('-a', '--assembly', nargs='+', type=str, required=True,
                        help='FASTA file(s) for assemblies')
    parser.add_argument('-k', '--k_refs', type=str, required=True,
                        help='GenBank file with reference loci')
    parser.add_argument('-g', '--allelic_typing', type=str, required=False,
                        help='SRST2-formatted FASTA file of allelic typing genes to include in '
                             'results')
    parser.add_argument('-o', '--out', type=str, required=False, default='./kaptive_results',
                        help='Output directory/prefix')
    parser.add_argument('-v', '--verbose', action='store_true',
                        help='Display detailed information about each assembly in stdout')
    parser.add_argument('-t', '--threads', type=int, required=False,
                        default=min(multiprocessing.cpu_count(), 4),
                        help='The number of threads to use for the BLAST searches')
    parser.add_argument('--no_seq_out', action='store_true',
                        help='Suppress output files of sequences matching locus')
    parser.add_argument('--no_table', action='store_true',
                        help='Suppress output of tab-delimited table')
    parser.add_argument('--no_json', action='store_true',
                        help='Suppress output of JSON file')
    parser.add_argument('--start_end_margin', type=int, required=False, default=10,
                        help='Missing bases at the ends of locus allowed in a perfect match.')
    parser.add_argument('--min_gene_cov', type=float, required=False, default=90.0,
                        help='minimum required %% coverage for genes')
    parser.add_argument('--min_gene_id', type=float, required=False, default=80.0,
                        help='minimum required %% identity for genes')
    parser.add_argument('--low_gene_id', type=float, required=False, default=95.0,
                        help='genes with a %% identity below this value will be flagged as low '
                             'identity')
    parser.add_argument('--min_assembly_piece', type=int, required=False, default=100,
                        help='minimum locus matching assembly piece to return')
    parser.add_argument('--gap_fill_size', type=int, required=False, default=100,
                        help='when separate parts of the assembly are found within this distance, '
                             'they will be merged')
    parser.add_argument('--locus_label', type=str, required=False,
                        default='automatically determined',
                        help='In the Genbank file, the source feature must have a note '
                             'identifying the locus name, starting with this label followed by '
                             'a colon (e.g. /note="K locus: K1")')


def check_for_blast():
    """Checks to make sure the required BLAST+ tools are available."""
    if not find_program('makeblastdb'):
        quit_with_error('could not find makeblastdb tool (part of BLAST+)')
    if not find_program('blastn'):
        quit_with_error('could not find blastn tool (part of BLAST+)')
    if not find_program('tblastn'):
        quit_with_error('could not find tblastn tool (part of BLAST+)')


def find_program(name):
    """Checks to see if a program exists."""
    process = subprocess.Popen(['which', name], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    out, err = process.communicate()
    return bool(out) and not bool(err)


def fix_paths(args):
    """
    Changes the paths given by the user to absolute paths, which are easier to work with later.
    Also creates the output directory, if necessary.
    """
    args.assembly = [os.path.abspath(x) for x in args.assembly]
    args.k_refs = os.path.abspath(args.k_refs)
    if args.out[-1] == '/':
        args.out += 'kaptive_results'
    args.out = os.path.abspath(args.out)
    out_dir = os.path.dirname(args.out)
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)


def make_temp_dir(args):
    """Makes the temporary directory, if necessary. Returns the temp directory path."""
    temp_dir_name = 'kaptive_temp_' + str(os.getpid()) + '_' + str(random.randint(0, 999999))
    temp_dir = os.path.join(os.path.dirname(args.out), temp_dir_name)
    if not os.path.exists(temp_dir):
        os.makedirs(temp_dir)
    return temp_dir


def clean_up(k_ref_seqs, gene_seqs, temp_dir):
    """
    Deletes the temporary FASTA files. If the temp directory is then empty, it is deleted too.
    """
    try:
        os.remove(k_ref_seqs)
    except OSError:
        pass
    try:
        os.remove(gene_seqs)
    except OSError:
        pass
    try:
        if not os.listdir(temp_dir):
            os.rmdir(temp_dir)
    except FileNotFoundError:
        pass


def parse_genbank(genbank, temp_dir, locus_label):
    """
    This function reads the input Genbank file and produces two temporary FASTA files: one with the
    loci nucleotide sequences and one with the gene sequences.
    It returns the file paths for these two FASTA files along with a dictionary that links genes to
    loci.
    """
    k_ref_genes = {}
    k_ref_seqs_filename = os.path.join(temp_dir, 'temp_k_ref_seqs.fasta')
    gene_seqs_filename = os.path.join(temp_dir, 'temp_gene_seqs.fasta')
    k_ref_seqs = open(k_ref_seqs_filename, 'wt')
    gene_seqs = open(gene_seqs_filename, 'wt')
    if locus_label == 'automatically determined':
        locus_label = find_locus_label(genbank)
    else:
        check_locus_label(genbank, locus_label)
    for record in SeqIO.parse(genbank, 'genbank'):
        k_locus_name = ''
        for feature in record.features:
            if feature.type == 'source' and 'note' in feature.qualifiers:
                for note in feature.qualifiers['note']:
                    if note.startswith(locus_label):
                        k_locus_name = get_locus_name_from_note(note, locus_label)
                    elif note.startswith('Extra genes'):
                        k_locus_name = note.replace(':', '').replace(' ', '_')
        if k_locus_name in k_ref_genes:
            quit_with_error('Duplicate reference locus name: ' + k_locus_name)
        k_ref_genes[k_locus_name] = []

        # Extra genes are only used for the gene search, not the nucleotide search.
        if not k_locus_name.startswith('Extra_genes'):
            k_ref_seqs.write('>' + k_locus_name + '\n')
            k_ref_seqs.write(add_line_breaks_to_sequence(str(record.seq), 60))

        gene_num = 1
        for feature in record.features:
            if feature.type == 'CDS':
                gene = Gene(k_locus_name, gene_num, feature, record.seq)
                k_ref_genes[k_locus_name].append(gene)
                gene_num += 1
                gene_seqs.write(gene.get_fasta())
    k_ref_seqs.close()
    gene_seqs.close()
    return k_ref_seqs_filename, gene_seqs_filename, k_ref_genes


def find_locus_label(genbank):
    """
    Automatically finds the label for the locus sequences. The Genbank file must have exactly one
    possible label that is present in a note qualifier in the source feature for every record. If
    not, Kaptive will quit with an error.
    """
    possible_locus_labels = set()
    for record in SeqIO.parse(genbank, 'genbank'):
        for feature in record.features:
            if feature.type == 'source' and 'note' in feature.qualifiers:
                for note in feature.qualifiers['note']:
                    if ':' in note:
                        possible_locus_labels.add(note.split(':')[0].strip())
                    if '=' in note:
                        possible_locus_labels.add(note.split('=')[0].strip())
    if not possible_locus_labels:
        quit_with_error('None of the records contain a valid locus label')
    available_locus_labels = possible_locus_labels.copy()
    for record in SeqIO.parse(genbank, 'genbank'):
        locus_labels = set()
        for feature in record.features:
            if feature.type == 'source' and 'note' in feature.qualifiers:
                for note in feature.qualifiers['note']:
                    if ':' in note:
                        locus_labels.add(note.split(':')[0].strip())
                    if '=' in note:
                        locus_labels.add(note.split('=')[0].strip())
        if any(x == 'Extra genes' for x in locus_labels):
            continue
        if not locus_labels:
            quit_with_error('no possible locus labels were found for ' + record.name)
        previous_labels = available_locus_labels.copy()
        available_locus_labels = available_locus_labels.intersection(locus_labels)
        if not available_locus_labels:
            error_message = record.name + ' does not have a locus label matching the previous ' \
                            'records\n'
            error_message += 'Previous record labels: ' + ', '.join(list(previous_labels)) + '\n'
            error_message += 'Labels in ' + record.name + ': ' + ', '.join(list(locus_labels))
            quit_with_error(error_message)
    if len(available_locus_labels) > 1:
        error_message = 'multiple possible locus labels were found: ' + \
                        ', '.join(list(available_locus_labels)) + '\n'
        error_message += 'Please use the --locus_label option to specify which to use'
        quit_with_error(error_message)
    return list(available_locus_labels)[0]


def check_locus_label(genbank, locus_label):
    """
    Makes sure that every record in the Genbank file contains a note in the source feature
    beginning with the given label.
    """
    for record in SeqIO.parse(genbank, 'genbank'):
        found_label = False
        for feature in record.features:
            if feature.type == 'source' and 'note' in feature.qualifiers:
                for note in feature.qualifiers['note']:
                    if note.startswith(locus_label):
                        k_locus_name = get_locus_name_from_note(note, locus_label)
                        if k_locus_name:
                            found_label = True
        if not found_label:
            error_message = record.name + ' is missing a locus label\n'
            error_message += 'The source feature must have a note qualifier beginning with "' + \
                             locus_label + ':" followed by the locus name'
            quit_with_error(error_message)


def get_locus_name_from_note(full_note, locus_label):
    """
    Extracts the part of the note following the label (and any colons, spaces or equals signs).
    """
    locus_name = full_note[len(locus_label):].strip()
    while locus_name.startswith(':') or locus_name.startswith(' ') or \
            locus_name.startswith('='):
        locus_name = locus_name[1:]
    return locus_name


def check_files_exist(filenames):
    """Checks to make sure each file in the given list exists."""
    for filename in filenames:
        if filename is not None:
            check_file_exists(filename)


def check_file_exists(filename):
    """Checks to make sure the single given file exists."""
    if not os.path.isfile(filename):
        quit_with_error('could not find ' + filename)


def check_assembly_format(filenames):
    """Tries to load each assembly and shows an error if it did not successfully load."""
    for assembly in filenames:
        fasta = load_fasta(assembly)
        if len(fasta) < 1:
            quit_with_error('invalid FASTA file: ' + assembly)
        for record in fasta:
            header, seq = record
            if len(seq) == 0:
                quit_with_error('invalid FASTA file (contains a zero-length sequence): ' +
                                assembly)


def quit_with_error(message):
    """Displays the given message and ends the program's execution."""
    print('Error:', message, file=sys.stderr)
    sys.exit(1)


def get_best_k_type_match(assembly, k_refs_fasta, k_refs, threads):
    """
    Searches for all known locus types in the given assembly and returns the best match.
    Best match is defined as the locus type for which the largest fraction of the locus has a BLAST
    hit to the assembly. In cases of a tie, the mean identity of the locus type BLAST hits are used
    to determine the best.
    """
    for k_ref in k_refs.values():
        k_ref.clear()
    blast_hits = get_blast_hits(assembly.fasta, k_refs_fasta, threads)

    for hit in blast_hits:
        if hit.qseqid not in k_refs:
            quit_with_error('BLAST hit (' + hit.qseqid + ') not found in locus references')
        k_refs[hit.qseqid].add_blast_hit(hit)
    best_k_ref = None
    best_cov = 0.0
    for k_ref in k_refs.values():
        cov = k_ref.get_coverage()
        if cov > best_cov:
            best_cov = cov
            best_k_ref = k_ref
        elif cov == best_cov and best_k_ref and \
                k_ref.get_mean_blast_hit_identity() > best_k_ref.get_mean_blast_hit_identity():
            best_k_ref = k_ref
    if best_k_ref is not None:
        best_k_ref.clean_up_blast_hits()
    return copy.copy(best_k_ref)


def type_gene_search(assembly_pieces_fasta, type_gene_names, args):
    if not type_gene_names or not args.allelic_typing:
        return {}
    if not assembly_pieces_fasta:
        return {x: None for x in type_gene_names}

    makeblastdb(assembly_pieces_fasta)
    all_gene_blast_hits = get_blast_hits(assembly_pieces_fasta, args.allelic_typing, args.threads,
                                         type_genes=True)
    clean_blast_db(assembly_pieces_fasta)

    # Filter out small hits.
    all_gene_blast_hits = [x for x in all_gene_blast_hits
                           if x.query_cov >= args.min_gene_cov and x.pident >= args.min_gene_id]

    type_gene_results = {}
    for gene_name in type_gene_names:
        blast_hits = sorted([x for x in all_gene_blast_hits if x.gene_name == gene_name],
                            reverse=True, key=lambda z: z.bitscore)
        if not blast_hits:
            hit = None
        else:
            perfect_match = None
            for hit in blast_hits:
                if hit.pident == 100.0 and hit.query_cov == 100.0:
                    perfect_match = hit
                    break
            if perfect_match:
                hit = perfect_match
                hit.result = str(perfect_match.allele_number)
            else:
                hit = blast_hits[0]
                hit.result = str(blast_hits[0].allele_number) + '*'
        type_gene_results[gene_name] = hit

    return type_gene_results


def find_assembly_pieces(assembly, k_locus, args):
    """
    This function uses the BLAST hits in the given locus type to find the corresponding pieces of
    the given assembly. It saves its results in the KLocus object.
    """
    if not k_locus.blast_hits:
        return
    assembly_pieces = [x.get_assembly_piece(assembly) for x in k_locus.blast_hits]
    merged_pieces = merge_assembly_pieces(assembly_pieces)
    length_filtered_pieces = [x for x in merged_pieces if x.get_length() >= args.min_assembly_piece]
    if not length_filtered_pieces:
        return
    k_locus.assembly_pieces = fill_assembly_piece_gaps(length_filtered_pieces, args.gap_fill_size)

    # Now check to see if the biggest assembly piece seems to capture the whole locus. If so, this
    # is an ideal match.
    biggest_piece = sorted(k_locus.assembly_pieces, key=lambda z: z.get_length(), reverse=True)[0]
    start = biggest_piece.earliest_hit_coordinate()
    end = biggest_piece.latest_hit_coordinate()
    if good_start_and_end(start, end, k_locus.get_length(), args.start_end_margin):
        k_locus.assembly_pieces = [biggest_piece]

    # If it isn't the ideal case, we still want to check if the start and end of the locus were
    # found in the same contig. If so, fill all gaps in between so we include the entire
    # intervening sequence.
    else:
        earliest, latest, same_contig_and_strand = k_locus.get_earliest_and_latest_pieces()
        k_start = earliest.earliest_hit_coordinate()
        k_end = latest.latest_hit_coordinate()
        if good_start_and_end(k_start, k_end, k_locus.get_length(), args.start_end_margin) and \
           same_contig_and_strand:
            gap_filling_piece = AssemblyPiece(assembly, earliest.contig_name, earliest.start,
                                              latest.end, earliest.strand)
            k_locus.assembly_pieces = merge_assembly_pieces(k_locus.assembly_pieces +
                                                            [gap_filling_piece])
    k_locus.identity = get_mean_identity(k_locus.assembly_pieces)


def protein_blast(assembly, k_locus, gene_seqs, args):
    """
    Conducts a BLAST search of all known locus proteins. Stores the results in the KLocus
    object.
    """
    hits = get_blast_hits(assembly.fasta, gene_seqs, args.threads, genes=True)
    hits = [x for x in hits if x.query_cov >= args.min_gene_cov and x.pident >= args.min_gene_id]

    best_hits = []
    for expected_gene in k_locus.gene_names:
        best_hit = get_best_hit_for_query(hits, expected_gene, k_locus)
        if best_hit is not None:
            best_hits.append(best_hit)
    best_hits = sorted(best_hits, key=lambda x: x.bitscore, reverse=True)
    for best_hit in best_hits:
        if best_hit in hits:
            hits = cull_conflicting_hits(best_hit, hits)

    expected_hits = []
    for expected_gene in k_locus.gene_names:
        best_hit = get_best_hit_for_query(hits, expected_gene, k_locus)
        if not best_hit:
            k_locus.missing_expected_genes.append(expected_gene)
        else:
            best_hit.over_identity_threshold = best_hit.pident >= args.low_gene_id
            expected_hits.append(best_hit)
            hits = [x for x in hits if x is not best_hit]
            hits = cull_conflicting_hits(best_hit, hits)
    other_hits = cull_all_conflicting_hits(hits)

    k_locus.expected_hits_inside_locus = [x for x in expected_hits
                                          if x.in_assembly_pieces(k_locus.assembly_pieces)]
    k_locus.expected_hits_outside_locus = [x for x in expected_hits
                                           if not x.in_assembly_pieces(k_locus.assembly_pieces)]
    k_locus.other_hits_inside_locus = [x for x in other_hits
                                       if x.in_assembly_pieces(k_locus.assembly_pieces)]
    k_locus.other_hits_outside_locus = [x for x in other_hits
                                        if not x.in_assembly_pieces(k_locus.assembly_pieces)]


def create_table_file(output_prefix, type_gene_names):
    """
    Creates the table file and writes a header line if necessary.
    If the file already exists and the header line is correct, then it does nothing (to allow
    multiple independent processes to append to the file).
    """
    table_path = output_prefix + '_table.txt'

    # If the table already exists, we don't need to do anything.
    if os.path.isfile(table_path):
        with open(table_path, 'r') as existing_table:
            first_line = existing_table.readline().strip()
            if first_line.startswith('Assembly\tBest match locus'):
                return

    headers = ['Assembly',
               'Best match locus',
               'Match confidence',
               'Problems',
               'Coverage',
               'Identity',
               'Length discrepancy',
               'Expected genes in locus',
               'Expected genes in locus, details',
               'Missing expected genes',
               'Other genes in locus',
               'Other genes in locus, details',
               'Expected genes outside locus',
               'Expected genes outside locus, details',
               'Other genes outside locus',
               'Other genes outside locus, details']

    if type_gene_names:
        headers += type_gene_names

    with open(table_path, 'w') as table:
        table.write('\t'.join(headers))
        table.write('\n')


def get_type_gene_names(type_genes_fasta):
    gene_names = []
    if type_genes_fasta:
        gene_names = set()
        with open(type_genes_fasta, 'rt') as type_genes_db:
            for line in type_genes_db:
                if line.startswith('>'):
                    try:
                        gene_names.add(line.split('>')[1].split('__')[1])
                    except IndexError:
                        quit_with_error(type_genes_fasta + ' not formatted as an SRST2 database '
                                                           'FASTA file')
        if not gene_names:
            quit_with_error(type_genes_fasta + ' not formatted as an SRST2 database FASTA file')
        gene_names = sorted(list(gene_names))
    return gene_names


def check_name_for_o1_o2(k_locus):
    """
    This function has special logic for dealing with the O1/O2 locus. If the wbbY and wbbZ genes
    are both found, then we call the locus O2 (instead of O1/O2). If neither are found, then we
    call the locus O1.
    """
    if not (k_locus.name == 'O1/O2v1' or k_locus.name == 'O1/O2v2'):
        return
    other_gene_names = [x.qseqid for x in k_locus.other_hits_outside_locus]
    both_present = ('Extra_genes_wbbY/wbbZ_01_wbbY' in other_gene_names and
                    'Extra_genes_wbbY/wbbZ_02_wbbZ' in other_gene_names)
    both_absent = ('Extra_genes_wbbY/wbbZ_01_wbbY' not in other_gene_names and
                   'Extra_genes_wbbY/wbbZ_02_wbbZ' not in other_gene_names)
    if both_present:
        k_locus.name = k_locus.name.replace('O1/O2', 'O1')
    elif both_absent:
        k_locus.name = k_locus.name.replace('O1/O2', 'O2')


def output(output_prefix, assembly, k_locus, args, type_gene_names, type_gene_results,
           json_list, output_table, output_json, all_gene_dict):
    """
    Writes a line to the output table describing all that we've learned about the given locus and
    writes to stdout as well.
    """
    uncertainty_chars = k_locus.get_match_uncertainty_chars()

    try:
        expected_in_locus_per = 100.0 * len(k_locus.expected_hits_inside_locus) / \
            len(k_locus.gene_names)
        expected_out_locus_per = 100.0 * len(k_locus.expected_hits_outside_locus) / \
            len(k_locus.gene_names)
        expected_genes_in_locus_str = str(len(k_locus.expected_hits_inside_locus)) + ' / ' + \
            str(len(k_locus.gene_names)) + ' (' + float_to_str(expected_in_locus_per) + '%)'
        expected_genes_out_locus_str = str(len(k_locus.expected_hits_outside_locus)) + ' / ' + \
            str(len(k_locus.gene_names)) + ' (' + float_to_str(expected_out_locus_per) + '%)'
        missing_per = 100.0 * len(k_locus.missing_expected_genes) / len(k_locus.gene_names)
        missing_genes_str = str(len(k_locus.missing_expected_genes)) + ' / ' + \
            str(len(k_locus.gene_names)) + ' (' + float_to_str(missing_per) + '%)'
    except ZeroDivisionError:
        expected_genes_in_locus_str, expected_genes_out_locus_str, missing_genes_str = '', '', ''

    output_to_stdout(assembly, k_locus, args.verbose, type_gene_names, type_gene_results,
                     uncertainty_chars, expected_genes_in_locus_str, expected_genes_out_locus_str,
                     missing_genes_str)
    if output_table:
        output_to_table(output_prefix, assembly, k_locus, type_gene_names, type_gene_results,
                        uncertainty_chars, expected_genes_in_locus_str,
                        expected_genes_out_locus_str)
    if output_json:
        add_to_json(assembly, k_locus, type_gene_names, type_gene_results, json_list,
                    uncertainty_chars, all_gene_dict)


def output_to_table(output_prefix, assembly, k_locus, type_gene_names, type_gene_results,
                    uncertainty_chars, expected_genes_in_locus_str, expected_genes_out_locus_str):
    line = [assembly.name,
            k_locus.name,
            k_locus.get_match_confidence(),
            uncertainty_chars,
            k_locus.get_coverage_string(),
            k_locus.get_identity_string(),
            k_locus.get_length_discrepancy_string(),
            expected_genes_in_locus_str,
            get_gene_info_string(k_locus.expected_hits_inside_locus),
            ';'.join(k_locus.missing_expected_genes),
            str(len(k_locus.other_hits_inside_locus)),
            get_gene_info_string(k_locus.other_hits_inside_locus),
            expected_genes_out_locus_str,
            get_gene_info_string(k_locus.expected_hits_outside_locus),
            str(len(k_locus.other_hits_outside_locus)),
            get_gene_info_string(k_locus.other_hits_outside_locus)]

    for gene_name in type_gene_names:
        hit = type_gene_results[gene_name]
        line.append('-' if not hit else hit.result)

    table_path = output_prefix + '_table.txt'
    table = open(table_path, 'at')
    table.write('\t'.join(line))
    table.write('\n')
    table.close()


def add_to_json(assembly, k_locus, type_gene_names, type_gene_results, json_list,
                uncertainty_chars, all_gene_dict):
    json_record = OrderedDict()
    json_record['Assembly name'] = assembly.name

    match_dict = OrderedDict()
    match_dict['Locus name'] = k_locus.name
    match_dict['Match confidence'] = k_locus.get_match_confidence()

    reference_dict = OrderedDict()
    reference_dict['Length'] = len(k_locus.seq)
    reference_dict['Sequence'] = k_locus.seq
    match_dict['Reference'] = reference_dict
    json_record['Best match'] = match_dict

    problems = OrderedDict()
    problems['Locus assembled in multiple pieces'] = str('?' in uncertainty_chars)
    problems['Missing genes in locus'] = str('-' in uncertainty_chars)
    problems['Extra genes in locus'] = str('+' in uncertainty_chars)
    problems['At least one low identity gene'] = str('*' in uncertainty_chars)
    json_record['Problems'] = problems

    blast_results = OrderedDict()
    blast_results['Coverage'] = k_locus.get_coverage_string()
    blast_results['Identity'] = k_locus.get_identity_string()
    blast_results['Length discrepancy'] = k_locus.get_length_discrepancy_string()
    assembly_pieces = []
    for i, piece in enumerate(k_locus.assembly_pieces):
        assembly_piece = OrderedDict()
        assembly_piece['Contig name'] = piece.contig_name
        assembly_piece['Contig start position'] = piece.start + 1
        assembly_piece['Contig end position'] = piece.end
        assembly_piece['Contig strand'] = piece.strand
        piece_seq = piece.get_sequence()
        assembly_piece['Length'] = len(piece_seq)
        assembly_piece['Sequence'] = piece_seq
        assembly_pieces.append(assembly_piece)
    blast_results['Locus assembly pieces'] = assembly_pieces
    json_record['blastn result'] = blast_results

    expected_genes_in_locus = {x.qseqid: x for x in k_locus.expected_hits_inside_locus}
    expected_hits_outside_locus = {x.qseqid: x for x in k_locus.expected_hits_outside_locus}
    other_hits_inside_locus = {x.qseqid: x for x in k_locus.other_hits_inside_locus}
    other_hits_outside_locus = {x.qseqid: x for x in k_locus.other_hits_outside_locus}

    k_locus_genes = []
    for gene in k_locus.genes:
        gene_dict = OrderedDict()
        gene_name = gene.full_name
        gene_dict['Name'] = gene_name
        if gene_name in expected_genes_in_locus:
            gene_dict['Result'] = 'Found in locus'
        elif gene_name in expected_hits_outside_locus:
            gene_dict['Result'] = 'Found outside locus'
        else:
            gene_dict['Result'] = 'Not found'
        gene_dict['Reference'] = gene.get_reference_info_json_dict()

        if gene_name in expected_genes_in_locus or gene_name in expected_hits_outside_locus:
            if gene_name in expected_genes_in_locus:
                hit = expected_genes_in_locus[gene_name]
            else:
                hit = expected_hits_outside_locus[gene_name]
            gene_dict['tblastn result'] = hit.get_blast_result_json_dict(assembly)
            gene_dict['Match confidence'] = hit.get_match_confidence()
        else:
            gene_dict['Match confidence'] = 'Not found'

        k_locus_genes.append(gene_dict)
    json_record['Locus genes'] = k_locus_genes

    extra_genes = OrderedDict()
    for gene_name, hit in other_hits_inside_locus.items():
        gene_dict = OrderedDict()
        gene = all_gene_dict[gene_name]
        gene_dict['Reference'] = gene.get_reference_info_json_dict()
        gene_dict['tblastn result'] = hit.get_blast_result_json_dict(assembly)
        extra_genes[gene_name] = gene_dict
    json_record['Other genes in locus'] = extra_genes

    other_genes = OrderedDict()
    for gene_name, hit in other_hits_outside_locus.items():
        gene_dict = OrderedDict()
        gene = all_gene_dict[gene_name]
        gene_dict['Reference'] = gene.get_reference_info_json_dict()
        gene_dict['tblastn result'] = hit.get_blast_result_json_dict(assembly)
        other_genes[gene_name] = gene_dict
    json_record['Other genes outside locus'] = other_genes

    if type_gene_names:
        allelic_typing = OrderedDict()
        for gene_name in type_gene_names:
            allelic_type = OrderedDict()
            if not type_gene_results[gene_name]:
                allelic_type['Allele'] = 'Not found'
            else:
                blast_hit = type_gene_results[gene_name]
                allele = blast_hit.result
                if allele.endswith('*'):
                    perfect_match = False
                    allele = allele[:-1]
                else:
                    perfect_match = True
                try:
                    allele = int(allele)
                except ValueError:
                    pass
                allelic_type['Allele'] = allele
                allelic_type['Perfect match'] = str(perfect_match)
                allelic_type['blastn result'] = blast_hit.get_blast_result_json_dict(assembly)
            allelic_typing[gene_name] = allelic_type
        json_record['Allelic_typing'] = allelic_typing

    json_list.append(json_record)


def write_json_file(output_prefix, json_list):
    json_filename = output_prefix + '.json'
    if not os.path.isfile(json_filename):
        with open(output_prefix + '.json', 'wt') as json_out:
            fcntl.flock(json_out, fcntl.LOCK_EX)
            json_out.write(json.dumps(json_list, indent=4))
            json_out.write('\n')
            fcntl.flock(json_out, fcntl.LOCK_UN)
    else:
        with open(output_prefix + '.json', 'r+t') as json_out:
            fcntl.flock(json_out, fcntl.LOCK_EX)
            file_data = json_out.read()
            try:
                existing_json_list = json.loads(file_data, object_pairs_hook=OrderedDict)
                json_list = existing_json_list + json_list
            except ValueError:
                pass
            json_out.seek(0)
            json_out.write(json.dumps(json_list, indent=4))
            json_out.write('\n')
            json_out.truncate()
            fcntl.flock(json_out, fcntl.LOCK_UN)


def output_to_stdout(assembly, k_locus, verbose, type_gene_names, type_gene_results,
                     uncertainty_chars, expected_genes_in_locus_str, expected_genes_out_locus_str,
                     missing_genes_str):
    if verbose:
        print()
        assembly_name_line = 'Assembly: ' + assembly.name
        print(assembly_name_line)
        print('-' * len(assembly_name_line))
        print('    Best match locus: ' + k_locus.name)
        print('    Match confidence: ' + k_locus.get_match_confidence())
        print('    Problems: ' + (uncertainty_chars if uncertainty_chars else 'None'))
        print('    Coverage: ' + k_locus.get_coverage_string())
        print('    Identity: ' + k_locus.get_identity_string())
        print('    Length discrepancy: ' + k_locus.get_length_discrepancy_string())
        print()
        print_assembly_pieces(k_locus.assembly_pieces)
        print_gene_hits('Expected genes in locus: ' + expected_genes_in_locus_str,
                        k_locus.expected_hits_inside_locus)
        print_gene_hits('Expected genes outside locus: ' + expected_genes_out_locus_str,
                        k_locus.expected_hits_outside_locus)
        print('    Missing expected genes: ' + missing_genes_str)
        for missing_gene in k_locus.missing_expected_genes:
            print('        ' + missing_gene)
        print()
        print_gene_hits('Other genes in locus: ' + str(len(k_locus.other_hits_inside_locus)),
                        k_locus.other_hits_inside_locus)
        print_gene_hits('Other genes outside locus: ' + str(len(k_locus.other_hits_outside_locus)),
                        k_locus.other_hits_outside_locus)

        for gene_name in type_gene_names:
            result = 'Not found' if not type_gene_results[gene_name] \
                else type_gene_results[gene_name].result
            print('    ' + gene_name + ' allele: ' + result)
        print()

    else:  # not verbose
        simple_output = assembly.name + ': ' + k_locus.name + uncertainty_chars
        for gene_name in type_gene_names:
            result = 'Not found' if not type_gene_results[gene_name] \
                else type_gene_results[gene_name].result
            simple_output += ', ' + gene_name + '=' + result
        print(simple_output)


def print_assembly_pieces(pieces):
    """This function prints assembly pieces nicely for verbose output."""
    print('    Locus assembly pieces:')
    if pieces:
        longest_header = max([len(x.get_nice_header()) for x in pieces])
        for piece in pieces:
            first_part = piece.get_nice_header()
            first_part = first_part.ljust(longest_header)
            second_part = piece.get_sequence_short()
            print('        ' + first_part + '  ' + second_part)
    print()


def print_gene_hits(title, hits):
    """This function prints gene hits nicely for verbose output."""
    print('    ' + title)
    if hits:
        longest_gene_name = max([len(x.qseqid) for x in hits])
        longest_contig_details = max([len(x.get_contig_details_string()) for x in hits])
        longest_coverage_details = max([len(x.get_coverage_details_string()) for x in hits])
        cov_space = max([x.query_cov for x in hits]) == 100.0
        id_space = max([x.pident for x in hits]) == 100.0
        spacing_1 = longest_gene_name + 2
        spacing_2 = spacing_1 + longest_contig_details + 2
        spacing_3 = spacing_2 + longest_coverage_details + 2
        for hit in hits:
            print('        ' + hit.get_aligned_string(spacing_1, spacing_2, spacing_3,
                                                      cov_space, id_space))
    print()


def float_to_str(float_in):
    """
    This function converts a float to a string in a special manner: if the float is an integer,
    the resulting string has no decimal point. Otherwise, one decimal point is used.
    """
    if float_in == int(float_in):
        return str(int(float_in))
    else:
        return '%.1f' % float_in


def get_blast_hits(database, query, threads, genes=False, type_genes=False):
    """Returns a list BlastHit objects for a search of the given query in the given database."""
    if genes:
        command = ['tblastn',
                   '-db_gencode', '11',  # bacterial translation table
                   '-seg', 'no']         # don't filter out low complexity regions
    else:
        command = ['blastn', '-task', 'blastn']
    command += ['-db', database, '-query', query, '-num_threads', str(threads), '-outfmt',
                '6 qseqid sseqid qstart qend sstart send evalue bitscore length pident qlen qseq '
                'sseq']
    process = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    out, err = process.communicate()
    out = convert_bytes_to_str(out)
    err = convert_bytes_to_str(err)
    if err:
        quit_with_error(command[0] + ' encountered an error:\n' + err)
    if process.returncode != 0:
        msg = command[0] + ' crashed!\n'

        # A known crash can occur with tblastn and recent versions of BLAST+ when multiple threads
        # are used. Check for this case and display an informative error message if so.
        version = get_blast_version(command[0])
        bad_version = (version == '2.4.0') or (version == '2.5.0') or (version == '2.6.0')
        if threads > 1 and bad_version:
            msg += '\nYou are using BLAST+ v' + version + ' which may crash when running with '
            msg += 'multiple threads.\n\n'
            msg += 'To avoid this issue, try one of the following:\n'
            msg += '  1) Use an unaffected version of BLAST+ (v2.3.0 or earlier should work)\n'
            msg += '  2) Run Kaptive with "--threads 1" (will probably be slower)\n'
        quit_with_error(msg)

    if genes:
        blast_hits = [GeneBlastHit(line) for line in line_iterator(out)]
    elif type_genes:
        blast_hits = [TypeGeneBlastHit(line) for line in line_iterator(out)]
    else:
        blast_hits = [BlastHit(line) for line in line_iterator(out)]
    return blast_hits


def get_blast_version(program):
    command = [program, '-version']
    process = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    out, err = process.communicate()
    out = convert_bytes_to_str(out)
    try:
        return out.split(': ')[1].split()[0].split('+')[0]
    except IndexError:
        return ''


def get_best_hit_for_query(blast_hits, query_name, k_locus):
    """
    Given a list of BlastHits, this function returns the best hit for the given query, based first
    on whether or not the hit is in the assembly pieces, then on bit score.
    It returns None if no BLAST hits match that query.
    """
    matching_hits = [x for x in blast_hits if x.qseqid == query_name]
    if matching_hits:
        return sorted(matching_hits,
                      key=lambda z: (z.in_assembly_pieces(k_locus.assembly_pieces), z.bitscore),
                      reverse=True)[0]
    else:
        return None


def cull_conflicting_hits(hit_to_keep, blast_hits):
    """
    This function returns a (potentially) reduced set of BLAST hits which excludes BLAST hits that
    overlap too much (same part of assembly) with the hit to keep.
    """
    return [x for x in blast_hits if not x.conflicts(hit_to_keep)]


def cull_all_conflicting_hits(blast_hits):
    """
    This function returns a (potentially) reduced set of BLAST hits where none of the remaining
    hits conflict.
    """
    blast_hits.sort(key=lambda x: x.bitscore, reverse=True)
    kept_hits = []
    while blast_hits:
        kept_hits.append(blast_hits.pop(0))
        blast_hits = cull_conflicting_hits(kept_hits[-1], blast_hits)
    return kept_hits


def merge_assembly_pieces(pieces):
    """
    Takes a list of AssemblyPiece objects and returns another list of AssemblyPiece objects where
    the overlapping pieces have been merged.
    """
    while True:
        merged_pieces = []
        merge_count = 0
        while pieces:
            merged_piece = pieces[0]
            unmerged = []
            for other_piece in pieces[1:]:
                combined = merged_piece.combine(other_piece)
                if not combined:
                    unmerged.append(other_piece)
                else:
                    merged_piece = combined
                    merge_count += 1
            merged_pieces.append(merged_piece)
            pieces = unmerged
        if merge_count == 0:
            break
        else:
            pieces = merged_pieces
    return merged_pieces


def fill_assembly_piece_gaps(pieces, max_gap_fill_size):
    """
    This function takes a list of assembly pieces, and if any of them are close enough to each
    other, the gap will be merged in.
    It assumes that all given pieces are from the same assembly.
    """
    pieces_by_contig_and_strand = {}
    fixed_pieces = []
    for piece in pieces:
        contig = piece.contig_name
        strand = piece.strand
        if (contig, strand) not in pieces_by_contig_and_strand:
            pieces_by_contig_and_strand[(contig, strand)] = []
        pieces_by_contig_and_strand[(contig, strand)].append(piece)
    for (contig, strand), pieces_in_contig_and_strand in pieces_by_contig_and_strand.items():
        gap_filling_pieces = []
        sorted_pieces = sorted(pieces_in_contig_and_strand, key=lambda x: x.start)
        max_end = sorted_pieces[0].end
        gaps = []
        for piece in sorted_pieces[1:]:
            if piece.start > max_end and piece.start - max_end <= max_gap_fill_size:
                gaps.append((max_end, piece.start))
            max_end = max(max_end, piece.end)
        assembly = sorted_pieces[0].assembly
        for gap in gaps:
            gap_filling_pieces.append(AssemblyPiece(assembly, contig, gap[0], gap[1], strand))
        before_merge = pieces_in_contig_and_strand + gap_filling_pieces
        filled_pieces = merge_assembly_pieces(before_merge)
        fixed_pieces += filled_pieces
    return fixed_pieces


def get_mean_identity(pieces):
    """Returns the mean identity (weighted by sequence length) for a list of assembly pieces."""
    identity_sum = 0.0
    length_sum = 0
    for piece in pieces:
        for hit in piece.blast_hits:
            length_sum += hit.length
            identity_sum += hit.length * hit.pident
    if identity_sum == 0.0:
        return 0.0
    else:
        return identity_sum / length_sum


def reverse_complement(seq):
    """Given a DNA sequences, this function returns the reverse complement sequence."""
    rev_comp = ''
    for i in reversed(range(len(seq))):
        rev_comp += complement_base(seq[i])
    return rev_comp


def complement_base(base):
    """Given a DNA base, this returns the complement."""
    forward = 'ATGCatgcRYSWKMryswkmBDHVbdhvNn.-?'
    reverse = 'TACGtacgYRSWMKyrswmkVHDBvhdbNn.-?N'
    return reverse[forward.find(base)]


def save_assembly_pieces_to_file(k_locus, assembly, output_prefix):
    """
    Creates a single FASTA file for all of the assembly pieces.
    Assumes all assembly pieces are from the same assembly.
    """
    if not k_locus.assembly_pieces:
        return None
    fasta_file_name = output_prefix + '_' + assembly.name + '.fasta'
    with open(fasta_file_name, 'w') as fasta_file:
        for piece in k_locus.assembly_pieces:
            fasta_file.write('>' + assembly.name + '_' + piece.get_header() + '\n')
            fasta_file.write(add_line_breaks_to_sequence(piece.get_sequence(), 60))
    return fasta_file_name


def add_line_breaks_to_sequence(sequence, length):
    """Wraps sequences to the defined length. All resulting sequences end in a line break."""
    seq_with_breaks = ''
    while len(sequence) > length:
        seq_with_breaks += sequence[:length] + '\n'
        sequence = sequence[length:]
    if sequence:
        seq_with_breaks += sequence
        seq_with_breaks += '\n'
    return seq_with_breaks


def line_iterator(string_with_line_breaks):
    """Iterates over a string containing line breaks, one line at a time."""
    prev_newline = -1
    while True:
        next_newline = string_with_line_breaks.find('\n', prev_newline + 1)
        if next_newline < 0:
            break
        yield string_with_line_breaks[prev_newline + 1:next_newline]
        prev_newline = next_newline


def load_k_locus_references(fasta, k_ref_genes):
    """Returns a dictionary of: key = locus name, value = KLocus object"""
    return {seq[0]: KLocus(seq[0], seq[1], k_ref_genes[seq[0]]) for seq in load_fasta(fasta)}


def load_fasta(filename):
    """Returns the names and sequences for the given fasta file."""
    fasta_seqs = []
    if get_compression_type(filename) == 'gz':
        open_func = gzip.open
    else:  # plain text
        open_func = open
    with open_func(filename, 'rt') as fasta_file:
        name = ''
        sequence = ''
        for line in fasta_file:
            line = line.strip()
            if not line:
                continue
            if line[0] == '>':  # Header line = start of new contig
                if name:
                    fasta_seqs.append((name.split()[0], sequence))
                    sequence = ''
                name = line[1:]
            else:
                sequence += line
        if name:
            fasta_seqs.append((name.split()[0], sequence))
    return fasta_seqs


def good_start_and_end(start, end, k_length, allowed_margin):
    """
    Checks whether the given start and end coordinates are within the accepted margin of error.
    """
    good_start = start <= allowed_margin
    good_end = end >= k_length - allowed_margin
    start_before_end = start < end
    return good_start and good_end and start_before_end


def get_gene_info_string(gene_hit_list):
    """Returns a single comma-delimited string summarising the gene hits in the given list."""
    return ';'.join([x.qseqid + ',' + str(x.pident) + '%' for x in gene_hit_list])


def is_contig_name_spades_format(contig_name):
    """
    Returns whether or not the contig name appears to be in the SPAdes/Velvet format.
    Example: NODE_5_length_150905_cov_4.42519
    """
    contig_name_parts = contig_name.split('_')
    return len(contig_name_parts) > 5 and contig_name_parts[0] == 'NODE' and \
        contig_name_parts[2] == 'length' and contig_name_parts[4] == 'cov'


def get_nice_contig_name(contig_name):
    """
    For a contig with a SPAdes/Velvet format, this function returns a simplified string that is
    just NODE_XX where XX is the contig number.
    For any other format, this function trims off everything following the first whitespace.
    """
    if is_contig_name_spades_format(contig_name):
        return 'NODE_' + contig_name.split('_')[1]
    else:
        return contig_name.split()[0]


class BlastHit(object):
    """
    Stores the BLAST hit output mostly verbatim. However, it does convert the BLAST ranges
    (1-based, inclusive end) to Python ranges (0-based, exclusive end).
    """
    def __init__(self, hit_string):
        parts = hit_string.split('\t')
        self.qseqid = parts[0]
        self.sseqid = parts[1]
        self.qstart = int(parts[2]) - 1
        self.qend = int(parts[3])
        self.sstart = int(parts[4])
        self.send = int(parts[5])
        if self.sstart <= self.send:
            self.strand = '+'
        else:
            self.sstart, self.send = self.send, self.sstart
            self.strand = '-'
        self.sstart -= 1
        self.evalue = float(parts[6])
        self.bitscore = float(parts[7])
        self.length = int(parts[8])
        self.pident = float(parts[9])
        self.query_cov = 100.0 * len(parts[11]) / float(parts[10])
        self.sseq = parts[12]

    def __repr__(self):
        return self.qseqid + ', ' + self.get_contig_details_string() + ', ' + \
               self.get_coverage_details_string() + ', ' + self.get_identity_details_string()

    def get_contig_details_string(self):
        """Returns a string describing the hit's range and strand in the contig."""
        return 'Contig: ' + get_nice_contig_name(self.sseqid) + ' (' + str(self.sstart) + '-' + \
               str(self.send) + ', ' + self.strand + ' strand)'

    def get_coverage_string(self):
        return '%.2f' % self.query_cov + '%'

    def get_coverage_details_string(self, extra_space=False):
        first_part = 'Cov: '
        second_part = self.get_coverage_string()
        if len(second_part) == 6 and extra_space:
            first_part += ' '
        return first_part + second_part

    def get_identity_string(self):
        return '%.2f' % self.pident + '%'

    def get_identity_details_string(self, extra_space=False):
        first_part = 'ID: '
        second_part = self.get_identity_string()
        if len(second_part) == 6 and extra_space:
            first_part += ' '
        return first_part + second_part

    def get_aligned_string(self, spacing_1, spacing_2, spacing_3, cov_space, id_space):
        """Returns a string describing the hit with spaced parts for alignment."""
        aligned_string = self.qseqid + '  '
        aligned_string = aligned_string.ljust(spacing_1)
        aligned_string += self.get_contig_details_string() + '  '
        aligned_string = aligned_string.ljust(spacing_2)
        aligned_string += self.get_coverage_details_string(cov_space) + '  '
        aligned_string = aligned_string.ljust(spacing_3)
        aligned_string += self.get_identity_details_string(id_space)
        return aligned_string

    def get_assembly_piece(self, assembly):
        """Returns the piece of the assembly which corresponds to this BLAST hit."""
        return AssemblyPiece(assembly, self.sseqid, self.sstart, self.send, self.strand, [self])

    def get_query_range(self):
        """Produces an IntRange object for the hit query."""
        return IntRange([(self.qstart, self.qend)])

    def in_assembly_pieces(self, assembly_pieces):
        """
        Returns True if the hit is in (or at least overlaps with) any of the given assembly pieces.
        """
        for piece in assembly_pieces:
            if piece.overlaps(self.sseqid, self.sstart, self.send):
                return True
        return False

    def get_blast_result_json_dict(self, assembly):
        blast_results = OrderedDict()
        blast_results['Coverage'] = self.get_coverage_string()
        blast_results['Identity'] = self.get_identity_string()
        blast_results['Contig name'] = self.sseqid
        blast_results['Contig start position'] = self.sstart
        blast_results['Contig end position'] = self.send
        blast_results['Contig strand'] = self.strand
        blast_results['Bit score'] = self.bitscore
        blast_results['E-value'] = self.evalue
        return blast_results


class GeneBlastHit(BlastHit):
    """This class adds a few gene-specific things to the BlastHit class."""
    def __init__(self, hit_string):
        BlastHit.__init__(self, hit_string)
        self.over_identity_threshold = False

    def conflicts(self, other):
        """
        Returns whether or not this hit conflicts with the other hit.
        A conflict is defined as the hits overlapping by 50% or more of the shortest hit's length.
        A hit is not considered to conflict with itself.
        """
        if self is other:
            return False
        if self.sseqid != other.sseqid:
            return False
        max_start = max(self.sstart, other.sstart)
        min_end = min(self.send, other.send)
        if max_start < min_end:
            overlap = min_end - max_start
        else:
            overlap = 0
        min_length = min(self.send - self.sstart, other.send - other.sstart)
        frac_overlap = overlap / min_length
        return frac_overlap > 0.5

    def get_blast_result_json_dict(self, assembly):
        blast_results = super(GeneBlastHit, self).get_blast_result_json_dict(assembly)
        nuc_seq = assembly.contigs[self.sseqid][self.sstart:self.send]
        if self.strand == '-':
            nuc_seq = reverse_complement(nuc_seq)
        blast_results['Nucleotide length'] = len(nuc_seq)
        blast_results['Protein length'] = len(self.sseq)
        blast_results['Nucleotide sequence'] = nuc_seq
        blast_results['Protein sequence'] = self.sseq
        return blast_results

    def get_match_confidence(self):
        cov = self.query_cov
        ident = self.pident
        if cov == 100.0 and ident >= 99.0:
            confidence = 'Very high'
        elif cov >= 99.0 and ident >= 95.0:
            confidence = 'High'
        elif cov >= 97.0 and ident >= 95.0:
            confidence = 'Good'
        elif cov >= 95.0 and ident >= 85.0:
            confidence = 'Low'
        else:
            confidence = 'None'
        return confidence


class TypeGeneBlastHit(BlastHit):
    """This class adds a couple type gene-specific things to the BlastHit class."""
    def __init__(self, hit_string):
        BlastHit.__init__(self, hit_string)
        try:
            name_parts = self.qseqid.split('__')
            self.gene_name = name_parts[1]
            self.allele_number = int(name_parts[2])
        except (IndexError, ValueError):
            self.gene_name = ''
            self.allele_number = 0

    def get_blast_result_json_dict(self, assembly):
        blast_results = OrderedDict()
        blast_results['Coverage'] = self.get_coverage_string()
        blast_results['Identity'] = self.get_identity_string()
        blast_results['Assembly piece name'] = self.sseqid
        blast_results['Assembly piece start position'] = self.sstart
        blast_results['Assembly piece end position'] = self.send
        blast_results['Assembly piece strand'] = self.strand
        blast_results['Bit score'] = self.bitscore
        blast_results['E-value'] = self.evalue
        blast_results['Length'] = len(self.sseq)
        blast_results['Sequence'] = self.sseq
        return blast_results


class KLocus(object):
    def __init__(self, name, seq, genes):
        self.name = name
        self.seq = seq
        self.genes = genes
        self.gene_names = [x.full_name for x in genes]
        self.blast_hits = []
        self.hit_ranges = IntRange()
        self.assembly_pieces = []
        self.identity = 0.0
        self.expected_hits_inside_locus = []
        self.missing_expected_genes = []
        self.expected_hits_outside_locus = []
        self.other_hits_inside_locus = []
        self.other_hits_outside_locus = []

    def __repr__(self):
        return 'Locus ' + self.name

    def get_length(self):
        """Returns the locus sequence length."""
        return len(self.seq)

    def add_blast_hit(self, hit):
        """Adds a BLAST hit and updates the hit ranges."""
        self.blast_hits.append(hit)
        self.hit_ranges.add_range(hit.qstart, hit.qend)

    def get_mean_blast_hit_identity(self):
        """Returns the mean identity (weighted by hit length) for all BLAST hits in the locus."""
        identity_sum = 0.0
        length_sum = 0
        for hit in self.blast_hits:
            length_sum += hit.length
            identity_sum += hit.length * hit.pident
        if identity_sum == 0.0:
            return 0.0
        else:
            return identity_sum / length_sum

    def clear(self):
        """
        Clears everything in the KLocus object relevant to a particular assembly - gets it ready
        for the next assembly.
        """
        self.blast_hits = []
        self.hit_ranges = IntRange()
        self.assembly_pieces = []
        self.identity = 0.0
        self.expected_hits_inside_locus = []
        self.missing_expected_genes = []
        self.expected_hits_outside_locus = []
        self.other_hits_inside_locus = []
        self.other_hits_outside_locus = []

    def get_coverage(self):
        """Returns the % of this locus which is covered by BLAST hits in the given assembly."""
        try:
            return 100.0 * self.hit_ranges.get_total_length() / len(self.seq)
        except ZeroDivisionError:
            return 0.0
    
    def get_coverage_string(self):
        return '%.2f' % self.get_coverage() + '%'

    def get_identity_string(self):
        return '%.2f' % self.identity + '%'
    
    def clean_up_blast_hits(self):
        """
        This function removes unnecessary BLAST hits from self.blast_hits.
        For each BLAST hit, we keep it if it offers new parts of the locus. If, on the other
        hand, it lies entirely within an existing hit (in locus positions), we ignore it. Since
        we first sort the BLAST hits longest to shortest, this strategy will prioritise long hits
        over short ones.
        """
        self.blast_hits.sort(key=lambda x: x.length, reverse=True)
        kept_hits = []
        k_range_so_far = IntRange()
        for hit in self.blast_hits:
            hit_range = hit.get_query_range()
            if not k_range_so_far.contains(hit_range):
                k_range_so_far.merge_in_range(hit_range)
                kept_hits.append(hit)
        self.blast_hits = kept_hits

    def get_match_uncertainty_chars(self):
        """
        Returns the character code which indicates uncertainty with how this locus was found in
        the current assembly.
        '?' means the locus was found in multiple discontinuous assembly pieces.
        '-' means that one or more expected genes were missing.
        '+' means that one or more additional genes were found in the locus assembly parts.
        '*' means that at least one of the expected genes in the locus is low identity.
        """
        uncertainty_chars = ''
        if len(self.assembly_pieces) > 1:
            uncertainty_chars += '?'
        if self.missing_expected_genes:
            uncertainty_chars += '-'
        if self.other_hits_inside_locus:
            uncertainty_chars += '+'
        if not all([x.over_identity_threshold for x in self.expected_hits_inside_locus]):
            uncertainty_chars += '*'
        return uncertainty_chars

    def get_length_discrepancy(self):
        """
        Returns an integer of the base discrepancy between the locus in the assembly and the
        reference locus sequence.
        E.g. if the assembly match was 5 bases shorter than the reference, this returns -5.
        This function only applies to cases where the locus was found in a single piece. In
        other cases, it returns None.
        """
        if len(self.assembly_pieces) != 1:
            return None
        only_piece = self.assembly_pieces[0]
        a_start = only_piece.start
        a_end = only_piece.end
        k_start = only_piece.earliest_hit_coordinate()
        k_end = only_piece.latest_hit_coordinate()
        expected_length = k_end - k_start
        actual_length = a_end - a_start
        return actual_length - expected_length

    def get_length_discrepancy_string(self):
        """
        Returns the length discrepancy, not as an integer but as a string with a sign and units.
        """
        length_discrepancy = self.get_length_discrepancy()
        if length_discrepancy is None:
            return 'n/a'
        length_discrepancy_string = str(length_discrepancy) + ' bp'
        if length_discrepancy > 0:
            length_discrepancy_string = '+' + length_discrepancy_string
        return length_discrepancy_string

    def get_earliest_and_latest_pieces(self):
        """
        Returns the AssemblyPiece with the earliest coordinate (closest to the locus start) and
        the AssemblyPiece with the latest coordinate (closest to the locus end)
        """
        earliest_piece = sorted(self.assembly_pieces, key=lambda x: x.earliest_hit_coordinate())[0]
        latest_piece = sorted(self.assembly_pieces, key=lambda x: x.latest_hit_coordinate())[-1]
        same_contig_and_strand = earliest_piece.contig_name == latest_piece.contig_name and \
            earliest_piece.strand == latest_piece.strand

        # Even though the pieces are on the same contig and strand, we still need to check whether
        # the earliest piece comes before the latest piece in that contig.
        if same_contig_and_strand:
            if earliest_piece.strand == '+' and earliest_piece.start > latest_piece.end:
                same_contig_and_strand = False
            elif earliest_piece.strand == '-' and earliest_piece.start < latest_piece.end:
                same_contig_and_strand = False
        return earliest_piece, latest_piece, same_contig_and_strand

    def get_match_confidence(self):
        """
        These confidence thresholds match those specified in the paper supp. text, with the
        addition of two new top-level categories: perfect and very high
        """
        single_piece = len(self.assembly_pieces) == 1
        cov = self.get_coverage()
        ident = self.identity
        missing = len(self.missing_expected_genes)
        extra = len(self.other_hits_inside_locus)
        if single_piece and cov == 100.0 and ident == 100.0 and missing == 0 and extra == 0 and \
                self.get_length_discrepancy() == 0:
            confidence = 'Perfect'
        elif single_piece and cov >= 99.0 and ident >= 95.0 and missing == 0 and extra == 0:
            confidence = 'Very high'
        elif single_piece and cov >= 99.0 and missing <= 3 and extra == 0:
            confidence = 'High'
        elif (single_piece or cov >= 95.0) and missing <= 3 and extra <= 1:
            confidence = 'Good'
        elif (single_piece or cov >= 90.0) and missing <= 3 and extra <= 2:
            confidence = 'Low'
        else:
            confidence = 'None'
        return confidence


class Assembly(object):
    def __init__(self, fasta_file):
        """Loads in an assembly and builds a BLAST database for it (if necessary)."""
        self.fasta = fasta_file
        self.name = fasta_file
        self.name = strip_extensions(fasta_file)
        self.contigs = {x[0]: x[1] for x in load_fasta(fasta_file)}  # key = name, value = sequence
        self.blast_db_already_present = self.blast_database_exists()
        if not self.blast_db_already_present:
            makeblastdb(self.fasta)

    def __del__(self):
        if not self.blast_db_already_present:
            clean_blast_db(self.fasta)

    def __repr__(self):
        return self.name

    def blast_database_exists(self):
        """Returns whether or not a BLAST database already exists for this assembly."""
        return os.path.isfile(self.fasta + '.nin') and os.path.isfile(self.fasta + '.nhr') and \
            os.path.isfile(self.fasta + '.nsq')


class AssemblyPiece(object):
    """
    This class describes a piece of an assembly: which contig the piece is on and what the range
    is.
    """
    def __init__(self, assembly, contig_name, contig_start, contig_end, strand, blast_hits=None):
        self.assembly = assembly
        self.contig_name = contig_name
        self.start = contig_start
        self.end = contig_end
        self.strand = strand
        if not blast_hits:
            blast_hits = []
        self.blast_hits = blast_hits

    def __repr__(self):
        return self.assembly.name + '_' + self.get_header()

    def get_header(self):
        """Returns a descriptive string for the FASTA header when saving this piece to file."""
        nice_contig_name = get_nice_contig_name(self.contig_name)
        return nice_contig_name + '_' + str(self.start + 1) + '_to_' + str(self.end) + \
            '_' + self.strand + '_strand'

    def get_nice_header(self):
        """Like get_header, but uses spaces/parentheses instead of underscores for readability."""
        nice_contig_name = get_nice_contig_name(self.contig_name)
        return nice_contig_name + ' (' + str(self.start + 1) + '-' + str(self.end) + \
            ', ' + self.strand + ' strand)'

    def get_bandage_range(self):
        """Returns the assembly piece in a Bandage path format."""
        if is_contig_name_spades_format(self.contig_name):
            name = self.contig_name.split('_')[1]
        else:
            name = self.contig_name.split()[0]
        return '(' + str(self.start + 1) + ') ' + name + '+ (' + str(self.end) + ')'

    def get_sequence(self):
        """Returns the DNA sequence for this piece of the assembly."""
        seq = self.assembly.contigs[self.contig_name][self.start:self.end]
        if self.strand == '+':
            return seq
        else:
            return reverse_complement(seq)

    def get_length(self):
        """Returns the sequence length for this piece."""
        return self.end - self.start

    def get_sequence_short(self):
        """Returns a shortened format of the sequence"""
        seq = self.get_sequence()
        length = len(seq)
        if len(seq) > 9:
            seq = seq[0:6] + '...' + seq[-6:]
        return seq + ' (' + str(length) + ' bp)'

    def combine(self, other):
        """
        If this assembly piece and the other can be combined, this function returns the combined
        piece. If they can't, it returns None.
        To be able to combine, pieces must be overlapping or directly adjacent and on the same
        strand.
        """
        if self.contig_name != other.contig_name or self.strand != other.strand:
            return None
        combined = IntRange([(self.start, self.end)])
        combined.add_range(other.start, other.end)
        if len(combined.ranges) == 1:
            new_start, new_end = combined.ranges[0]
            combined_hits = self.blast_hits + other.blast_hits
            return AssemblyPiece(self.assembly, self.contig_name, new_start, new_end, self.strand,
                                 combined_hits)
        else:
            return None

    def overlaps(self, contig_name, start, end):
        """Returns whether this assembly piece overlaps with the given parameters."""
        return self.contig_name == contig_name and self.start < end and start < self.end

    def earliest_hit_coordinate(self):
        """Returns the lowest query start coordinate in the BLAST hits."""
        if not self.blast_hits:
            return None
        return sorted([x.qstart for x in self.blast_hits])[0]

    def latest_hit_coordinate(self):
        """Returns the highest query end coordinate in the BLAST hits."""
        if not self.blast_hits:
            return None
        return sorted([x.qend for x in self.blast_hits])[-1]


class IntRange(object):
    """
    This class contains one or more integer ranges. Overlapping ranges will be merged together.
    It stores its ranges in a Python-like fashion where the last value in each range is
    exclusive.
    """
    def __init__(self, ranges=None):
        if not ranges:
            ranges = []
        self.ranges = []
        self.add_ranges(ranges)
        self.simplify()

    def __repr__(self):
        return str(self.ranges)

    def add_range(self, start, end):
        """Adds a single range."""
        self.add_ranges([(start, end)])

    def add_ranges(self, ranges):
        """Adds multiple ranges (list of tuples)."""
        self.ranges += ranges
        self.simplify()

    def merge_in_range(self, other):
        """Merges the other IntRange object into this one."""
        self.add_ranges(other.ranges)

    def get_total_length(self):
        """Returns the number of integers in the ranges."""
        return sum([x[1] - x[0] for x in self.ranges])

    def simplify(self):
        """Collapses overlapping ranges together."""
        fixed_ranges = []
        for int_range in self.ranges:
            if int_range[0] > int_range[1]:
                fixed_ranges.append((int_range[1], int_range[0]))
            elif int_range[0] < int_range[1]:
                fixed_ranges.append(int_range)
        starts_ends = [(x[0], 1) for x in fixed_ranges]
        starts_ends += [(x[1], -1) for x in fixed_ranges]
        starts_ends.sort(key=lambda z: z[0])
        current_sum = 0
        cumulative_sum = []
        for start_end in starts_ends:
            current_sum += start_end[1]
            cumulative_sum.append((start_end[0], current_sum))
        prev_depth = 0
        start = 0
        combined = []
        for pos, depth in cumulative_sum:
            if prev_depth == 0:
                start = pos
            elif depth == 0:
                combined.append((start, pos))
            prev_depth = depth
        self.ranges = combined

    def contains(self, other):
        """Returns True if the other IntRange is entirely contained within this IntRange."""
        for other_range in other.ranges:
            other_start, other_end = other_range
            contained = False
            for this_range in self.ranges:
                this_start, this_end = this_range
                if other_start >= this_start and other_end <= this_end:
                    contained = True
                    break
            if not contained:
                return False
        return True


class Gene(object):
    """This class prepares and stores a gene taken from the input Genbank file."""
    def __init__(self, k_locus_name, num, feature, k_locus_seq):
        self.k_locus_name = k_locus_name
        self.feature = feature
        gene_num_string = str(num).zfill(2)
        self.full_name = k_locus_name + '_' + gene_num_string
        if 'gene' in feature.qualifiers:
            self.gene_name = feature.qualifiers['gene'][0]
            self.full_name += '_' + self.gene_name
        else:
            self.gene_name = None
        if 'product' in feature.qualifiers:
            self.product = feature.qualifiers['product'][0]
        else:
            self.product = None
        if 'EC_number' in feature.qualifiers:
            self.ec_number = feature.qualifiers['EC_number'][0]
        else:
            self.ec_number = None
        self.nuc_seq = feature.extract(k_locus_seq)
        self.prot_seq = str(self.nuc_seq.translate(table=11))
        self.nuc_seq = str(self.nuc_seq)

    def get_fasta(self):
        """
        Returns the FASTA version of this gene: a header line followed by sequence lines (of
        protein sequence) ending in a line break.
        """
        return '>' + self.full_name + '\n' + \
               add_line_breaks_to_sequence(self.prot_seq, 60)

    def get_reference_info_json_dict(self):
        reference_dict = OrderedDict()
        if self.gene_name:
            reference_dict['Gene'] = self.gene_name
        if self.product:
            reference_dict['Product'] = self.product
        if self.ec_number:
            reference_dict['EC number'] = self.ec_number
        reference_dict['Nucleotide length'] = len(self.nuc_seq)
        reference_dict['Protein length'] = len(self.prot_seq)
        reference_dict['Nucleotide sequence'] = self.nuc_seq
        reference_dict['Protein sequence'] = self.prot_seq
        return reference_dict


def convert_bytes_to_str(bytes_or_str):
    """
    This function is for both Python2 and Python3. If the input is a str, it just returns that
    same str. If not, it assumes its bytes (and we're in Python3) and it returns it as a str.
    """
    if isinstance(bytes_or_str, str):
        return bytes_or_str
    else:
        return bytes_or_str.decode()


def makeblastdb(fasta):
    """
    If the FASTA file is not compressed, this just runs makeblastdb. If it is compressed,
    it runs gunzip and pipes into makeblastdb.
    """
    if ' ' in fasta:
        print('WARNING: spaces in file paths may not work in BLAST', file=sys.stderr)
    if get_compression_type(fasta) == 'gz':
        gunzip_command = ['gunzip', '-c', fasta]
        makeblastdb_command = ['makeblastdb', '-dbtype', 'nucl', '-in', '-', '-out', fasta,
                               '-title', fasta]
        gunzip = subprocess.Popen(gunzip_command, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        makeblastdb_process = subprocess.Popen(makeblastdb_command, stdin=gunzip.stdout,
                                               stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        gunzip.stdout.close()
        _, err = makeblastdb_process.communicate()
    else:  # plain text
        makeblastdb_command = ['makeblastdb', '-dbtype', 'nucl', '-in', fasta]
        makeblastdb_process = subprocess.Popen(makeblastdb_command, stdout=subprocess.PIPE,
                                               stderr=subprocess.PIPE)
        _, err = makeblastdb_process.communicate()
    if err:
        quit_with_error('makeblastdb encountered an error:\n' + convert_bytes_to_str(err))


def remove_if_exists(filename):
    try:
        os.remove(filename)
    except OSError:
        pass


def clean_blast_db(fasta):
    remove_if_exists(fasta + '.nsq')
    remove_if_exists(fasta + '.nhr')
    remove_if_exists(fasta + '.nin')


def get_compression_type(filename):
    """
    Attempts to guess the compression (if any) on a file using the first few bytes.
    http://stackoverflow.com/questions/13044562
    """
    magic_dict = {'gz': (b'\x1f', b'\x8b', b'\x08'),
                  'bz2': (b'\x42', b'\x5a', b'\x68'),
                  'zip': (b'\x50', b'\x4b', b'\x03', b'\x04')}
    max_len = max(len(x) for x in magic_dict)

    unknown_file = open(filename, 'rb')
    file_start = unknown_file.read(max_len)
    unknown_file.close()
    compression_type = 'plain'
    for filetype, magic_bytes in magic_dict.items():
        if file_start.startswith(magic_bytes):
            compression_type = filetype
    if compression_type == 'bz2':
        quit_with_error('cannot use bzip2 format - use gzip instead')
    if compression_type == 'zip':
        quit_with_error('cannot use zip format - use gzip instead')
    return compression_type


def strip_extensions(filename):
    """
    This function removes extensions from a file name. Examples:
      assembly.fasta -> assembly
      assembly.fa.gz -> assembly
      genome.assembly.fa.gz -> genome.assembly
    """
    name = os.path.basename(filename)
    if name.lower().endswith('.gz'):
        name = name[:-3]
    if name.lower().endswith('.fa'):
        name = name[:-3]
    elif name.lower().endswith('.fna'):
        name = name[:-4]
    elif name.lower().endswith('.fas'):
        name = name[:-4]
    elif name.lower().endswith('.fasta'):
        name = name[:-6]
    return name


if __name__ == '__main__':
    main()
