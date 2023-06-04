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
from tempfile import TemporaryDirectory
from pathlib import Path

from kaptive.log import warning
from kaptive.database import Database
from kaptive.assembly import Assembly, load_fasta


def run_align(args, db: Database):
    for fasta_file in args.assemblies:
        sample = load_fasta(fasta_file)
        if not sample:
            warning(f'No sequences found in {fasta_file.name}')
            continue
        sample = Assembly(fasta_file, sample)  # Create Assembly object

        if args.kaptive_refs and db.loci:
            sample.get_best_locus_match(db)

        if args.allelic_typing and db.type_gene_names:
            sample.type_gene_search(db)

        if sample.result is None:
            warning(f'No locus found for {sample.name}')
            continue
        else:
            sample.result.print_result()
            if not args.no_seq_out:
                sample.result.write_fasta()
            if not args.no_table:
                sample.result.add_to_table()
            if not args.no_json:
                sample.result.add_to_json()

        if not args.keep_db:
            for file in sample.blast_db:
                file.unlink(missing_ok=True)


def get_gene_info_string(gene_hit_list):
    """Returns a single comma-delimited string summarising the gene hits in the given list."""
    gene_hit_strings = []
    for gene in gene_hit_list:
        gene_hit_string = ""
        gene_hit_string += f'{gene.qseqid},{gene.pident}%'
        if gene.truncated_protein:
            gene_hit_string += f',truncated_protein({len(gene.prot_seq)}/{len(gene.gene.prot_seq)})'
        if gene.fragmented:
            gene_hit_string += f',fragmented({gene.get_coverage_string()})'
        if gene.insertion:
            gene_hit_string += ',insertion'
        if gene.deletion:
            gene_hit_string += ',deletion'
        if gene.edge_of_contig:
            gene_hit_string += ',edge_of_contig'
        gene_hit_strings.append(gene_hit_string)
    return ';'.join(gene_hit_strings)


def get_scores(hits: List[BlastHit], groups: dict, attr_list=['bitscore', 'pident', 'query_cov']) -> {}:
    """
    For a list of blast hits, calculates the total and mean scores per group.

    The key designates the group, the value should be a list of objects that a name attribute corresponding to the
    qseqid of the blast hit and a length attribute corresponding to the length of the query sequence.
    A group could be:
        {KL1: [Gene, Gene], KL2: [Gene, Gene]} for gene hits
         {KL1: [Locus], KL2: [Locus]} for locus hits
    """
    scores = {}
    for group, expected in groups.items():
        found = [i for i in hits if i.qseqid in [i.name for i in expected]]
        if found:
            scores[group] = {'expected': len(expected), 'found': len(found)}
            for attr in attr_list:
                scores[group][attr] = sum(getattr(hit, attr) for hit in found)
                scores[group][f'mean {attr} by expected hits'] = scores[group][attr] / len(expected)
                scores[group][f'mean {attr} by length'] = scores[group][attr] / sum(len(i) for i in found)
                scores[group][f'mean {attr} by found hits'] = scores[group][attr] / len(found)
    return scores


