"""
This module contains classes for interacting with bacterial genome assemblies and contigs and a pipeline
to type them.

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
from pathlib import Path
from statistics import stdev, mean
from itertools import groupby, chain
from subprocess import Popen, PIPE, DEVNULL

from Bio.Seq import Seq

from kaptive.typing import AssemblyTypingResult, ContigPiece, AssemblyGeneResult
from kaptive.misc import parse_fasta
from kaptive.intrange import merge_ranges, range_overlap
from kaptive.database import Database, Locus
from kaptive.alignment import Alignment, cull_conflicting_alignments, cull_all_conflicting_alignments
from kaptive.log import warning, log, quit_with_error


# Classes -------------------------------------------------------------------------------------------------------------
class AssemblyError(Exception):
    pass


class Assembly:
    def __init__(self, path: Path | None = None, name: str | None = None, contigs: dict[str: Contig] | None = None):
        self.path = path or Path()
        self.name = name or ''
        self.contigs = contigs or {}

    @classmethod
    def from_path(cls, path: Path, allow_plasmids: bool = True) -> Assembly:
        self = cls(path=path, name=path.name.strip('.gz').rsplit('.', 1)[0])
        for n, d, s in parse_fasta(path):
            self.contigs[c.name] = (c := Contig(self, n, d, Seq(s)))
        if not self.contigs:
            raise AssemblyError(f"No contigs found in {path}")
        return self

    def __repr__(self):
        return self.name

    def __hash__(self):
        return hash(self.name)

    def __len__(self):
        return self.path.stat().st_size

    def __iter__(self) -> tuple[str, Contig]:
        for contig_name, contig in self.contigs.items():
            yield contig_name, contig

    def __getitem__(self, contig_name):
        return self.contigs[contig_name]

    def __contains__(self, contig_name):
        return contig_name in self.contigs


class ContigError(Exception):
    pass


class Contig(object):
    """
    This class describes a contig in an assembly: the name, length, and sequence.
    """

    def __init__(self, assembly: Assembly | None = None, name: str | None = None, description: str | None = None,
                 sequence: Seq | None = Seq(''), replicon: str | None = None, topology: str | None = None):
        self.assembly = assembly
        self.name = name or ''
        self.description = description or ''
        self.sequence = sequence
        self.replicon = replicon or ''
        self.topology = topology or ''
        if self.description:
            self._get_replicon()
            self._get_topology()

    def __repr__(self):
        return self.name

    def __len__(self):
        return len(self.sequence)

    def __hash__(self):
        return hash(self.name)

    def _get_replicon(self):
        if any(i in self.name + self.description for i in {'plasmid', '__pl'}):
            self.replicon = "plasmid"
        elif any(i in self.name + self.description for i in {'chromosome', '__ch'}):
            self.replicon = "chromosome"
        else:
            self.replicon = "unknown"

    def _get_topology(self):
        if any(i in self.name + self.description for i in {'circular=true', 'complete', 'circular=Y'}):
            self.topology = "circular"
        elif any(i in self.name + self.description for i in {'circular=false', 'circular=N', 'linear'}):
            self.topology = "linear"
        else:
            self.topology = "unknown"


# Functions -----------------------------------------------------------------------------------------------------------
def get_scores(results: list[tuple[Locus: list[Alignment]]], score: str, weight: str,
               reverse_sort: bool = True) -> list[tuple[Locus, float]]:
    if weight == 'none':
        results = [(l, sum(getattr(x, score) for x in r)) for l, r in results]
    elif weight == "locus_length":
        results = [(l, sum(getattr(x, score) for x in r) / len(l)) for l, r in results]
    elif weight == "genes_expected":
        results = [(l, sum(getattr(x, score) for x in r) / len(l.genes)) for l, r in results]
    elif weight == "genes_found":
        results = [(l, sum(getattr(x, score) for x in r) / len(r)) for l, r in results]
    elif weight == "prop_genes_found":
        results = [(l, sum(getattr(x, score) for x in r) * (len(r) / len(l.genes))) for l, r in results]
    else:
        quit_with_error(f"Invalid weight {weight}")
    return sorted(results, key=lambda x: x[1], reverse=reverse_sort)


def score_stats(scores: list[tuple[Locus, float]]) -> list[tuple[Locus, float, float]]:
    """Returns locus, score and zscore"""
    _scores = {i[1] for i in scores}  # Get unique scores
    if len(_scores) == 1:
        return [(scores[0][0], scores[0][1], 0)]  # If only one score, return it with a zscore of 0
    m, sd = mean(_scores), stdev(_scores)  # Calculate mean and standard deviation
    return [(scores[0][0], scores[0][1], 0)] if sd == 0 else [  # If standard deviation is 0, return zscore of 0
        (scores[i][0], scores[i][1], (scores[i][1] - m) / sd) for i in range(len(scores))]


def typing_pipeline(
        assembly: Path, db: Database, threads: int | None = 1, min_cov: float | None = 0.5, min_zscore: float | None = 3,
        score_metric: str | None = 'AS', weight_metric: str | None = 'prop_genes_found', verbose: bool | None = False
) -> AssemblyTypingResult | None:
    """Typing function for assemblies. Aligns genes to the assembly and scores the results."""
    try:
        assembly = Assembly.from_path(assembly)
    except AssemblyError as e:
        warning(f'Error parsing {assembly}: {e}')
        return

    log(f'Typing {assembly}, min_cov={min_cov}, min_zscore={min_zscore}, score_metric={score_metric}, '
        f'weight_metric={weight_metric}', verbose)

    # assembly = Assembly.from_path(args.input[0])
    # threads = args.threads
    # min_cov = args.min_cov
    # min_zscore = args.min_zscore
    # score_metric = args.score_metric
    # weight_metric = args.weight_metric
    # verbose = args.verbose

    # ALIGN AND SORT GENES -------------------------------------------------------------------------------------------
    gene_alignments, for_scoring, extra_gene_alignments, is_element_alignments = [], {}, {}, []
    cmd = f"minimap2 -c -t {threads} {assembly.path} -"  # f"miniprot -T 11 -S -t {threads} {assembly.path} -"
    with Popen(cmd.split(), stdin=PIPE, stdout=PIPE, stderr=PIPE) as proc:  # Align genes to assembly
        stdout, stderr = proc.communicate(db.as_gene_fasta().encode())  # Send gene sequences to proc.stdin
        for a in stdout.splitlines():  # Read alignments from proc.stdout
            if not (a := Alignment.from_paf_line(a)).target_name in assembly.contigs:
                continue  # Skip alignment if it's not to a contig in the assembly, happens if plasmids aren't allowed
            if a.query_name in db.genes:  # If the alignment is to a gene
                gene_alignments.append(a)  # Add to list of gene alignments
                if a.percent_query_coverage >= min_cov:  # If the alignment has sufficient coverage
                    if (x := for_scoring.get(a.query_name, None)) and a.AS < x.AS:
                        continue  # If there's already an alignment to this gene, and it's better, skip this one
                    for_scoring[a.query_name] = a  # Add to list of alignments for scoring
            elif a.query_name in db.extra_genes:  # If the alignment is to an extra gene
                if (x := extra_gene_alignments.get(a.query_name, None)) and a.AS < x.AS:
                    continue  # If there's already an alignment to this gene, and it's better, skip this one
                extra_gene_alignments[a.query_name] = a
            else:
                warning(f'Alignment to unknown gene {a.query_name} in {assembly}')

    if not gene_alignments:  # If no gene alignments, skip this assembly
        warning(f'No alignments for {assembly}: {stderr.decode()}')
        return

    if not for_scoring:  # If no alignments, skip this assembly
        warning(f'No gene alignments sufficient for typing {assembly}')
        return

    # CALCULATE BEST MATCH AND SCORES --------------------------------------------------------------------------------
    for_scoring = [(db.loci[l], list(a)) for l, a in groupby(sorted(for_scoring.values(), key=lambda x: x.query_name.split('_')[0]
                                                                    ), key=lambda x: x.query_name.split('_')[0])]
    scores = {  # Get scores for each metric
        metric: {  # Get scores for each weight
            weight: score_stats(get_scores(for_scoring, metric, weight)) for weight in {
                'none', 'locus_length', 'genes_expected', 'genes_found', 'prop_genes_found', weight_metric}
        } for metric in {'percent_identity', 'matching_bases', 'AS', 'percent_query_coverage', score_metric}
    }
    best_match, score, zscore = scores[score_metric][weight_metric][0]  # Get best match, score, and zscore
    result = AssemblyTypingResult(  # Create result object
        db=db, best_match=best_match, assembly=assembly, score=score, zscore=zscore, scores=scores,
        min_zscore=min_zscore, min_cov = min_cov, score_metric=score_metric, weight_metric = weight_metric)

    # GET LOCUS PIECES PER CONTIG -------------------------------------------------------------------------------------
    for contig, alignments in groupby(sorted(chain(  # Group alignments by contig (target_name)
            gene_alignments, extra_gene_alignments.values(), is_element_alignments
    ), key=lambda x: x.target_name), key=lambda x: x.target_name):

        # SORT ALIGNMENTS ON CONTIG -----------------------------------------------------------------------------------
        contig = assembly.contigs[contig]  # type: Contig
        expected, unexpected = {}, []  # Init dicts to store alignments
        for a in alignments:  # Store best locus gene alignments in dict to get expected order later
            if a.query_name in best_match.genes:
                expected[a.query_name] = a  # Store best locus gene alignments in dict to get expected order later
            else:  # Store other alignments in a list for quick culling
                unexpected.append(a)  # Store other alignments in a list for quick culling
        unexpected = cull_all_conflicting_alignments(unexpected)  # Remove conflicting other alignments
        for a in expected.values():  # Remove other alignments overlapping best match gene alignments
            unexpected = list(cull_conflicting_alignments(a, unexpected))

        # GET PIECES FROM EXPECTED GENES ---------------------------------------------------------------------------
        contig_pieces = []
        if expected:  # If there are any expected gene alignments
            strand = "+" if max(  # Determine strand of the locus on the contig
                i.strand == best_match.genes[i.query_name].strand for i in expected.values()) else "-"
            last_pos, last_gene = max([(int(a.query_name.split("_")[1]), a) for a in expected.values()],
                                      key=lambda x: x[0])
            ranges = [  # Collect ranges of expected gene alignments
                (a.target_start, a.target_end) for a in expected.values()
            ] if last_pos != len(best_match.genes) else [  # If the last gene is not the last gene in the locus
                (a.target_start, a.target_end) for a in expected.values() if \
                (a.strand == "+" and a.target_end <= last_gene.target_end) or \
                (a.strand == "-" and a.target_start >= last_gene.target_start)
            ]  # Else only collect ranges that don't extend beyond the last gene
            for start, end in merge_ranges(ranges, tolerance=len(db.largest_locus)):
                contig_pieces.append(
                    ContigPiece(
                        contig, start=start, end=end, result=result, strand=strand,
                        sequence=contig.sequence[start:end] if strand == "+" else contig.sequence[
                                                                                  start:end].reverse_complement()
                    )
                )
        # SORT GENES ON CONTIG -------------------------------------------------------------------------------------
        previous_result, piece = None, None  # For determining neighbouring genes
        for a in sorted(chain(expected.values(), unexpected), key=lambda x: x.target_start):  # Add gene results
            # Determine gene type, gene, and piece edge threshold
            if a.query_name in expected:
                gene_type, gene, piece_edge = 'expected_genes', best_match.genes[a.query_name], 0
            elif (gene := db.extra_genes.get(a.query_name, None)):
                gene_type, piece_edge = 'extra_genes', 0
            elif (gene := db.is_elements.get(a.query_name, None)):
                gene_type, piece_edge = 'is_elements', 200  # IS elements are allowed to be 200 bp outside piece
            else:
                gene_type, gene, piece_edge = 'unexpected_genes', db.genes[a.query_name], 0
            # Determine if alignment is below threshold
            below_threshold = a.percent_identity < db.gene_threshold
            # Determine if gene outside locus (piece == None) or inside locus (piece != None)
            if contig_pieces:  # Get the piece the alignment is in if there are any pieces
                piece = next((p for p in contig_pieces if range_overlap((p.start - piece_edge, p.end + piece_edge),
                                                                        (a.target_start, a.target_end))), None)
            if not piece:  # Assert exceptions and additional logic
                if gene_type == 'is_elements' or below_threshold:
                    continue  # Ignore IS elements and below threshold genes outside the locus (orthologs)
            else:
                if gene_type == "extra_genes":
                    warning(
                        f"Extra gene {gene.name} found inside locus {best_match.name} on contig {contig.name}")
                    continue  # Ignore extra genes inside the locus
                if gene_type == "unexpected_genes" and not below_threshold:
                    gene_type = 'expected_genes'  # If unexpected gene is not below threshold, it is expected

            # Create gene result and add it to the result
            gene_result = AssemblyGeneResult(
                contig, gene, start=a.target_start, end=a.target_end, strand=a.strand,
                percent_identity=a.percent_identity, percent_coverage=a.percent_query_coverage,
                piece=piece, below_threshold=below_threshold, gene_type=gene_type,
                neighbour_left=previous_result,
                dna_seq=contig.sequence[a.target_start:a.target_end] if a.strand == "+" else \
                    contig.sequence[a.target_start:a.target_end].reverse_complement()
            )
            result.add_gene_result(gene_result)
            previous_result = gene_result

        result.pieces += contig_pieces  # Add the pieces to the result

    # Sort the pieces by the sum of the expected gene order to get the expected order of the pieces
    result.pieces.sort(key=lambda x: min(int(i.gene.name.split("_")[1]) for i in x.expected_genes))
    result.expected_genes_inside_locus.sort(key=lambda x: int(x.gene.name.split("_")[1]))
    result.expected_genes_outside_locus.sort(key=lambda x: int(x.gene.name.split("_")[1]))
    result.unexpected_genes_inside_locus.sort(key=lambda x: int(x.gene.name.split("_")[1]))
    result.unexpected_genes_outside_locus.sort(key=lambda x: int(x.gene.name.split("_")[1]))
    return result
