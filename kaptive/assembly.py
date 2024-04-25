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

from pathlib import Path
from statistics import stdev, mean
from itertools import groupby, chain
from subprocess import Popen, PIPE

from Bio.Seq import Seq

from kaptive.typing import TypingResult, LocusPiece, GeneResult
from kaptive.database import Database, Locus
from kaptive.alignment import Alignment, cull_conflicting_alignments, cull_all_conflicting_alignments
from kaptive.misc import parse_fasta, merge_ranges, range_overlap
from kaptive.log import warning, quit_with_error


# Classes -------------------------------------------------------------------------------------------------------------
class AssemblyError(Exception):
    pass


class Assembly:
    def __init__(self, path: Path | None = None, name: str | None = None, contigs: dict[str: Contig] | None = None):
        self.path = path or Path()
        self.name = name or ''
        self.contigs = contigs or {}

    @classmethod
    def from_path(cls, path: Path, **kwargs) -> Assembly:
        self = cls(path=path, name=path.name.strip('.gz').rsplit('.', 1)[0],
                   contigs={n: Contig(n, d, Seq(s)) for n, d, s in parse_fasta(path, **kwargs)})
        if not self.contigs:
            raise AssemblyError(f"No contigs found in {path}")
        return self

    def __repr__(self):
        return self.name

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

    def __init__(self, name: str | None = None, description: str | None = None,
                 sequence: Seq | None = Seq('')):
        self.name = name or ''
        self.description = description or ''
        self.sequence = sequence

    def __repr__(self):
        return self.name

    def __len__(self):
        return len(self.sequence)


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
        assembly: Path | Assembly, db: Database | Path, threads: int | None = 1, min_cov: float | None = 0.5,
        score_metric: str | None = 'AS', weight_metric: str | None = 'prop_genes_found',
        max_other_genes: float | None = 1, percent_expected_genes: float | None = 50,
        allow_below_threshold: bool | None = False, debug: bool | None = False,
        verbose: bool | None = False) -> TypingResult | None:
    if isinstance(db, Path):
        try:
            db = Database.from_genbank(db)
        except Exception as e:
            warning(f'Error parsing {db}: {e}')
            return None
    if isinstance(assembly, Path):
        try:
            assembly = Assembly.from_path(assembly, verbose=verbose)
        except AssemblyError as e:
            warning(f'Error parsing {assembly}: {e}')
            return None

    # ALIGN AND SORT GENES -------------------------------------------------------------------------------------------
    gene_alignments, for_scoring, extra_gene_alignments, is_element_alignments = [], {}, {}, []
    with Popen(cmd := f"minimap2 -c -t {threads} {assembly.path} -".split(), stdin=PIPE, stdout=PIPE,
               stderr=PIPE) as proc:  # Align genes to assembly
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

    if not any([gene_alignments, for_scoring]):  # If no alignments, skip this assembly
        warning(f'No gene alignments sufficient for typing {assembly}')
        return None

    # CALCULATE BEST MATCH AND SCORES --------------------------------------------------------------------------------
    for_scoring = [(
        db.loci[locus], list(a)) for locus, a in groupby(  # Group alignments by locus
        sorted(for_scoring.values(), key=lambda i: i.query_name.split('_')[0]), key=lambda i: i.query_name.split('_')[0]
    )]  # type: list[tuple[Locus, list[Alignment]]]
    weight_metrics = {'none', 'locus_length', 'genes_expected', 'genes_found', 'prop_genes_found',
                      weight_metric} if debug else {weight_metric}
    score_metrics = {'percent_identity', 'matching_bases', 'AS', 'percent_query_coverage', score_metric} if debug else {
        score_metric}
    scores = {s: {w: score_stats(get_scores(for_scoring, s, w)) for w in weight_metrics} for s in score_metrics}
    best_match, score, zscore = scores[score_metric][weight_metric][0]  # Get best match, score, and zscore
    # Create a TypingResult object to store the results
    result = TypingResult(
        assembly.name, db, best_match, score, zscore, scores=scores,
        confidence_args={
            'max_other_genes': max_other_genes, 'percent_expected_genes': percent_expected_genes,
            'allow_below_threshold': allow_below_threshold, 'gene_threshold': db.gene_threshold
        },
        scoring_args={'min_cov': min_cov, 'score_metric': score_metric, 'weight_metric': weight_metric}
    )
    with Popen(cmd, stdin=PIPE, stdout=PIPE, stderr=PIPE) as proc:  # Align the best match locus to the assembly
        locus_alignments = [Alignment.from_paf_line(
            i) for i in proc.communicate(best_match.as_fasta().encode())[0].decode().splitlines()]
    # GET ALIGNMENTS PER CONTIG ---------------------------------------------------------------------------------------
    for contig, alignments in groupby(sorted(chain(locus_alignments, gene_alignments, extra_gene_alignments.values()),
                                             key=lambda x: x.target_name), key=lambda x: x.target_name):
        # SORT ALIGNMENTS ON CONTIG -----------------------------------------------------------------------------------
        contig, locus, expected, other, pieces = assembly.contigs[contig], [], [], [], []
        [locus.append(a) if a.query_name == best_match.name else expected.append(
            a) if a.query_name in best_match.genes else other.append(a) for a in alignments]
        other = cull_all_conflicting_alignments(other)  # Remove conflicting other alignments
        for a in expected:  # Remove other alignments overlapping best match gene alignments
            other = list(cull_conflicting_alignments(a, other))
        # GET PIECES ---------------------------------------------------------------------------------------------------
        if locus and expected:  # If there are locus alignments and expected gene alignments
            pieces = [LocusPiece(contig.name, result, s, e) for s, e in merge_ranges(
                [(a.target_start, a.target_end) for a in locus], tolerance=len(db.largest_locus))]
        # SORT GENES ON CONTIG -------------------------------------------------------------------------------------
        previous_result, piece = None, None  # For determining neighbouring genes
        for a in sorted(chain(expected, other), key=lambda x: x.target_start):  # Add gene results
            # Determine gene type, gene, and piece edge threshold
            if (gene := best_match.genes.get(a.query_name, None)):  # If the alignment is to an expected gene
                gene_type, piece_edge = 'expected_genes', 0
            elif (gene := db.extra_genes.get(a.query_name, None)):  # If the alignment is to an extra gene
                gene_type, piece_edge = 'extra_genes', 0
            # elif (gene := db.is_elements.get(a.query_name, None)):
            #     gene_type, piece_edge = 'is_elements', 200  # IS elements are allowed to be 200 bp outside piece
            else:
                gene_type, gene, piece_edge = 'unexpected_genes', db.genes[a.query_name], 0
            piece = next((p for p in pieces if range_overlap(  # Get the piece the alignment
                (p.start - piece_edge, p.end + piece_edge), (a.target_start, a.target_end))), None) if pieces else None
            # if not piece and gene_type == 'is_elements':  # If the alignment is to an IS element
            #     continue  # Skip IS elements outside the locus
            # Test if gene is partial (would overlap with contig edge if alignment carried through the full gene)
            partial_gene = False
            if a.target_start <= a.query_start:  # If the alignment starts at the beginning of the contig
                partial_gene = True
            if a.target_end <= a.query_end:  # If the alignment ends at the end of the contig
                partial_gene = True
            if (a.target_length - a.target_end) <= (a.query_length - a.query_end):
                # If the gene had extended past the end of the contig if the alignment was longer
                partial_gene = True
                # TODO: Update the gene alignment to extend to the start/end of the contig
            gene_result = GeneResult(  # Create gene result and add it to the result
                contig.name, gene, result, piece, a.target_start, a.target_end, a.strand, previous_result,
                gene_type=gene_type, dna_seq=contig.sequence[a.target_start:a.target_end] if a.strand == "+" else \
                    contig.sequence[a.target_start:a.target_end].reverse_complement(), partial=partial_gene
            )  # Calculate gene frame (from alignment start) and translate to generate protein percent identity
            frame = next((i for i in [0, 1, 2] if (a.query_start + i) / 3 % 1 == 0), 0)
            gene_result.extract_translation(frame, table=11, to_stop=True)  # Extract translation and align to reference
            gene_result.below_threshold = gene_result.percent_identity < db.gene_threshold  # Check if below threshold
            if not piece and gene_result.below_threshold:  # If below threshold
                continue  # Skip this gene, homologue in another part of the genome
            result.add_gene_result(gene_result)  # Add the gene result to the result
            previous_result = gene_result  # Set the previous result to the current result
        # FINALISE CONTIG ---------------------------------------------------------------------------------------------
        for piece in pieces:  # Add sequences to pieces and add them to the result
            if not piece.expected_genes: # If the piece has no expected genes, skip it
                continue
            piece.strand = "+" if max(i.strand == i.gene.strand for i in piece.expected_genes) else "-"
            piece.sequence = contig.sequence[piece.start:piece.end] if piece.strand == "+" else \
                contig.sequence[piece.start:piece.end].reverse_complement()
            result.pieces.append(piece)  # Add the piece to the result
    # FINALISE RESULT -------------------------------------------------------------------------------------------------
    result.missing_genes = list(set(best_match.genes) - {
        i.gene.name for i in chain(result.expected_genes_inside_locus, result.expected_genes_outside_locus)})
    # Sort the pieces by the sum of the expected gene order to get the expected order of the pieces
    result.pieces.sort(key=lambda x: min(int(i.gene.name.split("_")[1]) for i in x.expected_genes))
    [l.sort(key=lambda x: int(x.gene.name.split("_")[1])) for l in (
        result.expected_genes_inside_locus, result.expected_genes_outside_locus, result.unexpected_genes_inside_locus,
        result.unexpected_genes_outside_locus)]
    result.missing_genes.sort(key=lambda x: int(x.split("_")[1]))
    return result
