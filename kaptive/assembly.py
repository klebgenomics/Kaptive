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
from itertools import chain
from json import loads
from subprocess import Popen, PIPE
from typing import TextIO, Pattern, Generator
from re import compile
from os import fstat

from Bio.Seq import Seq
from Bio.SeqIO.FastaIO import SimpleFastaParser
import numpy as np
np.seterr(divide='ignore', invalid='ignore')  # Ignore divide by zero and invalid value errors

from kaptive.typing import TypingResult, LocusPiece, GeneResult
from kaptive.database import Database, load_database
from kaptive.alignment import Alignment, iter_alns, group_alns, cull_filtered
from kaptive.misc import opener, merge_ranges, range_overlap, check_cpus, check_file
from kaptive.log import log, warning

# Constants -----------------------------------------------------------------------------------------------------------
_ASSEMBLY_FASTA_REGEX = compile(r'\.(fasta|fa|fna|ffn)(\.gz)?$')
# _ASSEMBLY_GRAPH_REGEX = compile(r'\.(gfa)(\.gz)?$')


# Classes -------------------------------------------------------------------------------------------------------------
class AssemblyError(Exception):
    pass


class Assembly:
    def __init__(self, path: Path | None = None, name: str | None = None, contigs: dict[str: Contig] | None = None):
        self.path = path or Path()
        self.name = name or path.name.strip('.gz').rsplit('.', 1)[0] if path else ''
        self.contigs = contigs or {}

    def __repr__(self):
        return self.name

    def __len__(self):
        return sum(len(i) for i in self.contigs.values())

    def seq(self, ctg: str, start: int, end: int, strand: str = "+") -> Seq:
        return self.contigs[ctg].seq[start:end] if strand == "+" else self.contigs[ctg].seq[
                                                                      start:end].reverse_complement()

    def map(self, stdin: str, threads: int, ) -> Generator[Alignment, None, None]:
        return iter_alns(Popen(f"minimap2 -c -t {threads} {self.path} -".split(), stdin=PIPE, stdout=PIPE,
                               stderr=PIPE).communicate(stdin.encode())[0].decode())


class ContigError(Exception):
    pass


class Contig(object):
    """
    This class describes a contig in an assembly: the name, length, and sequence.
    """

    def __init__(self, name: str | None = None, desc: str | None = None, seq: Seq | None = Seq('')):
        self.name = name or ''
        self.desc = desc or ''
        self.seq = seq
        # self.neighbours_L = neighbours_L or {}  # For assembly graphs
        # self.neighbours_R = neighbours_R or {}  # For assembly graphs

    def __repr__(self):
        return self.name

    def __len__(self):
        return len(self.seq)


# Functions -----------------------------------------------------------------------------------------------------------
def parse_assembly(file: Path | str, verbose: bool = False) -> Assembly | None:
    """Parse an assembly file and return an Assembly object"""
    if file := check_file(file):  # Check the file exists, warn if not (instead of quitting)
        if _ASSEMBLY_FASTA_REGEX.search(file.name):
            assembly = Assembly(file)
            try:
                with opener(file, mode='rt') as f:
                    for header, seq in SimpleFastaParser(f):
                        name, description = header.split(' ', 1) if ' ' in header else (header, '')
                        assembly.contigs[name] = Contig(name, description, Seq(seq))
            except Exception as e:
                return warning(f"Error parsing {file.name}\n{e}")
            log(f"Parsed {assembly} as FASTA", verbose=verbose)
            return assembly
        return warning(f"File extension must match {_ASSEMBLY_FASTA_REGEX.pattern}: {file.name}")


# class GFAError(Exception):
#     pass
#
#
# def parse_gfa(file: Path, **kwargs) -> Assembly:
#     """Parse a GFA file and return an Assembly object"""
#     links = []  # Store links as we need to ensure all ctgs are added first
#     with opener(file) as f, NamedTemporaryFile(mode='wt') as tmp:
#         for line in f:
#             if (parts := line.strip().split('\t')):  # See https://github.com/GFA-spec/GFA-spec for GFA format spec
#                 if parts[0] == 'S':
#                     tmp.write(f">{parts[1]}\n{parts[2]}\n")
#                 elif parts[0] == 'L':
#                     links.append((parts[1], parts[2], parts[3], parts[4]))
#         # tmp.flush()
#         try:
#             assembly = Assembly(file, aligner=mp.Aligner(tmp.name, **kwargs))
#         except Exception as e:
#             raise GFAError(f"Error parsing {file}: {e}")
#     for fr, to, from_orient, to_orient in links:
#         try:
#             from_ctg, to_ctg = assembly.contigs[fr], assembly.contigs[to]
#         except KeyError as e:
#             raise GFAError(f"Cannot link {fr} -> {to} in {file}: {e}")
#         from_dict = from_ctg.neighbours_R if from_orient == '+' else from_ctg.neighbours_L
#         to_dict = to_ctg.neighbours_L if to_orient == '+' else to_ctg.neighbours_R
#         from_dict[to_ctg.name] = to_ctg
#         to_dict[from_ctg.name] = from_ctg
#     return assembly


def parse_result(line: str, db: Database, regex: Pattern | None = None, samples: set[str] | None = None,
                 loci: set[str] | None = None) -> TypingResult | None:
    if regex and not regex.search(line):
        return None
    try:
        d = loads(line)
    except Exception as e:
        warning(f"Error parsing JSON line: {e}\n{line}")
        return None
    if samples and d['sample_name'] not in samples:
        return None
    if loci and d['best_match'] not in loci:
        return None
    try:
        return TypingResult.from_dict(d, db)
    except Exception as e:
        warning(f"Error converting JSON line to TypingResult: {e}")
        return None


def write_headers(tsv: TextIO | None = None, no_header: bool = False, scores: bool = False) -> None:
    """Write the headers to a file handle."""
    if tsv and not no_header and (tsv.name == '<stdout>' or fstat(tsv.fileno()).st_size == 0):
        if scores:
            tsv.write('Assembly\tLocus\tAS\tmlen\tblen\tq_len\tgenes_found\tgenes_expected\n')
        else:
            tsv.write(
                'Assembly\tBest match locus\tBest match type\tMatch confidence\tProblems\tIdentity\tCoverage\t'
                'Length discrepancy\tExpected genes in locus\tExpected genes in locus, details\t'
                'Missing expected genes\tOther genes in locus\tOther genes in locus, details\t'
                'Expected genes outside locus\tExpected genes outside locus, details\t'
                'Other genes outside locus\tOther genes outside locus, details\t'
                'Truncated genes, details\tExtra genes, details\n'
            )


def typing_pipeline(
        assembly: str | Path | Assembly, db: str | Path | Database, threads: None | int = 0,
        score_metric: int | None = 0, weight_metric: int | None = 3, min_cov: float | None = 50,
        n_best: int | None = 2, max_other_genes: int | None = 1, percent_expected_genes: float | None = 50,
        allow_below_threshold: bool | None = False, score_file: TextIO | None = None, verbose: bool | None = False
) -> TypingResult | None:
    """
    Performs *in silico* serotyping on a bacterial genome assembly using a database of known loci.
    :param assembly: Path to the assembly file or Assembly object
    :param db: Path to the database file or Database object
    :param threads: Number of threads to use for alignment
    :param score_metric: Alignment score to use: 0=AS, 1=mlen, 2=blen, 3=q_len
    :param weight_metric: Score weighting metric: 0=None, 1=Genes found, 2=Genes expected, 3=Prop genes, 4=blen, 5=q_len
    :param min_cov: Minimum coverage for a gene to be used for scoring
    :param n_best: Number of best loci from the 1st round of scoring to be fully aligned to the assembly
    :param max_other_genes: Max other genes to allow in the best locus to be considered Typeable
    :param percent_expected_genes: Percent of expected genes required to be considered Typeable
    :param allow_below_threshold: Allow genes below the threshold to be considered Typeable
    :param score_file: File handle to write the scores to, will not type the assembly if provided
    :param verbose: Print progress to stderr
    :return: TypingResult object or None
    """
    # CHECK ARGS -------------------------------------------------------------------------------------------------------
    threads = check_cpus(threads, verbose=verbose) if not threads else threads  # Get the number of threads if 0
    if not isinstance(db, Database) and not (db := load_database(db, verbose=verbose)):
        return None
    if not isinstance(assembly, Assembly) and not (assembly := parse_assembly(assembly, verbose=verbose)):
        return None

    # ALIGN GENES ------------------------------------------------------------------------------------------------------
    # TODO: Consider creating a Scores object to store the scores and methods for scoring/writing
    scores = np.zeros((len(db), 6))  # Init scores array with 6 columns: AS, mlen, blen, q_len, genes_found, genes_expected
    idx, alignments = {l: i for i, l in enumerate(db.loci)}, []  # Set up index for loci and list for alignments
    for q, alns in group_alns(assembly.map(db.format('ffn'), threads)):  # Group alignments by query gene (Alignment.q)
        if q.startswith("Extra"):
            alignments.append(max(alns, key=lambda x: x.mlen))  # Add the best alignment for extra genes
        else:
            alignments += (alns := list(alns))  # Add all alignments to the list
            # Use the best alignment for each gene for scoring, if the coverage is above the minimum
            if ((best := max(alns, key=lambda x: x.mlen)).blen / best.q_len) * 100 >= min_cov:
                scores[idx[q.split("_", 1)[0]]] += [best.tags['AS'], best.mlen, best.blen, best.q_len, 1, 0]
            # For each gene, add: AS, mlen, blen, q_len, genes_found (1), genes_expected (0 but will update later)

    if scores.max() == 0:  # If no gene alignments were found, return None so pipeline can continue
        return warning(f'No gene alignments sufficient for typing {assembly}')

    # SCORE LOCI -------------------------------------------------------------------------------------------------------
    scores[:, 5] = np.array([len(l.genes) for l in db.loci.values()])  # Add expected genes to the score matrix

    if score_file:  # If we are just scoring the assembly
        score_file.write(  # Write the scores to the file
            ''.join([f"{assembly}\t{k}\t" + '\t'.join(map(str, v)) + '\n' for k, v in zip(idx.keys(), scores)])
        )
        return log(f"Finished scoring {assembly}", verbose=verbose)  # Return without typing the assembly

    # Process the scores to get the best loci to fully align, this collapses the matrix to a 1D array
    if weight_metric:  # If we are using a weighted score
        if weight_metric == 1:
            scores = scores[:, score_metric] / scores[:, 4]  # Genes found
        elif weight_metric == 2:
            scores = scores[:, score_metric] / scores[:, 5]  # Genes expected
        elif weight_metric == 3:
            scores = scores[:, score_metric] * (scores[:, 4] / scores[:, 5])  # Prop genes
        elif weight_metric == 4:
            scores = scores[:, score_metric] / scores[:, 2]  # blen
        elif weight_metric == 5:
            scores = scores[:, score_metric] / scores[:, 3]  # q_len
    else:
        scores = scores[:, score_metric]  # Unweighted score

    best_loci = [db[int(i)] for i in np.argsort(scores)[::-1][:n_best]]  # Get the best loci to fully align
    scores, idx = np.zeros((len(best_loci), 4)), {l.name: i for i, l in enumerate(best_loci)}  # Init scores and index
    locus_alignments = {l.name: [] for l in best_loci}  # Init dict to store alignments for each locus
    for locus, alns in group_alns(assembly.map(''.join(i.format('fna') for i in best_loci), threads)):  # Group by locus
        for a in alns:  # For each alignment of the locus
            scores[idx[locus]] += [a.tags['AS'], a.mlen, a.blen, a.q_len]  # Add the alignment to the scores
            locus_alignments[locus].append(a)  # Add the alignment to the locus alignments
    best_match = best_loci[np.argmax(scores[:, score_metric])]  # Get the best match based on the highest score

    # RECONSTRUCT LOCUS ------------------------------------------------------------------------------------------------
    result = TypingResult(assembly.name, db, best_match)  # Create the result object
    pieces = {  # Init dict to store pieces for each contig
        ctg: [LocusPiece(ctg, result, s, e) for s, e in  # Create pieces for each merged contig range
              merge_ranges([(a.r_st, a.r_en) for a in alns], len(db.largest_locus))]  # Merge ranges by largest locus
        for ctg, alns in group_alns(locus_alignments[best_match.name], key='ctg')  # Group by contig
    }  # We can't add strand as the pieces may be merged from multiple alignments, we will determine from the genes

    # GET GENE RESULTS -------------------------------------------------------------------------------------------------
    for a in cull_filtered(lambda i: i.q in best_match.genes, alignments):  # For each non-overlapping gene alignment
        if (gene := best_match.genes.get(a.q)):  # Get gene reference from database and gene type
            gene_type = "expected_genes"
        elif (gene := db.extra_genes.get(a.q)):
            gene_type = "extra_genes"
        else:
            gene = db.genes.get(a.q)
            gene_type = "unexpected_genes"

        # Get Piece if gene range overlaps with a piece
        piece = next(filter(lambda p: range_overlap((p.start, p.end), (a.r_st, a.r_en)) > 0,
                            pieces.get(a.ctg, [])), None)
        # Create gene result and extract sequence from assembly
        gene_result = GeneResult(a.ctg, gene, result, piece, a.r_st, a.r_en, a.strand, gene_type=gene_type,
                                 partial=a.partial, dna_seq=assembly.seq(a.ctg, a.r_st, a.r_en, a.strand))
        # Evaluate the gene in protein space by comparing the translation to the reference gene
        gene_result.compare_translation(table=11, to_stop=True)  # This will also trigger the protein alignment
        gene_result.below_threshold = gene_result.percent_identity < db.gene_threshold  # Check if below threshold
        if not piece and gene_result.below_threshold:  # If below protein identity threshold
            continue  # Skip this gene, probably a homologue in another part of the genome
        result.add_gene_result(gene_result)  # Add the gene result to the result to get neighbouring genes
        # previous_result = gene_result  # Set the previous gene result to the current result

    # FINALISE PIECES --------------------------------------------------------------------------------------------------
    for ctg, pieces in pieces.items():  # Add sequences to pieces and add them to the result
        for piece in pieces:
            if piece.expected_genes:  # If the piece has expected genes
                piece.strand = "+" if max(i.strand == i.gene.strand for i in piece.expected_genes) else "-"
                # Piece strand is consensus of expected gene strands
                piece.sequence = assembly.seq(ctg, piece.start, piece.end, piece.strand)
                result.pieces.append(piece)  # Add the piece to the result

    # FINALISE RESULT -------------------------------------------------------------------------------------------------
    # Sort the pieces by the sum of the expected gene order to get the expected order of the pieces
    result.pieces.sort(key=lambda x: min(int(i.gene.name.split("_")[1]) for i in x.expected_genes))
    [l.sort(key=lambda x: int(x.gene.name.split("_")[1])) for l in (
        result.expected_genes_inside_locus, result.expected_genes_outside_locus, result.unexpected_genes_inside_locus,
        result.unexpected_genes_outside_locus)]
    result.missing_genes = sorted(list(set(best_match.genes) - {
        i.gene.name for i in chain(result.expected_genes_inside_locus, result.expected_genes_outside_locus)}),
                                  key=lambda x: int(x.split("_")[1]))
    result.get_confidence(allow_below_threshold, max_other_genes, percent_expected_genes)
    log(f"Finished typing {result}", verbose=verbose)
    return result
