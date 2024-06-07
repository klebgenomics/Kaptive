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

from Bio.Seq import Seq
from Bio.SeqIO.FastaIO import SimpleFastaParser
import numpy as np
np.seterr(divide='ignore', invalid='ignore')  # Ignore divide by zero and invalid value errors

from kaptive.typing import TypingResult, LocusPiece, GeneResult
from kaptive.database import Database, load_database
from kaptive.alignment import Alignment, iter_alns, group_alns, cull_filtered
from kaptive.misc import opener, merge_ranges, range_overlap, check_cpus, rtrn, check_file
from kaptive.log import warning, log

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
        return self.path.stat().st_size

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
def parse_assembly(file: Path | str, **kwargs) -> Assembly | None:
    """Parse an assembly file and return an Assembly object"""
    if not (file := check_file(file, warning)):  # Check the file exists, warn if not (instead of quitting)
        return None
    if _ASSEMBLY_FASTA_REGEX.search(file.name):
        try:
            assembly = Assembly(file)
            with opener(file, mode='rt') as f:
                for h, s in SimpleFastaParser(f):
                    n, d = h.split(' ', 1) if ' ' in h else (h, '')
                    assembly.contigs[n] = Contig(n, d, Seq(s))
            return rtrn(f"Parsed {assembly} as FASTA", return_obj=assembly, **kwargs)
        except Exception as e:
            return rtrn(f"Error parsing {file.name}\n{e}", warning)
    # elif _ASSEMBLY_GRAPH_REGEX.search(file.name):
    #     log(f"Parsing {file.name} as GFA", verbose=verbose)
    #     return parse_gfa(file, **kwargs)
    return rtrn(f"Unknown assembly format for {file.name}", warning)


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
        warning(f"Error converting JSON line to TypingResult: {e}\n{line}")
        return None


def write_headers(tsv: TextIO | None = None, no_header: bool = False, scores: bool = False) -> None:
    """Write the headers to a file handle."""
    if tsv:
        if tsv.name != '<stdout>' and tsv.tell() != 0:  # If file is path and not already written to
            no_header = True  # Headers already written, useful for running on HPC
        if not no_header:
            if scores:
                tsv.write('Assembly\tLocus\tAS\tmlen\tblen\tq_len\tgenes_found\tgenes_expected\n')
            else:
                tsv.write(
                    'Assembly\tBest match locus\tBest match type\tConfidence\tProblems\tIdentity\tCoverage\t'
                    'Length discrepancy\tExpected genes in locus\tExpected genes in locus, details\t'
                    'Missing expected genes\tOther genes in locus\tOther genes in locus, details\t'
                    'Expected genes outside locus\tExpected genes outside locus, details\t'
                    'Other genes outside locus\tOther genes outside locus, details\t'
                    'Truncated genes, details\tExtra genes\n'
                )


def typing_pipeline(
        assembly: str | Path | Assembly, db: str | Path | Database, threads: None | int = 0,
        score_metric: int | None = 0, weight_metric: int | None = 3, fallback_weight: int | None = 5,
        min_cov: float | None = 50, max_full: int | None = 1, max_other_genes: int | None = 1,
        percent_expected_genes: float | None = 50, allow_below_threshold: bool | None = False,
        score_file: TextIO | None = None, verbose: bool | None = False) -> TypingResult | None:
    """
    Performs *in silico* serotyping on a bacterial genome assembly using a database of known loci.
    It can be run as part of batch process (preloaded assembly and database) or as a standalone function.
    """
    # CHECK ARGS -------------------------------------------------------------------------------------------------------
    threads = check_cpus(threads) if not threads else threads  # Get the number of threads if 0
    if not isinstance(db, Database) and not (db := load_database(db, verbose=verbose)):
        return None
    if not isinstance(assembly, Assembly) and not (assembly := parse_assembly(assembly, verbose=verbose)):
        return None

    # assembly, db = parse_assembly('test/data/GCA_000584355.1_ASM58435v1_genomic.fna'), load_database('ab_k')
    # assembly, db = parse_assembly('test/kpsc/2018-01-1001_20.fasta'), load_database('kpsc_k')
    # threads, score_metric, weight_metric, max_full, min_cov, verbose = check_cpus(), 0, 3, 2, 50, True

    # ALIGN GENES ------------------------------------------------------------------------------------------------------
    scores = np.zeros((len(db), 6))  # Init scores array with 6 columns
    idx, alignments = {l: i for i, l in enumerate(db.loci)}, []  # Set up index for loci and list for alignments
    for q, alns in group_alns(assembly.map(db.format('ffn'), threads)):  # Group alignments by query gene (Alignment.q)
        if q.startswith("Extra"):
            alignments.append(max(alns, key=lambda x: x.mlen))  # Add the best alignment for extra genes
        else:
            alignments += (alns := list(alns))  # Add all alignments to the list
            # Use the best alignment for each gene for scoring, if the coverage is above the minimum
            if ((best := max(alns, key=lambda x: x.mlen)).blen / best.q_len) * 100 >= min_cov:
                scores[idx[q.split("_", 1)[0]]] += [best.tags['AS'], best.mlen, best.blen, best.q_len, 1, 0]

    if scores.max() == 0:
        return rtrn(f'No gene alignments sufficient for typing {assembly}', warning)

    # SCORE LOCI -------------------------------------------------------------------------------------------------------
    scores[:, 5] = np.array([len(l.genes) for l in db.loci.values()])  # Add expected genes to the score matrix

    if score_file:
        score_file.write(
            ''.join([f"{assembly}\t{k}\t" + '\t'.join(map(str, v)) + '\n' for k, v in zip(idx.keys(), scores)])
        )
        return rtrn(f"Finished scoring {assembly}", verbose=verbose)

    if (prop_genes := (scores[:, 4] / scores[:, 5])).max() * 100 < percent_expected_genes:
        warning(f"Max proportion of genes is < {percent_expected_genes}% for {assembly}, using fallback weight")
        weight_metric = fallback_weight

    # Matrix cols correspond to score_metric: 0 = AS, 1 = mlen, 2 = blen, 3 = q_len, 4 = genes found, 5 = genes expected
    # Weight arguments: 0 = None, 1 = Genes found, 2 = Genes expected, 3 = Prop genes, 4 = blen, 5 = q_len
    if weight_metric:  # If we are using a weighted score
        if weight_metric == 1:
            scores = scores[:, score_metric] / scores[:, 4]  # Genes found
        elif weight_metric == 2:
            scores = scores[:, score_metric] / scores[:, 5]  # Genes expected
        elif weight_metric == 3:
            scores = scores[:, score_metric] * prop_genes  # Prop genes
        elif weight_metric == 4:
            scores = scores[:, score_metric] / scores[:, 2]  # blen
        elif weight_metric == 5:
            scores = scores[:, score_metric] / scores[:, 3]  # q_len
    else:
        scores = scores[:, score_metric]  # Unweighted score

    best_loci = [db[int(i)] for i in np.argsort(scores)[::-1][:max_full]]
    scores, idx = np.zeros((len(best_loci), 4)), {l.name: i for i, l in enumerate(best_loci)}
    locus_alignments = {l.name: [] for l in best_loci}
    for locus, alns in group_alns(assembly.map(''.join(i.format('fna') for i in best_loci), threads)):
        for a in alns:
            scores[idx[locus]] += [a.tags['AS'], a.mlen, a.blen, a.q_len]
            locus_alignments[locus].append(a)
    best_match = best_loci[np.argmax(scores[:, score_metric])]  # Get the best match based on the highest score

    # RECONSTRUCT LOCUS ------------------------------------------------------------------------------------------------
    result = TypingResult(assembly.name, db, best_match)
    pieces = {  # Align best match seq to assembly, group alignments by ctg and create pieces
        ctg: [LocusPiece(ctg, result, s, e) for s, e in
              merge_ranges([(a.r_st, a.r_en) for a in alns], len(db.largest_locus))]
        for ctg, alns in group_alns(locus_alignments[best_match.name], key='ctg')
    }  # We can't add strand as the pieces may be merged from multiple alignments, we will determine from the genes

    # GET GENE RESULTS -------------------------------------------------------------------------------------------------
    for a in cull_filtered(lambda i: i.q in best_match.genes, alignments):  # For each gene alignment
        # Get gene reference from database and gene type
        if (gene := best_match.genes.get(a.q)):
            gene_type = "expected_genes"
        elif (gene := db.extra_genes.get(a.q)):
            gene_type = "extra_genes"
        else:
            gene = db.genes.get(a.q)
            gene_type = "unexpected_genes"

        # Get Piece if gene range overlaps with a piece
        piece = next(filter(lambda p: range_overlap((p.start, p.end), (a.r_st, a.r_en)) > 0,
                            pieces.get(a.ctg, [])), None)
        # Create gene result and add it to the result
        gene_result = GeneResult(  # Create gene result and add it to the result
            a.ctg, gene, result, piece, a.r_st, a.r_en, a.strand,
            gene_type=gene_type, partial=a.partial, dna_seq=assembly.seq(a.ctg, a.r_st, a.r_en, a.strand)
        )
        gene_result.compare_translation(table=11, to_stop=True)  # This will also trigger the protein alignment
        gene_result.below_threshold = gene_result.percent_identity < db.gene_threshold  # Check if below threshold
        if not piece and gene_result.below_threshold:  # If below threshold
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
    return rtrn(f"Finished typing {result}", return_obj=result, verbose=verbose)
