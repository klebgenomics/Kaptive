from pathlib import Path
from itertools import chain
from json import loads
from subprocess import Popen, PIPE
from typing import TextIO, Pattern, Generator, Union, Optional
from re import compile as re_compile
from os import fstat
from gzip import open as gzip_open
from warnings import warn

from Bio.Seq import Seq
from Bio.SeqIO.FastaIO import SimpleFastaParser
import numpy as np

np.seterr(divide='ignore', invalid='ignore')  # Ignore divide by zero and invalid value errors

from kaptive.typing import TypingResult, LocusPiece, GeneResult
from kaptive.database import Database, load_database
from kaptive.alignment import Alignment, group_alns, cull_filtered
from kaptive.utils import Resources, LiteralFile
from kaptive.interval import range_overlap, merge_ranges


# Constants -----------------------------------------------------------------------------------------------------------
# _ASSEMBLY_FASTA_REGEX = compile(r'\.(fasta|fa|fna|ffn)(\.gz|\.bz2|\.xz)?$')
_ASSEMBLY_FASTA_REGEX = re_compile(r'\.(fasta|fa|fna|ffn)(\.gz)?$')
_ASSEMBLY_HEADER = ('Assembly\tBest match locus\tBest match type\tMatch confidence\tProblems\tIdentity\tCoverage\t'
                    'Length discrepancy\tExpected genes in locus\tExpected genes in locus, details\t'
                    'Missing expected genes\tOther genes in locus\tOther genes in locus, details\t'
                    'Expected genes outside locus\tExpected genes outside locus, details\t'
                    'Other genes outside locus\tOther genes outside locus, details\t'
                    'Truncated genes, details\tExtra genes, details\n')
_SCORES_HEADER = 'Assembly\tLocus\tAS\tmlen\tblen\tq_len\tgenes_found\tgenes_expected\n'


# Exceptions and warnings ----------------------------------------------------------------------------------------------
class TyperWarning(UserWarning):
    pass


# Classes --------------------------------------------------------------------------------------------------------------
class Contig(object):
    """
    This class describes a contig in an assembly: the name, length, and sequence.
    """

    def __init__(self, name: str, desc: str, seq: Seq):
        self.name = name
        self.desc = desc
        self.seq = seq

    def __repr__(self):
        return self.name

    def __len__(self):
        return len(self.seq)


class Assembly:
    def __init__(self, path_: Path, name: str, contigs: dict[str, Contig] = None):
        self.path = path_
        self.name = name
        self.contigs = contigs or {}

    def __repr__(self):
        return self.name

    def __len__(self):
        return sum(len(i) for i in self.contigs.values())

    def seq(self, ctg: str, start: int, end: int, strand: str = "+") -> Seq:
        return self.contigs[ctg].seq[start:end] if strand == "+" else self.contigs[ctg].seq[
                                                                      start:end].reverse_complement()

    def map(self, query: str, threads: int, extra_args: str = '') -> Generator[Alignment, None, None]:
        cmd = "minimap2 -c " + (f"{extra_args} " if extra_args else '') + f'-t {threads} "{self.path}" -'
        stdout, stderr = Popen(cmd, stdin=PIPE, stdout=PIPE, stderr=PIPE, universal_newlines=True,
                               shell=True).communicate(query)
        for line in stdout.splitlines():
            yield Alignment.from_paf_line(line)


# Functions -----------------------------------------------------------------------------------------------------------
def parse_assembly(path: Union[str, Path]) -> Optional[Assembly]:
    """Parse an assembly file and return an Assembly object"""
    path = LiteralFile(path)
    if match := _ASSEMBLY_FASTA_REGEX.search(path.name):
        assembly = Assembly(path, path.name.rstrip(match.group()))
        try:
            with (gzip_open if path.suffix == '.gz' else open)(path, mode='rt') as f:
                for header, seq in SimpleFastaParser(f):
                    header = header.split(maxsplit=1)
                    name, description = header if len(header) == 2 else (header[0], '')
                    assembly.contigs[name] = Contig(name, description, Seq(seq))
        except Exception as e:
            raise RuntimeError(f"Error parsing {path}") from e
        return assembly
    raise RuntimeError(f"File extension must match {_ASSEMBLY_FASTA_REGEX.pattern}: {path}")


def parse_result(line: str, db: Database, regex: Pattern = None, samples: set[str] = None,
                 loci: set[str] = None) -> Optional[TypingResult]:
    if regex and not regex.search(line):
        return None
    try:
        d = loads(line)
    except Exception as e:
        raise RuntimeError(f"Error parsing JSON line") from e
    if samples and d['sample_name'] not in samples:
        return None
    if loci and d['best_match'] not in loci:
        return None
    try:
        return TypingResult.from_dict(d, db)
    except Exception as e:
        raise RuntimeError(f"Error converting JSON line to TypingResult") from e


def write_headers(tsv: TextIO = None, no_header: bool = False, scores: bool = False) -> int:
    """Write appropriate header to a file handle."""
    if tsv and not no_header and (tsv.name == '<stdout>' or fstat(tsv.fileno()).st_size == 0):
        return tsv.write(_SCORES_HEADER if scores else _ASSEMBLY_HEADER)


def typing_pipeline(
        assembly: Union[str, Path, Assembly], db: Union[str, Path, Database], threads: int = 1,
        score_metric: int = 0, weight_metric: int = 3, min_cov: float = 50, n_best: int = 2,
        max_other_genes: int = 1, percent_expected_genes: float = 50, allow_below_threshold: bool = False,
        score_file: TextIO = None) -> Optional[TypingResult]:
    """
    Performs *in silico* serotyping on a bacterial genome assembly using a database of known loci.
    :param assembly: Path to the assembly file or Assembly object
    :param db: Path to the database file or Database object
    :param threads: Number of threads to use for alignment
    :param score_metric:  score to use: 0=AS, 1=mlen, 2=blen, 3=q_len
    :param weight_metric: Score weighting metric: 0=None, 1=Genes found, 2=Genes expected, 3=Prop genes, 4=blen, 5=q_len
    :param min_cov: Minimum coverage for a gene to be used for scoring
    :param n_best: Number of top loci from the 1st round of scoring to be fully aligned to the assembly
    :param max_other_genes: Max other genes to allow in the best locus to be considered Typeable
    :param percent_expected_genes: Percent of expected genes required to be considered Typeable
    :param allow_below_threshold: Allow genes below the threshold to be considered Typeable
    :param score_file: File handle to write the scores to, will not type the assembly if provided
    :return: TypingResult object or None
    """
    # CHECK ARGS -------------------------------------------------------------------------------------------------------
    if not isinstance(db, Database) and not (db := load_database(db)):
        return None
    if not isinstance(assembly, Assembly) and not (assembly := parse_assembly(assembly)):
        return None

    # ALIGN GENES ------------------------------------------------------------------------------------------------------
    # Init scores array with 6 columns: AS, mlen, blen, q_len, genes_found, genes_expected
    scores, alignments = np.zeros((len(db), 6)), []
    # Group alignments by query gene (Alignment.q)
    for q, alns in group_alns(assembly.map(db.format('ffn'), threads)):
        if q.startswith("Extra_genes"):
            alignments.append(max(alns, key=lambda x: x.mlen))  # Add the best alignment for extra genes
        else:
            alignments.extend(alns := list(alns))  # Add all alignments to the list, convert generator to list too
            # Use the best alignment for each gene for scoring, if the coverage is above the minimum
            if ((best := max(alns, key=lambda x: x.mlen)).blen / best.q_len) * 100 >= min_cov:
                scores[db.loci[q.split('_', 1)[0]].index] += [
                    best.tags['AS'], best.mlen, best.blen, best.q_len, 1, 0]
            # For each gene, add: AS, mlen, blen, q_len, genes_found (1), genes_expected (0 but will update later)

    if scores.max() == 0:  # If no gene alignments were found, return None so pipeline can continue
        return warn(f'No gene alignments sufficient for typing {assembly}\n'
                       f'Have you used the appropriate database for your species?', TyperWarning)

    # SCORE LOCI -------------------------------------------------------------------------------------------------------
    scores[:, 5] = db.expected_gene_counts  # Add expected genes to the 6th column (0-based) score matrix

    if score_file:  # If we are just scoring the assembly
        score_file.write(  # Write the scores to the file
            ''.join([f"{assembly}\t{k}\t" + '\t'.join(map(str, v)) + '\n' for k, v in zip(db.loci.keys(), scores)])
        )
        return None

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

    best_loci = [db[int(i)] for i in
                 np.argsort(scores)[::-1][:min(n_best, len(scores))]]  # Get the best loci to fully align
    scores, idx = np.zeros((len(best_loci), 4)), {l.name: i for i, l in enumerate(best_loci)}  # Init scores and index
    locus_alignments = {l.name: [] for l in best_loci}  # Init dict to store alignments for each locus
    # Group alignments by locus
    for locus, alns in group_alns(assembly.map(''.join(i.format('fna') for i in best_loci), threads)):
        for a in alns:  # For each alignment of the locus
            scores[idx[locus]] += [a.tags['AS'], a.mlen, a.blen, a.q_len]  # Add alignment metrics to the scores
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
        if gene := best_match.genes.get(a.q):  # Get gene reference from database and gene type
            gene_type = "expected_genes"
        elif gene := db.extra_genes.get(a.q):
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
    result.pieces.sort(key=lambda x: min(i.gene.start for i in x.expected_genes))
    [l.sort(key=lambda x: gene.start) for l in (
        result.expected_genes_inside_locus, result.expected_genes_outside_locus, result.unexpected_genes_inside_locus,
        result.unexpected_genes_outside_locus)]
    result.missing_genes = list(set(best_match.genes) - {
        i.gene.name for i in chain(result.expected_genes_inside_locus, result.expected_genes_outside_locus)
    })
    result.get_confidence(allow_below_threshold, max_other_genes, percent_expected_genes)
    return result
