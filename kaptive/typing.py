"""
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

import re
from functools import cached_property
from json import dumps, loads
from typing import Generator
import warnings
from itertools import chain

from Bio.Seq import Seq
from dna_features_viewer import GraphicFeature, GraphicRecord

from kaptive.database import Database, Locus, Gene
from kaptive.log import warning


# Classes -------------------------------------------------------------------------------------------------------------
class TypingResultError(Exception):
    pass


class TypingResult:
    def __init__(self, db: Database | None = None, best_match: Locus | None = None,
                 score: float | None = 0, zscore: float | None = 0, pieces: list[LocusPiece] | None = None,
                 expected_genes_inside_locus: list[GeneResult] | None = None,
                 expected_genes_outside_locus: list[GeneResult] | None = None, missing_genes: list[str] | None = None,
                 unexpected_genes_inside_locus: list[GeneResult] | None = None,
                 unexpected_genes_outside_locus: list[GeneResult] | None = None,
                 extra_genes: list[GeneResult] | None = None, is_elements: list[GeneResult] | None = None,
                 truncated_genes: list[GeneResult] | None = None, scores: dict[str, dict[str, list[float]]] | None = None,
                 min_zscore: float | None = 0, min_cov: float | None = 0, gene_threshold: float | None = 0,
                 score_metric: str | None = None, weight_metric: str | None = None):
        self.db = db
        self.best_match = best_match
        self.score = score
        self.zscore = zscore
        self.pieces = pieces or []  # Pieces of locus reconstructed from alignments
        self.expected_genes_inside_locus = expected_genes_inside_locus or []  # Genes from best_match
        self.expected_genes_outside_locus = expected_genes_outside_locus or []  # Genes from best_match
        self.missing_genes = missing_genes or []  # Genes from best_match that were not found
        self.unexpected_genes_inside_locus = unexpected_genes_inside_locus or []  # Genes from other loci
        self.unexpected_genes_outside_locus = unexpected_genes_outside_locus or []  # Genes from other loci
        self.extra_genes = extra_genes or []  # in db.extra_genes, ALWAYS outside locus (gene_result.piece == None)
        self.truncated_genes = truncated_genes or []  # Genes that are truncated
        self.is_elements = is_elements or []  # in db.is_elements, ALWAYS inside locus (gene_result.piece != None)
        self.scores = scores or {}
        self.min_zscore = min_zscore
        self.min_cov = min_cov
        self.gene_threshold = gene_threshold
        self.score_metric = score_metric or ""
        self.weight_metric = weight_metric or ""

    def __len__(self):
        return sum(len(i) for i in self.pieces) if self.pieces else 0

    def add_gene_result(self, gene_result: GeneResult):
        if gene_result.neighbour_left:  # If gene_result.neighbour_left is not None, the gene is not the first gene
            gene_result.neighbour_left.neighbour_right = gene_result
        if gene_result.piece:  # If gene_result.piece is not None, the gene is inside the locus
            gene_result.piece.add_gene_result(gene_result)
            gene_type = f"{gene_result.gene_type}{'_inside_locus' if gene_result.gene_type.startswith(('expected', 'unexpected')) else ''}"
        else:  # If gene_result.piece is None, the gene is outside the locus
            gene_type = f"{gene_result.gene_type}{'_outside_locus' if gene_result.gene_type.startswith(('expected', 'unexpected')) else ''}"
        getattr(self, gene_type).append(gene_result)  # Add gene result to the appropriate list

    @cached_property  # Cache the percent identity so it is only calculated once
    def percent_identity(self) -> float:
        return (sum(i.percent_identity for i in x) / len(x)) if (x := self.expected_genes_inside_locus) else 0

    @cached_property  # Cache the percent coverage so it is only calculated once
    def percent_coverage(self) -> float:
        return sum(i.percent_coverage for i in x) / (y := len(x)) * (y / len(self.best_match.genes)) if (x := self.expected_genes_inside_locus) else 0

    @cached_property  # Cache the confidence so it is only calculated once
    def confidence(self) -> str:
        if self.zscore >= self.min_zscore:
            if len(self.pieces) == 1:
                if any(i.below_threshold for i in self.expected_genes_inside_locus):
                    return "Novel"
                else:
                    return "Medium" if self.missing_genes else "High"
            else:
                return "Low" if self.missing_genes else "Medium"
        else:
            return "Untypable"

    @cached_property  # Cache the phenotype so it is only calculated once
    def phenotype(self) -> str:
        gene_phenotypes = set()  # Init set to store gene phenotypes to be used as a key in the phenotypes dict
        for gene in chain(self.expected_genes_inside_locus, self.unexpected_genes_inside_locus,
                            self.expected_genes_outside_locus, self.unexpected_genes_outside_locus, self.extra_genes):
            gene.get_phenotype()  # Because neighbouring genes need to be calculated first to test for truncation
            if gene.phenotype == 'truncated':
                self.truncated_genes.append(gene)
            if gene.gene_type in {'expected_genes', 'extra_genes'}:  # The reported phenotype only considers expected
                gene_phenotypes.add((gene.gene.name, gene.phenotype))  # or extra genes
        # NOTE: The best_match.phenotypes MUST be sorted from largest to smallest gene set to make sure any sets with
        # extra genes are tested first.
        for gene_set, phenotype in self.best_match.phenotypes:  # Get the phenotype from the phenotypes dict
            if len(gene_set) == len(gene_phenotypes.intersection(gene_set)):
                return phenotype
        return self.best_match.phenotypes[0][1]  # If no phenotype is found, return the first (default) one

    @cached_property
    def problems(self) -> str:
        problems = '?' if len(self.pieces) > 1 else ''
        problems += '-' if self.missing_genes else ''
        problems += '+' if self.unexpected_genes_inside_locus else ''
        problems += '*' if any(i.below_threshold for i in self.expected_genes_inside_locus) else ''
        # problems += '!' if any(i.phenotype == 'truncated' for i in self.expected_genes_inside_locus) else ''
        return problems

    def as_fasta(self) -> str:
        """Returns a fasta-formatted nucleotide sequence of the locus with a newline character at the end."""
        return "".join(i.as_fasta() for i in self.pieces)

    def as_gene_fasta(self) -> str:
        """Returns a fasta-formatted nucleotide sequence of the locus genes with a newline character at the end."""
        return "".join(i.as_fasta() for i in self.pieces)

    def as_protein_fasta(self) -> str:
        """
        Returns a fasta-formatted protein sequence of the locus genes with a newline character at the end."""
        return "".join(i.as_protein_fasta() for i in self.pieces)

    def as_GraphicRecord(self) -> GraphicRecord:
        features, start = [], 0
        for piece in self.pieces:
            features.extend(piece.as_GraphicFeatures(start))
            start += len(piece)
        return GraphicRecord(sequence_length=self.__len__(), first_index=0, features=features,
                             sequence=[p.sequence for p in self.pieces])


class AssemblyTypingResult(TypingResult):
    def __init__(self, assembly: 'Assembly' | None = None, **kwargs):
        super().__init__(**kwargs)
        self.assembly = assembly

    # @classmethod
    # def from_dict(cls, result_dict: dict, db: Database):
    #     try:
    #         from kaptive.assembly import Assembly, Contig
    #         assembly = Assembly(name=result_dict['Assembly'], contigs=
    #         [assembly.contigs(Contig(name=contig)) for contig in {i['Contig'] for i in result_dict['Contig pieces']}]
    #         self = cls(db=db, assembly=assembly, best_match=db.loci[result_dict['Best match locus']],
    #                    missing_expected_genes=[db.genes[i] for i in result_dict['Missing expected genes']])
    #         for contig_piece in result_dict['Contig pieces']:
    #             contig_piece = ContigPiece.from_dict(contig_piece, self)
    #             for gene_result in contig_piece.gene_results.values():
    #                 self.gene_results[gene_result.gene.name] = gene_result
    #                 if gene_result.gene.name in self.best_match.genes:
    #                     self.expected_genes_inside_locus.append(gene_result)
    #                 else:
    #                     self.other_genes_inside_locus.append(gene_result)
    #             self.pieces.append(contig_piece)
    #         for gene_result in result_dict['Expected genes outside locus']:
    #             gene_result = AssemblyGeneResult.from_dict(gene_result, self)
    #             self.gene_results[gene_result.gene.name] = gene_result
    #             self.expected_genes_outside_locus.append(gene_result)
    #         for gene_result in result_dict['Other genes outside locus']:
    #             gene_result = AssemblyGeneResult.from_dict(gene_result, self)
    #             self.gene_results[gene_result.gene.name] = gene_result
    #             self.other_genes_outside_locus.append(gene_result)
    #     except ValueError:
    #         raise TypingResultError(f"Invalid dict: {result_dict}")
    #     return self

    def as_table(self, debug: bool = False) -> str:
        return '\t'.join(
            [
                self.assembly.name,
                self.best_match.name,
                self.phenotype,
                self.confidence,
                self.problems,
                f"{self.percent_identity:.2f}%",
                f"{self.percent_coverage:.2f}%",
                f"{self.__len__() - len(self.best_match)} bp" if len(self.pieces) == 1 else 'n/a',
                f"{(x := len(self.expected_genes_inside_locus))} / {(y := len(self.best_match.genes))} ({100 * x / y:.2f}%)",
                ';'.join(str(i) for i in x) if (x := self.expected_genes_inside_locus) else '',
                ';'.join(self.missing_genes),
                f"{len((x := self.unexpected_genes_inside_locus))}",
                ';'.join(str(i) for i in x) if x else '',
                f"{(x := len(self.expected_genes_outside_locus))} / {(y := len(self.best_match.genes))} ({100 * x / y:.2f}%)",
                ';'.join(str(i) for i in x) if (x := self.expected_genes_outside_locus) else '',
                f"{len(self.unexpected_genes_outside_locus)}",
                ';'.join(str(i) for i in x) if (x := self.unexpected_genes_outside_locus) else '',
                ';'.join([str(i) for i in self.truncated_genes]),
            ] + ([] if not debug else [  # Add debug columns if debug is True
                ';'.join([str(i) for i in self.extra_genes]),  # Extra genes
                # ';'.join([str(i) for i in self.is_elements]),
                str(len({i.contig.name for i in self.pieces})),  # Number of contigs
                str(len(self.pieces)),  # Number of contig pieces
                ';'.join(str(i) for i in self.pieces),  # Contig pieces
                f"{self.score:.1f}",  # Best score
                f"{self.zscore:.1f}",  # Best zscore
                ';'.join(
                    # Only print the top 2 scores per score metric and weight metric: score_weight:locus,score,zscore|...
                    f"{s}_{w}:{'|'.join(f'{l},{x:.4f},{y:.4f}' for l, x, y in self.scores[s][w][:2])}" for s in
                    self.scores for w in self.scores[s]
                ),
                f'min_zscore={self.min_zscore};min_cov={self.min_cov};gene_threshold={self.db.gene_threshold};'
                f'score_metric={self.score_metric};weight_metric={self.weight_metric}',
            ])
        ) + "\n"

    def as_json(self, **kwargs) -> str:
        """Returns a JSON string representation of the result with a newline character at the end."""
        return dumps({
            'Assembly': self.assembly.name,
            'Best match locus': self.best_match.name,
            'Best match type': self.phenotype,
            'Problems': self.problems,
            'Percent identity': f"{self.percent_identity:.2f}%",
            'Percent coverage': f"{self.percent_coverage:.2f}%",
            'Contig pieces': [i.as_dict() for i in self.pieces],
            # All genes inside the locus will be reported in Contig pieces, so only need missing and outside locus
            'Missing expected genes': self.missing_genes,
            'Expected genes outside locus': [i.as_dict() for i in self.expected_genes_outside_locus],
            'Unexpected genes outside locus': [i.as_dict() for i in self.unexpected_genes_outside_locus],
        }, **kwargs) + "\n"


class LocusPieceError(Exception):
    pass


class LocusPiece:
    def __init__(self, result: AssemblyTypingResult | None = None, start: int | None = 0, end: int | None = 0,
                 strand: str | None = None, expected_genes: list[GeneResult] | None = None,
                 unexpected_genes: list[GeneResult] | None = None,
                 is_elements: list[GeneResult] | None = None, sequence: Seq | None = None):
        self.result = result
        self.start = start or 0
        self.end = end or 0
        self.expected_genes = expected_genes or []  # Genes from best_match
        self.unexpected_genes = unexpected_genes or []  # Genes that were found from other loci
        self.is_elements = is_elements or []  # Specifically db.is_elements
        self.strand = strand or "unknown"
        self.sequence = sequence or Seq("")

    def __len__(self):
        return self.end - self.start

    def __iter__(self):
        return iter(self.expected_genes + self.unexpected_genes + self.is_elements)

    # def as_dict(self):
    #     self.extract_sequence()
    #     return {
    #         'Contig': self.contig.name,
    #         'Start': self.start + 1,  # Add 1 to start to make it 1-based
    #         'End': self.end,
    #         'Strand': self.strand,
    #         'Sequence': str(self.sequence),
    #         'Genes': [gene.as_dict() for gene in self]
    #     }

    def as_fasta(self) -> str:
        return f">{self.__repr__()}\n{self.sequence}\n"

    def add_gene_result(self, gene_result: AssemblyGeneResult):
        assert gene_result.piece == self  # Shouldn't be necessary, but just in case
        if gene_result.start < self.start:  # Update start and end if necessary
            self.start = gene_result.start
        if gene_result.end > self.end:
            self.end = gene_result.end
        getattr(self, gene_result.gene_type).append(gene_result)

    def as_GraphicFeatures(self, relative_start: int = 0) -> Generator[GraphicFeature, None, None]:
        start, end = relative_start, relative_start + len(self)
        yield GraphicFeature(start=start, end=end, strand=1, thickness=20,
                             color=('tab:blue', self.result.percent_identity / 100), label=self.__repr__())
        for gene in self:
            # Get relative gene start within piece
            gene_start = start + (gene.start - self.start) if gene.strand == "+" else end - (gene.end - self.start)
            gene_end = start + (gene.end - self.start) if gene.strand == "+" else end - (gene.start - self.start)
            if self.strand == "+":
                strand = gene.gene.strand if gene.strand == gene.gene.strand else gene.strand
            else:
                strand = gene.gene.strand if gene.strand != gene.gene.strand else gene.strand
            yield GraphicFeature(
                start=gene_start, end=gene_end, strand=0 if gene.phenotype == "truncated" else 1 if strand == "+" else -1,
                color=("green" if gene.gene_type == 'expected_genes' else "orange", gene.percent_identity / 100),
                linecolor='red' if gene.below_threshold else "yellow" if gene.phenotype == "truncated" else 'black',
                legend_text=gene.gene_type, label=gene.__repr__()
            )


class ContigPiece(LocusPiece):
    def __init__(self, contig: 'Contig', *args, **kwargs):
        """Represents a part of a contig derived from alignments within a certain proximity to each other."""
        super().__init__(*args, **kwargs)
        self.contig = contig

    def __repr__(self):
        return f"{self.contig.name}:{self.start}-{self.end}{self.strand}"

    # @classmethod
    # def from_dict(cls, contig_piece_dict: dict, result: AssemblyTypingResult, **kwargs):
    #     return cls(
    #         result=result, contig=result.assembly.contigs[contig_piece_dict['Contig']],
    #         start=contig_piece_dict['Start'] - 1, end=contig_piece_dict['End'], strand=contig_piece_dict['Strand'],
    #         sequence=Seq(contig_piece_dict['Sequence']),
    #         gene_results={i['Gene']: AssemblyGeneResult.from_dict(i, result) for i in contig_piece_dict['Genes']},
    #         **kwargs
    #     )


# class ReadsPiece(LocusPiece):
#     def __init__(self, read_group: 'ReadGroup' | None = None, reads: dict[str, list[Alignment]] | None = None,
#                  **kwargs):
#         super().__init__(**kwargs)
#         self.read_group = read_group
#         self.reads = reads or {}
#
#     def add_gene_result(self, gene_result: ReadsGeneResult):
#         gene_result.contig_piece = self
#         self.gene_results[gene_result.gene.name] = gene_result
#         self.start = min(g.relative_locus_start for g in self.gene_results.values())
#         self.end = max(g.relative_locus_end for g in self.gene_results.values())
#         for read, alignments in gene_result.reads.items():
#             self.reads[read] = self.reads.get(read, []) + alignments
#             self.alignments += alignments


class GeneResultError(Exception):
    pass


class GeneResult:
    """
    Class to store alignment results for a single gene in a locus for either a ReadResult or a AssemblyResult.
    """

    def __init__(self, gene: Gene, piece: LocusPiece | None = None, start: int | None = 0, end: int | None = 0,
                 strand: str | None = None, neighbour_left: GeneResult | None = None,
                 neighbour_right: GeneResult | None = None, dna_seq: Seq | None = None, protein_seq: Seq | None = None,
                 below_threshold: bool | None = False, phenotype: str | None = None, gene_type: str | None = None,
                 percent_identity: float | None = 0, percent_coverage: float | None = 0):
        self.gene = gene
        self.start = start
        self.end = end
        self.strand = strand
        self.piece = piece  # inside locus if not None
        self.neighbour_left = neighbour_left  # neighbour to the left of the gene
        self.neighbour_right = neighbour_right  # neighbour to the right of the gene
        self.dna_seq = dna_seq or Seq("")
        self.protein_seq = protein_seq or Seq("")
        self.below_threshold = below_threshold or False
        self.phenotype = phenotype or ""
        self.gene_type = gene_type or ""
        self.percent_identity = percent_identity
        self.percent_coverage = percent_coverage

    def __repr__(self):
        return self.gene.name

    def __len__(self):
        return self.end - self.start

    def extract_translation(self, table: int = 11, cds: bool = False, to_stop: bool = True, gap: str = '-',
                            stop_symbol: str = '*'):
        """
        Extracts the translation from the DNA sequence of the gene result. Implemented as a method so
        unnecessary translations are not performed.
        :param table: NCBI translation table number
        :param cds: if True, only translates the CDS
        :param to_stop: if True, stops translation at the first stop codon
        :param gap: gap character
        :param stop_symbol: stop codon character
        """
        if len(self.protein_seq) == 0:  # Only translate if the protein sequence is not already stored
            if len(self.dna_seq) == 0:
                raise GeneResultError(f'No DNA sequence for {self.__repr__()}')
            with warnings.catch_warnings(record=True) as w:
                self.protein_seq = self.dna_seq.translate(table=table, cds=cds, to_stop=to_stop, gap=gap,
                                                          stop_symbol=stop_symbol)
                # for i in w:
                #     warning(f"{i.message}: {self.__repr__()}")
            if len(self.protein_seq) == 0:
                warning(f'No protein sequence for {self.__repr__()}')

    def as_fasta(self) -> str:
        """Returns a fasta-formatted nucleotide sequence with a newline character at the end."""
        if len(self.dna_seq) == 0:
            warning(f'No DNA sequence for {self}')
            return ""
        return f'>{self.__repr__()}\n{self.dna_seq}\n'

    def as_protein_fasta(self) -> str:
        """Returns a fasta-formatted protein sequence with a newline character at the end."""
        try:  # Try to extract the translation or return an empty string, some genes (or IS) may be severely truncated
            self.extract_translation()
        except GeneResultError:
            warning(f'Could not translate {self}')
            return ""
        return f'>{self.__repr__()}\n{self.protein_seq}\n'


class AssemblyGeneResult(GeneResult):
    def __init__(self, contig: 'Contig', *args, **kwargs):
        """This class represents ONE alignment of a locus gene to an assembly."""
        super().__init__(*args, **kwargs)
        self.contig = contig

    # @classmethod
    # def from_dict(cls, gene_result_dict: dict, result: AssemblyTypingResult, **kwargs):
    #     if gene_result_dict['Gene'] in result.db.genes:
    #         gene = result.db.genes[gene_result_dict['Gene']]
    #     elif gene_result_dict['Gene'] in result.db.is_elements:
    #         gene = result.db.is_elements[gene_result_dict['Gene']]
    #     else:
    #         raise GeneResultError(f"Gene {gene_result_dict['Gene']} not found in database, did you pass --is-seqs?")
    #     return cls(
    #         gene=gene,
    #         alignment=Alignment(
    #             strand=gene_result_dict['Strand'], query_name=gene_result_dict['Gene'],
    #             target_name=gene_result_dict['Contig'], target_start=gene_result_dict['Start'] - 1,
    #             target_end=gene_result_dict['End'], percent_identity=gene_result_dict['Percent identity'],
    #             percent_query_coverage=gene_result_dict['Percent coverage']
    #         ),
    #         dna_seq=Seq(gene_result_dict['Sequence']),
    #         contig=result.assembly.contigs[gene_result_dict['Contig']],
    #         truncated_protein=gene_result_dict['Truncated protein'],
    #         below_threshold=gene_result_dict['Below threshold'],
    #         **kwargs
    #     )

    @cached_property
    def edge_of_contig(self) -> bool:
        return any((self.start == 0, self.end == len(self.contig)))

    def __str__(self) -> str:
        s = f'{self.gene.name},{self.percent_identity:.2f}%,{self.percent_coverage:.2f}%'
        s += f',truncated({(len(self.protein_seq) / len(self.gene.protein_seq)) * 100:.2f}%)' if self.phenotype == "truncated" else ""
        s += ",edge" if self.edge_of_contig else ""
        s += ",below_threshold" if self.below_threshold else ""
        # s += f",{self.phenotype}" if self.gene_type == "is_elements" else ""
        return s

    def get_phenotype(self, **kwargs):
        """
        Calculates the phenotype of the gene result, e.g. "present" / "truncated".
        IS elements will not be tested for presence/truncation but if they are next to
        any genes that are truncated.
        """
        # if self.gene_type == 'is_elements':  # Special phenotype for IS elements
        #     # Test to see if it interrupts any surrounding genes on a Piece
        #     if self.neighbour_left and self.neighbour_left.phenotype == "truncated":
        #         self.phenotype += f"{self.neighbour_left.gene.name}"
        #     if self.neighbour_right and self.neighbour_right.phenotype == "truncated":
        #         self.phenotype += ("," if self.phenotype else "") + f"{self.neighbour_right.gene.name}"
        # else:
        if self.percent_coverage == 100:  # If the gene is not fully covered
            self.extract_translation(**kwargs)  # Extract the translation
            self.gene.extract_translation(**kwargs)  # Extract the translation of the gene
            if len(self.protein_seq) < len(self.gene.protein_seq):  # If the protein is truncated
                self.phenotype = "truncated"
        # else:
        #     if not self.edge_of_contig:  # If the gene is not at the edge of the contig
        #         self.phenotype = "truncated"
        if not self.phenotype:
            self.phenotype = "present"

    def as_dict(self):
        return {
            'Gene': self.gene.name,
            'Contig': self.contig.name,
            'Start': self.start + 1,  # Add 1 to start to make it 1-based
            'End': self.end,
            'Strand': self.strand,
            'Percent identity': self.percent_identity,
            'Percent coverage': self.percent_coverage,
            'Phenotype': self.phenotype,
            'Edge of contig': self.edge_of_contig,
            'Below threshold': self.below_threshold,
            'Sequence': str(self.dna_seq),
        }

    # def as_feature(self) -> list[SeqFeature]:
    #     try:  # Try to extract the translation or return an empty string, some genes (or IS) may be severely truncated
    #         self.extract_translation()
    #     except GeneResultError:
    #         warning(f'Could not translate {self}')
    #     loc = FeatureLocation(self.start, self.end, strand=1 if self.strand == "+" else -1)
    #     return [
    #         SeqFeature(
    #             location=loc, type="CDS", id=self.gene.name, qualifiers={
    #                 "gene": [self.gene.gene_name], "locus_tag": [self.gene.name], "product": [self.gene.product],
    #                 "translation": [str(self.protein_seq)], "translation_table": [11], "codon_start": [1],
    #                 "note": [f"percent_identity={self.percent_identity:.2f}%, percent_coverage={self.percent_coverage:.2f}%"],
    #                 'inference': 'alignment:Kaptive', 'protein_id': f"gnl|Kaptive|{self.gene.name}"
    #             }
    #         ),
    #         SeqFeature(location=loc, type="gene", id=self.gene.name,
    #                    qualifiers={"gene": [self.gene.gene_name], "locus_tag": [self.gene.name]})
    #     ]


# class ReadsTypingResult(TypingResult):
#     """
#     Result class to hold results for a single locus
#     """
#
#     def __init__(self, db: Database | None = None, read_group: 'ReadGroup' | None = None,
#                  truncated_genes: list[ReadsGeneResult] | None = None,
#                  expected_genes_inside_locus: list[ReadsGeneResult] | None = None,
#                  expected_genes_outside_locus: list[ReadsGeneResult] | None = None,
#                  other_genes_inside_locus: list[ReadsGeneResult] | None = None,
#                  other_genes_outside_locus: list[ReadsGeneResult] | None = None,
#                  gene_results: dict[str: ReadsGeneResult] | None = None,
#                  reads: dict[str, list[Alignment]] | None = None, **kwargs
#                  ):
#
#         super().__init__(db, **kwargs)
#         self.read_group = read_group
#         self.expected_genes_inside_locus = expected_genes_inside_locus or []
#         self.expected_genes_outside_locus = expected_genes_outside_locus or []
#         self.other_genes_inside_locus = other_genes_inside_locus or []
#         self.other_genes_outside_locus = other_genes_outside_locus or []
#         self.gene_results = gene_results or {}
#         self.truncated_genes = truncated_genes or []
#         self.reads = reads or {}
#
#
# class ReadsGeneResult(GeneResult):
#     def __init__(self, alignments: list[Alignment] | None = None, truncated_protein: bool | None = None,
#                  below_threshold: bool | None = None, reads: set[str] | None = None,
#                  alignment_ranges: list[tuple[int, int]] | None = None, base_coverage: int | None = 0,
#                  relative_locus_start: int | None = 0, relative_locus_end: int | None = 0, **kwargs):
#         """
#         This class represents alignments of reads (queries) to the locus gene (target).
#         There are multiple alignments that may extend past the gene,
#         so coverage and identity are calculated a little differently.
#         The locus can be pieced together from the alignments.
#         """
#         super().__init__(**kwargs)
#         self.alignments = alignments or []
#         self.reads = reads or {}
#         self.truncated_protein = truncated_protein or False
#         self.below_threshold = below_threshold or False
#         self.alignment_ranges = alignment_ranges or []  # List of tuples of start and end positions of alignments
#         self.base_coverage = base_coverage  # Number of bases covered by alignments
#         self.relative_locus_start = relative_locus_start  # Start of gene relative to locus
#         self.relative_locus_end = relative_locus_end  # End of gene relative to locus
#
#     @classmethod
#     def from_alignments(cls, gene: Gene, alignments: list[Alignment], **kwargs):
#         """
#         Creates a ReadsGeneResult for a single gene from a list of alignments against the locus.
#         """
#         self = cls(gene=gene, **kwargs)
#         for alignment in sorted(alignments, key=lambda x: x.target_start):
#             self.alignment_ranges.append((alignment.target_start, alignment.target_end))  # Add alignment range
#             self.reads[alignment.query_name] = self.reads.get(alignment.query_name, []) + [alignment]  # Add read
#             self.alignments.append(alignment)  # Add alignment
#         self.alignment_ranges = merge_ranges(self.alignment_ranges)  # Merge overlapping alignment ranges
#         self.base_coverage += sum(end - start for start, end in self.alignment_ranges)
#         self.relative_locus_start = self.alignment_ranges[0][0] + self.gene.start  # Start of gene relative to locus
#         self.relative_locus_end = self.alignment_ranges[-1][1] + self.gene.end  # End of gene relative to locus
#         return self


def parse_results(json, db: Database, regex: re.Pattern | None, samples: list[str] | None, loci: list[str] | None
                  ) -> Generator[AssemblyTypingResult, None, None]:
    with json.open() as f:
        for json_line in f.readlines():
            if regex and not regex.search(json_line):
                continue
            try:
                result_dict = loads(json_line)
                if samples and result_dict['Assembly'] not in samples:
                    continue
                if loci and result_dict['Best match locus'] not in loci:
                    continue
                yield AssemblyTypingResult.from_dict(result_dict, db)
            except ValueError:
                raise f"Invalid JSON: {json_line}"
