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

from Bio.Seq import Seq
# from Bio.SeqFeature import SeqFeature, FeatureLocation
# from Bio.Graphics import GenomeDiagram
# from reportlab.lib.colors import black

from kaptive.intrange import merge_ranges
from kaptive.database import Database, Locus, Gene
from kaptive.alignment import Alignment, cull_conflicting_alignments, cull_all_conflicting_alignments
from kaptive.assembly import Assembly, Contig
from kaptive.log import warning


# Classes -------------------------------------------------------------------------------------------------------------
class ResultError(Exception):
    pass


class Result:
    def __init__(self, db: Database | None = None, best_match: Locus | None = None,
                 extra_genes: list[GeneResult] | None = None, gene_results: dict[str, GeneResult] | None = None,
                 truncated_genes: list[GeneResult] | None = None, missing_expected_genes: list[Gene] | None = None,
                 scores: dict[str, float] | None = None):
        self.db = db
        self.best_match = best_match
        self.missing_expected_genes = missing_expected_genes or []
        self.extra_genes = extra_genes or []
        self.gene_results = gene_results or {}
        self.truncated_genes = truncated_genes or []
        self.scores = scores or {}


class AssemblyResult(Result):
    def __init__(self, db: Database | None = None, assembly: 'Assembly' | None = None,
                 contig_pieces: list[ContigPiece] | None = None,
                 expected_genes_inside_locus: list[AssemblyGeneResult] | None = None,
                 expected_genes_outside_locus: list[AssemblyGeneResult] | None = None,
                 other_genes_inside_locus: list[AssemblyGeneResult] | None = None,
                 other_genes_outside_locus: list[AssemblyGeneResult] | None = None, **kwargs):

        super().__init__(db, **kwargs)
        self.assembly = assembly
        self.contig_pieces = contig_pieces or []
        self.expected_genes_inside_locus = expected_genes_inside_locus or []
        self.expected_genes_outside_locus = expected_genes_outside_locus or []
        self.other_genes_inside_locus = other_genes_inside_locus or []
        self.other_genes_outside_locus = other_genes_outside_locus or []

    @classmethod
    def from_best_locus(
            cls, db: Database, best_locus: Locus, assembly: 'Assembly', best_locus_alignments: list[Alignment],
            other_locus_alignments: list[Alignment | None], is_element_alignments: list[Alignment | None],
            gene_threshold: float = 95, gene_distance: int = 10000, **kwargs):
        self = cls(
            db=db, assembly=assembly, best_match=best_locus, missing_expected_genes=[
                v for k, v in best_locus.genes.items() if k not in [i.query_name for i in best_locus_alignments]],
            **kwargs
        )
        # Cull all overlapping alignments of genes which are not in the best locus
        other_locus_alignments = cull_all_conflicting_alignments(other_locus_alignments)
        for alignment in best_locus_alignments:  # Then cull other genes if they are overlapped by expected genes
            other_locus_alignments = list(cull_conflicting_alignments(alignment, other_locus_alignments))

        contig_results = {}  # Dict to store alignments on each contig
        for alignment in best_locus_alignments + other_locus_alignments:  # Iterate over all remaining alignments
            if alignment.target_name not in contig_results:
                contig_results[alignment.target_name] = []
            contig_results[alignment.target_name].append(alignment)

        for alignment in cull_all_conflicting_alignments(is_element_alignments):
            if alignment.target_name in contig_results:  # We only care about alignments that are in the locus
                contig_results[alignment.target_name].append(alignment)

        # Construct piece(s) of locus on each contig where the best locus gene alignments were found
        gene_pos = {g: n for n, g in enumerate(best_locus.genes, start=1)}  # Get expected positions of expected genes
        for contig, alignments in contig_results.items():  # Iterate over each contig
            # Collect the alignment ranges from each expected gene on the contig if they are in the expected order
            ranges, previous_alignment, contig = [], None, assembly[contig]
            for alignment in sorted(alignments, key=lambda x: x.target_start):  # Sort by start, can be on either strand
                if alignment.query_name in gene_pos:  # Test if the alignment is for an expected gene
                    if previous_alignment:  # Test if the position of alignment follows expected order
                        # Do this by checking if difference between positions is 1 or less
                        if abs(gene_pos[previous_alignment.query_name] - gene_pos[alignment.query_name]) <= 1:
                            # If so, add the previous alignment to the ranges
                            ranges.append((previous_alignment.target_start, previous_alignment.target_end))
                    previous_alignment = alignment  # Update previous alignment
            if previous_alignment:  # Add the last alignment to the ranges
                ranges.append((previous_alignment.target_start, previous_alignment.target_end))

            if ranges:  # Merge together ranges of expected genes in expected order with a tolerance from gene_distance
                contig_pieces = [  # Create a ContigPiece for each range
                    ContigPiece(assembly=assembly, contig=contig, start=start, end=end, result=self)
                    for start, end in merge_ranges(ranges, tolerance=gene_distance)]
            else:
                contig_pieces = []
            for alignment in alignments:  # Add gene alignments on this contig to the contig pieces and the result
                contig_piece, category = None, None
                below_threshold = alignment.percent_identity < gene_threshold
                gene, is_element = (db.is_elements[alignment.query_name], True) \
                    if alignment.query_name in db.is_elements else (db.genes[alignment.query_name], False)

                # Find the contig piece that the alignment belongs to
                contig_piece = next((p for p in contig_pieces if p.start <= alignment.target_start <= alignment.target_end <= p.end), None)
                if not contig_piece:  # OUTSIDE locus
                    if not is_element:  # Exclude IS elements outside the locus
                        category = 'expected_genes_outside_locus' if alignment.query_name in self.best_match.genes \
                            else 'other_genes_outside_locus'  # Add to the correct category
                else:  # INSIDE locus
                    category = 'expected_genes_inside_locus' if alignment.query_name in self.best_match.genes \
                        else 'other_genes_inside_locus'  # Add to the correct category

                if category:  # Gene result will be created and reported, then added to the correct category
                    self.gene_results[gene_result.gene.name] = (gene_result := AssemblyGeneResult(
                        gene, alignment, contig, below_threshold=below_threshold))

                    if gene_result.alignment.percent_query_coverage == 100 and not is_element:
                        gene_result.compare_translation()  # compare the translation to check for truncation

                    getattr(self, category).append(gene_result)  # Add to the correct category using attribute name
                    alignment.query_name.startswith('Extra') and self.extra_genes.append(gene_result)  # Add to extra
                    gene_result.truncated_protein and self.truncated_genes.append(gene_result)  # Add to truncated
                    contig_piece and contig_piece.add_gene_result(gene_result)  # Add gene result to the contig piece

            for contig_piece in contig_pieces:  # Add the contig pieces to the result
                contig_piece.guess_strand()  # Guess the strand of the contig piece once all alignments are added
                self.contig_pieces.append(contig_piece)
        return self

    @classmethod
    def from_dict(cls, result_dict: dict, db: Database):
        try:
            assembly = Assembly(name=result_dict['Assembly'])
            [assembly.add_contig(Contig(name=contig)) for contig in {i['Contig'] for i in result_dict['Contig pieces']}]
            self = cls(db=db, assembly=assembly, best_match=db.loci[result_dict['Best match locus']],
                       missing_expected_genes=[db.genes[i] for i in result_dict['Missing expected genes']])
            for contig_piece in result_dict['Contig pieces']:
                contig_piece = ContigPiece.from_dict(contig_piece, self)
                for gene_result in contig_piece.gene_results.values():
                    self.gene_results[gene_result.gene.name] = gene_result
                    if gene_result.gene.name in self.best_match.genes:
                        self.expected_genes_inside_locus.append(gene_result)
                    else:
                        self.other_genes_inside_locus.append(gene_result)
                self.contig_pieces.append(contig_piece)
            for gene_result in result_dict['Expected genes outside locus']:
                gene_result = AssemblyGeneResult.from_dict(gene_result, self)
                self.gene_results[gene_result.gene.name] = gene_result
                self.expected_genes_outside_locus.append(gene_result)
            for gene_result in result_dict['Other genes outside locus']:
                gene_result = AssemblyGeneResult.from_dict(gene_result, self)
                self.gene_results[gene_result.gene.name] = gene_result
                self.other_genes_outside_locus.append(gene_result)
        except ValueError:
            raise ResultError(f"Invalid dict: {result_dict}")
        return self

    def __len__(self):
        return sum([len(i) for i in self.contig_pieces]) if self.contig_pieces else 0

    @cached_property  # Cache the phenotype so it is only calculated once
    def phenotype(self) -> str:
        phenotype = self.best_match.type
        if phenotype == 'special logic':  # If there is special logic, apply it here
            g = {i.gene.gene_name for i in self.extra_genes}  # Names of all the extra genes that were found
            for locus_type, gene_names in self.best_match.special_logic.items():
                # Use set operations to match all extra loci, this way they don't need to be sorted
                if gene_names == g or (len(g) == 0 and len(gene_names) == 0):
                    phenotype = locus_type
                    break
        if self.truncated_genes:  # Check if CPS phenotype dependant genes are truncated
            for gene_result in self.truncated_genes:
                for core_gene in ["wcaJ", "wbaP"]:
                    if (x := gene_result.gene) and x.gene_name == core_gene:
                        phenotype += f" ({core_gene} truncated)"

        if not phenotype:  # Assert that the phenotype is not None
            raise ResultError(f"Phenotype is None: {self}")
        return phenotype

    @cached_property  # Cache the percent identity so it is only calculated once
    def percent_identity(self) -> float:
        return sum(i.alignment.percent_identity for i in self.expected_genes_inside_locus) / \
               len(self.expected_genes_inside_locus) if self.expected_genes_inside_locus else 0

    @cached_property  # Cache the percent coverage so it is only calculated once
    def percent_coverage(self) -> float:
        return sum(i.alignment.percent_query_coverage for i in self.expected_genes_inside_locus) / \
               len(self.expected_genes_inside_locus) if self.expected_genes_inside_locus else 0

    @cached_property  # Cache the problems string so it is only calculated once
    def problems(self) -> str:
        problems = "-" if self.missing_expected_genes else ""
        problems += "+" if self.other_genes_inside_locus else ""
        problems += "?" if len(self.contig_pieces) > 1 else ""
        problems += "*" if any(i.below_threshold for i in self.expected_genes_inside_locus) else ""
        problems += "!" if self.truncated_genes else ""
        return problems

    def as_list(self) -> list[str]:
        return [
            self.assembly.name,
            self.best_match.name,
            self.phenotype,
            self.problems,
            f"{self.percent_coverage:.2f}%",
            f"{self.percent_identity:.2f}%",
            f"{self.__len__() - len(self.best_match)} bp" if len(self.contig_pieces) == 1 else 'n/a',
            f"{(x := len(self.expected_genes_inside_locus))} / {(y := len(self.best_match.genes))} ({100 * x / y:.2f}%)",
            f"{';'.join(str(i) for i in self.expected_genes_inside_locus) if self.expected_genes_inside_locus else ''}",
            f"{';'.join(i.name for i in self.missing_expected_genes) if self.missing_expected_genes else ''}",
            f"{len(self.other_genes_inside_locus)}",
            f"{';'.join(str(i) for i in self.other_genes_inside_locus) if self.other_genes_inside_locus else ''}",
            f"{(x := len(self.expected_genes_outside_locus))} / {(y := len(self.best_match.genes))} ({100 * x / y:.2f}%)",
            f"{';'.join(str(i) for i in self.expected_genes_outside_locus) if self.expected_genes_outside_locus else ''}",
            f"{len(self.other_genes_outside_locus)}",
            f"{';'.join(str(i) for i in self.other_genes_outside_locus) if self.other_genes_outside_locus else ''}"
        ]

    def as_fasta(self) -> str:
        """Returns a fasta-formatted nucleotide sequence of the locus with a newline character at the end."""
        return "".join(i.as_fasta() for i in self.contig_pieces)

    def as_gene_fasta(self) -> str:
        """Returns a fasta-formatted nucleotide sequence of the locus genes with a newline character at the end."""
        return "".join(i.as_fasta() for i in self.gene_results.values())

    def as_protein_fasta(self) -> str:
        """
        Returns a fasta-formatted protein sequence of the locus genes with a newline character at the end."""
        return "".join(i.as_protein_fasta() for i in self.gene_results.values())

    def as_json(self, **kwargs) -> str:
        """Returns a JSON string representation of the result with a newline character at the end."""
        return dumps({
            'Assembly': self.assembly.name,
            'Best match locus': self.best_match.name,
            'Best match type': self.phenotype,
            'Problems': self.problems,
            'Percent identity': f"{self.percent_identity:.2f}%",
            'Percent coverage': f"{self.percent_coverage:.2f}%",
            'Contig pieces': [i.as_dict() for i in self.contig_pieces],
            # All genes inside the locus will be reported in Contig pieces, so only need missing and outside locus
            'Missing expected genes': [i.name for i in self.missing_expected_genes],
            'Expected genes outside locus': [i.as_dict() for i in self.expected_genes_outside_locus],
            'Other genes outside locus': [i.as_dict() for i in self.other_genes_inside_locus]
        }, **kwargs) + "\n"

    # def as_diagram(self):
    #     d = GenomeDiagram.Diagram(f"{self.assembly} {self.best_match}")
    #     t = d.new_track(track_level=1, name="test", greytrack=False, start=0, end=len(self.best_match),
    #                     scale_ticks=0)
    #     s = t.new_set()
    #     for gene_result in self.expected_genes_inside_locus + self.other_genes_inside_locus:
    #         s.add_feature(gene_result.gene.feature, sigil=get_gene_shape(gene_result), label_size=20,
    #                       label_angle=20, color=get_gene_colour(gene_result), label=True,
    #                       name=gene_result.gene.name, label_position="middle",
    #                       label_strand=1, arrowshaft_height=1.0, border=black)
    #     t.add_set(s)
    #     d.draw(format="linear", pagesize=(1800, 200), x=0, yt=0, yb=0, y=0, fragments=1, start=0, end=len(self.best_match))
    #     d.write(f"kaptive_{self.assembly.name}_{self.best_match}.svg", "SVG")


class ContigPieceError(Exception):
    pass


class ContigPiece:
    def __init__(self, contig: 'Contig' | None = None, strand: str | None = None, result: AssemblyResult | None = None,
                 alignments: list[Alignment] | None = None, assembly: 'Assembly' | None = None,
                 gene_results: dict[str, AssemblyGeneResult] | None = None, sequence: Seq | None = None,
                 start: int | None = None, end: int | None = None):
        """
        Represents a part of a contig derived from alignments within a certain proximity to each other.
        """
        self.result = result
        self.contig = contig
        self.assembly = assembly
        self.start = start or 0
        self.end = end or 0
        self.gene_results = gene_results or {}
        self.strand = strand
        self.alignments = alignments or []
        self.sequence = sequence or Seq("")

    @classmethod
    def from_dict(cls, contig_piece_dict: dict, result: AssemblyResult, **kwargs):
        return cls(
            result=result, contig=result.assembly.contigs[contig_piece_dict['Contig']],
            start=contig_piece_dict['Start'] - 1, end=contig_piece_dict['End'], strand=contig_piece_dict['Strand'],
            sequence=Seq(contig_piece_dict['Sequence']),
            gene_results={i['Gene']: AssemblyGeneResult.from_dict(i, result) for i in contig_piece_dict['Genes']},
            **kwargs
        )

    def __len__(self):
        return self.end - self.start

    def __repr__(self):
        return f"{self.contig.name}:{self.start}-{self.end}"

    def __iter__(self):
        return iter(self.gene_results.values())

    def guess_strand(self):
        """
        Guess the piece strand by comparing the consensus strand of gene references compared to the consensus strand
        of the gene alignments.
        """
        expected_strand = "+" if max(i.gene.strand for i in self.gene_results.values() if i.gene) == 1 else "-"
        alignment_strand = max(i.strand for i in self.alignments)
        if alignment_strand == expected_strand:
            self.strand = alignment_strand
        elif alignment_strand == "-" and expected_strand == "+":
            self.strand = "-"
        elif alignment_strand == "+" and expected_strand == "-":
            self.strand = "+"
        else:
            raise ContigPieceError(f"Strand mismatch: {self} has alignments on both strands")

    def add_gene_result(self, gene_result: AssemblyGeneResult):
        gene_result.contig_piece = self
        self.add_alignment(gene_result.alignment)
        self.gene_results[gene_result.gene.name] = gene_result

    def add_alignment(self, alignment: Alignment):
        # Update start and end if necessary
        if alignment.target_start < self.start:
            self.start = alignment.target_start
        if alignment.target_end > self.end:
            self.end = alignment.target_end
        self.alignments.append(alignment)

    def extract_sequence(self):
        if len(self.sequence) == 0:  # Only extract sequence if it is not already stored
            if self.strand == "-":
                self.sequence = self.contig.sequence[self.start:self.end].reverse_complement()
            else:
                self.sequence = self.contig.sequence[self.start:self.end]

    def as_dict(self):
        self.extract_sequence()
        return {
            'Contig': self.contig.name,
            'Start': self.start + 1,  # Add 1 to start to make it 1-based
            'End': self.end,
            'Strand': self.strand,
            'Sequence': str(self.sequence),
            'Genes': [gene.as_dict() for gene in self.gene_results.values()]
        }

    def as_fasta(self) -> str:
        self.extract_sequence()
        return f">{self.__repr__()}{self.strand}\n{self.sequence}\n"


class GeneResultError(Exception):
    pass


class GeneResult:
    """
    Class to store alignment results for a single gene in a locus for either a ReadResult or a AssemblyResult.
    """

    def __init__(self, gene: Gene, dna_seq: Seq | None = None, protein_seq: Seq | None = None):
        self.gene = gene
        self.dna_seq = dna_seq or Seq("")
        self.protein_seq = protein_seq or Seq("")


class AssemblyGeneResult(GeneResult):
    def __init__(self, gene: Gene, alignment: Alignment, contig: 'Contig',
                 contig_piece: ContigPiece | None = None,
                 truncated_protein: bool | None = None, below_threshold: bool | None = None, **kwargs):
        """
        This class represents ONE alignment of a locus gene to an assembly.
        """
        super().__init__(gene, **kwargs)
        self.alignment = alignment
        alignment.add_query_sequence(self.gene.dna_seq)
        self.contig = contig
        self.contig_piece = contig_piece
        self.edge_of_contig = alignment.target_start == 0 or alignment.target_end == alignment.target_length
        self.truncated_protein = truncated_protein or False
        self.below_threshold = below_threshold or False

    @classmethod
    def from_dict(cls, gene_result_dict: dict, result: AssemblyResult, **kwargs):
        if gene_result_dict['Gene'] in result.db.genes:
            gene = result.db.genes[gene_result_dict['Gene']]
        elif gene_result_dict['Gene'] in result.db.is_elements:
            gene = result.db.is_elements[gene_result_dict['Gene']]
        else:
            raise GeneResultError(f"Gene {gene_result_dict['Gene']} not found in database, did you pass --is-seqs?")
        return cls(
            gene=gene,
            alignment=Alignment(
                strand=gene_result_dict['Strand'], query_name=gene_result_dict['Gene'],
                target_name=gene_result_dict['Contig'], target_start=gene_result_dict['Start'] - 1,
                target_end=gene_result_dict['End'], percent_identity=gene_result_dict['Percent identity'],
                percent_query_coverage=gene_result_dict['Percent coverage']
            ),
            dna_seq=Seq(gene_result_dict['Sequence']),
            contig=result.assembly.contigs[gene_result_dict['Contig']],
            truncated_protein=gene_result_dict['Truncated protein'],
            below_threshold=gene_result_dict['Below threshold'],
            **kwargs
        )

    def __str__(self):
        return f'{self.gene.name},{self.alignment.percent_identity:.2f}%,{self.alignment.percent_query_coverage:.2f}%' \
               f'{",trunc" if self.truncated_protein else ""}' \
               f'{",edge" if self.edge_of_contig else ""}'

    def extract_dna_seq(self):
        """
        Extracts the DNA sequence from the alignment and stores it in self.dna_seq.
        """
        if len(self.dna_seq) == 0:  # Only extract sequence if it is not already stored
            if self.alignment and self.contig:
                self.dna_seq = self.contig.sequence[self.alignment.target_start:self.alignment.target_end]
                if self.alignment.strand == "-":
                    self.dna_seq = self.dna_seq.reverse_complement()
            else:
                raise GeneResultError(f'No DNA sequence for {self}')

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
            self.extract_dna_seq()
            self.protein_seq = self.dna_seq.translate(table=table, cds=cds, to_stop=to_stop, gap=gap,
                                                      stop_symbol=stop_symbol)
            if len(self.protein_seq) == 0:
                raise GeneResultError(f'No protein sequence for {self}')

    def compare_translation(self, **kwargs):
        """
        Convenience method to compare the protein sequence of the gene result to the protein sequence of the gene.
        """
        self.extract_translation(**kwargs)
        self.gene.extract_translation(**kwargs)
        if len(self.protein_seq) < len(self.gene.protein_seq):
            self.truncated_protein = True

    def as_dict(self):
        self.extract_dna_seq()
        return {
            'Gene': self.gene.name,
            'Contig': self.alignment.target_name,
            'Start': self.alignment.target_start + 1,  # Add 1 to start to make it 1-based
            'End': self.alignment.target_end,
            'Strand': self.alignment.strand,
            'Percent identity': self.alignment.percent_identity,
            'Percent coverage': self.alignment.percent_query_coverage,
            'Truncated protein': self.truncated_protein,
            'Edge of contig': self.edge_of_contig,
            'Below threshold': self.below_threshold,
            'Sequence': str(self.dna_seq),
        }

    def as_fasta(self) -> str:
        """Returns a fasta-formatted nucleotide sequence with a newline character at the end."""
        self.extract_dna_seq()
        return f'>{self.contig.assembly.name}_{self.gene.name}\n{self.dna_seq}\n'

    def as_protein_fasta(self) -> str:
        """Returns a fasta-formatted protein sequence with a newline character at the end."""
        try:  # Try to extract the translation or return an empty string, some genes (or IS) may be severely truncated
            self.extract_translation()
        except GeneResultError:
            warning(f'Could not translate {self}')
            return ""
        return f'>{self.contig.assembly.name}_{self.gene.name}\n{self.protein_seq}\n'

    # def as_feature(self) -> SeqFeature:
    #     try:  # Try to extract the translation or return an empty string, some genes (or IS) may be severely truncated
    #         self.extract_translation()
    #     except GeneResultError:
    #         warning(f'Could not translate {self}')
    #     return SeqFeature(
    #         type="CDS", id=self.gene.name,
    #         qualifiers={"gene": [self.gene.gene_name], "locus_tag": [self.gene.name], "product": [self.gene.product],
    #                     "translation": [str(self.protein_seq)]},
    #         location=FeatureLocation(self.alignment.target_start, self.alignment.target_end,
    #                                  strand=-1 if self.alignment.strand == "-" else 1)
    #     )


# class ReadsResult(Result):
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
#                  gene_results: dict[str: ReadsGeneResult] | None = None, **kwargs):
#
#         super().__init__(db, **kwargs)
#         self.read_group = read_group
#         self.expected_genes_inside_locus = expected_genes_inside_locus or []
#         self.expected_genes_outside_locus = expected_genes_outside_locus or []
#         self.other_genes_inside_locus = other_genes_inside_locus or []
#         self.other_genes_outside_locus = other_genes_outside_locus or []
#         self.gene_results = gene_results or {}
#         self.truncated_genes = truncated_genes or []
#
#     def add_best_locus(self, db: Database, best_locus: Locus, read_group: 'ReadGroup',
#                        best_locus_alignments: list[Alignment]):
#         self.best_match, self.read_group, self.db = best_locus, read_group, db
#
#         # Add gene results for each gene in the locus if reads overlap the gene range
#         for g in best_locus.genes.values():
#             gene_alignments = [i for i in best_locus_alignments if  i.target_start <= g.gene.feature.location.start <= i.target_end or
#                                i.target_start <= g.gene.feature.location.end <= i.target_end]
#             if gene_alignments:
#                 gene_result = ReadsGeneResult.from_alignments(g, gene_alignments)
#                 self.gene_results[g.name] = gene_result
#
#
# class ReadsGeneResult(GeneResult):
#     def __init__(self, alignments: list[Alignment] | None = None, truncated_protein: bool | None = None,
#                  below_threshold: bool | None = None, **kwargs):
#         """
#         This class represents alignments of reads overlapping a range representing the locus gene.
#         There are multiple alignments that may extend past the gene,
#         so coverage and identity are calculated a little differently.
#         """
#         super().__init__(**kwargs)
#         self.alignments = alignments
#         self.truncated_protein = truncated_protein or False
#         self.below_threshold = below_threshold or False
#
#     @classmethod
#     def from_alignments(cls, gene: Gene, alignments: list[Alignment], **kwargs):
#         """
#         Creates a ReadsGeneResult for a single gene from a list of alignments against the locus.
#         """
#         return cls(gene=gene, alignments=alignments, **kwargs)
#
#     def add_snps(self, vcf: str, all_csqs: bool):
#         """
#         Function to call SNPs in the locus using the snp_calling_pipeline function
#         """
#         # Get consensus fasta
#         # self.consensus_fasta = seq_to_dict(vcf2consensus(vcf, self.locus.get_fasta()))
#         # if len(self.consensus_fasta) > 1:
#         #     raise ValueError("More than one sequence in consensus fasta")
#
#         # Parse vcf output into VcfRecord objects
#         snps = [VcfRecord(record) for record in vcf.splitlines() if record and not record.startswith("#")]
#         # If SNP is part of same locus, add to snp dict
#         self.snps |= {snp.pos: snp for snp in snps if snp.chrom == self.locus.name}
#         # Add snps to genes
#         for gene in self.expected_genes:
#             gene.snps = [snp for snp in snps if snp.gene_id == gene.name]
#             if gene.snps:
#                 gene.mutation = get_mutation(gene.snps, all_mutations=all_csqs)
#
#     def add_is_elements(self, ise_paf: Alignments, ise_range_merge_tolerance: int = 0):
#         # Get reads that aligned to the best locus that are also in the ISE PAF
#         locus_reads = {}
#         for paf in self.alignments.targets[self.locus.name]:
#             if paf.read in ise_paf.queries:
#                 if paf.read in locus_reads:
#                     locus_reads[paf.read].append(paf)
#                 else:
#                     locus_reads[paf.read] = [paf]
#
#         if not locus_reads:
#             return  # Failsafe 1
#
#         # Get is_elements that have reads that aligned to the best locus
#         ise_queries_in_locus = {}
#         for read in locus_reads.keys():
#             for paf in ise_paf.queries[read]:
#                 if paf.target_name in ise_queries_in_locus:
#                     ise_queries_in_locus[paf.target_name].append(paf)
#                 else:
#                     ise_queries_in_locus[paf.target_name] = [paf]
#
#         if not ise_queries_in_locus:
#             return  # Failsafe 2
#
#         # Regions on the locus that are covered by alignments corresponding to is_elements alignments
#         locus_paf = [i for paf in locus_reads.values() for i in paf]
#         locus_ise_regions = target_ranges_covered_by_alignments(locus_paf, ise_range_merge_tolerance)[self.locus.name]
#
#         # Iterate over each region and determine which ISE is in the region
#         for n, region in enumerate(locus_ise_regions):
#             best_ise = None
#             best_overlap = 0
#             # Iterate over each ISE query
#             for ise, pafs in ise_queries_in_locus.items():
#                 # Get all alignments for the ISE that are in the locus
#                 locus_alignments_for_ise = [i for paf in pafs for i in locus_reads[paf.read]]
#
#                 ise_regions = target_ranges_covered_by_alignments(locus_alignments_for_ise, ise_range_merge_tolerance)[self.locus.name]
#
#                 for ise_region in ise_regions:
#                     overlap = range_overlap(ise_region, region)
#                     if overlap > best_overlap:
#                         best_overlap = overlap
#
#                         is_element = ISElement(ise, pafs, self, ise_region)
#                         is_element.locus_region = region
#                         if not best_ise or is_element.identity > best_ise.identity:
#                             best_ise = is_element
#
#             if best_ise:
#                 self.is_elements[f"ise_{n + 1}"] = best_ise
#                 for gene in self.expected_genes:
#                     if range_overlap(region, (int(gene.feature.location.start), int(gene.feature.location.end))):
#                         gene.is_elements.append(best_ise)
#                         best_ise.genes.append(gene)
#             else:
#                 warning(f"No ISE found for region {region} in {self.locus.name}")


# class ISElementResult:
#     """
#     Class that represents alignments of a single ISE corresponding to a region on the locus
#     """
#
#     def __init__(self, name: str, alignments: list[Alignment], result: ReadResult, alignment_range: tuple[int, int]):
#         self.result = result
#         self.locus = result.locus
#         self.db = result.locus.db
#         self.name = name
#         self.alignments = alignments
#         self.reads = set(paf.read for paf in alignments)
#         self.length = alignments[0].target_length
#         # self.alignment_ranges = target_ranges_covered_by_alignments(alignments)[name]
#         self.alignment_range = alignment_range
#         # self.overlapping_bases = sum(end - start for start, end in self.alignment_ranges)
#         self.overlapping_bases = self.alignment_range[1] - self.alignment_range[0]
#         self.coverage = self.overlapping_bases / self.length * 100
#         self.identity = weighted_identity(alignments, 0, self.length)
#         self.genes = []
#         self.locus_region = None
#
#     def __repr__(self):
#         return f'{self.name} {self.coverage:.2f}% {self.identity:.2f}%'


# def get_gene_colour(gene_result: AssemblyGeneResult):
#     if gene_result.alignment.percent_identity >= 99.5:
#         return "#7F407F"
#     elif gene_result.alignment.percent_identity >= 98:
#         return "#995C99"
#     elif gene_result.alignment.percent_identity >= 95:
#         return "#B27DB2"
#     elif gene_result.alignment.percent_identity >= 90:
#         return "#CCA3CC"
#     elif gene_result.alignment.percent_identity >= 80:
#         return "#D9C3D9"
#     else:
#         return "#E5E5E5"
#
#
# def get_gene_shape(gene_result: AssemblyGeneResult):
#     if gene_result.truncated_protein:
#         return "JAGGY"
#     elif gene_result.edge_of_contig:
#         return "BOX"
#     else:
#         return "BIGARROW"

def parse_results(json, db: Database, regex: re.Pattern | None, samples: list[str] | None, loci: list[str] | None
                  ) -> Generator[AssemblyResult, None, None]:
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
                yield AssemblyResult.from_dict(result_dict, db)
            except ValueError:
                raise f"Invalid JSON: {json_line}"
