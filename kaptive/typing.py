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

from itertools import chain
from warnings import catch_warnings
from json import dumps
from typing import TextIO
from io import TextIOBase
from os import PathLike, path

from Bio.Seq import Seq
from Bio.Align import PairwiseAligner
import matplotlib
matplotlib.use('Agg')  # Prevents the need for a display when plotting
from dna_features_viewer import GraphicFeature, GraphicRecord

from kaptive.database import Database, Locus, Gene
from kaptive.log import warning

# Constants -----------------------------------------------------------------------------------------------------------
_PROTEIN_ALIGNER = PairwiseAligner(scoring='blastp', mode='local')


# Classes -------------------------------------------------------------------------------------------------------------
class TypingResultError(Exception):
    pass


class TypingResult:
    """
    This is a class to store the results of a typing analysis for a single sample. It is designed to be flexible
    enough to store results from both read and assembly typing, and can be easily reconstructed from JSON to use
    with the `convert` utility. It should not store any information that is not directly related to the typing
    such as contig sequences or read alignments.
    """

    def __init__(
            self, sample_name: str | None, db: Database | None, best_match: Locus | None = None,
            pieces: list[LocusPiece] | None = None, expected_genes_inside_locus: list[GeneResult] | None = None,
            expected_genes_outside_locus: list[GeneResult] | None = None, missing_genes: list[str] | None = None,
            unexpected_genes_inside_locus: list[GeneResult] | None = None,
            unexpected_genes_outside_locus: list[GeneResult] | None = None,
            extra_genes: list[GeneResult] | None = None):
        self.sample_name = sample_name or ""
        self.db = db
        self.best_match = best_match
        self.pieces = pieces or []  # Pieces of locus reconstructed from alignments
        self.expected_genes_inside_locus = expected_genes_inside_locus or []  # Genes from best_match
        self.expected_genes_outside_locus = expected_genes_outside_locus or []  # Genes from best_match
        self.missing_genes = missing_genes or []  # Genes from best_match that were not found
        self.unexpected_genes_inside_locus = unexpected_genes_inside_locus or []  # Genes from other loci
        self.unexpected_genes_outside_locus = unexpected_genes_outside_locus or []  # Genes from other loci
        self.extra_genes = extra_genes or []  # in db.extra_genes, ALWAYS outside locus (gene_result.piece == None)
        # Properties to cache the values, these are protected to prevent accidental modification
        self._percent_identity = None
        self._percent_coverage = None
        self._phenotype = None
        self._problems = None
        self._confidence = None

    def __repr__(self):
        return f"{self.sample_name} {self.best_match.name}"

    def __len__(self):
        return sum(len(i) for i in self.pieces) if self.pieces else 0

    def __iter__(self):
        return chain(
            self.expected_genes_inside_locus, self.unexpected_genes_inside_locus,
            self.expected_genes_outside_locus, self.unexpected_genes_outside_locus, self.extra_genes)

    def add_gene_result(self, gene_result: GeneResult):
        if gene_result.piece:  # If gene_result.piece is not None, the gene is inside the locus
            gene_result.piece.add_gene_result(gene_result)
            gene_type = f"{gene_result.gene_type}{'_inside_locus' if gene_result.gene_type.startswith(('expected', 'unexpected')) else ''}"
        else:  # If gene_result.piece is None, the gene is outside the locus
            gene_type = f"{gene_result.gene_type}{'_outside_locus' if gene_result.gene_type.startswith(('expected', 'unexpected')) else ''}"
        getattr(self, gene_type).append(gene_result)  # Add gene result to the appropriate list

    @property
    def percent_identity(self) -> float:
        if self._percent_identity is None:
            self._percent_identity = (sum(i.percent_identity for i in x) / len(x)) if (
                x := self.expected_genes_inside_locus) else 0
        return self._percent_identity

    @property
    def percent_coverage(self) -> float:
        if self._percent_coverage is None:
            self._percent_coverage = sum(len(i) for i in x) / sum(len(i) for i in self.best_match.genes.values()) * 100 \
                if (x := self.expected_genes_inside_locus) else 0
        return self._percent_coverage

    @property
    def phenotype(self) -> str:
        if self._phenotype is None:
            gene_phenotypes = set()  # Init set to store gene phenotypes to be used as a key in the phenotypes dict
            for gene in self:
                if gene.gene_type in {'expected_genes', 'extra_genes'}:  # The reported phenotype only considers these
                    gene_phenotypes.add((gene.gene.name, gene.phenotype))  # or extra genes
            # NOTE: The best_match.phenotypes MUST be sorted from largest to smallest gene set to make sure any sets with
            # extra genes are tested first.
            self._phenotype = next(
                (p for gs, p in self.best_match.phenotypes if len(gs) == len(gene_phenotypes.intersection(gs))),
                self.best_match.type_label)  # If no phenotype is found, return the type label
        return self._phenotype

    @property
    def problems(self) -> str:
        if self._problems is None:
            self._problems = f'?{x}' if (x := len(self.pieces)) != 1 else ''
            self._problems += '-' if self.missing_genes else ''
            self._problems += '+' if self.unexpected_genes_inside_locus else ''
            self._problems += '*' if any(
                i.percent_coverage >= 90 and i.below_threshold for i in self.expected_genes_inside_locus) else ''
            self._problems += '!' if any(i.phenotype == "truncated" for i in self) else ''
        return self._problems

    @property
    def confidence(self) -> str:
        return self._confidence if self._confidence is not None else "Not calculated"

    def get_confidence(self, allow_below_threshold: bool, max_other_genes: int, percent_expected_genes: float):
        p = len(set(i.gene.name for i in self.expected_genes_inside_locus)) / len(self.best_match.genes) * 100
        other_genes = len(set(i.gene.name for i in self.unexpected_genes_inside_locus if not i.phenotype == "truncated"))
        if not allow_below_threshold and "*" in self.problems:
            self._confidence = "Untypeable"
        else:
            if len(self.pieces) == 1 and not self.missing_genes and not other_genes:
                self._confidence = "Typeable"
            elif other_genes <= max_other_genes and p >= percent_expected_genes:
                self._confidence = "Typeable"
            else:
                self._confidence = "Untypeable"

    @classmethod
    def from_dict(cls, d: dict, db: Database) -> TypingResult:
        if not (best_match := db.loci.get(d['best_match'])):
            raise TypingResultError(f"Best match {d['best_match']} not found in database")
        self = TypingResult(sample_name=d['sample_name'], db=db, best_match=best_match, missing_genes=d['missing_genes'])
        # Set the cached properties
        self._percent_identity = float(d['percent_identity'])
        self._percent_coverage = float(d['percent_coverage'])
        self._phenotype = d['phenotype']
        self._problems = d['problems']
        self._confidence = d['confidence']
        # Add the pieces and create the gene results
        self.pieces = [LocusPiece.from_dict(i, result=self) for i in d['pieces']]
        pieces = {i.__repr__(): i for i in self.pieces}
        gene_results = {}  # This was previously a dict comp, but we need to check the gene is in the database, see #31
        for r in chain(d['expected_genes_inside_locus'], d['unexpected_genes_inside_locus'],
                       d['expected_genes_outside_locus'], d['unexpected_genes_outside_locus'],
                       d['extra_genes']):
            if not (gene := db.genes.get(r['gene'])) and not (gene := db.extra_genes.get(r['gene'])):
                raise TypingResultError(f"Gene {r['gene']} not found in database")
            x = GeneResult.from_dict(r, result=self, piece=pieces.get(r['piece']), gene=gene)
            gene_results[x.__repr__()] = x

        for gene_result in gene_results.values():
            self.add_gene_result(gene_result)
        return self

    def format(self, format_spec) -> str | GraphicRecord | dict:
        if format_spec == 'tsv':
            return '\t'.join(
                [
                    self.sample_name, self.best_match.name, self.phenotype, self.confidence, self.problems,
                    f"{self.percent_identity:.2f}%", f"{self.percent_coverage:.2f}%",
                    f"{self.__len__() - len(self.best_match)} bp" if len(self.pieces) == 1 else 'n/a',
                    f"{(x := len({i.gene.gene_name for i in self.expected_genes_inside_locus}))} / {(y := len(self.best_match.genes))} ({100 * x / y:.2f}%)",
                    ';'.join(str(i) for i in x) if (x := self.expected_genes_inside_locus) else '',
                    ';'.join(self.missing_genes), f"{len(x := self.unexpected_genes_inside_locus)}",
                    ';'.join(str(i) for i in x) if x else '',
                    f"{len(x := self.expected_genes_outside_locus)} / {(y := len(self.best_match.genes))} ({100 * len(x) / y:.2f}%)",
                    ';'.join(str(i) for i in x) if x else '',
                    f"{len(x := self.unexpected_genes_outside_locus)}",
                    ';'.join(str(i) for i in x) if x else '',
                    ';'.join(str(i) for i in filter(lambda z: z.phenotype == "truncated", self)),
                    ';'.join([str(i) for i in self.extra_genes])
                ]
            ) + "\n"
        if format_spec == 'fna':  # Return the nucleotide sequence of the locus
            return "".join([i.format(format_spec) for i in self.pieces])
        if format_spec in {'faa', 'ffn'}:  # Return the protein or nucleotide sequence of gene results
            return "".join([i.format(format_spec) for i in self])
        if format_spec in {'png', 'svg'}:
            features, start = [], 0
            for piece in self.pieces if self.pieces[0].strand == "+" else reversed(self.pieces):
                features.extend(piece.format(format_spec, start))
                start += len(piece)
            return GraphicRecord(sequence_length=self.__len__(), first_index=0, features=features,
                                 feature_level_height=1.5, sequence=[p.sequence for p in self.pieces])
        if format_spec == 'json':
            return dumps(
                {
                    'sample_name': self.sample_name, 'best_match': self.best_match.name, 'confidence': self.confidence,
                    'phenotype': self.phenotype, 'problems': self.problems, 'percent_identity': str(self.percent_identity),
                    'percent_coverage': str(self.percent_coverage), 'missing_genes': self.missing_genes
                } | {
                    attr: [i.format(format_spec) for i in getattr(self, attr)] for attr in {
                        'pieces', 'expected_genes_inside_locus', 'unexpected_genes_inside_locus',
                        'expected_genes_outside_locus', 'unexpected_genes_outside_locus', 'extra_genes'
                    }
                }) + "\n"
        raise ValueError(f"Unknown format specifier {format_spec}")

    def write(self,
              tsv: TextIO | None = None,
              json: TextIO | None = None,
              fna: str | PathLike | TextIO | None = None,
              ffn: str | PathLike | TextIO | None = None,
              faa: str | PathLike | TextIO | None = None,
              plot: str | PathLike | None = None,
              plot_fmt: str = 'png'):
        """Write the typing result to files or file handles."""
        [f.write(self.format(fmt)) for f, fmt in [(tsv, 'tsv'), (json, 'json')] if isinstance(f, TextIOBase)]
        for f, fmt in [(fna, 'fna'), (ffn, 'ffn'), (faa, 'faa')]:
            if f:
                if isinstance(f, TextIOBase):
                    f.write(self.format(fmt))
                elif isinstance(f, PathLike) or isinstance(f, str):
                    with open(path.join(f,  f'{self.sample_name}_kaptive_results.{fmt}'), 'wt') as handle:
                        handle.write(self.format(fmt))
        if plot:
            ax = self.format(plot_fmt).plot(figure_width=18)[0]  # type: 'matplotlib.axes.Axes'
            ax.set_title(f"{self.sample_name} {self.best_match} ({self.phenotype}) - {self.confidence}")
            ax.figure.savefig(path.join(plot, f'{self.sample_name}_kaptive_results.{plot_fmt}'), bbox_inches='tight')
            ax.figure.clear()  # TODO: Check if this is necessary


class LocusPieceError(Exception):
    pass


class LocusPiece:
    def __init__(self, id: str | None = None, result: TypingResult | None = None, start: int | None = 0,
                 end: int | None = 0, strand: str | None = None, sequence: Seq | None = None,
                 expected_genes: list[GeneResult] | None = None, unexpected_genes: list[GeneResult] | None = None,
                 extra_genes: list[GeneResult] | None = None):
        self.id = id or ''  # TODO: rename as seq_id for clarity, actual id is self.__repr__()
        self.result = result
        self.start = start
        self.end = end
        self.strand = strand or "unknown"
        self.sequence = sequence or Seq("")
        self.expected_genes = expected_genes or []  # Genes from best_match
        self.unexpected_genes = unexpected_genes or []  # Genes that were found from other loci
        self.extra_genes = extra_genes or []  # Genes that were found outside the locus

    def __len__(self):
        return self.end - self.start

    def __iter__(self):
        return chain(self.expected_genes, self.unexpected_genes, self.extra_genes)

    def __str__(self):
        return self.id

    def __repr__(self):
        return f"{self.id}:{self.start}-{self.end}{self.strand}"

    @classmethod
    def from_dict(cls, d: dict, **kwargs) -> LocusPiece:
        return cls(id=d['id'], start=int(d['start']), end=int(d['end']), strand=d['strand'],
                   sequence=Seq(d['sequence']), **kwargs)

    def format(self, format_spec, relative_start: int = 0) -> str | dict | list[GraphicFeature]:
        if format_spec == 'fna':
            return f">{self.result.sample_name}|{self.id}:{self.start}-{self.end}{self.strand}\n{self.sequence}\n"
        if format_spec == 'json':
            return {'id': self.id, 'start': str(self.start), 'end': str(self.end), 'strand': self.strand,
                    'sequence': str(self.sequence)}
        if format_spec in {'png', 'svg'}:
            return [GraphicFeature(start=relative_start, end=relative_start + len(self), strand=1, thickness=30,
                                   color='#762a83', label=str(self), linewidth=0)] + [
                gene.format(format_spec, gene.start - self.start + relative_start) for gene in self]
        raise ValueError(f"Unknown format specifier {format_spec}")

    def add_gene_result(self, gene_result: GeneResult):
        if gene_result.start < self.start:  # Update start and end if necessary
            self.start = gene_result.start
        if gene_result.end > self.end:
            self.end = gene_result.end
        getattr(self, gene_result.gene_type).append(gene_result)


class GeneResultError(Exception):
    pass


class GeneResult:
    """
    Class to store alignment results for a single gene in a locus for either a ReadResult or a AssemblyResult.
    """

    def __init__(self, id: str, gene: Gene, result: TypingResult | None = None,
                 piece: LocusPiece | None = None, start: int | None = 0, end: int | None = 0, strand: str | None = None,
                 dna_seq: Seq | None = Seq(""), protein_seq: Seq | None = Seq(""), below_threshold: bool | None = False,
                 phenotype: str | None = "present", gene_type: str | None = None, partial: bool | None = False,
                 percent_identity: float | None = 0, percent_coverage: float | None = 0):
        self.id = id or ''  # TODO: rename as seq_id for clarity, actual id is self.__repr__()
        self.gene = gene
        self.result = result  # TODO: replace with sample_name (only TypingResult attribute used, needed for fa headers)
        self.start = start
        self.end = end
        self.strand = strand
        self.partial = partial
        self.piece = piece  # inside locus if not None
        self.dna_seq = dna_seq
        self.protein_seq = protein_seq
        self.below_threshold = below_threshold
        self.phenotype = phenotype
        self.gene_type = gene_type or ""
        self.percent_identity = percent_identity
        self.percent_coverage = percent_coverage

    def __repr__(self):
        return f"{self.gene.name} {self.id}:{self.start}-{self.end}{self.strand}"

    def __len__(self):
        return self.end - self.start

    def __str__(self) -> str:
        s = f'{self.gene.name},{self.percent_identity:.2f}%,{self.percent_coverage:.2f}%'
        s += ",partial" if self.partial else ""
        s += ',truncated' if self.phenotype == "truncated" else ""
        s += ",below_id_threshold" if self.below_threshold else ""
        return s

    @classmethod
    def from_dict(cls, d: dict, **kwargs) -> GeneResult:
        return cls(
            id=d['id'], start=int(d['start']), end=int(d['end']), strand=d['strand'], dna_seq=Seq(d['dna_seq']),
            protein_seq=Seq(d['protein_seq']), below_threshold=True if d['below_threshold'] == 'True' else False,
            phenotype=d['phenotype'], gene_type=d['gene_type'], partial=True if d['partial'] == 'True' else False,
            percent_identity=float(d['percent_identity']), percent_coverage=float(d['percent_coverage']), **kwargs
        )

    def format(self, format_spec, relative_start: int = 0) -> str | dict | GraphicFeature:
        if format_spec == 'ffn':
            if len(self.dna_seq) == 0:
                warning(f'No DNA sequence for {self}')
                return ""
            return (f'>{self.gene.name} {self.result.sample_name}|{self.id}:{self.start}-{self.end}{self.strand}\n'
                    f'{self.dna_seq}\n')
        if format_spec == 'faa':
            if len(self.protein_seq) == 0:
                warning(f'No protein sequence for {self.__repr__()}')
                return ""
            return (f'>{self.gene.name} {self.result.sample_name}|{self.id}:{self.start}-{self.end}{self.strand}\n'
                    f'{self.protein_seq}\n')
        if format_spec == 'json':
            return {
                'id': self.id, 'start': str(self.start), 'end': str(self.end), 'strand': self.strand,
                'dna_seq': str(self.dna_seq), 'protein_seq': str(self.protein_seq), 'partial': str(self.partial),
                'below_threshold': str(self.below_threshold), 'phenotype': self.phenotype, 'gene_type': self.gene_type,
                'percent_identity': str(self.percent_identity), 'percent_coverage': str(self.percent_coverage),
                'gene': self.gene.name, 'piece': self.piece.__repr__() if self.piece else '',
            }
        if format_spec in {'png', 'svg'}:
            strand = self.gene.strand if self.strand == self.gene.strand else self.strand
            return GraphicFeature(
                start=relative_start, end=relative_start + len(self), legend_text=self.gene_type, label=str(self),
                strand=0 if self.phenotype == "truncated" or self.partial else 1 if strand == "+" else -1,
                thickness=40, linewidth=3,
                color=('#762a83' if self.gene_type == 'expected_genes' else "orange", self.percent_identity / 100),
                linecolor='red' if self.below_threshold else "yellow" if self.phenotype == "truncated" else 'black'
            )
        raise ValueError(f"Unknown format specifier {format_spec}")

    def compare_translation(self, truncation_tolerance: float = 95, **kwargs):
        """
        Extracts the translation from the DNA sequence of the gene result.
        Will also extract the translation from the gene if it is not already stored.
        """
        self.gene.extract_translation(**kwargs)  # Extract the translation from the gene if it is not already stored
        if len(self.dna_seq) == 0:  # If the DNA sequence is empty, raise an error
            raise GeneResultError(f'No DNA sequence for {self.__repr__()}')
        with catch_warnings(record=True) as w:  # Catch Biopython warnings
            protein_seqs = [self.dna_seq[i:].translate(**kwargs) for i in range(3)]  # Translate in all 3 frames
        frame, self.protein_seq = max(enumerate(protein_seqs), key=lambda x: len(x[1]))  # Get the longest translation
        self.start += frame  # Update the start position to the frame with the longest translation
        if len(self.protein_seq) <= 1:  # If the protein sequence is still empty, raise a warning
            warning(f'No protein sequence for {self.__repr__()}')
        elif len(self.gene.protein_seq) > 1:  # If both sequences are not empty
            if alignments := _PROTEIN_ALIGNER.align(self.gene.protein_seq, self.protein_seq):  # Align the sequences
                alignment = max(alignments, key=lambda x: x.score)  # Get the best alignment
                self.percent_identity = alignment.counts().identities / alignment.length * 100
                self.percent_coverage = (len(self.protein_seq) / len(self.gene.protein_seq)) * 100
                if (self.percent_coverage < truncation_tolerance and  # If the coverage is less than the tolerance
                        # not partial, and not unexpected gene outside locus
                        not self.partial and not (self.gene_type == "unexpected_genes" and not self.piece)):
                    self.phenotype = "truncated"  # Set the phenotype to truncated
            else:
                warning(f'Error aligning {self.__repr__()}')

