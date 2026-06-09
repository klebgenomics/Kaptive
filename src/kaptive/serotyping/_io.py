from abc import ABC, abstractmethod
from enum import IntEnum, IntFlag, auto
from dataclasses import dataclass
from typing import Self

import numpy as np

from ._pipeline import SerotypingResult, GeneState


# Enums ----------------------------------------------------------------------------------------------------------------
class SerotypingProblem(IntFlag):
    NONE = 0
    FRAGMENTED = auto()  # Locus is broken up over pieces across the same/multiple contig(s)
    UNEXPECTED_GENES = auto()  # Unexpected genes from other loci inside the locus boundary
    MISSING_GENES = auto()  # Expected genes from the best match locus missing inside the locus boundary
    NOVEL_GENES = auto()  # Genes inside the locus boundary that fall below the pairwise amino acid identity threshold
    TRUNCATED_GENES = auto()  # Genes inside the locus boundary that do not form complete amino acid sequences

    @classmethod
    def evaluate(cls, result: SerotypingResult) -> int:
        problems = cls.NONE
        if len(result.locus_pieces) > 1:
            problems |= cls.FRAGMENTED
        if np.any(result.genes.is_inside & result.genes.is_expected):
            problems |= cls.UNEXPECTED_GENES


        # if np.any(genes.is_expected & (genes.states == GeneState.TRUNCATED)):
        #     problems |= self.TRUNCATED_GENES
        #
        # if np.any(active_extra):
        #     problems |= self.UNEXPECTED_GENES
        #
        # # Missing genes: Are any expected clusters absent from the active expected hits?
        # locus_slice = self._db.locus_gene_slices[best_locus_idx]
        # expected_clusters = set(self._db.gene_cluster_ids[locus_slice])
        # found_expected = set(self._db.gene_cluster_ids[genes.gene_indices[active_expected]])
        # if not expected_clusters.issubset(found_expected):
        #     problems |= self.MISSING_GENES
        #
        # # Novel genes: Anything physically inside the locus that falls below the pident threshold
        # novel_mask = genes.is_inside & (genes.states == GeneState.NOVEL)
        # if np.any(novel_mask):
        #     problems |= self.NOVEL_GENES
        return problems


# Classes --------------------------------------------------------------------------------------------------------------
@dataclass(slots=True, frozen=True)
class ConfidenceEvaluator:
    """
    Evaluates the final SerotyperResult to determine if the match is trustworthy.
    """
    max_other_genes: int = 1
    min_completeness: float = 0.5  # Replaced prop_expected_genes for SoA clarity
    allow_below_threshold: bool = False
    min_zscore: float | None = None

    def __call__(self, result: SerotypingResult) -> bool:
        pass


@dataclass(slots=True, frozen=True)
class ReportRow(ABC):
    """Base class for representing an in silico serotyping report for a single sample in a TSV row.
    Allows for the auto-documentation of columns using the class docstrings under the "Attributes" section."""

    @classmethod
    def header(cls) -> str:
        """Returns the TSV header string dynamically generated from the dataclass fields."""
        return "\t".join(cls.__annotations__.keys()) + "\n"

    def __str__(self) -> str:
        """Formats the row data into a tab-separated string."""
        return "\t".join(self.__annotations__.values()) + "\n"

    @classmethod
    @abstractmethod
    def from_result(cls, result: SerotypingResult, evaluator: ConfidenceEvaluator) -> Self:
        ...


@dataclass(slots=True, frozen=True)
class KaptiveRow(ReportRow):
    """Represents an in silico serotyping report for a single sample in a TSV row. in the original Kaptive format.

    Attributes:
        Assembly (str): The name of the input assembly, taken from the assembly filename.
        Best_match_locus (str): The locus type which most closely matches the assembly.
        Best_match_type (str): The predicted serotype/phenotype of the assembly.
        Match_confidence (str): A categorical measure of locus call quality.
        Problems (str): Characters indicating issues with the locus match.
        Identity (str): Weighted percent identity of the best matching locus to the assembly.
        Coverage (str): Weighted percent coverage of the best matching locus in the assembly.
        Length_discrepancy (str): If the locus was found in a single piece,
            this is the difference between the locus length and the assembly length.
        Expected_genes_in_locus (str): A fraction indicating how many of the genes in the best matching locus were
            found in the locus part of the assembly.
        Expected_genes_in_locus_details (str): Gene names for the expected genes found in the locus part of the assembly.
        Missing_expected_genes (str): A string listing the gene names of expected genes that were not found.
        Other_genes_in_locus (str): The number of unexpected genes (genes from loci other than the best match) which
            were found in the locus part of the assembly.
        Other_genes_in_locus_details (str): Gene names for the other genes found in the locus part of the assembly.
        Expected_genes_outside_locus (str): A fraction indicating how many of the expected genes which were found in
            the assembly but not in the locus part of the assembly (usually zero).
        Expected_genes_outside_locus_details (str):  Gene names for the expected genes found outside the locus part of
            the assembly.
        Other_genes_outside_locus (str): The number of unexpected genes (genes from loci other than the best match)
            which were found outside the locus part of the assembly.
        Other_genes_outside_locus_details (str): Gene names for the other genes found outside the locus part of the
            assembly.
        Truncated_genes_details (str): Gene names for the truncated genes found in the assembly.
        Extra_genes_details (str): Gene names for the extra genes found in the assembly.

    Note:
        Numbers beside gene names indicate the percent identity and percent coverage of the gene in the assembly.

    Warning:
        You may sometimes see two copies of the same gene in the `Expected genes in locus, details` column.
        These represent (likely) parts of the same gene which have usually been split over contigs. In
        `kaptive v3.0.0` onwards, we adopted this behavior to allow users to see where locus splitting has occurred,
        and determine the total percent identity of a gene that has been split.
    """
    Assembly: str
    Best_match_locus: str
    Best_match_type: str
    Match_confidence: str
    Problems: str
    Identity: str
    Coverage: str
    Length_discrepancy: str
    Expected_genes_in_locus: str
    Expected_genes_in_locus_details: str
    Missing_expected_genes: str
    Other_genes_in_locus: str
    Other_genes_in_locus_details: str
    Expected_genes_outside_locus: str
    Expected_genes_outside_locus_details: str
    Other_genes_outside_locus: str
    Other_genes_outside_locus_details: str
    Truncated_genes_details: str
    Extra_genes_details: str

    @staticmethod
    def format_gene(gene):
        gene_string = [gene.name, f'{gene.percent_identity:.2f}%', f'{gene.percent_coverage:.2f}%']
        if gene.partial:
            gene_string.append(f'partial')
        # s += ',truncated' if gene.phenotype == "truncated" else ""
        # s += ",below_id_threshold" if gene.below_threshold else ""
        return ','.join(gene_string)

    @classmethod
    def from_result(cls, result: SerotypingResult, evaluator: ConfidenceEvaluator):
        ...


@dataclass(slots=True, frozen=True)
class Pha4geRow(ReportRow):
    """Represents an in silico serotyping report for a single sample in a TSV row. in the original Kaptive format.

    Attributes:
        sample (str): Lorem ipsum
        genotyping_method (str): Lorem ipsum
        genotyping_schema_taxon (str): Lorem ipsum
        genotyping_database_name (str): Lorem ipsum
        genotyping_database_version (str): Lorem ipsum
        genotyping_schema_name (str): Lorem ipsum
        genotyping_software_name (str): Lorem ipsum
        genotyping_software_version (str): Lorem ipsum
        genotype (str): Lorem ipsum
        genotype_confidence_value (str): Lorem ipsum
        genotype_predicted_phenotype (str): Lorem ipsum
    """
    sample: str
    genotyping_method: str
    genotyping_schema_taxon: str
    genotyping_database_name: str
    genotyping_database_version: str
    genotyping_schema_name: str
    genotyping_software_name: str
    genotyping_software_version: str
    genotype: str
    genotype_confidence_value: str
    genotype_predicted_phenotype: str

    @classmethod
    def from_result(cls, result: SerotypingResult, evaluator: ConfidenceEvaluator):
        ...
