r"""I/O formatting and TSV report generation for in silico serotyping results.

This module provides abstract base class [`ReportRow`][kaptive.serotyping.io.ReportRow] and concrete implementation
dataclasses [`KaptiveRow`][kaptive.serotyping.io.KaptiveRow] and [`Pha4geRow`][kaptive.serotyping.io.Pha4geRow]
for exporting [`SerotypingResult`][kaptive.serotyping.models.SerotypingResult] objects into tab-separated value (TSV)
report files adhering to original Kaptive or standard PHA4GE reporting formats.
"""

from abc import ABC, abstractmethod
from dataclasses import dataclass, fields
from typing import Self

import numpy as np

from kaptive.serotyping.models import GeneState, SerotypingProblem, SerotypingResult


@dataclass(slots=True, frozen=True)
class ReportRow(ABC):
    r"""Abstract base class for tabular *in silico* serotyping report rows.

    Provides a uniform interface and binary serialization methods (`__bytes__` and `header`) for converting
    [`SerotypingResult`][kaptive.serotyping.models.SerotypingResult] instances into tab-separated (TSV) outputs.
    Attributes documented in subclass docstrings correspond directly to TSV report column headers.
    """

    @classmethod
    def header(cls) -> bytes:
        r"""Generate the TSV header row as UTF-8 encoded bytes.

        Returns:
            bytes: Tab-separated column header line ending with a newline (`b"\n"`).
        """
        return ("\t".join(f.name for f in fields(cls)) + "\n").encode("utf-8")

    def __bytes__(self) -> bytes:
        r"""Serialize the report row fields into a tab-separated binary TSV row.

        Returns:
            bytes: Tab-separated field values ending with a newline (`b"\n"`).
        """
        return b"\t".join(getattr(self, f.name) for f in fields(self)) + b"\n"

    @classmethod
    @abstractmethod
    def from_result(cls, result: SerotypingResult) -> Self:
        r"""Construct a report row instance from a serotyping result.

        Args:
            result (SerotypingResult): The serotyping analysis result to format.
                See [`SerotypingResult`][kaptive.serotyping.models.SerotypingResult].

        Returns:
            Self: Instance of the concrete [`ReportRow`][kaptive.serotyping.io.ReportRow] subclass
                populated with result fields.
        """
        ...


@dataclass(slots=True, frozen=True)
class KaptiveRow(ReportRow):
    r"""Report row representation matching the classic Kaptive TSV output format.

    Encapsulates all summary statistics, locus match calls, problem flags, gene details, and coverage metrics for a
    single genome assembly in tab-separated binary format compatible with traditional Kaptive output parsers.

    Attributes:
        Kaptive_version (bytes): The version of Kaptive used to perform serotyping.
        Database_name (bytes): Name of the reference database used for serotyping.
        Database_version (bytes): Version of the reference database used.
        Assembly (bytes): Identifier/filename of the analyzed genome assembly.
        Best_match_locus (bytes): Best matching reference locus type identifier.
        Best_match_type (bytes): Predicted serotype/phenotype call for the genome.
        Match_confidence (bytes): Confidence classification (`b"Typeable"` or `b"Untypeable"`).
        Problems (bytes): Symbolic character flags representing
            [`SerotypingProblem`][kaptive.serotyping.models.SerotypingProblem]
            locus match issues (`?`, `+`, `-`, `*`, `!`).
        Identity (bytes): Mean percentage amino acid identity across intact expected locus genes.
        Coverage (bytes): Percentage coverage of the best matching reference locus by assembly contigs.
        Length_discrepancy (bytes): Difference in base pairs between assembly locus length and reference locus length
            (or `"n/a"`).
        Expected_genes_in_locus (bytes): Count and fraction of expected locus genes found inside locus boundary.
        Expected_genes_in_locus_details (bytes): Detailed identity and coverage specs for expected genes inside locus.
        Missing_expected_genes (bytes): Semicolon-separated names of expected genes not found.
        Other_genes_in_locus (bytes): Count of unexpected genes from other loci found inside locus boundary.
        Other_genes_in_locus_details (bytes): Detailed specs for unexpected genes inside locus.
        Expected_genes_outside_locus (bytes): Count and fraction of expected locus genes found outside locus boundary.
        Expected_genes_outside_locus_details (bytes): Detailed specs for expected genes found outside locus.
        Other_genes_outside_locus (bytes): Count of unexpected genes found outside locus boundary.
        Other_genes_outside_locus_details (bytes): Detailed specs for unexpected genes found outside locus.
        Truncated_genes_details (bytes): Detailed specs for truncated or partial genes.
        Extra_genes_details (bytes): Detailed specs for allowed extra database genes.

    Note:
        Numbers beside gene names indicate percentage identity and percentage coverage of the gene in the genome.

    Warning:
        You may sometimes see two copies of the same gene in the `Expected_genes_in_locus_details` column.
        These represent parts of the same gene split over contig boundaries.
    """

    Kaptive_version: bytes
    Database_name: bytes
    Database_version: bytes
    Assembly: bytes
    Best_match_locus: bytes
    Best_match_type: bytes
    Match_confidence: bytes
    Problems: bytes
    Identity: bytes
    Coverage: bytes
    Length_discrepancy: bytes
    Expected_genes_in_locus: bytes
    Expected_genes_in_locus_details: bytes
    Missing_expected_genes: bytes
    Other_genes_in_locus: bytes
    Other_genes_in_locus_details: bytes
    Expected_genes_outside_locus: bytes
    Expected_genes_outside_locus_details: bytes
    Other_genes_outside_locus: bytes
    Other_genes_outside_locus_details: bytes
    Truncated_genes_details: bytes
    Extra_genes_details: bytes

    @classmethod
    def header(cls) -> bytes:
        r"""Generate backwards-compatible column header bytes for classic Kaptive reports.

        Replaces internal field name underscores with spaces and `_details` with `, details` to maintain exact
        compatibility with legacy Kaptive TSV headers.

        Returns:
            bytes: Tab-separated legacy header line ending with a newline (`b"\n"`).
        """
        headers = [
            f.name.encode("utf-8")
            .replace(b"_details", b", details")
            .replace(b"_", b" ")
            for f in fields(cls)
        ]
        return b"\t".join(headers) + b"\n"

    @classmethod
    def from_result(cls, result: SerotypingResult) -> "KaptiveRow":
        r"""Construct a classic [`KaptiveRow`][kaptive.serotyping.io.KaptiveRow] from a serotyping result.

        Calculates gene counts, percentage coverages, identity metrics, and problem symbol codes, formatting all fields
        into UTF-8 encoded bytes for backwards-compatible TSV output.

        Args:
            result (SerotypingResult): The serotyping call result.
                See [`SerotypingResult`][kaptive.serotyping.models.SerotypingResult].

        Returns:
            KaptiveRow: Formatted report row object.
        """
        hits = result.gene_hits
        states = result.gene_states

        # High-performance boolean masks
        in_loc = hits.is_inside
        out_loc = ~hits.is_inside
        exp = hits.is_expected
        extra = hits.is_extra
        unexp = ~exp & ~extra

        def _format_genes(mask: np.ndarray) -> bytes:
            r"""Helper to rapidly construct the details string for a specific masked subset of genes.

            Args:
                mask (np.ndarray): Boolean mask selecting the subset of gene hits to format.

            Returns:
                bytes: Semicolon-separated details string for the selected genes.
            """
            indices = np.where(mask)[0]
            if indices.size == 0:
                return b""

            details = []
            for i in indices:
                gene_name = result.gene_seqs.ids[i].encode("utf-8")
                parts = [
                    gene_name,
                    b"%.2f%%" % result.protein_identities[i],
                    b"%.2f%%" % result.gene_hits.coverages[i],
                ]

                if states[i] == GeneState.PARTIAL.value:
                    parts.append(b"partial")
                elif states[i] == GeneState.TRUNCATED.value:
                    parts.append(b"truncated")
                elif states[i] == GeneState.NOVEL.value:
                    parts.append(b"below_id_threshold")

                details.append(b",".join(parts))
            return b";".join(details)

        # Expected Inside
        mask_exp_in = in_loc & exp
        n_exp_in = len(np.unique(result.gene_hits.gene_indices[mask_exp_in]))

        # Expected Outside
        mask_exp_out = out_loc & exp
        n_exp_out = len(np.unique(result.gene_hits.gene_indices[mask_exp_out]))

        expected_total = n_exp_in + n_exp_out + len(result.missing_expected_genes)

        in_comp = (n_exp_in / expected_total * 100.0) if expected_total > 0 else 0.0
        exp_in_str = (
            b"%d / %d (%.2f%%)" % (n_exp_in, expected_total, in_comp)
            if expected_total
            else b"0 / 0 (0.00%)"
        )

        out_comp = (n_exp_out / expected_total * 100.0) if expected_total > 0 else 0.0
        exp_out_str = (
            b"%d / %d (%.2f%%)" % (n_exp_out, expected_total, out_comp)
            if expected_total
            else b"0 / 0 (0.00%)"
        )

        # Other counts
        n_unexp_in = len(np.unique(result.gene_hits.gene_indices[in_loc & unexp]))
        n_unexp_out = len(np.unique(result.gene_hits.gene_indices[out_loc & unexp]))

        return cls(
            Kaptive_version=result.kaptive_version.encode(),
            Database_name=result.database_name.encode(),
            Database_version=result.database_version.encode(),
            Assembly=result.genome.encode(),
            Best_match_locus=result.best_locus_name.encode(),
            Best_match_type=result.phenotype.encode(),
            Match_confidence=b"Typeable" if result.typeable else b"Untypeable",
            Problems=result.problems.to_symbols(),
            Identity=b"%.2f%%" % result.percent_identity,
            Coverage=b"%.2f%%" % result.percent_coverage,
            Length_discrepancy=b"n/a"
            if (
                result.length_discrepancy is None or np.isnan(result.length_discrepancy)
            )
            else b"%d" % int(result.length_discrepancy),
            Expected_genes_in_locus=exp_in_str,
            Expected_genes_in_locus_details=_format_genes(mask_exp_in),
            Missing_expected_genes=b";".join(g.encode("utf-8") for g in result.missing_expected_genes),
            Other_genes_in_locus=b"%d" % n_unexp_in,
            Other_genes_in_locus_details=_format_genes(in_loc & unexp),
            Expected_genes_outside_locus=exp_out_str,
            Expected_genes_outside_locus_details=_format_genes(mask_exp_out),
            Other_genes_outside_locus=b"%d" % n_unexp_out,
            Other_genes_outside_locus_details=_format_genes(out_loc & unexp),
            Truncated_genes_details=_format_genes(
                (states == GeneState.TRUNCATED.value)
                | (states == GeneState.PARTIAL.value)
            ),
            Extra_genes_details=_format_genes(extra),
        )


@dataclass(slots=True, frozen=True, kw_only=True)
class Pha4geRow(ReportRow):
    r"""Report row representation adhering to Public Health Alliance for Genomic Epidemiology (PHA4GE) standards.

    Encapsulates sample metadata, taxonomy, software versioning, genotype calls, and confidence values in tab-separated
    binary format standardized for public health surveillance data exchange. 
    
    For more information on the rationale and specifics of the PHA4GE genotyping specification, 
    please see: https://github.com/pha4ge/genotyping-specification

    Attributes:
        sample (bytes): Sample identifier taken from genome assembly filename.
        genotyping_method (bytes): Genotyping methodology string (default `b"In silico serotyping"`).
        genotyping_schema_taxon (bytes): NCBITaxon formatted organism species string and taxon ID.
        genotyping_database_name (bytes): Name of reference database used for serotyping.
        genotyping_database_version (bytes): Version of reference database used.
        genotyping_schema_name (bytes): Schema name (default `b"Kaptive"`).
        genotyping_software_name (bytes): Software name (default `b"Kaptive"`).
        genotyping_software_version (bytes): Kaptive software version used for analysis.
        genotype (bytes): Best matching locus type identifier call.
        genotype_predicted_phenotype (bytes): Predicted surface antigen phenotype/serotype string.
        genotype_confidence_value (bytes): Confidence assessment (`b"Typeable"` or `b"Untypeable"`).
        genotyping_details (bytes): Human-readable descriptions of any locus match problems detected.
        genotyping_method_url (bytes): Repository URL for methodology documentation.
    """

    sample: bytes
    genotyping_method: bytes = b"In silico serotyping"
    genotyping_schema_taxon: bytes
    genotyping_database_name: bytes
    genotyping_database_version: bytes
    genotyping_schema_name: bytes = b"Kaptive"
    genotyping_software_name: bytes = b"Kaptive"
    genotyping_software_version: bytes
    genotype: bytes
    genotype_predicted_phenotype: bytes
    genotype_confidence_value: bytes
    genotyping_details: bytes
    genotyping_method_url: bytes = b"https://github.com/klebgenomics/Kaptive"

    @classmethod
    def from_result(cls, result: SerotypingResult) -> "Pha4geRow":
        r"""Construct a standardized [`Pha4geRow`][kaptive.serotyping.io.Pha4geRow] from a serotyping result.

        Transforms numeric taxon IDs and problem flags into human-readable PHA4GE-compliant strings and binary bytes.

        Args:
            result (SerotypingResult): The serotyping call result.
                See [`SerotypingResult`][kaptive.serotyping.models.SerotypingResult].

        Returns:
            Pha4geRow: Formatted PHA4GE report row object.
        """
        # Transform problem enum into human-readable string
        if result.problems:
            detail_parts = []
            if SerotypingProblem.TRUNCATED_GENES in result.problems:
                detail_parts.append(b"truncated gene/s in locus")
            if SerotypingProblem.NOVEL_GENES in result.problems:
                detail_parts.append(b"low identity gene/s")
            if SerotypingProblem.FRAGMENTED in result.problems:
                detail_parts.append(b"match broken into %d pieces" % len(result.locus_pieces))
            if SerotypingProblem.MISSING_GENES in result.problems:
                detail_parts.append(b"missing expected gene/s")
            if SerotypingProblem.UNEXPECTED_GENES in result.problems:
                detail_parts.append(b"unexpected gene/s in locus")
            details = b"Best locus match: %b. Problems: %b" % (result.best_locus_name.encode(), b", ".join(detail_parts))
        else:
            details = b"Best locus match: %b." % result.best_locus_name.encode()

        return cls(
            sample=result.genome.encode(),
            genotyping_schema_taxon=b"%s [NCBITaxon:%d]"
            % (result.database_organism.encode(), result.database_taxon),
            genotyping_database_name=result.database_name.encode(),
            genotyping_database_version=result.database_version.encode(),
            genotyping_software_version=result.kaptive_version.encode(),
            genotype=result.best_locus_name.encode(),
            genotype_confidence_value=b"Typeable" if result.typeable else b"Untypeable",
            genotype_predicted_phenotype=result.phenotype.encode(),
            genotyping_details=details,
        )
