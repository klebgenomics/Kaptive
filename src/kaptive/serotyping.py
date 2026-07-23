import rammappy
from abc import ABC, abstractmethod
from collections.abc import Iterable
from dataclasses import dataclass, fields
from enum import IntEnum, IntFlag, auto
from pathlib import Path
from typing import Any, ClassVar, Self
from itertools import chain
import numba as nb

import numpy as np
import numpy.typing as npt
from numba import njit
from rammappy.align import Aligner, Index

from kaptive._version import __version__
from kaptive.core.alignment import Alignments
from kaptive.core.collections import BatchedContainer
from kaptive.core.genome import GenomeAssembly
from kaptive.core.interval import Intervals
from kaptive.core.pairwise import PairwiseAligner
from kaptive.core.seq import Sequences
from kaptive.db import Database
from kaptive.compare import LocusData


# Enums ----------------------------------------------------------------------------------------------------------------
class GeneState(IntEnum):
    """Mutually exclusive states for locus genes found in a genome assembly

    Attributes:
        NORMAL (int): The gene was found intact as expected
        PARTIAL (int): The gene was broken up over a contig edge
        TRUNCATED (int): The gene does not form a complete amino acid sequence
        NOVEL (int): The gene translation diverges significantly from the closest reference
    """

    NORMAL = 0
    PARTIAL = 1
    TRUNCATED = 2
    NOVEL = 3


class SerotypingProblem(IntFlag):
    """Symbolic problems with the serotype call - mostly for reporting

    Attributes:
        NONE (int): No problems with the serotype call.
        FRAGMENTED (int): Locus is broken up over pieces across the same/multiple contig(s). Symbol: `?`
        UNEXPECTED_GENES (int): Unexpected genes from other loci inside the locus boundary. Symbol: `+`
        MISSING_GENES (int): Expected genes from the best match locus missing inside the locus boundary. Symbol: `-`
        NOVEL_GENES (int): Genes inside the locus boundary that fall below the pairwise amino acid identity threshold. Symbol: `*`
        TRUNCATED_GENES (int): Genes inside the locus boundary that do not form complete amino acid sequences. Symbol: `!`
        SYMBOLS (tuple[bytes, ...]): Precomputed lookup table mapping flag combinations to their respective symbols.
    """

    NONE = 0
    FRAGMENTED = auto()
    UNEXPECTED_GENES = auto()
    MISSING_GENES = auto()
    NOVEL_GENES = auto()
    TRUNCATED_GENES = auto()

    SYMBOLS: ClassVar[tuple[bytes, ...]]

    def to_symbols(self) -> bytes:
        """Renders the problem code for the Kaptive TSV output"""
        return self.SYMBOLS[self.value]


_serotyping_flags = (
    (SerotypingProblem.FRAGMENTED.value, b"?"),
    (SerotypingProblem.UNEXPECTED_GENES.value, b"+"),
    (SerotypingProblem.MISSING_GENES.value, b"-"),
    (SerotypingProblem.NOVEL_GENES.value, b"*"),
    (SerotypingProblem.TRUNCATED_GENES.value, b"!"),
)
SerotypingProblem.SYMBOLS = tuple(
    b"".join(sym for flag, sym in _serotyping_flags if i & flag)
    for i in range(1 << max(SerotypingProblem).value.bit_length())
)


# Classes --------------------------------------------------------------------------------------------------------------
@dataclass(slots=True, frozen=True)
class GeneHits(BatchedContainer[Any, 'GeneHits']):
    """A high-performance SoA container for the final classified gene alignments.

    Encapsulating these arrays prevents `SerotypingResult` from being cluttered
    with 11 parallel arrays and enables synchronized vectorized filtering and
    dynamic property calculations.
    """

    gene_indices: npt.NDArray[np.int32]  # dtype=np.int32 (Global DB index)
    q_starts: npt.NDArray[np.int32]  # dtype=np.int32
    q_ends: npt.NDArray[np.int32]  # dtype=np.int32
    t_indices: npt.NDArray[np.uint32]  # dtype=np.uint32 (GenomeAssembly contig index)
    t_starts: npt.NDArray[np.int32]  # dtype=np.int32
    t_ends: npt.NDArray[np.int32]  # dtype=np.int32
    strands: npt.NDArray[np.int8]  # dtype=np.int8
    is_expected: npt.NDArray[np.bool_]  # dtype=bool
    is_inside: npt.NDArray[np.bool_]  # dtype=bool
    is_extra: npt.NDArray[np.bool_]  # dtype=bool
    expected_positions: npt.NDArray[np.int32]
    expected_strands: npt.NDArray[np.int8]
    gene_ids: tuple[str, ...]
    cluster_names: tuple[str, ...]
    product_descriptions: tuple[str, ...]
    coverages: npt.NDArray[np.float32]

    @classmethod
    def empty(cls) -> "GeneHits":
        """Creates an empty GeneHits container."""
        return cls(
            np.empty(0, dtype=np.int32),
            np.empty(0, dtype=np.int32),
            np.empty(0, dtype=np.int32),
            np.empty(0, dtype=np.uint32),
            np.empty(0, dtype=np.int32),
            np.empty(0, dtype=np.int32),
            np.empty(0, dtype=np.int8),
            np.empty(0, dtype=bool),
            np.empty(0, dtype=bool),
            np.empty(0, dtype=bool),
            np.empty(0, dtype=np.int32),
            np.empty(0, dtype=np.int8),
            (),
            (),
            (),
            np.empty(0, dtype=np.float32),
        )

    @classmethod
    def concat(cls, batches: Iterable['GeneHits']) -> 'GeneHits':
        batches = list(batches)
        if not batches:
            return cls.empty()
        return cls(
            gene_indices=np.concatenate([b.gene_indices for b in batches]),
            q_starts=np.concatenate([b.q_starts for b in batches]),
            q_ends=np.concatenate([b.q_ends for b in batches]),
            t_indices=np.concatenate([b.t_indices for b in batches]),
            t_starts=np.concatenate([b.t_starts for b in batches]),
            t_ends=np.concatenate([b.t_ends for b in batches]),
            strands=np.concatenate([b.strands for b in batches]),
            is_expected=np.concatenate([b.is_expected for b in batches]),
            is_inside=np.concatenate([b.is_inside for b in batches]),
            is_extra=np.concatenate([b.is_extra for b in batches]),
            expected_positions=np.concatenate([b.expected_positions for b in batches]),
            expected_strands=np.concatenate([b.expected_strands for b in batches]),
            gene_ids=tuple(x for b in batches for x in b.gene_ids),
            cluster_names=tuple(x for b in batches for x in b.cluster_names),
            product_descriptions=tuple(x for b in batches for x in b.product_descriptions),
            coverages=np.concatenate([b.coverages for b in batches])
        )

    def __len__(self) -> int:
        return len(self.gene_indices)

    def __getitem__(self, item) -> "GeneHits":
        """Allows vectorised slicing and boolean masking of all arrays simultaneously."""
        return GeneHits(
            gene_indices=self.gene_indices[item],
            q_starts=self.q_starts[item],
            q_ends=self.q_ends[item],
            t_indices=self.t_indices[item],
            t_starts=self.t_starts[item],
            t_ends=self.t_ends[item],
            strands=self.strands[item],
            is_expected=self.is_expected[item],
            is_inside=self.is_inside[item],
            is_extra=self.is_extra[item],
            expected_positions=self.expected_positions[item],
            expected_strands=self.expected_strands[item],
            gene_ids=tuple(np.array(self.gene_ids, dtype=object)[item]),
            cluster_names=tuple(np.array(self.cluster_names, dtype=object)[item]),
            product_descriptions=tuple(
                np.array(self.product_descriptions, dtype=object)[item]
            ),
            coverages=self.coverages[item],
        )

    @property
    def frames(self):
        return (-self.q_starts) % 3

    @property
    def query_lengths(self) -> npt.NDArray[np.int32]:
        """Dynamic calculation of alignment lengths on the query assembly."""
        return self.q_ends - self.q_starts

    @property
    def target_lengths(self) -> npt.NDArray[np.int32]:
        """Dynamic calculation of alignment lengths on the database targets."""
        return self.t_ends - self.t_starts

    @property
    def q_intervals(self) -> Intervals:
        return Intervals(self.q_starts, self.q_ends, self.strands)

    @property
    def t_intervals(self) -> Intervals:
        return Intervals(self.t_starts, self.t_ends, self.strands)

    @classmethod
    def from_dict(cls, data: dict) -> "GeneHits":
        return cls(
            gene_indices=np.array(data["gene_indices"], dtype=np.int32),
            q_starts=np.array(data["q_starts"], dtype=np.int32),
            q_ends=np.array(data["q_ends"], dtype=np.int32),
            t_indices=np.array(data["t_indices"], dtype=np.uint32),
            t_starts=np.array(data["t_starts"], dtype=np.int32),
            t_ends=np.array(data["t_ends"], dtype=np.int32),
            strands=np.array(data["strands"], dtype=np.int8),
            is_expected=np.array(data["is_expected"], dtype=bool),
            is_inside=np.array(data["is_inside"], dtype=bool),
            is_extra=np.array(data["is_extra"], dtype=bool),
            expected_positions=np.array(
                data.get("expected_positions", []), dtype=np.int32
            ),
            expected_strands=np.array(data.get("expected_strands", []), dtype=np.int8),
            gene_ids=tuple(data.get("gene_ids", [])),
            cluster_names=tuple(data.get("cluster_names", [])),
            product_descriptions=tuple(data.get("product_descriptions", [])),
            coverages=np.array(data.get("coverages", []), dtype=np.float32),
        )

    def to_dict(self) -> dict:
        """Returns a dictionary containing the SoA arrays for orjson serialization."""
        return {
            k: getattr(self, k)
            for k in (
                "gene_indices",
                "q_starts",
                "q_ends",
                "t_indices",
                "t_starts",
                "t_ends",
                "strands",
                "is_expected",
                "is_inside",
                "is_extra",
                "expected_positions",
                "expected_strands",
                "gene_ids",
                "cluster_names",
                "product_descriptions",
                "coverages",
            )
        }


@dataclass(slots=True, frozen=True)
class LocusPieces(BatchedContainer[Any, 'LocusPieces']):
    """
    A high-performance SoA container for the bounding coordinates of the locus pieces.
    """

    ctg_indices: npt.NDArray[np.uint32]
    starts: npt.NDArray[np.int32]
    ends: npt.NDArray[np.int32]
    strands: npt.NDArray[np.int8]

    def __len__(self) -> int:
        return len(self.ctg_indices)

    def __getitem__(self, item: int | slice | npt.NDArray | list) -> 'Any | LocusPieces':
        if isinstance(item, (int, np.integer)):
            raise NotImplementedError("Single item access not implemented for LocusPieces")
        return LocusPieces(
            ctg_indices=self.ctg_indices[item],
            starts=self.starts[item],
            ends=self.ends[item],
            strands=self.strands[item]
        )

    @classmethod
    def concat(cls, batches: Iterable['LocusPieces']) -> 'LocusPieces':
        batches = list(batches)
        if not batches:
            return cls.empty()
        return cls(
            ctg_indices=np.concatenate([b.ctg_indices for b in batches]),
            starts=np.concatenate([b.starts for b in batches]),
            ends=np.concatenate([b.ends for b in batches]),
            strands=np.concatenate([b.strands for b in batches])
        )

    @classmethod
    def empty(cls) -> "LocusPieces":
        return cls(
            np.empty(0, dtype=np.uint32),
            np.empty(0, dtype=np.int32),
            np.empty(0, dtype=np.int32),
            np.empty(0, dtype=np.int8),
        )

    @classmethod
    def from_dict(cls, data: dict) -> "LocusPieces":
        return cls(
            ctg_indices=np.array(data["ctg_indices"], dtype=np.uint32),
            starts=np.array(data["starts"], dtype=np.int32),
            ends=np.array(data["ends"], dtype=np.int32),
            strands=np.array(data["strands"], dtype=np.int8),
        )

    def to_dict(self) -> dict:
        """Returns a dictionary containing the SoA arrays for orjson serialization."""
        return {
            k: getattr(self, k) for k in ("ctg_indices", "starts", "ends", "strands")
        }


@dataclass(slots=True, frozen=True)
class SerotypingResult:
    """
    Efficient, immutable container for a result generated by the Serotyper representing an _in silico_ serotyping call.

    Designed to be lightweight for JSON serialisation and SQL database storage, but contian enough information
    to re-create and inspect the result in memory.

    Several fields are included for the sake of tabular outputs
    that are not calculated by the Row formatters in the `kaptive.serotyping.io` module because they require
    information from `kaptive.database.Database` instances that may be unavailable, such as length discrepancy.

    This dataclass houses several SoA dataclasses for efficient downstream computation.
    """

    kaptive_version: str
    database_name: str
    database_version: str
    database_organism: str
    database_taxon: int
    genome: str
    best_locus_idx: int
    best_locus_name: str
    best_locus_score: float
    best_locus_completeness: float  # Proportion of genes found
    locus_pieces: LocusPieces
    length_discrepancy: float
    locus_seqs: Sequences
    gene_hits: GeneHits
    gene_states: npt.NDArray[np.int8]
    gene_seqs: Sequences
    translations: Sequences
    percent_identity: float
    percent_coverage: float
    protein_identities: npt.NDArray[np.float32]
    phenotype: str
    typeable: bool
    missing_expected_genes: tuple[str, ...]

    @property
    def problems(self):
        p = SerotypingProblem.NONE
        if len(self.locus_pieces) > 1:
            p |= SerotypingProblem.FRAGMENTED
        # Unexpected genes from other loci (not expected and not an allowed extra gene)
        if np.any(
            self.gene_hits.is_inside
            & ~self.gene_hits.is_expected
            & ~self.gene_hits.is_extra
        ):
            p |= SerotypingProblem.UNEXPECTED_GENES
        # Missing genes: Overall completeness is not 100%, or expected genes are found but outside the locus boundary
        if self.best_locus_completeness < 1.0 or np.any(
            ~self.gene_hits.is_inside & self.gene_hits.is_expected
        ):
            p |= SerotypingProblem.MISSING_GENES
        # Novel genes: inside boundary but below identity threshold
        if np.any(
            self.gene_hits.is_inside & (self.gene_states == GeneState.NOVEL.value)
        ):
            p |= SerotypingProblem.NOVEL_GENES
        # Truncated genes: inside boundary but partial or truncated
        if np.any(
            self.gene_hits.is_inside
            & (
                (self.gene_states == GeneState.TRUNCATED.value)
                | (self.gene_states == GeneState.PARTIAL.value)
            )
        ):
            p |= SerotypingProblem.TRUNCATED_GENES
        return p

    @classmethod
    def from_dict(cls, data: dict) -> "SerotypingResult":
        return cls(
            kaptive_version=data["kaptive_version"],
            database_name=data["database_name"],
            database_version=data["database_version"],
            database_organism=data["database_organism"],
            database_taxon=data["database_taxon"],
            genome=data["genome"],
            best_locus_idx=data["best_locus_idx"],
            best_locus_name=data["best_locus_name"],
            best_locus_score=data["best_locus_score"],
            best_locus_completeness=data["best_locus_completeness"],
            length_discrepancy=data["length_discrepancy"],
            locus_pieces=LocusPieces.from_dict(data["locus_pieces"]),
            gene_hits=GeneHits.from_dict(data["gene_hits"]),
            gene_states=np.array(data["gene_states"], dtype=np.int8),
            percent_identity=data["percent_identity"],
            percent_coverage=data["percent_coverage"],
            phenotype=data["phenotype"],
            typeable=data["typeable"],
            missing_expected_genes=tuple(data.get("missing_expected_genes", [])),
            locus_seqs=Sequences.from_dict(data["locus_seqs"]),
            gene_seqs=Sequences.from_dict(data["gene_seqs"]),
            translations=Sequences.from_dict(data["translations"]),
            protein_identities=np.array(data["protein_identities"], dtype=np.float32)
        )

    def to_locus_data(self) -> "LocusData":
        from kaptive.compare import LocusData

        mask = self.gene_hits.is_inside & ~self.gene_hits.is_extra
        
        return LocusData(
            proteins=self.translations[mask],
            name=self.genome,
            backbone=self.gene_hits.t_intervals[mask],
            pieces=self.locus_pieces,
            gene_ctg_indices=self.gene_hits.t_indices[mask],
        )

    def to_dict(self) -> dict:
        """Converts the result into a lightweight, JSON-safe dictionary."""
        return {
            "kaptive_version": self.kaptive_version,
            "database_name": self.database_name,
            "database_version": self.database_version,
            "database_organism": self.database_organism,
            "database_taxon": self.database_taxon,
            "genome": self.genome,
            "best_locus_idx": self.best_locus_idx,
            "best_locus_name": self.best_locus_name,
            "best_locus_score": self.best_locus_score,
            "best_locus_completeness": self.best_locus_completeness,
            "length_discrepancy": self.length_discrepancy,
            "percent_identity": self.percent_identity,
            "percent_coverage": self.percent_coverage,
            "phenotype": self.phenotype,
            "typeable": self.typeable,
            "missing_expected_genes": self.missing_expected_genes,
            "problems": int(self.problems),
            "locus_pieces": self.locus_pieces.to_dict(),
            "gene_hits": self.gene_hits.to_dict(),
            "gene_states": self.gene_states,
            "protein_identities": self.protein_identities,
            "locus_seqs": self.locus_seqs.to_dict(),
            "gene_seqs": self.gene_seqs.to_dict(),
            "translations": self.translations.to_dict()
        }


class Serotyper:
    """
    Performs _in silico_ serotyping on bacterial genome assemblies.
    """
    def __init__(
        self,
        db: Database,
        max_other_genes: int = 1, 
        min_completeness: float = 0.5, 
        allow_below_threshold: bool = False,
        preset: rammappy.Preset | None = None,
        scoring_metric: str = "scores",
        min_gene_coverage: float = 0.20,
    ):
        """Serotyper object that performs _in silico_ serotyping on bacterial genome assemblies.
        """
        self._db: Database = db
        self.max_other_genes: int = max_other_genes
        self.min_completeness: float = min_completeness
        self.allow_below_threshold: bool = allow_below_threshold
        self.preset: rammappy.Preset | None = preset
        self.scoring_metric: str = scoring_metric
        self.min_gene_coverage: float = min_gene_coverage
        self._protein_aligner = PairwiseAligner()
        
        # Count expected genes per locus for weighting
        self._expected_genes_per_locus = np.zeros(len(self._db.loci), dtype=np.float32)
        np.add.at(self._expected_genes_per_locus, self._db.gene_locus_indices[~self._db.extra_genes], 1.0)
        self._expected_genes_per_locus = np.maximum(self._expected_genes_per_locus, 1.0)
        
        # Prepare metadata for alignment iterator parsing
        self._gene_seqs = [
            (str(i).encode(), bytes(self._db.genes.seqs[self._db.genes.offsets[i] : self._db.genes.offsets[i] + self._db.genes.lengths[i]]))
            for i in range(len(self._db.genes))
        ]
        self._gene_meta = [(str(i), self._db.genes.lengths[i]) for i in range(len(self._db.genes))]


    def __call__(self, genome: GenomeAssembly | str | Path) -> SerotypingResult | None:
        genome = GenomeAssembly.ensure(genome)

        contig_index = genome.get_rammappy_index()
        aligner = Aligner(index=contig_index, preset=None, do_cigar=True, do_cs=False, do_md=False)
        opts = aligner.options
        opts.filtering.best_n = 50000
        opts.filtering.pri_ratio = 0.0
        aligner.options = opts

        gene_aln_batch = aligner.map_batch(self._gene_seqs)
        gene_alns = Alignments.from_mapping_iterators(self._gene_meta, gene_aln_batch)

        # Calculate total coverage per gene across all alignments for reporting
        q_indices = gene_alns.q_names.astype(np.int32)
        q_lengths = gene_alns.q_aln_lens
        total_q_covs = np.zeros(len(self._db.genes), dtype=np.float32)
        np.add.at(total_q_covs, q_indices, q_lengths)
        total_q_covs /= self._db.genes.lengths

        # Scoring phase ------------------------------------------------------------------------------------------------
        # Calculate alignment query coverage per alignment
        q_covs = gene_alns.q_covs
        valid_cov_mask = q_covs >= self.min_gene_coverage

        valid_alns = gene_alns[valid_cov_mask]
        valid_q_covs = q_covs[valid_cov_mask]
        valid_gene_indices = valid_alns.q_names.astype(np.int32)

        # Sort to find the best alignment per gene (highest q_cov, tie-breaker: highest score)
        order = np.lexsort((-valid_alns.scores, -valid_q_covs, valid_gene_indices))
        valid_alns = valid_alns[order]
        valid_gene_indices = valid_gene_indices[order]
        valid_q_covs = valid_q_covs[order]
        
        # Take the first (best) alignment for each gene
        _, unique_indices = np.unique(valid_gene_indices, return_index=True)
        best_alns = valid_alns[unique_indices]
        best_gene_indices = valid_gene_indices[unique_indices]
        best_q_covs = valid_q_covs[unique_indices]

        valid_locus_indices = self._db.gene_locus_indices[best_gene_indices]
        valid_not_extra = ~self._db.extra_genes[best_gene_indices]

        # Calculate locus scores by summing best q_covs of valid expected genes
        locus_scores = np.zeros(len(self._db.loci), dtype=np.float64)
        np.add.at(locus_scores, valid_locus_indices[valid_not_extra], best_q_covs[valid_not_extra])

        # Calculate completeness count per locus (unique expected genes matched per locus)
        locus_counts = np.zeros(len(self._db.loci), dtype=np.float32)
        matched_expected_genes = best_gene_indices[valid_not_extra]
        np.add.at(locus_counts, self._db.gene_locus_indices[matched_expected_genes], 1.0)

        locus_completeness = locus_counts / self._expected_genes_per_locus
        final_locus_scores = locus_scores * (locus_completeness ** 3)

        self._last_scores = final_locus_scores.copy()
        self._last_completeness = locus_completeness.copy()
        
        best_locus_idx = np.argmax(final_locus_scores)
        best_locus_name = self._db.loci.ids[best_locus_idx]
        best_locus_completeness = locus_completeness[best_locus_idx]

        # Reconstruction phase -----------------------------------------------------------------------------------------
        valid_alns = gene_alns
        
        # Cull alignments, prioritizing genes belonging to the best match locus
        valid_indices = valid_alns.q_names.astype(np.int32)
        priority_mask = self._db.gene_locus_indices[valid_indices] == best_locus_idx
        
        culled_alns = valid_alns.cull_overlaps(by_query=False, priority_mask=priority_mask, max_overlap_fraction=0.1)
        
        # Re-extract arrays for the culled batch
        culled_gene_indices = culled_alns.q_names.astype(np.int32)
        t_indices = np.array([genome.id_map[n] for n in culled_alns.t_names], dtype=np.uint32)
        # Cluster intervals by contig using max_locus_length as tolerance
        culled_intervals = culled_alns.to_intervals(by_query=False)
        piece_ids = culled_intervals.cluster_spatial(tolerance=self._db.max_locus_length, group_by=t_indices)
        
        # Identify expected genes and the clusters they fall into
        is_expected = (self._db.gene_locus_indices[culled_gene_indices] == best_locus_idx) & ~self._db.extra_genes[
            culled_gene_indices
        ]
        valid_cluster_ids = np.unique(piece_ids[is_expected])
        is_extra = self._db.extra_genes[culled_gene_indices]

        # Calculate gene alignment coverages based on total coverage across all non-overlapping fragments
        coverages = np.clip(total_q_covs[culled_gene_indices] * 100.0, 0.0, 100.0)

        # Identify the primary hit for each expected gene for bounding box calculation
        primary_expected = np.zeros(len(culled_alns), dtype=bool)
        is_expected_hits = np.where(is_expected)[0]
        if len(is_expected_hits) > 0:
            exp_gene_indices = culled_gene_indices[is_expected_hits]
            exp_scores = culled_alns.scores[is_expected_hits]
            order = np.lexsort((-exp_scores, exp_gene_indices))
            sorted_exp_gene_indices = exp_gene_indices[order]
            _, unique_indices = np.unique(sorted_exp_gene_indices, return_index=True)
            best_hits = is_expected_hits[order[unique_indices]]
            primary_expected[best_hits] = True

        # Construct bounding locus pieces from valid pieces
        l_ctg_indices, l_starts, l_ends, l_strands = [], [], [], []
        l_expected_means = []
        for c_id in valid_cluster_ids:
            piece_mask = piece_ids == c_id
            
            piece_primary = piece_mask & primary_expected
            if np.any(piece_primary):
                ctg_idx = t_indices[piece_mask][0]
                l_ctg_indices.append(ctg_idx)
                start = np.min(culled_intervals.starts[piece_primary])
                end = np.max(culled_intervals.ends[piece_primary])
                l_starts.append(start)
                l_ends.append(end)
                
                exp_genes = culled_gene_indices[piece_primary]
                l_expected_means.append(np.mean(self._db.gene_positions[exp_genes]))

                exp_strands = self._db.gene_intervals.strands[exp_genes]
                found_strands = culled_alns.strands[piece_primary]
                if np.sum(found_strands * exp_strands) < 0:
                    l_strands.append(-1)
                else:
                    l_strands.append(1)

        # Recompute is_inside using the strict expected boundaries
        is_inside = np.zeros(len(culled_alns), dtype=bool)
        for ctg_idx, start, end in zip(l_ctg_indices, l_starts, l_ends):
            on_ctg = t_indices == ctg_idx
            # An alignment is inside if it overlaps the bounding region
            inside_this = on_ctg & (culled_intervals.starts <= end) & (culled_intervals.ends >= start)
            is_inside |= inside_this

        # Sort pieces by expected mean position
        piece_order = np.argsort(l_expected_means)

        locus_pieces = LocusPieces(
            ctg_indices=np.array(l_ctg_indices, dtype=np.uint32)[piece_order],
            starts=np.array(l_starts, dtype=np.int32)[piece_order],
            ends=np.array(l_ends, dtype=np.int32)[piece_order],
            strands=np.array(l_strands, dtype=np.int8)[piece_order],
        )

        # Identify missing expected genes
        expected_genes_mask = (self._db.gene_locus_indices == best_locus_idx) & ~self._db.extra_genes
        expected_gene_indices = np.where(expected_genes_mask)[0]
        # Which ones did we find inside the locus?
        found_expected_gene_indices = culled_gene_indices[is_expected & is_inside]
        missing_indices = np.setdiff1d(expected_gene_indices, found_expected_gene_indices, assume_unique=True)
        missing_expected_genes = tuple(self._db.genes.ids[i] for i in missing_indices)
        
        # Calculate actual completeness based on reconstructed locus
        actual_locus_completeness = 1.0 - (len(missing_indices) / len(expected_gene_indices)) if len(expected_gene_indices) > 0 else 1.0


        gene_hits = GeneHits(
            gene_indices=culled_gene_indices,
            q_starts=culled_alns.q_starts,
            q_ends=culled_alns.q_ends,
            t_indices=t_indices,
            t_starts=culled_alns.t_starts,
            t_ends=culled_alns.t_ends,
            strands=culled_alns.strands,
            is_expected=is_expected,
            is_inside=is_inside,
            is_extra=is_extra,
            expected_positions=self._db.gene_positions[culled_gene_indices].astype(np.int32),
            expected_strands=self._db.gene_intervals.strands[culled_gene_indices],
            gene_ids=tuple(self._db.genes.ids[i] for i in culled_gene_indices),
            cluster_names=tuple(self._db.cluster_keys[self._db.gene_cluster_ids[i]] for i in culled_gene_indices),
            product_descriptions=tuple(
                self._db.description_keys[self._db.gene_description_ids[i]] for i in culled_gene_indices
            ),
            coverages=coverages,
        )

        # Locus extraction phase ---------------------------------------------------------------------------------------
        if len(locus_pieces) > 0:  # Extract locus sequences using the batched SoA locus pieces
            locus_seqs = genome.contigs.extract(
                locus_pieces.ctg_indices,
                locus_pieces.starts,
                locus_pieces.ends,
                locus_pieces.strands,
            )
        else:
            locus_seqs = Sequences.empty()

        # Calculate coverage and length discrepancy
        assem_len = np.sum(locus_pieces.ends - locus_pieces.starts)
        ref_len = self._db.loci.lengths[best_locus_idx]
        pcov = float(min(100.0, (assem_len / ref_len) * 100.0)) if ref_len > 0 else 0.0
        if len(locus_pieces) == 1:
            length_discrepancy = float(assem_len - ref_len)
        else:
            length_discrepancy = float("nan")

        # Gene state phase ---------------------------------------------------------------------------------------------
        gene_seqs = genome.contigs.extract_intervals(  # Extract gene nucleotides from their contigs
            gene_hits.t_indices,
            gene_hits.t_intervals,
            new_ids=tuple(self._db.genes.ids[i] for i in gene_hits.gene_indices),
        )
        # Translate nucleotides to amino acids, compensating for the reading frames of the alignments
        # Truncate at the first stop codon to match Old Kaptive's behavior, which prevents frameshifts 
        # from pulling down the identity score of the valid upstream alignment.
        prot_seqs = gene_seqs.translate(frames=gene_hits.frames, to_stop=True)
        
        # Initialize states
        gene_states = np.full(len(gene_hits), GeneState.NORMAL.value, dtype=np.int8)
        is_partial = culled_alns.is_partial
        db_gene_lengths = self._db.genes.lengths[gene_hits.gene_indices]
        
        # In Old Kaptive, truncation was calculated based on the *translated protein* length 
        # (up to the first stop codon) divided by the reference protein length.
        prot_covs = (prot_seqs.lengths * 3.0) / db_gene_lengths
        
        # We must also update the coverage to reflect the protein coverage, so it prints correctly in the output
        gene_hits.coverages[:] = np.clip(prot_covs * 100.0, 0.0, 100.0)
        
        # A partial gene colliding with a contig edge is excluded from being truncated
        is_truncated = (~is_partial) & (prot_covs < 0.90)
        gene_states[is_partial] = GeneState.PARTIAL.value
        gene_states[is_truncated] = GeneState.TRUNCATED.value
        prot_alns = self._protein_aligner(prot_seqs, self._db.translations[gene_hits.gene_indices])
        prot_idents = prot_alns.pidents.astype(np.float32)

        # Drop genes outside the locus that fall below the identity threshold,
        # mirroring Old Kaptive's behavior of ignoring spurious homologies
        is_spurious = (~gene_hits.is_inside) & (prot_idents < self._db.metadata.id_threshold)
        if np.any(is_spurious):
            keep_mask = ~is_spurious
            gene_hits = gene_hits[keep_mask]
            gene_seqs = gene_seqs[keep_mask]
            prot_seqs = prot_seqs[keep_mask]
            gene_states = gene_states[keep_mask]
            prot_idents = prot_idents[keep_mask]

        # Normal genes that fall below the identity threshold are considered NOVEL
        below_threshold = (gene_states == GeneState.NORMAL.value) & (prot_idents < self._db.metadata.id_threshold)
        gene_states[below_threshold] = GeneState.NOVEL.value
        valid_pidents = prot_idents[gene_states == GeneState.NORMAL.value]
        pident = float(np.mean(valid_pidents)) if valid_pidents.size > 0 else 0.0

        # Phenotype Evaluation phase -----------------------------------------------------------------------------------
        base_phenotype = self._db.serotypes[best_locus_idx]
        phenotypes = self._db.phenotypes

        if len(phenotypes) > 0:
            # A cluster is considered 'active' if it's found NORMAL or PARTIAL
            q_active = np.zeros(len(self._db.cluster_keys), dtype=bool)
            is_active = (gene_states == GeneState.NORMAL.value) | (gene_states == GeneState.PARTIAL.value)
            if np.any(is_active):
                active_clusters = self._db.gene_cluster_ids[gene_hits.gene_indices[is_active]]
                q_active[active_clusters] = True

            # Vectorized rule evaluation
            locus_match = phenotypes.locus_masks[:, best_locus_idx]
            extra_match = (phenotypes.extra_masks & ~q_active).sum(axis=1) == 0
            inactive_match = (phenotypes.inactive_masks & q_active).sum(axis=1) == 0

            if np.any(valid_mask := locus_match & extra_match & inactive_match):
                valid_indices = np.where(valid_mask)[0]
                is_suffix = phenotypes.as_suffix[valid_indices]

                if len(replacements := valid_indices[~is_suffix]) > 0:
                    best_rep_idx = replacements[np.argmax(phenotypes.priorities[replacements])]
                    base_phenotype = phenotypes.ids[best_rep_idx]

                if len(suffixes := valid_indices[is_suffix]) > 0:
                    sorted_suffixes = suffixes[np.argsort(-phenotypes.priorities[suffixes])]
                    suffix_strs = [phenotypes.ids[i] for i in sorted_suffixes]
                    base_phenotype = f"{base_phenotype}{''.join(suffix_strs)}"

        # Confidence evaluation phase ----------------------------------------------------------------------------------
        typeable = True
        if actual_locus_completeness < self.min_completeness:
            typeable = False
            
        # Truncated unexpected genes do not count towards the limit, mirroring Old Kaptive
        is_unexpected = gene_hits.is_inside & ~gene_hits.is_expected & ~gene_hits.is_extra
        is_not_truncated = gene_states != GeneState.TRUNCATED.value
        unexpected_count = np.count_nonzero(is_unexpected & is_not_truncated)
        if unexpected_count > self.max_other_genes:
            typeable = False

        # 3. Check for any genes falling below the identity threshold
        if not self.allow_below_threshold:
            if np.any(gene_hits.is_inside & (gene_states == GeneState.NOVEL.value)):
                typeable = False

        # Return result object -----------------------------------------------------------------------------------------
        return SerotypingResult(
            kaptive_version=__version__,
            database_name=self._db.metadata.name,
            database_version=self._db.metadata.version,
            database_organism=self._db.metadata.organism,
            database_taxon=self._db.metadata.taxon,
            genome=genome.id,
            best_locus_idx=best_locus_idx,
            best_locus_name=best_locus_name,
            best_locus_score=locus_scores[best_locus_idx],
            best_locus_completeness=actual_locus_completeness,
            length_discrepancy=length_discrepancy,
            gene_hits=gene_hits,
            gene_states=gene_states,
            locus_pieces=locus_pieces,
            locus_seqs=locus_seqs,
            gene_seqs=gene_seqs,
            translations=prot_seqs,
            percent_identity=pident,
            percent_coverage=pcov,
            protein_identities=prot_idents,
            phenotype=base_phenotype,
            typeable=typeable,
            missing_expected_genes=missing_expected_genes,
        )


@dataclass(slots=True, frozen=True)
class ReportRow(ABC):
    """Base class for representing an in silico serotyping report for a single sample in a TSV row.
    Allows for the auto-documentation of columns using the class docstrings under the "Attributes" section."""

    @classmethod
    def header(cls) -> "bytes":
        """Returns the TSV header as UTF-8 encoded bytes for fast binary I/O."""
        return ("\t".join(f.name for f in fields(cls)) + "\n").encode("utf-8")

    def __bytes__(self) -> bytes:
        """Formats the row data into a tab-separated string."""
        return b"\t".join(getattr(self, f.name) for f in fields(self)) + b"\n"

    @classmethod
    @abstractmethod
    def from_result( cls, result: SerotypingResult) -> "Self": ...


@dataclass(slots=True, frozen=True)
class KaptiveRow(ReportRow):
    """Represents an _in silico_ serotyping report for a single sample in a TSV row in the original Kaptive format.

    Attributes:
        Kaptive_version (bytes): The version of Kaptive used to perform _in silico_ serotyping.
        Database_name (bytes): The name of the database used for _in silico_ serotyping.
        Database_version (bytes): The version of the database used for _in silico_ serotyping.
        Assembly (bytes): The name of the sample, taken from the gneome assembly filename.
        Best_match_locus (bytes): The locus type which most closely matches the genome.
        Best_match_type (bytes): The predicted serotype/phenotype of the genome.
        Match_confidence (bytes): A binary measure of whether the serotyping call was "Typeable" or "Untypeable".
        Problems (bytes): Characters representing `SerotypingProblem`s  indicating issues with the locus match.
        Identity (bytes): Weighted percent identity of the best matching locus to the genome.
        Coverage (bytes): Weighted percent coverage of the best matching locus in the genome.
        Length_discrepancy (bytes): If the locus was found in a single piece,
            this is the difference between the locus length and the genome length.
        Expected_genes_in_locus (bytes): A fraction indicating how many of the genes in the best matching locus were
            found in the locus part of the genome.
        Expected_genes_in_locus_details (bytes): Gene names for the expected genes found in the locus part of the genome.
        Missing_expected_genes (bytes): A string listing the gene names of expected genes that were not found.
        Other_genes_in_locus (bytes): The number of unexpected genes (genes from loci other than the best match) which
            were found in the locus part of the genome.
        Other_genes_in_locus_details (bytes): Gene names for the other genes found in the locus part of the genome.
        Expected_genes_outside_locus (bytes): A fraction indicating how many of the expected genes which were found in
            the genome but not in the locus part of the genome (usually zero).
        Expected_genes_outside_locus_details (bytes):  Gene names for the expected genes found outside the locus part of
            the genome.
        Other_genes_outside_locus (bytes): The number of unexpected genes (genes from loci other than the best match)
            which were found outside the locus part of the genome.
        Other_genes_outside_locus_details (bytes): Gene names for the other genes found outside the locus part of the
            genome.
        Truncated_genes_details (bytes): Gene names for the truncated genes found in the genome.
        Extra_genes_details (bytes): Gene names for the extra genes found in the genome.

    Note:
        Numbers beside gene names indicate the percent identity and percent coverage of the gene in the genome.

    Warning:
        You may sometimes see two copies of the same gene in the `Expected genes in locus, details` column.
        These represent (likely) parts of the same gene which have usually been split over contigs. In
        `kaptive v3.0.0` onwards, we adopted this behavior to allow users to see where locus splitting has occurred,
        and determine the total percent identity of a gene that has been split.
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
    def header(cls) -> "bytes":
        """Returns the TSV header formatted for backwards compatibility."""
        headers = [
            f.name.encode("utf-8")
            .replace(b"_details", b", details")
            .replace(b"_", b" ")
            for f in fields(cls)
        ]
        return b"\t".join(headers) + b"\n"

    @classmethod
    def from_result(cls, result: SerotypingResult):
        hits = result.gene_hits
        states = result.gene_states

        # High-performance boolean masks
        in_loc = hits.is_inside
        out_loc = ~hits.is_inside
        exp = hits.is_expected
        extra = hits.is_extra
        unexp = ~exp & ~extra

        def _format_genes(mask: np.ndarray) -> bytes:
            """Helper to rapidly construct the details string for a specific masked subset of genes."""
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


@dataclass(slots=True, frozen=True)
class Pha4geRow(ReportRow):
    """Represents an in silico serotyping report for a single sample in a TSV row in the PHA4GE format.

    Attributes:
        sample (bytes): The name of the sample, taken from the assembly filename.
        genotyping_method (bytes): The genotyping method (in silico serotyping)
        genotyping_schema_taxon (bytes): The taxonomy of the database used for _in silico_ serotyping.
        genotyping_database_name (bytes): The name of the database used for _in silico_ serotyping.
        genotyping_database_version (bytes): The version of the database used for _in silico_ serotyping.
        genotyping_schema_name (bytes): The name of the database used for _in silico_ serotyping.
        genotyping_software_name (bytes): The name of the software (Kaptive)
        genotyping_software_version (bytes): The version of Kaptive used to perform _in silico_ serotyping.
        genotype (bytes): The locus type which most closely matches the genome.
        genotype_confidence_value (bytes): A binary measure of whether the serotyping call was "Typeable" or "Untypeable".
        genotype_predicted_phenotype (bytes): The predicted serotype/phenotype of the genome.
        genotyping_details (bytes): The details of the serotyping result.
    """

    sample: bytes
    genotyping_schema_taxon: bytes
    genotyping_database_name: bytes
    genotyping_database_version: bytes
    genotyping_software_version: bytes
    genotype: bytes
    genotype_confidence_value: bytes
    genotype_predicted_phenotype: bytes
    genotyping_details: bytes
    genotyping_method: bytes = b"In silico serotyping"
    genotyping_schema_name: bytes = b"Kaptive"
    genotyping_software_name: bytes = b"Kaptive"
    genotyping_method_url: bytes = b"https://github.com/klebgenomics/Kaptive"

    @classmethod
    def from_result(cls, result: SerotypingResult):

        # Transform problem enum into human-readable string
        if result.problems:
            detail_parts = []
            if SerotypingProblem.TRUNCATED_GENES in result.problems:
                detail_parts.append(b'truncated gene(s) in locus')
            if SerotypingProblem.NOVEL_GENES in result.problems:
                detail_parts.append(b'gene(s) below identity threshold')
            if SerotypingProblem.FRAGMENTED in result.problems:
                detail_parts.append(b'locus fragmented into %d pieces' % len(result.locus_pieces))
            if SerotypingProblem.MISSING_GENES in result.problems:
                detail_parts.append(b'missing expected gene(s)')
            if SerotypingProblem.UNEXPECTED_GENES in result.problems:
                detail_parts.append(b"unexpected gene(s) in locus")
            details = b'Problems: %b' % b', '.join(detail_parts)
        else:
            details = b''

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
            genotyping_details=details
        )
