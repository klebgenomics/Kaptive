r"""Data models and container classes for serotyping analysis.

This module provides enumerations for gene classification states and serotyping quality problems,
along with high-performance Structure-of-Arrays (SoA) containers and data models for storing
gene alignment hits, locus fragments, and complete serotyping analysis results.

Classes:
    [`GeneState`][kaptive.serotyping.models.GeneState]: Enumeration of mutually exclusive gene classification states.
    [`SerotypingProblem`][kaptive.serotyping.models.SerotypingProblem]: Bitflag enumeration of quality problems.
    [`GeneHits`][kaptive.serotyping.models.GeneHits]: SoA container storing classified gene alignments.
    [`LocusPieces`][kaptive.serotyping.models.LocusPieces]: SoA container storing bounding coordinates.
    [`SerotypingResult`][kaptive.serotyping.models.SerotypingResult]: Complete immutable serotyping result record.
"""

from collections.abc import Iterable
from dataclasses import dataclass, field
from enum import IntEnum, IntFlag, auto
from typing import TYPE_CHECKING, Any, ClassVar, Self

import numpy as np
import numpy.typing as npt

from kaptive.core.collections import BatchedContainer
from kaptive.core.interval import Intervals
from kaptive.core.seq import Sequences

if TYPE_CHECKING:
    from kaptive.compare import LocusData


class GeneState(IntEnum):
    r"""Mutually exclusive states for locus genes found in a genome assembly.

    Attributes:
        NORMAL (int): The gene was found intact as expected.
        PARTIAL (int): The gene was broken up over a contig edge.
        TRUNCATED (int): The gene does not form a complete amino acid sequence.
        NOVEL (int): The gene translation diverges significantly from the closest reference.
    """

    NORMAL = 0
    PARTIAL = 1
    TRUNCATED = 2
    NOVEL = 3


class SerotypingProblem(IntFlag):
    r"""Symbolic problems with the serotype call used for report formatting.

    Bitflag values represent distinct issues detected during locus assembly analysis and
    can be combined bitwise.

    Attributes:
        NONE (int): No problems detected in the serotype call.
        FRAGMENTED (int): Locus is broken up into multiple pieces across contigs (Symbol: `?`).
        UNEXPECTED_GENES (int): Unexpected genes from non-target loci present inside locus boundary (Symbol: `+`).
        MISSING_GENES (int): Expected genes from target locus missing inside locus boundary (Symbol: `-`).
        NOVEL_GENES (int): Genes inside locus boundary falling below identity threshold (Symbol: `*`).
        TRUNCATED_GENES (int): Genes inside locus boundary that are truncated or partial (Symbol: `!`).
        SYMBOLS (ClassVar[tuple[bytes, ...]]): Precomputed lookup table mapping integer bitflag combinations
            to symbol byte strings.
    """

    NONE = 0
    FRAGMENTED = auto()
    UNEXPECTED_GENES = auto()
    MISSING_GENES = auto()
    NOVEL_GENES = auto()
    TRUNCATED_GENES = auto()

    SYMBOLS: ClassVar[tuple[bytes, ...]]

    def to_symbols(self) -> bytes:
        r"""Render the bitflag combination into formatted symbol bytes for TSV reporting.

        Returns:
            bytes: ASCII byte string containing concatenation of active problem symbols.
        """
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
class GeneHits(BatchedContainer[Any, "GeneHits"]):
    r"""A high-performance SoA container for classified gene alignments.

    Encapsulates parallel NumPy arrays and metadata tuples for gene alignments, enabling
    synchronized vectorised filtering and dynamic interval calculations. Inherits from
    [`BatchedContainer`][kaptive.core.collections.BatchedContainer].

    Attributes:
        gene_indices (npt.NDArray[np.int32]): Global database gene indices.
        q_starts (npt.NDArray[np.int32]): Alignment start positions on query contigs (0-indexed).
        q_ends (npt.NDArray[np.int32]): Alignment end positions on query contigs (0-indexed).
        t_indices (npt.NDArray[np.uint32]): Target contig indices in genome assembly.
        t_starts (npt.NDArray[np.int32]): Alignment start positions on database reference genes.
        t_ends (npt.NDArray[np.int32]): Alignment end positions on database reference genes.
        strands (npt.NDArray[np.int8]): Alignment strand orientations (+1 or -1).
        is_expected (npt.NDArray[np.bool_]): Boolean mask indicating expected locus genes.
        is_inside (npt.NDArray[np.bool_]): Boolean mask indicating hits within locus boundaries.
        is_extra (npt.NDArray[np.bool_]): Boolean mask indicating extra allowed genes.
        expected_positions (npt.NDArray[np.int32]): Expected relative gene order positions.
        expected_strands (npt.NDArray[np.int8]): Expected strand orientations (+1 or -1).
        gene_ids (npt.NDArray[np.bytes_]): 1D byte string array (`S32`) of gene identifier strings.
        cluster_names (npt.NDArray[np.bytes_]): 1D byte string array (`S10`) of gene cluster or family names.
        product_descriptions (npt.NDArray[np.bytes_]): 1D byte string array (`S64`) of functional gene product annotations.
        coverages (npt.NDArray[np.float32]): Gene alignment coverage proportions.
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
    gene_ids: npt.NDArray[np.bytes_]
    cluster_names: npt.NDArray[np.bytes_]
    product_descriptions: npt.NDArray[np.bytes_]
    coverages: npt.NDArray[np.float32]

    def __post_init__(self) -> None:
        for field_name, dtype in (
            ("gene_ids", "S32"),
            ("cluster_names", "S10"),
            ("product_descriptions", "S64"),
        ):
            val = getattr(self, field_name)
            if not isinstance(val, np.ndarray) or val.dtype.kind not in ("S", "a"):
                if isinstance(val, np.ndarray) and val.dtype.kind == "U":
                    encoded = [x.encode("utf-8") for x in val.flat]
                    arr = np.array(encoded, dtype=dtype).reshape(val.shape)
                elif isinstance(val, (list, tuple)):
                    encoded = [x.encode("utf-8") if isinstance(x, str) else x for x in val]
                    arr = np.array(encoded, dtype=dtype)
                else:
                    arr = np.asarray(val, dtype=dtype)
                object.__setattr__(self, field_name, arr)

    @classmethod
    def empty(cls) -> "GeneHits":
        r"""Create an empty `GeneHits` container with zero-length arrays and empty tuples.

        Returns:
            GeneHits: An empty [`GeneHits`][kaptive.serotyping.models.GeneHits] instance.
        """
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
            np.empty(0, dtype="S32"),
            np.empty(0, dtype="S10"),
            np.empty(0, dtype="S64"),
            np.empty(0, dtype=np.float32),
        )

    @classmethod
    def concat(cls, batches: Iterable[Self]) -> Self:  # type: ignore
        r"""Concatenate multiple `GeneHits` batches into a single container.

        Args:
            batches (Iterable[GeneHits]): An iterable of [`GeneHits`][kaptive.serotyping.models.GeneHits] instances.

        Returns:
            GeneHits: Combined [`GeneHits`][kaptive.serotyping.models.GeneHits] container.
        """
        batches_list = list(batches)
        if not batches_list:
            return cls.empty()  # type: ignore
        return cls(
            gene_indices=np.concatenate([b.gene_indices for b in batches_list]),
            q_starts=np.concatenate([b.q_starts for b in batches_list]),
            q_ends=np.concatenate([b.q_ends for b in batches_list]),
            t_indices=np.concatenate([b.t_indices for b in batches_list]),
            t_starts=np.concatenate([b.t_starts for b in batches_list]),
            t_ends=np.concatenate([b.t_ends for b in batches_list]),
            strands=np.concatenate([b.strands for b in batches_list]),
            is_expected=np.concatenate([b.is_expected for b in batches_list]),
            is_inside=np.concatenate([b.is_inside for b in batches_list]),
            is_extra=np.concatenate([b.is_extra for b in batches_list]),
            expected_positions=np.concatenate([b.expected_positions for b in batches_list]),
            expected_strands=np.concatenate([b.expected_strands for b in batches_list]),
            gene_ids=np.concatenate([b.gene_ids for b in batches_list]),
            cluster_names=np.concatenate([b.cluster_names for b in batches_list]),
            product_descriptions=np.concatenate([b.product_descriptions for b in batches_list]),
            coverages=np.concatenate([b.coverages for b in batches_list]),
        )

    def __len__(self) -> int:
        r"""Return total number of gene hit alignments in container.

        Returns:
            int: Number of elements in parallel arrays.
        """
        return len(self.gene_indices)

    def __getitem__(self, item: Any) -> "GeneHits":
        r"""Slice or boolean-mask all parallel array fields simultaneously.

        Args:
            item (Any): Slice, integer array, or boolean mask.

        Returns:
            GeneHits: A sliced [`GeneHits`][kaptive.serotyping.models.GeneHits] instance.
        """
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
            gene_ids=self.gene_ids[item],
            cluster_names=self.cluster_names[item],
            product_descriptions=self.product_descriptions[item],
            coverages=self.coverages[item],
        )

    @property
    def frames(self) -> npt.NDArray[np.int32]:
        r"""Calculate reading frame offsets for query alignments.

        Returns:
            npt.NDArray[np.int32]: Reading frame offsets calculated as `(-q_starts) % 3`.
        """
        return (-self.q_starts) % 3

    @property
    def query_lengths(self) -> npt.NDArray[np.int32]:
        r"""Calculate alignment spans on query assembly contigs.

        Returns:
            npt.NDArray[np.int32]: Alignment spans calculated as `q_ends - q_starts`.
        """
        return self.q_ends - self.q_starts

    @property
    def target_lengths(self) -> npt.NDArray[np.int32]:
        r"""Calculate alignment spans on database target references.

        Returns:
            npt.NDArray[np.int32]: Alignment spans calculated as `t_ends - t_starts`.
        """
        return self.t_ends - self.t_starts

    @property
    def q_intervals(self) -> Intervals:
        r"""Construct query genomic intervals container.

        Returns:
            Intervals: An [`Intervals`][kaptive.core.interval.Intervals] object wrapper for `q_starts`, `q_ends`,
                and `strands`.
        """
        return Intervals(self.q_starts, self.q_ends, self.strands)

    @property
    def t_intervals(self) -> Intervals:
        r"""Construct database target intervals container.

        Returns:
            Intervals: An [`Intervals`][kaptive.core.interval.Intervals] object wrapper for `t_starts`, `t_ends`,
                and `strands`.
        """
        return Intervals(self.t_starts, self.t_ends, self.strands)

    @classmethod
    def from_dict(cls, data: dict[str, Any]) -> "GeneHits":
        r"""Reconstruct a `GeneHits` container from a deserialized dictionary.

        Args:
            data (dict[str, Any]): Dictionary containing array data lists and tuple metadata.

        Returns:
            GeneHits: Reconstructed [`GeneHits`][kaptive.serotyping.models.GeneHits] instance.
        """

        def _to_bytes_array(val: Any, dtype: str) -> npt.NDArray[np.bytes_]:
            if val is None or len(val) == 0:
                return np.empty(0, dtype=dtype)
            if isinstance(val, np.ndarray) and val.dtype.kind in ("S", "a"):
                return val.astype(dtype)
            encoded = [x.encode("utf-8") if isinstance(x, str) else x for x in val]
            return np.array(encoded, dtype=dtype)

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
            expected_positions=np.array(data.get("expected_positions", []), dtype=np.int32),
            expected_strands=np.array(data.get("expected_strands", []), dtype=np.int8),
            gene_ids=_to_bytes_array(data.get("gene_ids", []), "S32"),
            cluster_names=_to_bytes_array(data.get("cluster_names", []), "S10"),
            product_descriptions=_to_bytes_array(data.get("product_descriptions", []), "S64"),
            coverages=np.array(data.get("coverages", []), dtype=np.float32),
        )

    def to_dict(self) -> dict[str, Any]:
        r"""Convert SoA array fields to a dictionary for JSON serialization.

        Returns:
            dict[str, Any]: Dictionary mapping field names to NumPy arrays and metadata tuples.
        """
        d = {
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
                "coverages",
            )
        }
        d["gene_ids"] = np.char.decode(self.gene_ids, "utf-8").tolist()
        d["cluster_names"] = np.char.decode(self.cluster_names, "utf-8").tolist()
        d["product_descriptions"] = np.char.decode(self.product_descriptions, "utf-8").tolist()
        return d


@dataclass(slots=True, frozen=True)
class LocusPieces(BatchedContainer[Any, "LocusPieces"]):
    r"""A high-performance SoA container for bounding coordinates of locus fragments.

    Stores contig indices, coordinate spans, and strand directions for locus pieces when a locus
    is fragmented across multiple contigs.
    Inherits from [`BatchedContainer`][kaptive.core.collections.BatchedContainer].

    Attributes:
        ctg_indices (npt.NDArray[np.uint32]): Target contig indices in assembly.
        starts (npt.NDArray[np.int32]): Locus fragment start coordinates (0-indexed).
        ends (npt.NDArray[np.int32]): Locus fragment end coordinates (0-indexed).
        strands (npt.NDArray[np.int8]): Locus fragment strand orientations (+1 or -1).
    """

    ctg_indices: npt.NDArray[np.uint32]
    starts: npt.NDArray[np.int32]
    ends: npt.NDArray[np.int32]
    strands: npt.NDArray[np.int8]

    def __len__(self) -> int:
        r"""Return total number of locus pieces in container.

        Returns:
            int: Number of fragment elements.
        """
        return len(self.ctg_indices)

    def __getitem__(self, item: int | slice | npt.NDArray[Any] | list[int]) -> "Any | LocusPieces":
        r"""Slice or array-mask all parallel fields of locus pieces simultaneously.

        Args:
            item (int | slice | npt.NDArray | list): Slice range, boolean mask, or index list.

        Returns:
            LocusPieces: A sliced [`LocusPieces`][kaptive.serotyping.models.LocusPieces] instance.

        Raises:
            NotImplementedError: If single integer key access is attempted.
        """
        if isinstance(item, (int, np.integer)):
            raise NotImplementedError("Single item access not implemented for LocusPieces")
        return LocusPieces(
            ctg_indices=self.ctg_indices[item],
            starts=self.starts[item],
            ends=self.ends[item],
            strands=self.strands[item],
        )

    @classmethod
    def concat(cls, batches: Iterable[Self]) -> Self:  # type: ignore
        r"""Concatenate multiple `LocusPieces` batches into a single container.

        Args:
            batches (Iterable[LocusPieces]): An iterable of
                [`LocusPieces`][kaptive.serotyping.models.LocusPieces] instances.

        Returns:
            LocusPieces: Combined [`LocusPieces`][kaptive.serotyping.models.LocusPieces] container.
        """
        batches_list = list(batches)
        if not batches_list:
            return cls.empty()  # type: ignore
        return cls(
            ctg_indices=np.concatenate([b.ctg_indices for b in batches_list]),
            starts=np.concatenate([b.starts for b in batches_list]),
            ends=np.concatenate([b.ends for b in batches_list]),
            strands=np.concatenate([b.strands for b in batches_list]),
        )

    @classmethod
    def empty(cls) -> "LocusPieces":
        r"""Create an empty `LocusPieces` container with zero-length arrays.

        Returns:
            LocusPieces: An empty [`LocusPieces`][kaptive.serotyping.models.LocusPieces] instance.
        """
        return cls(
            np.empty(0, dtype=np.uint32),
            np.empty(0, dtype=np.int32),
            np.empty(0, dtype=np.int32),
            np.empty(0, dtype=np.int8),
        )

    @classmethod
    def from_dict(cls, data: dict[str, Any]) -> "LocusPieces":
        r"""Reconstruct a `LocusPieces` container from a deserialized dictionary.

        Args:
            data (dict[str, Any]): Dictionary containing coordinate array lists.

        Returns:
            LocusPieces: Reconstructed [`LocusPieces`][kaptive.serotyping.models.LocusPieces] instance.
        """
        return cls(
            ctg_indices=np.array(data["ctg_indices"], dtype=np.uint32),
            starts=np.array(data["starts"], dtype=np.int32),
            ends=np.array(data["ends"], dtype=np.int32),
            strands=np.array(data["strands"], dtype=np.int8),
        )

    def to_dict(self) -> dict[str, Any]:
        r"""Convert array fields to a dictionary for JSON serialization.

        Returns:
            dict[str, Any]: Dictionary mapping field names to NumPy array data.
        """
        return {k: getattr(self, k) for k in ("ctg_indices", "starts", "ends", "strands")}


@dataclass(slots=True, frozen=True)
class SerotypingResult:
    r"""Efficient, immutable container representing an *in silico* serotyping call.

    Designed to be lightweight for JSON serialization and database storage while retaining full
    information needed to inspect and reconstruct alignment details. Houses nested SoA containers
    ([`LocusPieces`][kaptive.serotyping.models.LocusPieces] and [`GeneHits`][kaptive.serotyping.models.GeneHits])
    and sequence objects ([`Sequences`][kaptive.core.seq.Sequences]) for downstream processing.

    Attributes:
        kaptive_version (str): Version of Kaptive software that produced result.
        database_name (str): Name of target locus reference database.
        database_version (str): Version tag of reference database.
        database_organism (str): Target organism description in database.
        database_taxon (int): NCBI taxonomy ID of database.
        genome (str): Sample genome assembly identifier or filename.
        best_locus_idx (int): Index of best-matching locus in database.
        best_locus_name (str): Identifier name of best-matching locus.
        best_locus_score (float): Alignment score for best-matching locus.
        best_locus_completeness (float): Proportion of expected genes found in locus (0.0 to 1.0).
        locus_pieces (LocusPieces): Locus piece bounding coordinates container.
        length_discrepancy (float): Length discrepancy relative to reference locus.
        locus_seqs (Sequences): Sequences of identified locus region fragments.
        gene_hits (GeneHits): High-performance SoA container for gene alignment hits.
        gene_states (npt.NDArray[np.int8]): Gene classification state array matching
            [`GeneState`][kaptive.serotyping.models.GeneState] values.
        gene_seqs (Sequences): Extracted nucleotide sequences of locus genes.
        translations (Sequences): Translated amino acid sequences of locus genes.
        percent_identity (float): Overall nucleotide identity percentage across locus.
        percent_coverage (float): Overall reference coverage percentage.
        protein_identities (npt.NDArray[np.float32]): Per-gene protein identity percentages.
        phenotype (str): Inferred serotype phenotype description.
        typeable (bool): Flag indicating if confidence criteria for serotype call were met.
        missing_expected_genes (tuple[str, ...]): Identifiers of missing expected locus genes.
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
    problems: SerotypingProblem = field(init=False)

    def __post_init__(self) -> None:
        p = SerotypingProblem.NONE
        if len(self.locus_pieces) > 1:
            p |= SerotypingProblem.FRAGMENTED
        # Unexpected genes from other loci (not expected and not an allowed extra gene)
        if np.any(self.gene_hits.is_inside & ~self.gene_hits.is_expected & ~self.gene_hits.is_extra):
            p |= SerotypingProblem.UNEXPECTED_GENES
        # Missing genes: Overall completeness is not 100%, or expected genes are found but outside the locus boundary
        if self.best_locus_completeness < 1.0 or np.any(~self.gene_hits.is_inside & self.gene_hits.is_expected):
            p |= SerotypingProblem.MISSING_GENES
        # Novel genes: inside boundary but below identity threshold
        if np.any(self.gene_hits.is_inside & (self.gene_states == GeneState.NOVEL.value)):
            p |= SerotypingProblem.NOVEL_GENES
        # Truncated genes: inside boundary but partial or truncated
        if np.any(
            self.gene_hits.is_inside
            & ((self.gene_states == GeneState.TRUNCATED.value) | (self.gene_states == GeneState.PARTIAL.value))
        ):
            p |= SerotypingProblem.TRUNCATED_GENES

        object.__setattr__(self, "problems", p)

    @classmethod
    def from_dict(cls, data: dict[str, Any]) -> "SerotypingResult":
        r"""Reconstruct a `SerotypingResult` instance from a deserialized dictionary.

        Args:
            data (dict[str, Any]): Dictionary containing serialized fields and nested sub-dictionaries.

        Returns:
            SerotypingResult: Reconstructed [`SerotypingResult`][kaptive.serotyping.models.SerotypingResult] instance.
        """
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
            protein_identities=np.array(data["protein_identities"], dtype=np.float32),
        )

    def to_locus_data(self) -> "LocusData":
        r"""Convert result into a `LocusData` container for comparative multi-locus visualization.

        Extracts translations, locus backbone intervals, locus pieces, contig indices, gene states,
        and functional product descriptions for non-extra inside genes.

        Returns:
            LocusData: A [`LocusData`][kaptive.compare.LocusData] instance for multi-locus alignment plotting.
        """
        from kaptive.compare import LocusData

        mask = self.gene_hits.is_inside & ~self.gene_hits.is_extra
        descriptions = np.asarray(
            np.char.decode(self.gene_hits.product_descriptions[mask], "utf-8"),
            dtype=object,
        )

        return LocusData(
            proteins=self.translations[mask],  # type: ignore
            name=self.genome,
            backbone=self.gene_hits.t_intervals[mask],  # type: ignore
            pieces=self.locus_pieces,
            gene_ctg_indices=self.gene_hits.t_indices[mask],
            gene_states=self.gene_states[mask],
            gene_descriptions=descriptions,
        )

    def to_dict(self) -> dict[str, Any]:
        r"""Convert serotyping result into a dictionary suitable for JSON serialization.

        Returns:
            dict[str, Any]: Lightweight dictionary containing primitive types, lists, and nested dictionaries.
        """
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
            "problems": self.problems,
            "locus_pieces": self.locus_pieces.to_dict(),
            "gene_hits": self.gene_hits.to_dict(),
            "gene_states": self.gene_states,
            "protein_identities": self.protein_identities,
            "locus_seqs": self.locus_seqs.to_dict(),
            "gene_seqs": self.gene_seqs.to_dict(),
            "translations": self.translations.to_dict(),
        }
