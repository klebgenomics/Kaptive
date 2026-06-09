from pathlib import Path
from enum import IntEnum
from dataclasses import dataclass
from tempfile import NamedTemporaryFile
from threading import local as ThreadLocal
from concurrent.futures import ThreadPoolExecutor, as_completed

try:
    from os import process_cpu_count as cpu_count
except ImportError:
    from os import cpu_count

import numpy as np
import numpy.typing as npt
from numba import njit

from mappy import Aligner, ThreadBuffer

from kaptive.db import Database
from kaptive.core.alignment import Alignments
from kaptive.core.interval import Intervals
from kaptive.core.genome import GenomeAssembly
from kaptive.core.seq import Sequences, SeqRecord
from kaptive.core.pairwise import PairwiseAligner, PairwiseAlignments
from kaptive._version import __version__


# Enums ----------------------------------------------------------------------------------------------------------------
class GeneState(IntEnum):
    NORMAL = 0
    PARTIAL = 1
    TRUNCATED = 2
    NOVEL = 3


# Classes --------------------------------------------------------------------------------------------------------------
@dataclass(slots=True, frozen=True)
class GeneHits:
    """
    A high-performance SoA container for the final classified gene alignments.
    
    Encapsulating these arrays prevents `SerotypingResult` from being cluttered 
    with 11 parallel arrays and enables synchronized vectorized filtering and 
    dynamic property calculations.
    """
    gene_indices: npt.NDArray[np.int32]  # dtype=np.int32 (Global DB index)
    q_starts: npt.NDArray[np.int32]      # dtype=np.int32
    q_ends: npt.NDArray[np.int32]        # dtype=np.int32
    t_indices: npt.NDArray[np.uint32]    # dtype=np.uint32 (Contig index in GenomeAssembly)
    t_starts: npt.NDArray[np.int32]      # dtype=np.int32
    t_ends: npt.NDArray[np.int32]        # dtype=np.int32
    strands: npt.NDArray[np.int8]        # dtype=np.int8
    is_expected: npt.NDArray[np.bool_]   # dtype=bool
    is_inside: npt.NDArray[np.bool_]     # dtype=bool
    is_extra: npt.NDArray[np.bool_]      # dtype=bool

    @classmethod
    def empty(cls) -> 'GeneHits':
        """Creates an empty GeneHits container."""
        return cls(
            np.empty(0, dtype=np.int32), np.empty(0, dtype=np.int32),
            np.empty(0, dtype=np.int32), np.empty(0, dtype=np.uint32),
            np.empty(0, dtype=np.int32), np.empty(0, dtype=np.int32),
            np.empty(0, dtype=np.int8), np.empty(0, dtype=bool),
            np.empty(0, dtype=bool), np.empty(0, dtype=bool)
        )

    def __len__(self) -> int:
        return len(self.gene_indices)

    def __getitem__(self, item) -> 'GeneHits':
        """Allows vectorized slicing and boolean masking of all arrays simultaneously."""
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
            is_extra=self.is_extra[item]
        )

    @property
    def frames(self):
        return (3 - (self.t_starts % 3)) % 3
        
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


@dataclass(slots=True, frozen=True)
class SerotypingResult:
    kaptive_version: str
    database_name: str
    database_version: str
    genome: str
    best_locus_idx: int
    best_locus_name: str
    best_locus_score: float
    best_locus_completeness: float
    genes: GeneHits
    states: npt.NDArray[np.int8]
    locus_pieces: dict[str, Intervals]
    gene_seqs: Sequences
    protein_seqs: Sequences
    percent_identity: float
    percent_coverage: float
    protein_alignments: PairwiseAlignments
    phenotype: str


class Serotyper:
    """
    Performs _in silico_ serotyping on bacterial genome assemblies.
    """
    __slots__ = (
        '_db', '_executor', '_gene_aligner', '_thread_local', '_max_workers', '_protein_aligner', '_indexing_threads',
        '_gene_entropy_weights', '_expected_gene_counts'
    )
    def __init__(self, db: Database, max_workers: int | None = None, indexing_threads: int | None = None):
        self._db: Database = db
        self._max_workers = max_workers
        self._executor = None
        self._thread_local = ThreadLocal()
        self._protein_aligner = PairwiseAligner()
        self._indexing_threads = indexing_threads or cpu_count()
        cluster_counts = np.bincount(db.gene_cluster_ids)
        # safe_counts = np.maximum(cluster_counts, 1)  # Prevent division by zero
        # idf_weights = np.log2(len(db.loci) / safe_counts) + 1.0
        self._gene_entropy_weights: npt.NDArray[np.float64] = (1.0 / cluster_counts)[db.gene_cluster_ids]
        self._expected_gene_counts = np.array([len(i) for i in self._db.locus_intervals], dtype=np.uint8)
        # Initialise mappy index by writing the genes to a temporary fasta
        with NamedTemporaryFile(mode="wb", suffix=".fasta") as tmp:
            # We use the gene index for the header to act as array pointers later
            tmp.write(db.genes.to_fasta(use_indices=True))
            tmp_name = tmp.name
            self._gene_aligner = Aligner(fn_idx_in=tmp_name, n_threads=self._indexing_threads)
        Path(tmp_name).unlink(missing_ok=True)

    def __enter__(self):
        self._executor = ThreadPoolExecutor(max_workers=self._max_workers)
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        if self._executor is not None:
            self._executor.shutdown(wait=True)
            self._executor = None

    def __del__(self):
        if self._executor is not None:
            self._executor.shutdown(wait=False)

    @property
    def executor(self) -> ThreadPoolExecutor:
        if self._executor is None:
            self._executor = ThreadPoolExecutor(max_workers=self._max_workers)
        return self._executor

    def _process_contig(self, ctg: SeqRecord) -> Alignments | None:
        if not hasattr(self._thread_local, "buf"):
            self._thread_local.buf = ThreadBuffer()
        if gene_alns := list(self._gene_aligner.map(ctg.seq.decode('ascii'), buf=self._thread_local.buf)):
            return Alignments.from_mappy(ctg.id, len(ctg), gene_alns)
        return None

    def __call__(self, genome: GenomeAssembly | str | Path) -> SerotypingResult | None:
        genome = GenomeAssembly.ensure(genome)

        # Alignment phase ----------------------------------------------------------------------------------------------
        futures = [self.executor.submit(self._process_contig, i) for i in genome]
        if not (gene_alns := [i for future in as_completed(futures) if (i := future.result())]):
            return None
        else: # Concatenate all alignment-per-contig into a single batch for vectorization
            gene_alns = Alignments.concat(gene_alns).swap_sides()

        # Scoring phase ------------------------------------------------------------------------------------------------
        best_gene_alns = gene_alns.best()  # We take the best alignment per gene
        # Convert stringified FASTA headers from mappy directly back to global gene indices
        best_gene_aln_indices = best_gene_alns.q_names.astype(np.int32)
        best_gene_aln_scores = best_gene_alns.query_weighted_scores #* self._gene_entropy_weights[best_gene_aln_indices]
        # Filter out extra genes using the database's pre-computed boolean mask
        core_mask = ~self._db.extra_genes[best_gene_aln_indices]
        # Map each alignment to its corresponding locus index instantly
        best_aln_locus_indices = self._db.gene_locus_indices[best_gene_aln_indices]
        # Single-pass, zero-allocation accumulation of counts and core scores via Numba
        locus_counts, locus_core_scores = _accumulate_locus_metrics(
            best_aln_locus_indices, best_gene_aln_scores, core_mask, len(self._db.loci)
        )
        locus_completeness = locus_counts / self._expected_gene_counts
        locus_scores = locus_core_scores * locus_completeness
        best_locus_idx = np.argmax(locus_scores)
        best_locus_name = self._db.loci.ids[best_locus_idx]
        best_locus_completeness = locus_completeness[best_locus_idx]
        best_locus_score = locus_scores[best_locus_idx]

        # Reconstruction phase -----------------------------------------------------------------------------------------
        # Cull alignments, prioritizing genes belonging to the best match locus
        gene_indices = gene_alns.q_names.astype(np.int32)
        priority_mask = self._db.gene_locus_indices[gene_indices] == best_locus_idx
        culled_alns = gene_alns.cull_overlaps(by_query=False, priority_mask=priority_mask)
        # Re-extract arrays for the culled batch
        culled_gene_indices = culled_alns.q_names.astype(np.int32)
        t_indices = np.array([genome.id_map[n] for n in culled_alns.t_names], dtype=np.uint32)
        # Cluster intervals by contig using max_locus_length as tolerance
        culled_intervals = culled_alns.to_intervals(by_query=False)
        cluster_ids = culled_intervals.cluster_spatial(tolerance=self._db.max_locus_length, group_by=t_indices)
        # Identify expected genes and the clusters they fall into
        is_expected = (self._db.gene_locus_indices[culled_gene_indices] == best_locus_idx) & ~self._db.extra_genes[culled_gene_indices]
        valid_cluster_ids = np.unique(cluster_ids[is_expected])
        is_inside = np.isin(cluster_ids, valid_cluster_ids)
        is_extra = self._db.extra_genes[culled_gene_indices]

        # Construct bounding locus pieces per contig
        locus_pieces = {}
        for ctg_idx in np.unique(t_indices[is_inside]):
            ctg_mask = is_inside & (t_indices == ctg_idx)
            ctg_name = genome.contigs.ids[ctg_idx]
            # Use the heavily optimized SoA merge instead of manual looping and appending
            merged_intervals = culled_intervals[ctg_mask].sort().merge(tolerance=self._db.max_locus_length)
            # Force unstranded (0) for locus pieces as per original logic
            locus_pieces[ctg_name] = Intervals(merged_intervals.starts, merged_intervals.ends,
                                               np.zeros(len(merged_intervals), dtype=np.int8))
        genes = GeneHits(
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
        )

        # Gene state phase ---------------------------------------------------------------------------------------------
        # Extract gene nucleotides from their contigs
        gene_seqs = genome.contigs.extract_intervals(genes.t_indices, genes.t_intervals)
        # Translate nucleotides to amino acids, compensating for the reading frames of the alignments
        protein_seqs = gene_seqs.translate(frames=genes.frames)
        # Initialize states
        states = np.full(len(genes), GeneState.NORMAL.value, dtype=np.int8)
        is_partial = culled_alns.is_partial
        db_gene_lengths = self._db.genes.lengths[genes.gene_indices]
        # A partial gene colliding with a contig edge is excluded from being truncated
        is_truncated = (~is_partial) & (genes.query_lengths < (db_gene_lengths * 0.90))
        states[is_partial] = GeneState.PARTIAL.value
        states[is_truncated] = GeneState.TRUNCATED.value
        protein_alns = self._protein_aligner(protein_seqs, self._db.translations[genes.gene_indices])
        # Normal genes that fall below the identity threshold are considered NOVEL
        below_threshold = (states == GeneState.NORMAL.value) & (protein_alns.pidents < self._db.metadata.id_threshold)
        states[below_threshold] = GeneState.NOVEL.value

        valid_pidents = protein_alns.pidents[states == GeneState.NORMAL.value]
        percent_identity = float(np.mean(valid_pidents)) if valid_pidents.size > 0 else 0.0
        percent_coverage = float(best_locus_completeness * 100.0)

        return SerotypingResult(
            kaptive_version=__version__,
            database_name=self._db.metadata.name,
            database_version=self._db.metadata.version,
            genome=genome.id,
            best_locus_idx=best_locus_idx,
            best_locus_name=best_locus_name,
            best_locus_score=best_locus_score,
            best_locus_completeness=best_locus_completeness,
            genes=genes,
            states=states,
            locus_pieces=locus_pieces,
            gene_seqs=gene_seqs,
            protein_seqs=protein_seqs,
            percent_identity=percent_identity,
            percent_coverage=percent_coverage,
            protein_alignments=protein_alns,
            phenotype="Unknown"
        )


# Kernels --------------------------------------------------------------------------------------------------------------
@njit(cache=True, nogil=True)
def _accumulate_locus_metrics(locus_indices: npt.NDArray[np.int32], scores: npt.NDArray[np.float64],
                              core_mask: npt.NDArray[np.bool_], n_loci: int) -> tuple[npt.NDArray[np.int32], npt.NDArray[np.float64]]:
    """Single-pass fused kernel to calculate locus counts and core scores simultaneously without memory allocation."""
    counts = np.zeros(n_loci, dtype=np.int32)
    core_sums = np.zeros(n_loci, dtype=np.float64)

    for i in range(len(locus_indices)):
        locus_idx = locus_indices[i]
        counts[locus_idx] += 1
        if core_mask[i]:
            core_sums[locus_idx] += scores[i]

    return counts, core_sums