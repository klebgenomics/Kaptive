from pathlib import Path
from enum import IntEnum, IntFlag, auto
from dataclasses import dataclass, replace
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
from kaptive.core.seq import Sequences
from kaptive.core.pairwise import PairwiseAligner, PairwiseAlignments


# Enums ----------------------------------------------------------------------------------------------------------------
class SerotypingProblem(IntFlag):
    NONE = 0
    FRAGMENTED = auto()
    EXTRA_GENES = auto()
    MISSING_GENES = auto()
    NOVEL_GENES = auto()
    TRUNCATED_GENES = auto()


class GeneState(IntEnum):
    NORMAL = 0
    PARTIAL = 1
    TRUNCATED = 2
    BELOW_THRESHOLD = 3


# Classes --------------------------------------------------------------------------------------------------------------
@dataclass(slots=True, frozen=True)
class GeneHitBatch:
    """
    SoA container for the final classified gene alignments.
    Replaces the legacy object-oriented lists.
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
    states: npt.NDArray[np.int8]         # dtype=np.int8

    def __len__(self) -> int:
        return len(self.gene_indices)

    def extract_sequences(self, genome: 'GenomeAssembly', db: Database) -> tuple[Sequences, npt.NDArray[np.int32], 'GeneHitBatch']:
        """
        Extracts the full expected gene sequences by projecting the local alignment 
        boundaries outwards to capture the start and stop codons.
        Returns the sequence batch, the corrected reading frames, and a coordinate-updated batch.
        """
        q_lengths = db.genes.lengths[self.gene_indices].astype(np.int32)
        ctg_lengths = genome.contigs.lengths[self.t_indices].astype(np.int32)

        proj_left = np.where(self.strands == 1, self.q_starts, q_lengths - self.q_ends)
        proj_right = np.where(self.strands == 1, q_lengths - self.q_ends, self.q_starts)

        # Centralize coordinate projection and clipping inside Intervals
        intervals = Intervals(self.t_starts, self.t_ends, self.strands, self.gene_indices)
        intervals = intervals.expand(proj_left, proj_right, clip_lengths=ctg_lengths)

        actual_proj_left = self.t_starts - intervals.starts
        actual_proj_right = intervals.ends - self.t_ends

        # Calculate the new query-relative start positions to determine the reading frame
        new_q_starts = np.where(self.strands == 1, 
                                self.q_starts - actual_proj_left, 
                                self.q_starts - actual_proj_right)
                                
        new_q_ends = np.where(self.strands == 1,
                              self.q_ends + actual_proj_right,
                              self.q_ends + actual_proj_left)

        frames = (3 - (new_q_starts % 3)) % 3

        ids = tuple(db.genes.ids[x] for x in self.gene_indices)
        seqs = genome.contigs.extract_intervals(self.t_indices, intervals, new_ids=ids)
        
        updated_batch = replace(
            self,
            q_starts=new_q_starts.astype(np.int32),
            q_ends=new_q_ends.astype(np.int32),
            t_starts=intervals.starts,
            t_ends=intervals.ends
        )
        return seqs, frames.astype(np.int32), updated_batch


@dataclass(slots=True, frozen=True)
class SerotyperResult:
    database_name: str
    database_version: str
    genome: str
    best_match_locus: str
    best_match_score: float
    best_match_zscore: float
    genes: GeneHitBatch
    locus_boundaries: dict[str, Intervals]
    metrics: 'LocusScoreBatch'
    gene_seqs: Sequences
    protein_seqs: Sequences
    percent_identity: float
    percent_coverage: float
    protein_alignments: PairwiseAlignments
    phenotype: str
    problems: SerotypingProblem
    is_typeable: bool


@dataclass(frozen=True, slots=True)
class LocusScoreBatch:
    """
    SoA container for locus scoring metrics.
    Allows vectorised operations and mixed dtypes while keeping a clean namespace.
    """
    scores: npt.NDArray[np.float64]           # dtype=np.float64
    z_scores: npt.NDArray[np.float64]         # dtype=np.float64
    completeness: npt.NDArray[np.float64]     # dtype=np.float64
    synteny_weights: npt.NDArray[np.float64]  # dtype=np.float64
    contig_counts: npt.NDArray[np.int32]    # dtype=np.int32

    @classmethod
    def empty(cls, n_loci: int) -> 'LocusScoreBatch':
        return cls(
            scores=np.zeros(n_loci, dtype=np.float64),
            z_scores=np.zeros(n_loci, dtype=np.float64),
            completeness=np.zeros(n_loci, dtype=np.float64),
            synteny_weights=np.zeros(n_loci, dtype=np.float64),
            contig_counts=np.zeros(n_loci, dtype=np.int32)
        )


@dataclass(slots=True, frozen=True)
class ConfidenceEvaluator:
    """
    Evaluates the final SerotyperResult to determine if the match is trustworthy.
    """
    max_other_genes: int = 1
    min_completeness: float = 0.5  # Replaced prop_expected_genes for SoA clarity
    allow_below_threshold: bool = False
    min_zscore: float | None = None
    
    def __call__(self, result: SerotyperResult) -> bool:
        genes = result.genes
        best_idx = int(np.argmax(result.metrics.scores))
        
        # 1. Completeness (Proportion of expected genes found)
        if result.metrics.completeness[best_idx] < self.min_completeness:
            return False
            
        # 2. Novel genes (physically inside the locus that fall below the pident threshold)
        novel_count = np.sum(genes.is_inside & (genes.states == GeneState.BELOW_THRESHOLD))
        if novel_count > self.max_other_genes:
            return False
            
        # 3. Below threshold expected genes
        if not self.allow_below_threshold and np.any(genes.is_expected & (genes.states == GeneState.BELOW_THRESHOLD)):
            return False
                
        # 4. Z-score check (Safeguarded against small databases where max possible Z < min_zscore)
        if self.min_zscore is not None and len(result.metrics.scores) > 2 and result.best_match_zscore < self.min_zscore:
            return False
            
        return True


class Serotyper:
    """
    Performs _in silico_ serotyping on bacterial genome assemblies.
    """
    __slots__ = (
        '_db', '_executor', '_gene_aligner', '_thread_local',
        '_max_workers', '_protein_aligner', '_indexing_threads', '_confidence_evaluator', '_gene_weights'
    )
    def __init__(self, db: Database, confidence_evaluator: ConfidenceEvaluator | None = None,
                 max_workers: int | None = None, indexing_threads: int | None = None):
        self._db: Database = db
        self._max_workers = max_workers
        self._executor = None
        self._thread_local = ThreadLocal()
        self._protein_aligner = PairwiseAligner()
        self._confidence_evaluator = confidence_evaluator or ConfidenceEvaluator()
        
        # Calculate custom weights for locus-defining genes
        self._gene_weights = np.ones(len(self._db.genes), dtype=np.float64)
        for i, name in enumerate(self._db.genes.ids):
            if 'wzx' in name or 'wzy' in name:
                self._gene_weights[i] = 100.0
        # Initialise mappy index by writing the genes to a temporary fasta
        with NamedTemporaryFile(mode="wb", suffix=".fasta", delete=False) as tmp:
            # We use the gene index for the header to act as array pointers later
            tmp.write(b''.join(b">%i\n%b\n" % (i, rec.seq) for i, rec in enumerate(db.genes)))
            tmp_name = tmp.name
        # Pass the closed file to mappy, then delete it immediately
        self._gene_aligner = Aligner(fn_idx_in=tmp_name, preset="sr", best_n=10000, n_threads=indexing_threads or cpu_count())
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

    def _process_contig_chunk(self, ctg_id: str, ctg_len: int, seq: str, offset: int) -> Alignments | None:
        if not hasattr(self._thread_local, "buf"):
            self._thread_local.buf = ThreadBuffer()

        if alns := list(self._gene_aligner.map(seq, buf=self._thread_local.buf)):
            # The length of the query string given to from_mappy is ctg_len
            batch = Alignments.from_mappy(ctg_id, ctg_len, alns)
            return batch.shift_query(offset)

        return None

    def __call__(self, genome: GenomeAssembly | str | Path) -> SerotyperResult | None:
        if not isinstance(genome, GenomeAssembly):
            genome = GenomeAssembly.from_file(genome)

        # Align contigs to genes in chunks to massively parallelize large chromosomes
        CHUNK_SIZE = 250_000
        OVERLAP = 50_000
        futures = []
        for ctg in genome:
            seq_str = ctg.seq.decode('ascii')
            ctg_len = len(ctg.seq)
            
            if ctg_len <= CHUNK_SIZE:
                futures.append(self.executor.submit(self._process_contig_chunk, ctg.id, ctg_len, seq_str, 0))
            else:
                for i in range(0, ctg_len, CHUNK_SIZE - OVERLAP):
                    chunk = seq_str[i:i + CHUNK_SIZE]
                    futures.append(self.executor.submit(self._process_contig_chunk, ctg.id, ctg_len, chunk, i))

        if not (alignments := [i for future in as_completed(futures) if (i := future.result())]):
            return None  # Comment out during testing

        # Concatenate all alignment-per-contig into a single batch for vectorization
        alignments = Alignments.concat(alignments).swap_sides()

        # TODO: Implement graph stitching logic
        # Use available graph for resolving alignments split over contigs
        # if genome.edges:  # If the user provided a GFA, the genome will have edges.
        #     graph = genome.as_graph()
        #     alignments = graph.stitch_alignments(alignments)

        # Convert stringified FASTA headers from mappy directly back to global gene indices
        # Map each alignment to its corresponding locus index instantly
        aln_locus_indices = self._db.gene_locus_indices[alignments.q_names.astype(np.int32)]

        # Init score arrays
        metrics = LocusScoreBatch.empty(len(self._db.loci))
        match_bonuses = np.zeros(len(self._db.loci), dtype=np.float64)

        # Cull overlaps globally per-locus
        locus_hits = alignments.cull_overlaps(max_overlap_fraction=0.1, group_by=aln_locus_indices)
        
        if len(locus_hits) > 0:
            culled_gene_indices = locus_hits.q_names.astype(np.int32)
            culled_locus_indices = self._db.gene_locus_indices[culled_gene_indices]
            
            # Completeness calculation
            global_gene_scores = np.maximum(0, np.bincount(culled_gene_indices, weights=locus_hits.scores, minlength=len(self._db.genes)))
            global_gene_scores = global_gene_scores * self._gene_weights
            gene_coverage = np.clip(global_gene_scores / np.maximum(1, self._db.genes.lengths), 0.0, 1.0)
            
            _, contig_indices = np.unique(locus_hits.t_names, return_inverse=True)
            expected_pos_all = self._db.gene_positions[culled_gene_indices]
            
            sort_order = np.lexsort((locus_hits.t_starts, contig_indices, culled_locus_indices))
            locus_starts = np.array([s.start for s in self._db.locus_gene_slices], dtype=np.int32)
            locus_ends = np.array([s.stop for s in self._db.locus_gene_slices], dtype=np.int32)
            
            _metric_kernel(
                locus_starts, locus_ends, self._db.extra_genes,
                gene_coverage, global_gene_scores,
                culled_locus_indices[sort_order],
                contig_indices[sort_order].astype(np.int32),
                expected_pos_all[sort_order],
                len(self._db.loci),
                metrics.completeness, match_bonuses, metrics.synteny_weights, metrics.contig_counts
            )
            
        metrics.scores[:] = (metrics.completeness * metrics.synteny_weights) + match_bonuses

        # Calculate Z-scores vectorially, guarding against divide-by-zero
        std_score = np.std(metrics.scores)
        if std_score > 0:
            metrics.z_scores[:] = (metrics.scores - np.mean(metrics.scores)) / std_score
        else:
            metrics.z_scores[:] = 0.0

        # Locus Reconstruction & Classification
        best_locus_idx = int(np.argmax(metrics.scores))
        best_match_name = self._db.loci.ids[best_locus_idx]

        # Isolate and define the boundaries of the best match
        best_mask = (culled_locus_indices == best_locus_idx) if len(locus_hits) > 0 else np.zeros(0, dtype=bool)
        best_hits = locus_hits.filter(best_mask)

        locus_boundaries: dict[str, Intervals] = {}
        for ctg, batch in best_hits.split(by_query=False):
            itv = batch.expand().to_intervals(by_query=False).sort().merge(tolerance=self._db.max_locus_length)
            locus_boundaries[ctg] = itv

        # Calculate locus coverage as the physical length of the locus boundaries on the contigs
        # divided by the expected total length of the reference locus
        total_boundary_length = sum(np.sum(batch.lengths) for batch in locus_boundaries.values())
        expected_locus_length = self._db.loci.lengths[best_locus_idx]
        real_coverage = min(100.0, (total_boundary_length / expected_locus_length) * 100.0) if expected_locus_length > 0 else 0.0

        # Global Culling
        # Prioritize genes that are part of the best locus, or are known extra genes.
        # This prevents a slightly longer allele from another locus from displacing the correct locus's gene.
        tmp_gene_indices = alignments.q_names.astype(np.int32)
        tmp_locus_indices = self._db.gene_locus_indices[tmp_gene_indices]
        tmp_is_extra = self._db.extra_genes[tmp_gene_indices]
        priority_mask = (tmp_locus_indices == best_locus_idx) | tmp_is_extra

        global_culled = alignments.cull_overlaps(max_overlap_fraction=0.1, priority_mask=priority_mask)
        culled_gene_indices = global_culled.q_names.astype(np.int32)
        culled_locus_indices = self._db.gene_locus_indices[culled_gene_indices]

        # Vectorized Classification
        is_expected = (culled_locus_indices == best_locus_idx)
        is_inside = np.zeros(len(global_culled), dtype=bool)
        
        # Convert once for spatial queries
        global_itv = global_culled.to_intervals(by_query=False)

        for ctg, boundary_batch in locus_boundaries.items():
            if not np.any(ctg_mask := (global_culled.t_names == ctg)):
                continue

            # Query overlaps effortlessly via Intervals!
            ctg_itv = global_itv[ctg_mask]
            is_inside[ctg_mask] = ctg_itv.overlaps_with(boundary_batch)

        culled_t_indices = np.array([genome.id_map[n] for n in global_culled.t_names], dtype=np.uint32)

        genes = GeneHitBatch(
            gene_indices=culled_gene_indices,
            q_starts=global_culled.q_starts,
            q_ends=global_culled.q_ends,
            t_indices=culled_t_indices,
            t_starts=global_culled.t_starts,
            t_ends=global_culled.t_ends,
            strands=global_culled.strands,
            is_expected=is_expected,
            is_inside=is_inside,
            is_extra=self._db.extra_genes[culled_gene_indices],
            states=np.full(len(culled_gene_indices), GeneState.NORMAL, dtype=np.int8)
        )

        gene_seqs, frames, genes = genes.extract_sequences(genome, self._db)
        protein_seqs = gene_seqs.translate(frames)
        protein_alignments = self._protein_aligner(protein_seqs, self._db.translations[genes.gene_indices])

        # Calculate mutually exclusive states
        is_below_threshold = protein_alignments.pidents < self._db.metadata.id_threshold
        is_truncated = protein_seqs.internal_stops

        # Check if the expanded alignment hangs off the edge of the contig
        ctg_lengths = genome.contigs.lengths[genes.t_indices]
        q_lengths = self._db.genes.lengths[genes.gene_indices]
        is_partial = ((genes.t_starts == 0) & (genes.q_starts > 0)) | \
                     ((genes.t_ends == ctg_lengths) & (genes.q_ends < q_lengths))

        # Apply states hierarchically. 
        # As you noted: if a gene is partial, it is impossible to reliably assess if it is 
        # truly truncated or below threshold. PARTIAL safely overrides them!
        genes.states[is_below_threshold] = GeneState.BELOW_THRESHOLD
        genes.states[is_truncated] = GeneState.TRUNCATED
        genes.states[is_partial] = GeneState.PARTIAL
        
        # Exclude partial genes from mean_identity so assembly fragmentation doesn't tank the score
        identity_mask = genes.is_expected & genes.is_inside & (genes.states != GeneState.PARTIAL)
        if not np.any(identity_mask):
            # Fallback to including them if literally every single expected gene is partial
            identity_mask = genes.is_expected & genes.is_inside
            
        if np.any(identity_mask):
            mean_identity = float(np.mean(protein_alignments.pidents[identity_mask]))
        else:
            mean_identity = 0.0

        # Evaluate Phenotypes
        # 1. Identify which gene clusters are present and "active"
        # Expected genes can be NORMAL or PARTIAL.
        # Extra genes MUST be strictly NORMAL to be applied to logic.
        active_expected = genes.is_expected & (genes.states < GeneState.TRUNCATED)
        active_extra = genes.is_extra & (genes.states == GeneState.NORMAL)
        active_mask = active_expected | active_extra
        
        active_cluster_ids = np.unique(self._db.gene_cluster_ids[genes.gene_indices[active_mask]])
        active_clusters = {self._db.cluster_keys[cid] for cid in active_cluster_ids}

        # 2. Find matching phenotypes
        valid_phenotypes = []
        for pheno in self._db.phenotypes:
            if pheno.loci and best_match_name not in pheno.loci:
                continue
            if not pheno.extra_genes.issubset(active_clusters):
                continue
            if not pheno.inactive_genes.isdisjoint(active_clusters):
                continue
            valid_phenotypes.append(pheno)

        # 3. Resolve final phenotype
        if valid_phenotypes:
            # Sort by priority (lower number = higher priority) then alphabetically by ID
            valid_phenotypes.sort(key=lambda p: (p.priority, p.id))
            final_phenotype = valid_phenotypes[0].id
        else:
            # Fallback to the locus's default serotype, or the locus name itself if missing
            final_phenotype = self._db.serotypes[best_locus_idx] or best_match_name
            
        # 4. Populate Serotyping Problems
        problems = SerotypingProblem.NONE
        if metrics.contig_counts[best_locus_idx] > 1:
            problems |= SerotypingProblem.FRAGMENTED
            
        if np.any(genes.is_expected & (genes.states == GeneState.TRUNCATED)):
            problems |= SerotypingProblem.TRUNCATED_GENES
            
        if np.any(active_extra):
            problems |= SerotypingProblem.EXTRA_GENES
            
        # Missing genes: Are any expected clusters absent from the active expected hits?
        locus_slice = self._db.locus_gene_slices[best_locus_idx]
        expected_clusters = set(self._db.gene_cluster_ids[locus_slice])
        found_expected = set(self._db.gene_cluster_ids[genes.gene_indices[active_expected]])
        if not expected_clusters.issubset(found_expected):
            problems |= SerotypingProblem.MISSING_GENES
            
        # Novel genes: Anything physically inside the locus that falls below the pident threshold
        novel_mask = genes.is_inside & (genes.states == GeneState.BELOW_THRESHOLD)
        if np.any(novel_mask):
            problems |= SerotypingProblem.NOVEL_GENES

        result = SerotyperResult(
            genome=genome.id,
            database_name=self._db.metadata.name,
            database_version=self._db.metadata.version,
            best_match_locus=best_match_name,
            best_match_score=metrics.scores.max(),
            best_match_zscore=float(metrics.z_scores[best_locus_idx]),
            genes=genes,
            locus_boundaries=locus_boundaries,
            metrics=metrics,
            gene_seqs=gene_seqs,
            protein_seqs=protein_seqs,
            percent_identity=mean_identity,
            percent_coverage=real_coverage,
            protein_alignments=protein_alignments,
            phenotype=final_phenotype,
            problems=problems,
            is_typeable=False  # Placeholder
        )
        
        is_typeable = self._confidence_evaluator(result)
        return replace(result, is_typeable=is_typeable)


# Kernels --------------------------------------------------------------------------------------------------------------
@njit(cache=True, nogil=True)
def _lis_kernel(arr: np.ndarray) -> int:
    """
    Finds the length of the longest strictly increasing biological sequence.
    Highly optimized O(N log N) implementation using Numba and pre-allocation.
    """
    n = len(arr)
    if n == 0:
        return 0

    # Pre-allocate array to avoid dynamic resizing overhead
    tails = np.empty(n, dtype=arr.dtype)
    tails[0] = arr[0]
    length = 1

    for i in range(1, n):
        x = arr[i]
        if x > tails[length - 1]:
            tails[length] = x
            length += 1
        else:
            idx = np.searchsorted(tails[:length], x)
            tails[idx] = x

    return length


@njit(cache=True, nogil=True)
def _metric_kernel(
    locus_starts: np.ndarray, locus_ends: np.ndarray, extra_genes: np.ndarray,
    gene_coverage: np.ndarray, global_gene_scores: np.ndarray,
    locus_indices: np.ndarray, contig_indices: np.ndarray, expected_pos: np.ndarray,
    n_loci: int, out_completeness: np.ndarray, out_bonuses: np.ndarray, 
    out_synteny: np.ndarray, out_contigs: np.ndarray
):
    # 1. Completeness and Bonuses
    for i in range(n_loci):
        s = locus_starts[i]
        e = locus_ends[i]
        
        if e <= s or extra_genes[s]:
            continue
            
        total_expected = e - s
        if total_expected > 0:
            out_completeness[i] = np.sum(gene_coverage[s:e]) / total_expected
        out_bonuses[i] = np.sum(global_gene_scores[s:e]) / 10_000_000.0

    # 2. Synteny and Contig counts
    n_hits = len(locus_indices)
    if n_hits == 0:
        return

    curr_locus = locus_indices[0]
    curr_contig = contig_indices[0]
    ctg_start_idx = 0
    in_order_count = 0
    unique_contigs = 1
    total_locus_hits = 0

    for i in range(1, n_hits + 1):
        if i == n_hits or locus_indices[i] != curr_locus or contig_indices[i] != curr_contig:
            arr = expected_pos[ctg_start_idx:i]
            
            lis = _lis_kernel(arr)
            lds = _lis_kernel(arr[::-1])
            in_order_count += max(lis, lds)
            total_locus_hits += len(arr)
            
            if i < n_hits:
                if locus_indices[i] != curr_locus:
                    if total_locus_hits > 0:
                        out_synteny[curr_locus] = in_order_count / total_locus_hits
                    out_contigs[curr_locus] = unique_contigs
                    
                    curr_locus = locus_indices[i]
                    in_order_count = 0
                    unique_contigs = 1
                    total_locus_hits = 0
                else:
                    unique_contigs += 1
                
                curr_contig = contig_indices[i]
                ctg_start_idx = i

    if total_locus_hits > 0:
        out_synteny[curr_locus] = in_order_count / total_locus_hits
    out_contigs[curr_locus] = unique_contigs