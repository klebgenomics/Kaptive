from pathlib import Path
from concurrent.futures import ThreadPoolExecutor, as_completed
from dataclasses import dataclass
from enum import IntEnum

try:
    from os import process_cpu_count as cpu_count
except ImportError:
    from os import cpu_count

import numpy as np
import numpy.typing as npt

from pyfgs import GeneFinder, Model

from kaptive.db import Database
from kaptive.core.pairwise import PairwiseAligner, RandstrobeIndex
from kaptive.core.genome import GenomeAssembly
from kaptive.core.seq import Sequences
from kaptive.core.interval import Intervals


# Classes --------------------------------------------------------------------------------------------------------------
@dataclass(frozen=True, slots=True)
class GeneClusters:
    """Structure-of-Arrays container for annotated genes grouped into spatial clusters.

    This batch holds information about genes (ORFs) from the query genome that have successfully
    matched against the database and passed size filtering, grouped into distinct spatial clusters
    on their parent contigs.

    Attributes:
        cluster_ids (npt.NDArray[np.int32]): An array where each element identifies the spatial
            cluster to which the corresponding gene belongs.
        query_indices (npt.NDArray[np.uint32]): An array of indices identifying the original ORF
            in the global `query_proteins` batch.
        target_indices (npt.NDArray[np.uint32]): An array of indices identifying the best-matching
            target gene in the database (`target_proteins` batch).
        contig_indices (npt.NDArray[np.int32]): An array of integer indices mapping back to the
            original contig string IDs.
        intervals (Intervals): An `Intervals` representing the physical coordinates of
            each gene on its parent contig.
        seed_scores (npt.NDArray[np.uint32]): An array of randstrobe intersection counts
             that seeded the match.
        seed_coverage (npt.NDArray[np.float64]): An array of dynamic proportional scores representing
            the fraction of the expected database randstrobes found in the query ORF.
        pidents (npt.NDArray[np.float64]): An array of percent identity values from the pairwise
            alignment between the query ORF and the target database gene.
    """
    spatial_ids: npt.NDArray[np.int32]
    query_indices: npt.NDArray[np.uint32]
    target_indices: npt.NDArray[np.uint32]
    contig_indices: npt.NDArray[np.int32]
    intervals: Intervals
    seed_scores: npt.NDArray[np.uint32]
    seed_coverage: npt.NDArray[np.float64]
    pidents: npt.NDArray[np.float64]

    def __len__(self) -> int:
        """Returns the total number of annotated genes in this batch."""
        return len(self.spatial_ids)

    def filter(self, mask: np.ndarray) -> 'GeneClusters':
        """Returns a new batch containing only records where the mask is True."""
        return GeneClusters(
            spatial_ids=self.spatial_ids[mask],
            query_indices=self.query_indices[mask],
            target_indices=self.target_indices[mask],
            contig_indices=self.contig_indices[mask],
            intervals=self.intervals[mask],
            seed_scores=self.seed_scores[mask],
            seed_coverage=self.seed_coverage[mask],
            pidents=self.pidents[mask]
        )


@dataclass(frozen=True, slots=True)
class Homologs:
    """Structure-of-Arrays container for high-scoring homologs found outside the primary locus."""
    query_indices: npt.NDArray[np.uint32]
    target_indices: npt.NDArray[np.uint32]
    contig_indices: npt.NDArray[np.int32]
    intervals: Intervals
    seed_scores: npt.NDArray[np.uint32]
    seed_coverage: npt.NDArray[np.float64]

    def __len__(self) -> int:
        return len(self.query_indices)

    @classmethod
    def empty(cls) -> 'Homologs':
        """Creates an empty Homologs collection."""
        return cls(
            np.empty(0, dtype=np.uint32),
            np.empty(0, dtype=np.uint32),
            np.empty(0, dtype=np.int32),
            Intervals.empty(),
            np.empty(0, dtype=np.uint32),
            np.empty(0, dtype=np.float64)
        )


@dataclass(frozen=True, slots=True)
class LocusComparisons:
    """Structure-of-Arrays container for vectorized locus comparisons."""
    statuses: npt.NDArray[np.int8]
    intersect_genes: npt.NDArray[np.float64]
    completeness: npt.NDArray[np.float64]
    allele_scores: npt.NDArray[np.float64]
    missing_genes: npt.NDArray[np.int32]
    added_genes: npt.NDArray[np.int32]
    
    def __len__(self) -> int:
        return len(self.statuses)

    @property
    def best_locus_idx(self) -> int:
        if len(self) == 0:
            return -1
        # Primary key: status, Secondary: allele_scores, Tertiary: intersect_genes, Quaternary: completeness
        order = np.lexsort((self.completeness, self.intersect_genes, self.allele_scores, self.statuses))
        return int(order[-1])


@dataclass(frozen=True, slots=True)
class AnnotationResult:
    """A container for the final results of annotating a genome against a database.

    Attributes:
        genome_id (str): The identifier of the annotated genome assembly.
        clusters (GeneClusters): The detailed, clustered gene matches found in the genome.
        match_status (LocusStatus): An enum indicating the overall quality and type of the best locus match.
        match_locus (str): The identifier of the reference locus from the database that best matches
            the findings in the query genome.
        metrics (LocusComparisons): Vectorized comparison metrics for the final match.
        novel_orfs (Sequences): A batch of unannotated ORFs physically located within the bounds of the locus.
    """
    genome_id: str
    clusters: GeneClusters
    match_status: 'LocusStatus'
    match_locus: str
    metrics: LocusComparisons
    novel_orfs: Sequences
    extra_homologs: Homologs
    
    
class LocusStatus(IntEnum):
    """Enumeration of possible match states when comparing a query locus against a reference locus.

    Values are ordered roughly by confidence/completeness, from lowest (NOVEL) to highest
    (GENE_ORDER_MATCH).

    Attributes:
        NOVEL (0): The found genes do not form a recognizable match to any reference locus.
        GENE_SET_MATCH (1): Contains exactly the expected set of core genes for a reference locus,
            but they may be out of order, or the locus is split across multiple contigs.
        GENE_ORDER_MATCH (2): A perfect match: contains exactly the expected set of core genes,
            in the correct relative order (or perfectly reversed), and contained on a single contig.
    """
    NOVEL = 0
    GENE_SET_MATCH = 1
    GENE_ORDER_MATCH = 2


class Annotator:
    """The main orchestration engine for identifying and classifying antigen loci in a genome.

    This class manages the pipeline of predicting ORFs, searching them against a specialized
    database using a fast `RandstrobeIndex`, performing detailed pairwise alignments to confirm
    hits, clustering valid hits physically on contigs, and ultimately determining the `LocusStatus`
    against known reference loci.

    It utilizes a `ThreadPoolExecutor` for parallelizing the initial ORF prediction step across
    contigs.
    """
    __slots__ = ('_db', '_max_workers', '_executor', '_translations', '_orf_tolerance', '_enforce_strand',
                 '_min_cluster_size', '_min_strobe_coverage', '_db_strobe_counts', '_gene_finder', '_index', 
                 '_aligner', '_expected_matrix', '_expected_lists', '_extra_cluster_mask')
    _extra_cluster_mask: npt.NDArray[np.bool_]

    def __init__(
            self, db: Database, orf_tolerance: int = 4, enforce_strand: bool = False,
            min_cluster_size: int = 1, k: int=6, s: int=3, w_max: int=4, 
            min_strobe_coverage: float = 0.02, max_workers: int | None = None,
    ):
        """Initializes the Annotator with a target database and configuration parameters.

        Args:
            db (Database): The pre-compiled Kaptive database object.
            orf_tolerance (int, optional): The maximum number of unannotated ORFs allowed
                between two valid gene hits on the same contig to group them into the same
                cluster. Defaults to 4.
            enforce_strand (bool, optional): If True, hits must be on the same strand to
                be grouped into a cluster. Defaults to False.
            min_cluster_size (int, optional): The minimum number of genes a cluster must contain
                to be retained for locus evaluation. Defaults to 3.
            k (int, optional): The k-mer size used for building the target `RandstrobeIndex`. Defaults to 6.
            s (int, optional): The s-mer (syncmer) size used for the `RandstrobeIndex`. Defaults to 3.
            w_max (int, optional): The maximum window size for the second strobe in the `RandstrobeIndex`. Defaults to 4.
            min_strobe_coverage (float, optional): The minimum proportion of a database gene's expected randstrobes
                that must be found to confirm a hit. Defaults to 0.02 (2%).
            max_workers (int | None, optional): The number of threads to use for ORF prediction. If None,
                uses the system CPU count. Defaults to None.
        """
        self._db: Database = db
        self._max_workers = max_workers
        self._executor = None
        self._translations = db.translations
        self._orf_tolerance = orf_tolerance
        self._enforce_strand = enforce_strand
        self._min_cluster_size = min_cluster_size
        self._min_strobe_coverage = min_strobe_coverage
        self._gene_finder = GeneFinder(Model.Complete, whole_genome=False)
        self._index = RandstrobeIndex.build(self._translations, k=k, s=s, w_max=w_max, sort_by_hash=False)
        
        # Pre-compute the exact maximum number of randstrobes for every database gene to allow dynamic normalization
        self._db_strobe_counts = np.maximum(1, np.bincount(self._index.seq_indices, minlength=len(self._translations)))
        self._aligner = PairwiseAligner()
        
        # Pre-format database arrays for vectorized locus comparisons
        n_loci = len(self._db.loci)
        n_clusters = len(self._db.cluster_keys)
        
        self._expected_matrix = np.zeros((n_loci, n_clusters), dtype=bool)
        self._expected_lists = []
        
        # Mask of universally recognized "extra" genes (so they don't penalize core set matches)
        self._extra_cluster_mask = np.zeros(n_clusters, dtype=bool)
        self._extra_cluster_mask[self._db.gene_cluster_ids[self._db.extra_genes]] = True
        
        for i in range(n_loci):
            locus_slice = self._db.locus_gene_slices[i]
            core_mask = ~self._db.extra_genes[locus_slice]
            expected_core = self._db.gene_cluster_ids[locus_slice][core_mask]
            
            self._expected_matrix[i, expected_core] = True
            self._expected_lists.append(expected_core.tolist())

    def __enter__(self):
        """Context manager entry point. Initializes the ThreadPoolExecutor."""
        self._executor = ThreadPoolExecutor(max_workers=self._max_workers)
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        """Context manager exit point. Ensures the ThreadPoolExecutor is properly shut down."""
        if self._executor is not None:
            self._executor.shutdown(wait=True)
            self._executor = None

    def __del__(self):
        """Destructor to ensure executor shutdown if not used as a context manager."""
        if self._executor is not None:
            self._executor.shutdown(wait=False)

    @property
    def executor(self) -> ThreadPoolExecutor:
        """Provides access to the ThreadPoolExecutor, initializing it lazily if necessary."""
        if self._executor is None:
            self._executor = ThreadPoolExecutor(max_workers=self._max_workers)
        return self._executor

    def _predict_orfs(self, contig_idx: int, contig_seq: bytes) -> tuple[int, Intervals, list[bytes]] | None:
        """ThreadPool helper for predicting ORFs in a single contig.

        Args:
            contig_idx (int): The integer index of the contig.
            contig_seq (bytes): The raw nucleotide sequence of the contig.

        Returns:
            tuple[int, Intervals, list[bytes]] | None: The contig index, the spatial intervals of 
                the predicted genes, and a list of translated amino acid byte strings. Returns None if no ORFs found.
        """
        if orfs := self._gene_finder.find_genes(contig_seq):
            starts = np.array([g.start for g in orfs], dtype=np.int32)
            ends = np.array([g.end for g in orfs], dtype=np.int32)
            strands = np.array([g.strand for g in orfs], dtype=np.int8)
            translations = [g.translation() for g in orfs]
            
            intervals = Intervals(starts, ends, strands)
            return contig_idx, intervals, translations
        return None

    def _evaluate_locus_status(self, clusters: GeneClusters, allele_scores: npt.NDArray[np.float64]) -> tuple[LocusStatus, str, LocusComparisons]:
        """Evaluates the found gene clusters against the database to determine the best locus match.

        This method implements the core logic for identifying loci:
        
        1.  Extracts the set of unique database cluster IDs (representing gene families) found in the query.
        2.  Iterates through all known loci in the database.
        3.  For each locus, extracts the set of expected *core* genes (ignoring 'extra' genes).
        4.  Compares the found set to the expected set to determine missing and added genes.
        5.  Assigns a `LocusStatus` based on the comparison (e.g., perfect set match, single deletion).
        6.  If it's a perfect set match and all genes are on one contig, it checks the physical order
            to see if it qualifies as a `GENE_ORDER_MATCH`.
        7.  Maintains the best match found across all loci based on status rank and Jaccard similarity.

        Args:
            clusters (GeneClusters): The filtered and clustered gene matches found in the genome.

        Returns:
            tuple[LocusStatus, str]: A tuple containing the determined status enum and the name
                of the best-matching reference locus.
        """
        n_loci = len(self._db.loci)

        found_clusters = self._db.gene_cluster_ids[clusters.target_indices]
        
        # 1. Create a boolean mask of found clusters
        found_mask = np.zeros(len(self._db.cluster_keys), dtype=bool)
        found_mask[found_clusters] = True
        
        # Ignore universally recognized extra genes when validating core additions/deletions
        found_core_mask = found_mask & ~self._extra_cluster_mask
        
        # 2. Vectorized Set Operations
        missing_counts = (self._expected_matrix & ~found_core_mask).sum(axis=1, dtype=np.int32)
        added_counts = (~self._expected_matrix & found_core_mask).sum(axis=1, dtype=np.int32)
        
        intersect_counts = (self._expected_matrix & found_core_mask).sum(axis=1, dtype=np.float64)
        expected_counts = self._expected_matrix.sum(axis=1, dtype=np.float64)
        
        completeness = np.zeros(n_loci, dtype=np.float64)
        valid_expected = expected_counts > 0
        completeness[valid_expected] = intersect_counts[valid_expected] / expected_counts[valid_expected]

        # 3. Vectorized Status Assignment
        statuses = np.full(n_loci, LocusStatus.NOVEL.value, dtype=np.int8)
        exact_set_matches = (missing_counts <= 0) & valid_expected
        statuses[exact_set_matches] = LocusStatus.GENE_SET_MATCH.value

        # 4. Determine Ordering status for exact set matches
        if is_single_contig := len(np.unique(clusters.contig_indices)) == 1:
            order = np.argsort(clusters.intervals.starts)
            ordered_found_list = found_clusters[order].tolist()
            
            for idx in np.where(exact_set_matches)[0]:
                exp_list = self._expected_lists[idx]
                if len(ordered_found_list) == len(exp_list):
                    if ordered_found_list == exp_list or ordered_found_list == exp_list[::-1]:
                        statuses[idx] = LocusStatus.GENE_ORDER_MATCH.value

        metrics = LocusComparisons(
            statuses=statuses,
            intersect_genes=intersect_counts,
            completeness=completeness,
            allele_scores=allele_scores,
            missing_genes=missing_counts,
            added_genes=added_counts
        )

        if (best_idx := metrics.best_locus_idx) >= 0 and valid_expected[best_idx]:
            best_status = LocusStatus(statuses[best_idx])
            match_locus = self._db.loci.ids[best_idx]
        else:
            best_status = LocusStatus.NOVEL
            match_locus = "Unknown"

        return best_status, match_locus, metrics

    def __call__(self, genome: GenomeAssembly | str | Path) -> AnnotationResult | None:
        """Executes the full annotation pipeline on a given genome assembly.

        Args:
            genome (GenomeAssembly | str | Path): The genome to annotate.

        Returns:
            AnnotationResult | None: The final annotation result, or None if no valid loci were found.
        """
        genome = GenomeAssembly.ensure(genome)

        futures = [self.executor.submit(self._predict_orfs, genome.id_map[c.id], c.seq) for c in genome.contigs]
        results = [i for future in as_completed(futures) if (i := future.result())]
        
        if not results:
            return None

        # Consolidate results into flat SoA arrays
        intervals_list, contigs_list, translations_list = [], [], []
        for ctg_idx, itv, trans in results:
            intervals_list.append(itv)
            contigs_list.append(np.full(len(itv), ctg_idx, dtype=np.int32))
            translations_list.extend(trans)

        global_intervals = Intervals.concat(intervals_list)
        global_contigs = np.concatenate(contigs_list)
        
        # Generate informative IDs once globally (e.g., "contig_name_1", "contig_name_2")
        global_ids = tuple(f"{genome.contigs.ids[c]}_{i+1}" for c, i in zip(global_contigs, global_intervals.original_indices))
        translations = Sequences.from_bytes(translations_list, ids=global_ids)

        # Ensure sort_by_hash=True so the genome index can be the target of the search
        translation_index = RandstrobeIndex.build(translations, k=self._index.k, s=self._index.s,
                                                  w_max=self._index.w_max, sort_by_hash=True)

        # Search the genome ORFs using the database index
        # This yields the best genome ORF for EVERY allele in the database
        seeds = translation_index.top_hits(self._index, min_score=1)

        # 1. Dynamic Thresholding: Normalize scores by expected gene lengths to remove size bias
        seed_coverage = seeds.scores / self._db_strobe_counts[seeds.query_indices]
        
        # Combine proportional coverage with an absolute minimum floor (e.g. 3 strobes) to prevent noise on tiny genes
        valid_mask = (seed_coverage >= self._min_strobe_coverage) & (seeds.scores >= 3)
        valid_seeds = seeds.filter(valid_mask)
        valid_coverage = seed_coverage[valid_mask]

        if len(valid_seeds) == 0:
            return None

        # 2. Calculate unbiased allele scores for all loci using ALL valid seeds
        parent_loci = self._db.gene_locus_indices[valid_seeds.query_indices]
        allele_scores = np.bincount(parent_loci, weights=valid_coverage, minlength=len(self._db.loci)).astype(np.float64)

        # 3. Reduce to the single best database hit per genome ORF
        # By sorting on the negative float coverage, we dynamically pick the allele that best fits the genome ORF!
        order = np.lexsort((-valid_coverage, valid_seeds.target_indices))
        _, best_idx = np.unique(valid_seeds.target_indices[order], return_index=True)
        best_indices = order[best_idx]
        
        best_mask = np.zeros(len(valid_seeds), dtype=bool)
        best_mask[best_indices] = True
        confident_seeds = valid_seeds.filter(best_mask)
        confident_coverage = valid_coverage[best_mask]

        # Extract batches based on seeds
        hit_genome_indices = confident_seeds.target_indices
        hit_db_indices = confident_seeds.query_indices

        hit_intervals = global_intervals[hit_genome_indices]
        hit_contigs = global_contigs[hit_genome_indices]

        # First pass: Cluster the confident seed hits sequentially before alignment
        raw_cluster_ids = hit_intervals.cluster_sequential(
            tolerance=self._orf_tolerance, 
            group_by=hit_contigs, 
            enforce_strand=self._enforce_strand
        )

        # Filter Clusters by Minimum Size
        unique_ids, counts = np.unique(raw_cluster_ids, return_counts=True)
        large_clusters = unique_ids[counts >= self._min_cluster_size]
        size_mask = np.isin(raw_cluster_ids, large_clusters)
        
        if not np.any(size_mask):
            return None

        clustered_seeds = confident_seeds.filter(size_mask)
        clustered_coverage = confident_coverage[size_mask]
        clustered_genome_indices = hit_genome_indices[size_mask]
        clustered_db_indices = hit_db_indices[size_mask]

        # Second pass: Pairwise alignment strictly on the clustered hits
        aligned_db_genes = self._translations[clustered_db_indices]
        aligned_genome_orfs = translations[clustered_genome_indices]
        alignments = self._aligner(aligned_db_genes, aligned_genome_orfs, clustered_seeds)

        global_batch = GeneClusters(
            spatial_ids=raw_cluster_ids[size_mask],
            query_indices=clustered_genome_indices,
            target_indices=clustered_db_indices,
            contig_indices=hit_contigs[size_mask],
            intervals=hit_intervals[size_mask],
            seed_scores=clustered_seeds.scores,
            seed_coverage=clustered_coverage,
            pidents=alignments.pidents
        )
        
        # 6. Evaluate the combined global batch to support loci split across contigs
        best_status, best_locus, best_metrics = self._evaluate_locus_status(global_batch, allele_scores)
        
        # Filter the global batch to ONLY include hits that belong to the winning locus
        best_idx = best_metrics.best_locus_idx
        if best_idx >= 0:
            found_clusters = self._db.gene_cluster_ids[global_batch.target_indices]
            is_expected = self._expected_matrix[best_idx, found_clusters]
            is_extra = self._extra_cluster_mask[found_clusters]
            belongs_mask = is_expected | is_extra
            
            # Resolve Duplicates Dynamically: If the cluster grabbed multiple ORFs hitting the same 
            # gene family (e.g. two `wza`s), keep ONLY the ORF with the highest seed coverage!
            belongs_idx = np.where(belongs_mask)[0]
            
            order = np.argsort(-global_batch.seed_coverage[belongs_idx])
            sorted_belongs_idx = belongs_idx[order]
            
            _, unique_idx = np.unique(found_clusters[sorted_belongs_idx], return_index=True)
            
            final_locus_mask = np.zeros(len(global_batch), dtype=bool)
            final_locus_mask[sorted_belongs_idx[unique_idx]] = True
            
            best_cluster_batch = global_batch.filter(final_locus_mask)
        else:
            best_cluster_batch = global_batch
            
        # 7. Extract Novel/Unannotated ORFs physically located within the locus boundaries (per contig)
        novel_orfs_list = []
        for ctg in np.unique(best_cluster_batch.contig_indices):
            ctg_mask = best_cluster_batch.contig_indices == ctg
            ctg_batch = best_cluster_batch.filter(ctg_mask)
            bounds = ctg_batch.intervals.envelope
            
            if bounds is None:
                continue

            # Find all global ORFs on this specific contig that fall within the BGC bounds
            in_bounds_mask = (
                (global_contigs == ctg) & 
                (global_intervals.starts >= bounds.start) & 
                (global_intervals.ends <= bounds.end)
            )
            
            # Mathematically subtract the ORFs we successfully annotated to leave only the novel ones
            all_locus_orf_indices = np.where(in_bounds_mask)[0]
            novel_indices = np.setdiff1d(all_locus_orf_indices, ctg_batch.query_indices)
            if len(novel_indices) > 0:
                novel_orfs_list.append(translations[novel_indices])
                
        if novel_orfs_list:
            novel_orfs = Sequences.concat(novel_orfs_list)
        else:
            novel_orfs = Sequences.empty()

        # Extract Extra Homologs (High-scoring hits outside the primary BGC)
        # We use np.isin to instantly find confident seeds that were NOT part of the winning cluster
        is_in_locus = np.isin(confident_seeds.target_indices, best_cluster_batch.query_indices)
        outside_mask = ~is_in_locus
        outside_seeds = confident_seeds.filter(outside_mask)

        extra_homologs = Homologs(
            query_indices=outside_seeds.target_indices,
            target_indices=outside_seeds.query_indices,
            contig_indices=hit_contigs[outside_mask],
            intervals=hit_intervals[outside_mask],
            seed_scores=outside_seeds.scores,
            seed_coverage=confident_coverage[outside_mask]
        )

        return AnnotationResult(genome.id, best_cluster_batch, best_status, best_locus, best_metrics, novel_orfs,
                                extra_homologs)


def test():
    self = Annotator(Database.load('kpsc_k'))  # For testing, to be removed
    genome = GenomeAssembly.from_file(Path('..') / 'kpsc' / 'T7-472.fasta.gz')  # For testing, to be removed
    # This has the known locus KL63
    result = self(genome)