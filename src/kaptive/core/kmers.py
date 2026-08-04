r"""K-mer indexing, FracMinHash, Randstrobe sketch construction, and alignment seed management.

This module provides sequence sketch indices ([`FracMinHashIndex`][kaptive.core.kmers.FracMinHashIndex],
[`RandstrobeIndex`][kaptive.core.kmers.RandstrobeIndex]), alignment seed batch containers
([`Seeds`][kaptive.core.kmers.Seeds]), and Numba JIT parallel search kernels.
"""

from __future__ import annotations

from collections.abc import Iterable
from dataclasses import dataclass
from functools import cache
from typing import Any, NamedTuple, Self

import numpy as np
import numpy.typing as npt
from numba import njit, prange

from kaptive.core.collections import BatchedContainer
from kaptive.core.interval import Intervals
from kaptive.core.seq import Sequences

# Constants ------------------------------------------------------------------------------------------------------------
MINHASH_DTYPE = np.dtype(
    [
        ("hash", np.uint64),
        ("seq_idx", np.uint32),
        ("pos1", np.uint32),
    ]
)

RANDSTROBE_DTYPE = np.dtype(
    [
        ("hash", np.uint64),
        ("seq_idx", np.uint32),
        ("pos1", np.uint32),
        ("pos2", np.uint32),
    ]
)


# Classes --------------------------------------------------------------------------------------------------------------
class Seed(NamedTuple):
    r"""Alignment seed representing a potential matching region between query and target.

    Attributes:
        query_index (int): Index of the query sequence.
        target_index (int): Index of the target sequence.
        score (int): Match score (number of shared k-mer hashes or randstrobes).
        offset (int): Diagonal offset calculated as `query_pos - target_pos`.
    """

    query_index: int
    target_index: int
    score: int
    offset: int


@dataclass(frozen=True, slots=True)
class Seeds(BatchedContainer[Seed, "Seeds"]):
    r"""Structure-of-Arrays (SoA) batch container for alignment seeds.

    Stores query indices, target indices, scores, and diagonal offsets as 1D NumPy arrays.

    Attributes:
        query_indices (npt.NDArray[np.uint32]): 1D array of query sequence indices.
        target_indices (npt.NDArray[np.uint32]): 1D array of target sequence indices.
        scores (npt.NDArray[np.uint32]): 1D array of alignment match scores.
        offsets (npt.NDArray[np.int32]): 1D array of diagonal offset values.
    """

    query_indices: npt.NDArray[np.uint32]
    target_indices: npt.NDArray[np.uint32]
    scores: npt.NDArray[np.uint32]
    offsets: npt.NDArray[np.int32]

    def __len__(self) -> int:
        r"""Return the number of seeds in the batch.

        Returns:
            int: Total seed count.
        """
        return len(self.query_indices)

    def __getitem__(self, item: Any) -> Seed | Seeds:
        r"""Access seeds by integer index, slice, or boolean/integer NumPy array.

        Args:
            item (Any): Integer index, slice, or mask array.

        Returns:
            Seed | Seeds: A single [`Seed`][kaptive.core.kmers.Seed] or a
                sliced [`Seeds`][kaptive.core.kmers.Seeds] collection.

        Raises:
            IndexError: If integer index is out of bounds.
        """
        if isinstance(item, (int, np.integer)):
            if item < 0:
                item += len(self)
            if item < 0 or item >= len(self):
                raise IndexError("Batch index out of range")
            return Seed(
                int(self.query_indices[item]),
                int(self.target_indices[item]),
                int(self.scores[item]),
                int(self.offsets[item]),
            )

        if isinstance(item, slice):
            indices = np.arange(len(self))[item]
        else:
            item_arr = np.asarray(item)
            indices = np.nonzero(item_arr)[0] if item_arr.dtype.kind == "b" else item_arr

        return Seeds(
            self.query_indices[indices],
            self.target_indices[indices],
            self.scores[indices],
            self.offsets[indices],
        )

    @classmethod
    def concat(cls, batches: Iterable[Self]) -> Self:  # type: ignore
        r"""Concatenate multiple [`Seeds`][kaptive.core.kmers.Seeds] collections into a single batch.

        Args:
            batches (Iterable[Seeds]): Iterable of [`Seeds`][kaptive.core.kmers.Seeds] objects.

        Returns:
            Seeds: Concatenated [`Seeds`][kaptive.core.kmers.Seeds] collection.
        """
        batches_list = list(batches)
        if not batches_list:
            return cls.empty()  # type: ignore
        return cls(
            np.concatenate([b.query_indices for b in batches_list]),
            np.concatenate([b.target_indices for b in batches_list]),
            np.concatenate([b.scores for b in batches_list]),
            np.concatenate([b.offsets for b in batches_list]),
        )

    @classmethod
    def empty(cls) -> Seeds:
        r"""Create an empty [`Seeds`][kaptive.core.kmers.Seeds] collection.

        Returns:
            Seeds: Empty [`Seeds`][kaptive.core.kmers.Seeds] collection.
        """
        return cls(
            np.empty(0, dtype=np.uint32),
            np.empty(0, dtype=np.uint32),
            np.empty(0, dtype=np.uint32),
            np.empty(0, dtype=np.int32),
        )

    def filter(self, mask: npt.NDArray[np.bool_]) -> Seeds:
        r"""Return a new Seeds collection containing only records where the mask is True.

        Args:
            mask (npt.NDArray[np.bool_]): 1D boolean mask array.

        Returns:
            Seeds: Filtered [`Seeds`][kaptive.core.kmers.Seeds] collection.
        """
        return Seeds(
            self.query_indices[mask],
            self.target_indices[mask],
            self.scores[mask],
            self.offsets[mask],
        )

    def to_intervals(self, query_lengths: npt.NDArray[np.int32]) -> Intervals:
        r"""Convert seeds into target coordinates [`Intervals`][kaptive.core.interval.Intervals].

        Args:
            query_lengths (npt.NDArray[np.int32]): Query sequence lengths array.

        Returns:
            Intervals: Spatial [`Intervals`][kaptive.core.interval.Intervals] collection on target sequences.
        """
        t_starts = -self.offsets
        q_lens = query_lengths[self.query_indices]
        t_ends = t_starts + q_lens

        return Intervals(
            starts=t_starts,
            ends=t_ends,
            strands=np.ones(len(self), dtype=np.int8),
            original_indices=np.arange(len(self), dtype=np.int32),
        )

    def cull_overlaps(
        self,
        query_lengths: npt.NDArray[np.int32],
        max_overlap_fraction: float = 0.1,
        priority_mask: npt.NDArray[np.bool_] | None = None,
    ) -> Seeds:
        r"""Greedily cull seeds that overlap significantly on the target sequence.

        Args:
            query_lengths (npt.NDArray[np.int32]): Query sequence lengths array.
            max_overlap_fraction (float): Maximum allowable overlap fraction. Defaults to 0.1.
            priority_mask (npt.NDArray[np.bool_] | None): Optional boolean mask for priority score boost.

        Returns:
            Seeds: Filtered, non-overlapping [`Seeds`][kaptive.core.kmers.Seeds] collection.
        """
        n = len(self)
        if n == 0:
            return self

        # Sort order: priority_mask (True first), then scores (descending)
        if priority_mask is None:
            priority_mask = np.zeros(n, dtype=np.bool_)

        order = np.lexsort((-self.scores, ~priority_mask)).astype(np.int32)

        # We pass target_indices as the 'names' array to only cull overlaps on the same contig
        intervals = self.to_intervals(query_lengths)
        kept_mask = intervals.cull_overlaps(
            order=order,
            max_overlap_fraction=max_overlap_fraction,
            group_by=self.target_indices.astype(np.int32),
        )
        return self.filter(kept_mask)

    def top_hits(self, min_score: int = 1) -> Seeds:
        r"""Reduce the batch to only the highest-scoring seed for each query.

        Args:
            min_score (int): Minimum match score threshold. Defaults to 1.

        Returns:
            Seeds: Filtered [`Seeds`][kaptive.core.kmers.Seeds] collection containing top hit per query.
        """
        if len(self) == 0:
            return self

        # np.lexsort prioritizes the last array.
        # Primary key: query_indices (ascending). Secondary: scores (descending via ~)
        order = np.lexsort((~self.scores, self.query_indices))

        # Find the first occurrence of each query (which corresponds to its highest score)
        _, unique_idx = np.unique(self.query_indices[order], return_index=True)

        # Map back to the original array indices and sort to preserve the original query order
        best_idx = order[unique_idx]
        best_idx.sort()

        best_batch = Seeds(
            self.query_indices[best_idx],
            self.target_indices[best_idx],
            self.scores[best_idx],
            self.offsets[best_idx],
        )

        if min_score > 0:
            return best_batch.filter(best_batch.scores >= min_score)

        return best_batch

    def extract_sequences(self, queries: Sequences, targets: Sequences) -> tuple[Sequences, Sequences]:
        r"""Extract parallel batches of query and target sequences mapped by this seed batch.

        Args:
            queries (Sequences): Source collection of query sequences.
            targets (Sequences): Source collection of target sequences.

        Returns:
            tuple[Sequences, Sequences]: Paired query and target sequence collections.
        """
        return queries[self.query_indices], targets[self.target_indices]  # type: ignore


@dataclass(frozen=True, slots=True, kw_only=True)
class BaseKmerIndex:
    r"""Abstract base class for high-performance k-mer based sequence indices using Array of Structs (AoS).

    Attributes:
        records (npt.NDArray): Structured NumPy array holding index records.
        n_seqs (int): Total number of indexed sequences.
        is_sorted (bool): True if records are sorted by hash.
        k (int): K-mer length parameter.
    """

    records: npt.NDArray
    n_seqs: int = 0
    is_sorted: bool = False
    k: int = 10

    def __len__(self) -> int:
        r"""Return the number of records in the index.

        Returns:
            int: Record count.
        """
        return len(self.records)

    @classmethod
    def empty(cls) -> BaseKmerIndex:
        r"""Create an empty BaseKmerIndex.

        Raises:
            NotImplementedError: Must be implemented by subclasses.
        """
        raise NotImplementedError

    @classmethod
    def build(cls, batch: Sequences, **kwargs: Any) -> BaseKmerIndex:
        r"""Build a k-mer index from a sequence collection.

        Args:
            batch (Sequences): Sequence collection to index.
            **kwargs (Any): Subclass-specific options.

        Raises:
            NotImplementedError: Must be implemented by subclasses.
        """
        raise NotImplementedError

    def _build_queries(self, queries: Sequences) -> BaseKmerIndex:
        r"""Build query index with matching parameters.

        Args:
            queries (Sequences): Query sequence collection.

        Raises:
            NotImplementedError: Must be implemented by subclasses.
        """
        raise NotImplementedError

    def _prep_queries(self, queries: BaseKmerIndex | Sequences) -> BaseKmerIndex:
        r"""Prepare and validate query index against target index state.

        Args:
            queries (BaseKmerIndex | Sequences): Input query index or raw sequences.

        Returns:
            BaseKmerIndex: Compiled query index.

        Raises:
            ValueError: If target index is unsorted or query index is sorted.
        """
        if not self.is_sorted:
            raise ValueError("Target index must be sorted by hash for binary search. Build it with sort_by_hash=True.")

        if isinstance(queries, Sequences):
            queries = self._build_queries(queries)

        if queries.is_sorted:
            raise ValueError("Query index must NOT be sorted by hash. Build it with sort_by_hash=False.")

        return queries

    def top_hits(self, queries: BaseKmerIndex | Sequences, min_score: int = 1) -> Seeds:
        r"""Find the single best-matching target sequence for each query sequence.

        Args:
            queries (BaseKmerIndex | Sequences): Query sequence index or raw sequences.
            min_score (int): Minimum match score threshold. Defaults to 1.

        Returns:
            Seeds: Matched [`Seeds`][kaptive.core.kmers.Seeds] collection.
        """
        if len(queries) == 0 or len(self) == 0:
            return Seeds.empty()

        queries_idx = self._prep_queries(queries)
        q_offsets = _compute_query_offsets(queries_idx.records, queries_idx.n_seqs)

        seeds = Seeds(
            *_intersect_top_hit_kernel(queries_idx.records, q_offsets, queries_idx.n_seqs, self.records, self.n_seqs)
        )

        if min_score > 0:
            return seeds.filter(seeds.scores >= min_score)

        return seeds


@dataclass(frozen=True, slots=True, kw_only=True)
class FracMinHashIndex(BaseKmerIndex):
    r"""Specialized index for fast nucleotide sequence comparisons using FracMinHash.

    Attributes:
        scaled (int): FracMinHash scale factor (e.g. 100).
        canonical (bool): True if canonical k-mers (min of fwd/rev) are hashed.
        bits_per_char (int): Alphabet bits per character (2 for DNA).
        lut (npt.NDArray[np.uint8] | None): Alphabet character lookup table.
    """

    scaled: int = 100
    canonical: bool = True
    bits_per_char: int = 2
    lut: npt.NDArray[np.uint8] | None = None

    @classmethod
    def empty(cls) -> FracMinHashIndex:
        r"""Create an empty FracMinHashIndex.

        Returns:
            FracMinHashIndex: Empty index instance.
        """
        return cls(
            records=np.empty(0, dtype=MINHASH_DTYPE),
            n_seqs=0,
            is_sorted=False,
            k=21,
            scaled=100,
            canonical=True,
            bits_per_char=2,
            lut=None,
        )

    @classmethod
    def build(
        cls,
        batch: Sequences,
        k: int = 21,
        scaled: int = 100,
        canonical: bool = True,
        seed: int = 42,
        sort_by_hash: bool = False,
        lut: npt.NDArray[np.uint8] | None = None,
        bits_per_char: int = 2,
        **kwargs: Any,
    ) -> FracMinHashIndex:
        r"""Build a FracMinHashIndex from nucleotide sequence batch.

        Args:
            batch (Sequences): Input sequence collection.
            k (int): K-mer length. Defaults to 21.
            scaled (int): Sampling scale factor. Defaults to 100.
            canonical (bool): True if using canonical k-mers. Defaults to True.
            seed (int): Seed for hashing. Defaults to 42.
            sort_by_hash (bool): True to sort output records by hash. Defaults to False.
            lut (npt.NDArray[np.uint8] | None): Optional DNA lookup table.
            bits_per_char (int): Bits per character (2 for DNA).
            **kwargs (Any): Additional parameters.

        Returns:
            FracMinHashIndex: Built [`FracMinHashIndex`][kaptive.core.kmers.FracMinHashIndex] instance.
        """
        if len(batch) == 0:
            return cls.empty()

        n_seqs = len(batch)
        kernel_lut = lut if lut is not None else _dna_lut()

        # Pass 1: Count fracminhashes per sequence to allocate exact memory footprint
        counts = _count_fracminhash_kernel(
            batch.seqs, batch.offsets, batch.lengths, kernel_lut, k, scaled, canonical, bits_per_char
        )

        if (total_hashes := np.sum(counts)) == 0:
            return cls.empty()

        # Generate exact write offsets via cumulative sum
        out_offsets = np.empty(n_seqs, dtype=np.uint32)
        current_offset = 0
        for i in range(n_seqs):
            out_offsets[i] = current_offset
            current_offset += counts[i]

        out_records = np.empty(total_hashes, dtype=MINHASH_DTYPE)

        # Pass 2: Populate the fracminhashes in parallel
        _populate_fracminhash_kernel(
            batch.seqs,
            batch.offsets,
            batch.lengths,
            kernel_lut,
            k,
            scaled,
            canonical,
            bits_per_char,
            out_offsets,
            out_records,
        )

        if sort_by_hash:
            out_records = _radix_sort_records(out_records)

        return cls(
            records=out_records,
            n_seqs=n_seqs,
            is_sorted=sort_by_hash,
            k=k,
            scaled=scaled,
            canonical=canonical,
            bits_per_char=bits_per_char,
            lut=lut,
        )

    def _build_queries(self, queries: Sequences) -> FracMinHashIndex:
        r"""Build query FracMinHashIndex from sequences matching current settings.

        Args:
            queries (Sequences): Input query sequences.

        Returns:
            FracMinHashIndex: Unsorted query index.
        """
        return self.build(
            queries,
            k=self.k,
            scaled=self.scaled,
            canonical=self.canonical,
            sort_by_hash=False,
            lut=self.lut,
            bits_per_char=self.bits_per_char,
        )

    def to_sorted(self) -> FracMinHashIndex:
        r"""Return a new FracMinHashIndex with records sorted by hash.

        Returns:
            FracMinHashIndex: Sorted [`FracMinHashIndex`][kaptive.core.kmers.FracMinHashIndex].
        """
        if self.is_sorted:
            return self
        return self.__class__(
            records=_radix_sort_records(self.records),
            n_seqs=self.n_seqs,
            is_sorted=True,
            k=self.k,
            scaled=self.scaled,
            canonical=self.canonical,
            bits_per_char=self.bits_per_char,
            lut=self.lut,
        )


@dataclass(frozen=True, slots=True, kw_only=True)
class RandstrobeIndex(BaseKmerIndex):
    r"""Specialized index for fast amino-acid sequence comparisons using syncmer-linked randstrobes.

    Attributes:
        s (int): Sub-k-mer size for open syncmer evaluation. Defaults to 5.
        w_min (int): Minimum syncmer offset window bound. Defaults to 1.
        w_max (int): Maximum syncmer offset window bound. Defaults to 5.
        lut (npt.NDArray[np.uint8] | None): Optional amino acid lookup table.
    """

    s: int = 5
    w_min: int = 1
    w_max: int = 5
    lut: npt.NDArray[np.uint8] | None = None

    @classmethod
    def empty(cls) -> RandstrobeIndex:
        r"""Create an empty RandstrobeIndex.

        Returns:
            RandstrobeIndex: Empty index instance.
        """
        return cls(
            records=np.empty(0, dtype=RANDSTROBE_DTYPE),
            n_seqs=0,
            is_sorted=False,
            k=10,
            s=5,
            w_min=1,
            w_max=5,
            lut=None,
        )

    @classmethod
    def build(
        cls,
        batch: Sequences,
        k: int = 10,
        s: int = 5,
        w_min: int = 1,
        w_max: int = 5,
        canonical: bool = True,
        seed: int = 42,
        sort_by_hash: bool = False,
        lut: npt.NDArray[np.uint8] | None = None,
        **kwargs: Any,
    ) -> RandstrobeIndex:
        r"""Build a RandstrobeIndex from sequence batch.

        Args:
            batch (Sequences): Input sequence collection.
            k (int): K-mer length. Defaults to 10.
            s (int): Sub-k-mer size. Defaults to 5.
            w_min (int): Min syncmer window bound. Defaults to 1.
            w_max (int): Max syncmer window bound. Defaults to 5.
            canonical (bool): Canonical option. Defaults to True.
            seed (int): Random seed. Defaults to 42.
            sort_by_hash (bool): True to sort records by hash. Defaults to False.
            lut (npt.NDArray[np.uint8] | None): Optional alphabet lookup table.
            **kwargs (Any): Additional parameters.

        Returns:
            RandstrobeIndex: Constructed [`RandstrobeIndex`][kaptive.core.kmers.RandstrobeIndex].

        Raises:
            ValueError: If `s >= k`.
        """
        if s >= k:
            raise ValueError("Sub-k-mer size (s) must be strictly less than k-mer size (k).")

        if len(batch) == 0:
            return cls.empty()

        # Pass 1: Count randstrobes per sequence to allocate exact memory footprint
        kernel_lut = lut if lut is not None else _mmseqs12_lut()
        counts = _count_randstrobes_kernel(batch.seqs, batch.offsets, batch.lengths, kernel_lut, k, s, w_min)
        if (total_randstrobes := np.sum(counts)) == 0:
            return cls.empty()

        # Generate exact write offsets via cumulative sum
        out_offsets = np.empty(len(counts), dtype=np.uint32)
        current_offset = 0
        for i in range(len(counts)):
            out_offsets[i] = current_offset
            current_offset += counts[i]

        # Pre-allocate output arrays
        out_records = np.empty(total_randstrobes, dtype=RANDSTROBE_DTYPE)

        # Pass 2: Populate the randstrobes in parallel
        _populate_randstrobes_kernel(
            batch.seqs, batch.offsets, batch.lengths, kernel_lut, k, s, w_min, w_max, out_offsets, out_records
        )

        if sort_by_hash:
            out_records = _radix_sort_records(out_records)

        return cls(
            records=out_records,
            n_seqs=len(batch),
            is_sorted=sort_by_hash,
            k=k,
            s=s,
            w_min=w_min,
            w_max=w_max,
            lut=lut,
        )

    def _build_queries(self, queries: Sequences) -> RandstrobeIndex:
        r"""Build query RandstrobeIndex from sequences matching current settings.

        Args:
            queries (Sequences): Input query sequences.

        Returns:
            RandstrobeIndex: Unsorted query index.
        """
        return self.build(
            queries, k=self.k, s=self.s, w_min=self.w_min, w_max=self.w_max, sort_by_hash=False, lut=self.lut
        )


# Functions ------------------------------------------------------------------------------------------------------------
@cache
def _mmseqs12_lut(fill_value: int = 12) -> npt.NDArray[np.uint8]:
    r"""Generate a 256-byte array for O(1) MMSEQS12 alphabet mapping.

    Args:
        fill_value (int): Default integer index for unknown characters. Defaults to 12.

    Returns:
        npt.NDArray[np.uint8]: Read-only 256-element byte lookup table array.
    """
    mapping = {
        b"A": 0,
        b"S": 0,
        b"T": 0,
        b"L": 1,
        b"M": 1,
        b"I": 2,
        b"V": 2,
        b"K": 3,
        b"R": 3,
        b"E": 4,
        b"Q": 4,
        b"N": 5,
        b"D": 5,
        b"F": 6,
        b"Y": 6,
        b"C": 7,
        b"G": 8,
        b"H": 9,
        b"P": 10,
        b"W": 11,
    }
    # Map unknowns (e.g. X) to 12
    lut = np.full(256, fill_value, dtype=np.uint8)
    for source_byte, target_int in mapping.items():
        lut[source_byte[0]] = target_int
        lut[source_byte.lower()[0]] = target_int
    lut.flags.writeable = False
    return lut


@cache
def _dna_lut(fill_value: int = 4) -> npt.NDArray[np.uint8]:
    r"""Generate a 256-byte array for O(1) DNA alphabet mapping.

    Args:
        fill_value (int): Default integer index for non-ACGT characters. Defaults to 4.

    Returns:
        npt.NDArray[np.uint8]: Read-only 256-element byte lookup table array.
    """
    lut = np.full(256, fill_value, dtype=np.uint8)
    mapping = {b"A": 0, b"C": 1, b"T": 2, b"G": 3}
    for source_byte, target_int in mapping.items():
        lut[source_byte[0]] = target_int
        lut[source_byte.lower()[0]] = target_int
    lut.flags.writeable = False
    return lut


@cache
def _aa_lut(fill_value: int = 22) -> npt.NDArray[np.uint8]:
    r"""Generate a 256-byte array for O(1) Amino Acid alphabet mapping.

    Args:
        fill_value (int): Default integer index for non-standard characters. Defaults to 22.

    Returns:
        npt.NDArray[np.uint8]: Read-only 256-element byte lookup table array.
    """
    lut = np.full(256, fill_value, dtype=np.uint8)
    # 20 standard AAs + U (Selenocysteine) + O (Pyrrolysine) = 22 valid chars
    mapping = {
        b"A": 0,
        b"C": 1,
        b"D": 2,
        b"E": 3,
        b"F": 4,
        b"G": 5,
        b"H": 6,
        b"I": 7,
        b"K": 8,
        b"L": 9,
        b"M": 10,
        b"N": 11,
        b"P": 12,
        b"Q": 13,
        b"R": 14,
        b"S": 15,
        b"T": 16,
        b"V": 17,
        b"W": 18,
        b"Y": 19,
        b"U": 20,
        b"O": 21,
    }
    for source_byte, target_int in mapping.items():
        lut[source_byte[0]] = target_int
        lut[source_byte.lower()[0]] = target_int
    lut.flags.writeable = False
    return lut


# Kernels --------------------------------------------------------------------------------------------------------------
@njit(cache=True, nogil=True, inline="always")
def _splitmix64(x: np.uint64) -> np.uint64:
    r"""Fast, invertible 64-bit integer mixing function for randomizing hashes.

    Args:
        x (np.uint64): 64-bit unsigned integer value.

    Returns:
        np.uint64: Mixed 64-bit unsigned hash value.
    """
    x = np.uint64(x + np.uint64(0x9E3779B97F4A7C15))
    x = (x ^ (x >> np.uint64(30))) * np.uint64(0xBF58476D1CE4E5B9)
    x = (x ^ (x >> np.uint64(27))) * np.uint64(0x94D049BB133111EB)
    return x ^ (x >> np.uint64(31))


@njit(cache=True, nogil=True)
def _radix_sort_records(records: npt.NDArray) -> npt.NDArray:
    r"""In-place radix sort of records by the 'hash' field.

    Args:
        records (npt.NDArray): Structured 1D array containing a `'hash'` uint64 field.

    Returns:
        npt.NDArray: Sorted array of records.
    """
    n = len(records)
    src = records.copy()
    dst = np.empty_like(records)

    hist = np.zeros(256, dtype=np.uint32)
    pos = np.zeros(256, dtype=np.uint32)

    for p in range(8):
        shift = np.uint64(p * 8)

        hist.fill(0)
        for i in range(n):
            byte_val = (src[i]["hash"] >> shift) & np.uint64(0xFF)
            hist[byte_val] += 1

        current_pos = np.uint32(0)
        for i in range(256):
            pos[i] = current_pos
            current_pos += hist[i]

        for i in range(n):
            byte_val = (src[i]["hash"] >> shift) & np.uint64(0xFF)
            write_idx = pos[byte_val]
            dst[write_idx] = src[i]
            pos[byte_val] += 1

        tmp = src
        src = dst
        dst = tmp

    return src


@njit(parallel=True, cache=True, nogil=True)
def _count_fracminhash_kernel(
    seqs: npt.NDArray[np.uint8],
    offsets: npt.NDArray[np.int32],
    lengths: npt.NDArray[np.int32],
    lut: npt.NDArray[np.uint8],
    k: int,
    scaled: int,
    canonical: bool,
    bits_per_char: int,
) -> npt.NDArray[np.uint32]:
    r"""Count FracMinHash occurrences per sequence in parallel.

    Args:
        seqs (npt.NDArray[np.uint8]): Concatenated sequence bytes.
        offsets (npt.NDArray[np.int32]): Sequence start offsets.
        lengths (npt.NDArray[np.int32]): Sequence lengths.
        lut (npt.NDArray[np.uint8]): Alphabet character mapping table.
        k (int): K-mer length.
        scaled (int): Sampling scale factor.
        canonical (bool): True if computing canonical hashes.
        bits_per_char (int): Alphabet bits per character.

    Returns:
        npt.NDArray[np.uint32]: 1D array of FracMinHash counts per sequence.
    """
    n_seqs = len(offsets)
    counts = np.zeros(n_seqs, dtype=np.uint32)
    mask = (np.uint64(1) << np.uint64(bits_per_char * k)) - np.uint64(1)
    max_val = np.uint64(1) << np.uint64(bits_per_char)
    threshold = ~np.uint64(0) // np.uint64(scaled)

    for idx in prange(n_seqs):  # type: ignore
        seq_len = lengths[idx]
        if seq_len < k:
            continue

        seq_start = offsets[idx]
        k_val_fwd = np.uint64(0)
        k_val_rev = np.uint64(0)
        valid_len = 0
        n_hashes = 0

        for j in range(seq_len):
            char = seqs[seq_start + j]
            val = np.uint64(lut[char])

            if val < max_val:
                k_val_fwd = ((k_val_fwd << np.uint64(bits_per_char)) & mask) | val
                valid_len += 1

                if canonical:
                    comp = val ^ np.uint64(2)
                    k_val_rev = (k_val_rev >> np.uint64(bits_per_char)) | (comp << np.uint64(bits_per_char * (k - 1)))

                if valid_len >= k:
                    if canonical:
                        h_val = min(k_val_fwd, k_val_rev)
                    else:
                        h_val = k_val_fwd

                    h = _splitmix64(h_val)
                    if h <= threshold:
                        n_hashes += 1
            else:
                valid_len = 0
                k_val_fwd = np.uint64(0)
                k_val_rev = np.uint64(0)

        counts[idx] = n_hashes
    return counts


@njit(parallel=True, cache=True, nogil=True)
def _populate_fracminhash_kernel(
    seqs: npt.NDArray[np.uint8],
    offsets: npt.NDArray[np.int32],
    lengths: npt.NDArray[np.int32],
    lut: npt.NDArray[np.uint8],
    k: int,
    scaled: int,
    canonical: bool,
    bits_per_char: int,
    out_offsets: npt.NDArray[np.uint32],
    out_records: npt.NDArray,
) -> None:
    r"""Populate FracMinHash records array in parallel.

    Args:
        seqs (npt.NDArray[np.uint8]): Concatenated sequence bytes.
        offsets (npt.NDArray[np.int32]): Sequence start offsets.
        lengths (npt.NDArray[np.int32]): Sequence lengths.
        lut (npt.NDArray[np.uint8]): Alphabet mapping table.
        k (int): K-mer length.
        scaled (int): Sampling scale factor.
        canonical (bool): True if computing canonical hashes.
        bits_per_char (int): Alphabet bits per character.
        out_offsets (npt.NDArray[np.uint32]): Output array write offsets.
        out_records (npt.NDArray): Output records array to populate.
    """
    n_seqs = len(offsets)
    mask = (np.uint64(1) << np.uint64(bits_per_char * k)) - np.uint64(1)
    max_val = np.uint64(1) << np.uint64(bits_per_char)
    threshold = ~np.uint64(0) // np.uint64(scaled)

    for idx in prange(n_seqs):  # type: ignore
        seq_len = lengths[idx]
        if seq_len < k:
            continue

        seq_start = offsets[idx]
        out_s = out_offsets[idx]

        k_val_fwd = np.uint64(0)
        k_val_rev = np.uint64(0)
        valid_len = 0
        written = 0

        for j in range(seq_len):
            char = seqs[seq_start + j]
            val = np.uint64(lut[char])

            if val < max_val:
                k_val_fwd = ((k_val_fwd << np.uint64(bits_per_char)) & mask) | val
                valid_len += 1

                if canonical:
                    comp = val ^ np.uint64(2)
                    k_val_rev = (k_val_rev >> np.uint64(bits_per_char)) | (comp << np.uint64(bits_per_char * (k - 1)))

                if valid_len >= k:
                    pos = j - k + 1
                    if canonical:
                        h_val = min(k_val_fwd, k_val_rev)
                    else:
                        h_val = k_val_fwd

                    h = _splitmix64(h_val)
                    if h <= threshold:
                        out_records[out_s + written]["hash"] = h
                        out_records[out_s + written]["seq_idx"] = idx
                        out_records[out_s + written]["pos1"] = pos
                        written += 1
            else:
                valid_len = 0
                k_val_fwd = np.uint64(0)
                k_val_rev = np.uint64(0)


@njit(parallel=True, cache=True, nogil=True)
def _pack_kernel(
    read_offsets: npt.NDArray[np.uint32],
    write_offsets: npt.NDArray[np.uint32],
    counts: npt.NDArray[np.uint32],
    in_records: npt.NDArray,
    out_records: npt.NDArray,
) -> None:
    r"""Pack contiguous record segments into output buffer in parallel.

    Args:
        read_offsets (npt.NDArray[np.uint32]): Read source offsets array.
        write_offsets (npt.NDArray[np.uint32]): Write destination offsets array.
        counts (npt.NDArray[np.uint32]): Item counts array.
        in_records (npt.NDArray): Input records array.
        out_records (npt.NDArray): Output records array.
    """
    for idx in prange(len(counts)):  # type: ignore
        c = counts[idx]
        if c > 0:
            r = read_offsets[idx]
            w = write_offsets[idx]
            out_records[w : w + c] = in_records[r : r + c]


@njit(parallel=True, cache=True, nogil=True)
def _count_randstrobes_kernel(
    seqs: npt.NDArray[np.uint8],
    offsets: npt.NDArray[np.int32],
    lengths: npt.NDArray[np.int32],
    lut: npt.NDArray[np.uint8],
    k: int,
    s: int,
    w_min: int,
) -> npt.NDArray[np.uint32]:
    r"""Count randstrobes per sequence in parallel.

    Args:
        seqs (npt.NDArray[np.uint8]): Concatenated sequence bytes.
        offsets (npt.NDArray[np.int32]): Sequence start offsets.
        lengths (npt.NDArray[np.int32]): Sequence lengths.
        lut (npt.NDArray[np.uint8]): Alphabet mapping table.
        k (int): K-mer length.
        s (int): Sub-k-mer window size.
        w_min (int): Minimum window offset.

    Returns:
        npt.NDArray[np.uint32]: 1D array of randstrobe counts.
    """
    n_seqs = len(offsets)
    counts = np.zeros(n_seqs, dtype=np.uint32)

    for idx in prange(n_seqs):  # type: ignore
        seq_len = lengths[idx]
        if seq_len < k:
            continue

        seq_start = offsets[idx]
        n_syncmers = 0

        # Scan for random open syncmers
        for i in range(seq_len - k + 1):
            min_s_hash = ~np.uint64(0)
            min_s_idx = -1

            # Check all s-mers within this k-mer window
            for j in range(k - s + 1):
                s_val = np.uint64(0)
                for c in range(s):
                    char = seqs[seq_start + i + j + c]
                    s_val = s_val * np.uint64(12) + np.uint64(lut[char])

                h = _splitmix64(s_val)
                if h < min_s_hash:
                    min_s_hash = h
                    min_s_idx = j

            # Open Syncmer condition: minimum s-mer is at the first or last position
            if min_s_idx == 0 or min_s_idx == (k - s):
                n_syncmers += 1

        # A sequence yields exactly (n_syncmers - w_min) randstrobes (if positive)
        if n_syncmers > w_min:
            counts[idx] = n_syncmers - w_min

    return counts


@njit(parallel=True, cache=True, nogil=True)
def _populate_randstrobes_kernel(
    seqs: npt.NDArray[np.uint8],
    offsets: npt.NDArray[np.int32],
    lengths: npt.NDArray[np.int32],
    lut: npt.NDArray[np.uint8],
    k: int,
    s: int,
    w_min: int,
    w_max: int,
    out_offsets: npt.NDArray[np.uint32],
    out_records: npt.NDArray,
) -> None:
    r"""Populate randstrobe records array in parallel.

    Args:
        seqs (npt.NDArray[np.uint8]): Concatenated sequence bytes.
        offsets (npt.NDArray[np.int32]): Sequence start offsets.
        lengths (npt.NDArray[np.int32]): Sequence lengths.
        lut (npt.NDArray[np.uint8]): Alphabet mapping table.
        k (int): K-mer length.
        s (int): Sub-k-mer window size.
        w_min (int): Minimum window offset bound.
        w_max (int): Maximum window offset bound.
        out_offsets (npt.NDArray[np.uint32]): Output array write offsets.
        out_records (npt.NDArray): Output records array to populate.
    """
    n_seqs = len(offsets)

    for idx in prange(n_seqs):  # type: ignore
        seq_len = lengths[idx]
        if seq_len < k:
            continue

        seq_start = offsets[idx]
        max_syncmers = seq_len - k + 1

        # Thread-local allocations to store localized syncmer data safely within prange
        syncmer_pos = np.empty(max_syncmers, dtype=np.uint32)
        syncmer_hash = np.empty(max_syncmers, dtype=np.uint64)
        n_syncmers = 0

        # 1. Identify and cache the syncmers
        for i in range(seq_len - k + 1):
            min_s_hash = ~np.uint64(0)
            min_s_idx = -1

            for j in range(k - s + 1):
                s_val = np.uint64(0)
                for c in range(s):
                    char = seqs[seq_start + i + j + c]
                    s_val = s_val * np.uint64(12) + np.uint64(lut[char])

                h = _splitmix64(s_val)
                if h < min_s_hash:
                    min_s_hash = h
                    min_s_idx = j

            if min_s_idx == 0 or min_s_idx == (k - s):
                syncmer_pos[n_syncmers] = i

                # Hash the full k-mer window to use for randstrobe linkage
                k_val = np.uint64(0)
                for c in range(k):
                    char = seqs[seq_start + i + c]
                    k_val = k_val * np.uint64(12) + np.uint64(lut[char])

                syncmer_hash[n_syncmers] = _splitmix64(k_val)
                n_syncmers += 1

        # 2. Build Order-2 Randstrobes from the cached syncmers
        out_idx = out_offsets[idx]

        for i in range(n_syncmers - w_min):
            h1 = syncmer_hash[i]
            min_randstrobe_hash = ~np.uint64(0)
            best_j = -1

            end_j = i + w_max + 1
            if end_j > n_syncmers:
                end_j = n_syncmers

            for j in range(i + w_min, end_j):
                h2 = syncmer_hash[j]

                # Combine hashes asymmetrically (XOR of hashed values)
                combined = _splitmix64(h1 ^ _splitmix64(h2))
                if combined < min_randstrobe_hash:
                    min_randstrobe_hash = combined
                    best_j = j

            if best_j != -1:
                out_records[out_idx]["hash"] = min_randstrobe_hash
                out_records[out_idx]["seq_idx"] = idx
                out_records[out_idx]["pos1"] = syncmer_pos[i]
                out_records[out_idx]["pos2"] = syncmer_pos[best_j]
                out_idx += 1


@njit(cache=True, nogil=True)
def _compute_query_offsets(q_records: npt.NDArray, num_queries: int) -> npt.NDArray[np.uint32]:
    r"""Find boundary indices for each query sequence in flat hash array.

    Args:
        q_records (npt.NDArray): Query records array.
        num_queries (int): Total query count.

    Returns:
        npt.NDArray[np.uint32]: Query boundary offsets array.
    """
    q_offsets = np.zeros(num_queries + 1, dtype=np.uint32)
    if len(q_records) == 0:
        return q_offsets

    current_query = 0
    for i in range(len(q_records)):
        q = q_records[i]["seq_idx"]
        while current_query < q:
            current_query += 1
            q_offsets[current_query] = i

    while current_query < num_queries:
        current_query += 1
        q_offsets[current_query] = len(q_records)

    return q_offsets


@njit(inline="always")
def _tally_single_query(
    start: int,
    end: int,
    q_records: npt.NDArray,
    t_records: npt.NDArray,
    tally: npt.NDArray[np.uint32],
    anchors: npt.NDArray[np.int32],
) -> None:
    r"""Tally matching hashes for a single query sequence against target records.

    Args:
        start (int): Start index in query records array.
        end (int): End index in query records array.
        q_records (npt.NDArray): Query records array.
        t_records (npt.NDArray): Target records array.
        tally (npt.NDArray[np.uint32]): Match count tally array.
        anchors (npt.NDArray[np.int32]): Diagonal offset anchors array.
    """
    for i in range(start, end):
        q_rec = q_records[i]
        h = q_rec["hash"]
        q_p = q_rec["pos1"]

        # Manual binary search (searchsorted_left)
        low, high = 0, len(t_records)
        while low < high:
            mid = (low + high) // 2
            if t_records[mid]["hash"] < h:
                low = mid + 1
            else:
                high = mid

        curr = low

        # Walk forward to handle all collisions (handles out-of-bounds seamlessly)
        while curr < len(t_records) and t_records[curr]["hash"] == h:
            t_rec = t_records[curr]
            t_id = t_rec["seq_idx"]
            tally[t_id] += 1
            if tally[t_id] == 1:
                anchors[t_id] = np.int32(q_p) - np.int32(t_rec["pos1"])
            curr += 1


@njit(parallel=True, cache=True, nogil=True)
def _intersect_top_hit_kernel(
    q_records: npt.NDArray,
    q_offsets: npt.NDArray[np.uint32],
    num_queries: int,
    t_records: npt.NDArray,
    num_targets: int,
) -> tuple[npt.NDArray[np.uint32], npt.NDArray[np.uint32], npt.NDArray[np.uint32], npt.NDArray[np.int32]]:
    r"""Numba parallel kernel finding top-scoring target hit per query.

    Args:
        q_records (npt.NDArray): Query records array.
        q_offsets (npt.NDArray[np.uint32]): Query boundary offsets array.
        num_queries (int): Number of queries.
        t_records (npt.NDArray): Target records array.
        num_targets (int): Number of targets.

    Returns:
        tuple[npt.NDArray[np.uint32], npt.NDArray[np.uint32], npt.NDArray[np.uint32], npt.NDArray[np.int32]]:
            Tuple of (query_indices, target_indices, scores, offsets).
    """
    best_t = np.zeros(num_queries, dtype=np.uint32)
    max_s = np.zeros(num_queries, dtype=np.uint32)
    offsets = np.zeros(num_queries, dtype=np.int32)

    for q_idx in prange(num_queries):  # type: ignore
        start, end = q_offsets[q_idx], q_offsets[q_idx + 1]
        if start == end:
            continue

        # Thread-local allocations (Fits cleanly in L1 Cache)
        tally = np.zeros(num_targets, dtype=np.uint32)
        anchors = np.zeros(num_targets, dtype=np.int32)

        # The inline directive injects the binary search code directly here
        _tally_single_query(start, end, q_records, t_records, tally, anchors)

        # Argmax Search
        best_score, best_target = 0, 0
        for t in range(num_targets):
            if tally[t] > best_score:
                best_score = tally[t]
                best_target = t

        best_t[q_idx] = best_target
        max_s[q_idx] = best_score
        offsets[q_idx] = anchors[best_target]

    # Return dense arrays, returning q_idx explicitly for easy zipping
    q_idx_arr = np.arange(num_queries, dtype=np.uint32)
    return q_idx_arr, best_t, max_s, offsets
