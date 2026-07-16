from collections.abc import Iterable
from dataclasses import dataclass
from enum import IntEnum
from typing import NamedTuple, Self

import numpy as np
import numpy.typing as npt
from numba import njit

from kaptive.core.collections import BatchedContainer, RaggedArrayContainer
from kaptive.core.interval import Intervals, Strand


# Classes --------------------------------------------------------------------------------------------------------------
class CigarOp(IntEnum):
    """BAM CIGAR operation encodings as a Python Enum.

    This class provides a standardized, integer-based representation for CIGAR (Concise Idiosyncratic
    Gapped Alignment Report) operations, which describe how an alignment is constructed from pieces
    of the query and target sequences. Using an `IntEnum` allows for both readable access (e.g., `CigarOp.M`)
    and efficient integer-based comparisons in performance-critical code.

    The values correspond to the official BAM specification.

    Attributes:
        M (0): Alignment match (can be a sequence match or mismatch).
        I (1): Insertion to the reference.
        D (2): Deletion from the reference.
        N (3): Skipped region from the reference (e.g., intron).
        S (4): Soft clipping (clipped sequences present in the sequence record).
        H (5): Hard clipping (clipped sequences NOT present in the sequence record).
        P (6): Padding (silent deletion from a padded reference).
        EQ (7): Sequence match (explicitly a match).
        X (8): Sequence mismatch (explicitly a mismatch).
        B (9): Backwards compatibility operation.
    """
    M = 0  # Alignment match (can be a sequence match or mismatch)
    I = 1  # Insertion to the reference
    D = 2  # Deletion from the reference
    N = 3  # Skipped region from the reference
    S = 4  # Soft clipping (clipped sequences present in SEQ)
    H = 5  # Hard clipping (clipped sequences NOT present in SEQ)
    P = 6  # Padding (silent deletion from padded reference)
    EQ = 7  # Sequence match
    X = 8  # Sequence mismatch
    B = 9  # Backwards compatibility

    @property
    def char(self) -> str:
        """Returns the single-character string representation of the CIGAR operation (e.g., 'M', 'I')."""
        return "MIDNSHP=XB"[self.value]


@dataclass(frozen=True, slots=True)
class Cigars(RaggedArrayContainer[npt.NDArray[np.uint32], 'Cigars']):
    """A high-performance, batched container for CIGAR data using a flat memory layout.

    Instead of storing CIGAR strings or lists of tuples for each alignment, this class concatenates all
    CIGAR operations into a single, large NumPy array (`data`). This "ragged array" is managed by `offsets`
    and `lengths` arrays, which define the slice of the `data` array corresponding to each individual alignment's
    CIGAR sequence.

    This Structure-of-Arrays (SoA) approach offers several advantages:

    1.  **Memory Efficiency:** It avoids the overhead of millions of small Python list/tuple objects.
    2.  **Cache Locality:** All CIGAR data is contiguous in memory, improving access speed.
    3.  **Vectorization:** While CIGAR operations are inherently sequential per-alignment, batch-level
        operations like swapping insertions/deletions can be vectorized.

    Each CIGAR operation is encoded into a single 32-bit unsigned integer, following the BAM specification:

    - The upper 28 bits store the length of the operation.
    - The lower 4 bits store the operation type (corresponding to `CigarOp` values).

    Attributes:
        data (npt.NDArray[np.uint32]): A 1D array containing all concatenated, BAM-encoded CIGAR operations.
        offsets (npt.NDArray[np.int32]): A 1D array where `offsets[i]` gives the starting index in `data`
            for the i-th alignment's CIGAR sequence.
        lengths (npt.NDArray[np.int32]): A 1D array where `lengths[i]` gives the number of CIGAR operations
            for the i-th alignment.
    """
    data: npt.NDArray[np.uint32]
    offsets: npt.NDArray[np.int32]
    lengths: npt.NDArray[np.int32]

    def __len__(self) -> int:
        """Returns the number of CIGAR sequences in the batch."""
        return len(self.offsets)

    def __getitem__(self, item) -> 'npt.NDArray[np.uint32] | Cigars':
        """Accesses CIGAR data by index, slice, or boolean mask.

        - If `item` is an integer, it returns a NumPy array of the encoded CIGAR operations for that single alignment.
        - If `item` is a slice or mask, it returns a new, smaller `Cigars` containing only the selected CIGARs.
          This is an expensive operation as it requires reconstructing the flat `data` array.

        Args:
            item: An integer index, slice, or boolean NumPy array.

        Returns:
            npt.NDArray[np.uint32] | Cigars: A single CIGAR array or a new batch.
        """
        if isinstance(item, (int, np.integer)):
            if item < 0:
                item += len(self)
            if item < 0 or item >= len(self):
                raise IndexError("Batch index out of range")
            s, l = self.offsets[item], self.lengths[item]
            return self.data[s:s + l]

        if isinstance(item, slice):
            indices = np.arange(len(self))[item]
        else:
            item_arr = np.asarray(item)
            indices = np.nonzero(item_arr)[0] if item_arr.dtype.kind == 'b' else item_arr

        if len(indices) == 0:
            return self.empty()

        new_lengths = self.lengths[indices]
        new_offsets = np.zeros(len(new_lengths), dtype=np.int32)
        if len(new_lengths) > 1:
            np.cumsum(new_lengths[:-1], out=new_offsets[1:])

        extracted = np.concatenate(
            [self.data[self.offsets[i]:self.offsets[i] + self.lengths[i]] for i in indices])
        return Cigars(extracted, new_offsets, new_lengths)

    @classmethod
    def empty(cls) -> 'Cigars':
        """Creates an empty Cigars."""
        return cls(np.empty(0, dtype=np.uint32),
                   np.empty(0, dtype=np.int32),
                   np.empty(0, dtype=np.int32))

    @classmethod
    def concat(cls, batches: Iterable['Cigars']) -> 'Cigars':
        """Concatenates multiple Cigars objects into a single, larger batch.

        Args:
            batches (Iterable[Cigars]): A list or iterable of batches to concatenate.

        Returns:
            Cigars: A new, combined batch.
        """
        batches = list(batches)
        if not batches:
            return cls.empty()
        lengths = np.concatenate([b.lengths for b in batches])
        offsets = np.zeros(len(lengths), dtype=np.int32)
        if len(lengths) > 1:
            np.cumsum(lengths[:-1], out=offsets[1:])
        return cls(np.concatenate([b.data for b in batches]), offsets, lengths)

    def swap_sides(self) -> 'Cigars':
        """Returns a new Cigars with Insertion (I) and Deletion (D) operations swapped.

        This is used when swapping the query and target roles of an alignment. An insertion relative
        to the original target becomes a deletion relative to the new target (the original query),
        and vice-versa. The operation is performed by a fast, Numba-jitted kernel that manipulates
        the CIGAR encoding at the bit level.

        Returns:
            Cigars: A new batch with I and D operations swapped.
        """
        return Cigars(_swap_cigar_kernel(self.data), self.offsets, self.lengths)

    @classmethod
    def from_lists(cls, cigar_lists: list[npt.NDArray[np.uint32]]) -> 'Cigars':
        """Constructs a Cigars from a list of individual CIGAR NumPy arrays.

        This factory method handles the core logic of flattening the list of arrays into the single
        `data` array and calculating the corresponding `offsets` and `lengths`.

        Args:
            cigar_lists (list[npt.NDArray[np.uint32]]): A list where each element is a NumPy array
                of encoded CIGAR operations for one alignment.

        Returns:
            Cigars: The newly constructed batch.
        """
        if not cigar_lists:
            return cls.empty()
        lengths = np.array([len(c) for c in cigar_lists], dtype=np.int32)
        offsets = np.zeros(len(lengths), dtype=np.int32)
        if len(lengths) > 1:
            np.cumsum(lengths[:-1], out=offsets[1:])
        return cls(np.concatenate(cigar_lists), offsets, lengths)


class Alignment(NamedTuple):
    """A lightweight, read-only view of a single alignment record's data.

    This `NamedTuple` provides a convenient, object-like interface for accessing the data of a single
    alignment within an `Alignments`. It is created on-demand when a single item is requested
    from the batch (e.g., `my_batch[0]`). It holds references to the data but does not duplicate it,
    making it efficient for inspection and debugging.

    Attributes:
        idx (int): The original index of the record within its source `Alignments`.
        q_name (str): The name of the query sequence.
        q_length (int): The total length of the query sequence.
        q_start (int): The start position of the alignment on the query sequence (0-based, inclusive).
        q_end (int): The end position of the alignment on the query sequence (0-based, exclusive).
        t_name (str): The name of the target (reference) sequence.
        t_length (int): The total length of the target sequence.
        t_start (int): The start position of the alignment on the target sequence (0-based, inclusive).
        t_end (int): The end position of the alignment on the target sequence (0-based, exclusive).
        strand (Strand): The alignment orientation (`Strand.FORWARD` or `Strand.REVERSE`).
        length (int): The total length of the alignment block.
        match (int): The number of matching bases in the alignment.
        mismatch (int): The number of mismatches (from the 'NM' tag in a SAM/BAM file).
        quality (int): The mapping quality (MAPQ) score, indicating confidence in the alignment's location.
        cigar (npt.NDArray[np.uint32]): The NumPy array of BAM-encoded CIGAR operations for this alignment.
    """
    idx: int
    q_name: str
    q_length: int
    q_start: int
    q_end: int
    t_name: str
    t_length: int
    t_start: int
    t_end: int
    strand: Strand
    length: int
    match: int
    mismatch: int
    score: int
    quality: int
    cigar: npt.NDArray[np.uint32]
    is_primary: bool
    is_supplementary: bool
    is_spliced: bool
    divergence: float
    cs: bytes | None
    md: bytes | None

    from rammappy.align import Mapping


@dataclass(frozen=True, slots=True)
class Alignments(BatchedContainer[Alignment, 'Alignments']):
    """A high-performance, vectorized batch of alignment records.

    This class stores all data for a large collection of alignments in a Structure-of-Arrays (SoA)
    layout using NumPy arrays. This design is central to the performance of the library, enabling:

    -   **Vectorized Filtering:** Applying filters (e.g., by score, length, or location) to hundreds of
        thousands of alignments simultaneously using boolean masks.
    -   **Efficient Transformations:** Modifying all alignment coordinates at once (e.g., `shift_query`, `swap_sides`).
    -   **Reduced Memory Footprint:** Minimizing the overhead associated with large numbers of Python objects.

    The class is immutable. All methods that filter or modify the data return a new `Alignments` instance.

    Attributes:
        q_names (npt.NDArray[np.object_]): Array of query sequence names.
        q_lengths (npt.NDArray[np.int32]): Array of total query sequence lengths.
        q_starts (npt.NDArray[np.int32]): Array of alignment start positions on the query.
        q_ends (npt.NDArray[np.int32]): Array of alignment end positions on the query.
        t_names (npt.NDArray[np.object_]): Array of target sequence names.
        t_lengths (npt.NDArray[np.int32]): Array of total target sequence lengths.
        t_starts (npt.NDArray[np.int32]): Array of alignment start positions on the target.
        t_ends (npt.NDArray[np.int32]): Array of alignment end positions on the target.
        strands (npt.NDArray[np.int8]): Array of strand orientations.
        lengths (npt.NDArray[np.int32]): Array of alignment block lengths.
        matches (npt.NDArray[np.int32]): Array of the number of matching bases.
        mismatches (npt.NDArray[np.int32]): Array of the number of mismatches.
        scores (npt.NDArray[np.int32]): Array of alignment scores.
        qualities (npt.NDArray[np.uint8]): Array of mapping qualities (MAPQ).
        cigars (Cigars): The batched CIGAR data associated with these alignments.
    """
    q_names: npt.NDArray[np.object_]
    q_lengths: npt.NDArray[np.int32]
    q_starts: npt.NDArray[np.int32]
    q_ends: npt.NDArray[np.int32]
    t_names: npt.NDArray[np.object_]
    t_lengths: npt.NDArray[np.int32]
    t_starts: npt.NDArray[np.int32]
    t_ends: npt.NDArray[np.int32]
    strands: npt.NDArray[np.int8]
    lengths: npt.NDArray[np.int32]
    matches: npt.NDArray[np.int32]
    mismatches: npt.NDArray[np.int32]
    scores: npt.NDArray[np.int32]
    qualities: npt.NDArray[np.uint8]
    cigars: Cigars
    is_primary: npt.NDArray[np.bool_]
    is_supplementary: npt.NDArray[np.bool_]
    is_spliced: npt.NDArray[np.bool_]
    divergence: npt.NDArray[np.float64]
    cs: npt.NDArray[np.object_]
    md: npt.NDArray[np.object_]

    @property
    def q_aln_lens(self):
        return self.q_ends - self.q_starts

    @property
    def t_aln_lens(self):
        return self.t_ends - self.t_starts

    @property
    def q_covs(self):
        return np.divide(
            self.q_aln_lens, self.q_lengths,
            out=np.zeros_like(self.q_lengths, dtype=np.float64),
            where=self.q_lengths > 0
        )

    @property
    def t_covs(self):
        return np.divide(
            self.t_aln_lens, self.t_lengths,
            out=np.zeros_like(self.t_lengths, dtype=np.float64),
            where=self.t_lengths > 0
        )

    @property
    def q_weighted_scores(self):
        return self.scores * self.q_covs

    @property
    def t_weighted_scores(self):
        return self.scores * self.t_covs

    @property
    def minus_edit_distances(self):
        return -self.mismatches

    @property
    def q_weighted_minus_edit_distances(self):
        return self.minus_edit_distances * self.q_covs

    @property
    def minus_divergences(self):
        return -self.divergence

    @property
    def q_weighted_minus_divergences(self):
        return self.minus_divergences * self.q_covs

    @property
    def q_weighted_matches(self):
        return self.matches * self.q_covs

    @property
    def identities(self):
        return np.divide(
            self.matches, self.lengths,
            out=np.zeros_like(self.lengths, dtype=np.float64),
            where=self.lengths > 0
        )

    @property
    def q_weighted_identities(self):
        return self.identities * self.q_covs

    def __len__(self) -> int:
        """Returns the number of alignments in the batch."""
        return len(self.q_starts)

    @classmethod
    def from_mapping_iterators(cls, queries: list[tuple[str, int]], iterators: Iterable['rammappy.align.MappingIterator']) -> Self:
        """A factory method to create an Alignments from mapping iterators."""
        qn, ql, qs, qe, tn, tl, ts, te, st, bl, ml, nm, sc, mq = [], [], [], [], [], [], [], [], [], [], [], [], [], []
        ip, isu, isp, div, cs_list, md_list = [], [], [], [], [], []
        cigar_lists = []

        for (q_name, q_length), it in zip(queries, iterators):
            for h in it:
                qn.append(q_name)
                ql.append(q_length)
                qs.append(h.query_start)
                qe.append(h.query_end)
                tn.append(h.target_name.decode('ascii'))
                tl.append(h.target_len)
                ts.append(h.target_start)
                te.append(h.target_end)
                st.append(1 if 'Forward' in repr(h.strand) else -1)
                bl.append(h.block_len)
                ml.append(h.matches)
                nm.append(h.edit_distance)
                sc.append(h.score)
                mq.append(h.mapq)
                
                ip.append(h.is_primary)
                isu.append(h.is_supplementary)
                isp.append(h.is_spliced)
                div.append(h.divergence)
                cs_list.append(h.cs)
                md_list.append(h.md)
                
                cigar_bytes = h.cigar
                if cigar_bytes:
                    cigar_lists.append(parse_cigar_string(cigar_bytes))
                else:
                    cigar_lists.append(np.empty(0, dtype=np.uint32))

        if not qn:
            return cls.empty()

        return cls(
            q_names=np.array(qn, dtype=object),
            q_lengths=np.array(ql, dtype=np.int32),
            q_starts=np.array(qs, dtype=np.int32),
            q_ends=np.array(qe, dtype=np.int32),
            t_names=np.array(tn, dtype=object),
            t_lengths=np.array(tl, dtype=np.int32),
            t_starts=np.array(ts, dtype=np.int32),
            t_ends=np.array(te, dtype=np.int32),
            strands=np.array(st, dtype=np.int8),
            lengths=np.array(bl, dtype=np.int32),
            matches=np.array(ml, dtype=np.int32),
            mismatches=np.array(nm, dtype=np.int32),
            scores=np.array(sc, dtype=np.int32),
            qualities=np.array(mq, dtype=np.uint8),
            cigars=Cigars.from_lists(cigar_lists),
            is_primary=np.array(ip, dtype=bool),
            is_supplementary=np.array(isu, dtype=bool),
            is_spliced=np.array(isp, dtype=bool),
            divergence=np.array(div, dtype=np.float64),
            cs=np.array(cs_list, dtype=object),
            md=np.array(md_list, dtype=object)
        )

    @classmethod
    def concat(cls, batches: Iterable['Alignments']) -> Self:
        """Concatenates multiple Alignments objects into a single, larger batch.

        Args:
            batches (Iterable[Alignments]): A list or iterable of batches to concatenate.

        Returns:
            Alignments: A new, combined batch.
        """
        batches = list(batches)
        if not batches:
            raise ValueError("Cannot concatenate an empty iterable of batches")

        kwargs = {}
        for field_name in cls.__dataclass_fields__:
            if field_name == 'cigars':
                kwargs[field_name] = Cigars.concat([b.cigars for b in batches])
                continue
            first_val = getattr(batches[0], field_name)
            if isinstance(first_val, np.ndarray):
                kwargs[field_name] = np.concatenate([getattr(b, field_name) for b in batches])
            else:
                if any(getattr(b, field_name) != first_val for b in batches):
                    raise ValueError(f"Cannot concatenate batches with mismatched '{field_name}' values")
                kwargs[field_name] = first_val

        return cls(**kwargs)

    def __getitem__(self, item) -> 'Alignment | Alignments':
        """Accesses alignment records by index, slice, or boolean mask.

        - If `item` is an integer, it returns a single, lightweight `AlignmentRecord` view.
        - If `item` is a slice or a boolean mask, it returns a new `Alignments` containing
          only the selected records.

        Args:
            item: An integer index, slice, or boolean NumPy array.

        Returns:
            Alignment | Alignments: A single record view or a new, filtered batch.
        """
        if isinstance(item, (int, np.integer)):
            if item < 0:
                item += len(self)
            if item < 0 or item >= len(self):
                raise IndexError("Batch index out of range")
            return Alignment(
                idx=item, q_name=self.q_names[item], q_length=self.q_lengths[item],
                q_start=self.q_starts[item], q_end=self.q_ends[item], t_name=self.t_names[item],
                t_length=self.t_lengths[item], t_start=self.t_starts[item], t_end=self.t_ends[item],
                strand=Strand(self.strands[item]), length=self.lengths[item], match=self.matches[item],
                mismatch=self.mismatches[item], score=self.scores[item], quality=self.qualities[item], cigar=self.cigars[item],
                is_primary=self.is_primary[item], is_supplementary=self.is_supplementary[item],
                is_spliced=self.is_spliced[item], divergence=self.divergence[item],
                cs=self.cs[item], md=self.md[item]
            )
            
        return Alignments(
            q_names=self.q_names[item], q_lengths=self.q_lengths[item], q_starts=self.q_starts[item],
            q_ends=self.q_ends[item], t_names=self.t_names[item], t_lengths=self.t_lengths[item],
            t_starts=self.t_starts[item], t_ends=self.t_ends[item], strands=self.strands[item],
            lengths=self.lengths[item], matches=self.matches[item], mismatches=self.mismatches[item],
            scores=self.scores[item], qualities=self.qualities[item], cigars=self.cigars[item],
            is_primary=self.is_primary[item], is_supplementary=self.is_supplementary[item],
            is_spliced=self.is_spliced[item], divergence=self.divergence[item],
            cs=self.cs[item], md=self.md[item]
        )

    def best(self, by_query: bool = True) -> 'Alignments':
        """Returns Alignments with the best alignment per query or target."""
        if (n := len(self)) == 0:
            return self

        names = self.q_names if by_query else self.t_names

        _, name_ints = np.unique(names, return_inverse=True)

        # np.lexsort sorts by the last key first.
        # Primary: group by name
        # Secondary: highest weighted score
        # Tertiary: most matching bases (tie-breaker 1)
        # Quaternary: highest MAPQ (tie-breaker 2)
        order = np.lexsort((-self.qualities, -self.matches, -self.scores, name_ints))

        # Extract the first occurrence of each name using a highly optimized O(N) boundary mask
        name_sorted = name_ints[order]
        first_occurrence_mask = np.empty(n, dtype=bool)
        first_occurrence_mask[0] = True
        first_occurrence_mask[1:] = name_sorted[1:] != name_sorted[:-1]
        
        best_indices = order[first_occurrence_mask]

        # Sort indices to maintain the original relative order from the batch
        best_indices.sort()

        return self[best_indices]

    def cull_overlaps(self, max_overlap_fraction: float = 0.1, group_by: np.ndarray | None = None,
                      priority_mask: np.ndarray | None = None, by_query: bool = True) -> "Alignments":
        """Greedily culls alignments that overlap significantly with higher-scoring alignments.

        This method is used to resolve redundant or ambiguous mappings. It sorts all alignments by score
        (and optionally by a priority mask) in descending order. It then iterates through them, keeping an
        alignment only if it does not overlap "too much" with any already-kept, higher-scoring alignment.
        The core logic is implemented in a fast Numba-jitted kernel.

        Args:
            max_overlap_fraction (float, optional): The maximum allowable overlap as a fraction of an
                alignment's own length. If the overlap with a better alignment exceeds this, the alignment
                is culled. Defaults to 0.1.
            group_by (np.ndarray, optional): An integer array for grouping. Overlap checks will only be
                performed between alignments within the same group. This is useful for culling overlaps
                per-contig or per-gene. Defaults to None (all alignments in one group).
            priority_mask (np.ndarray, optional): A boolean mask to give certain alignments an artificial
                score boost, ensuring they are considered first in the greedy selection, regardless of their
                actual alignment score. Defaults to None.
            by_query (bool, optional): If True, culls based on overlaps in the query coordinates (the assembly).
                If False, culls based on the target coordinates (the database). Defaults to True, ensuring
                a physical piece of DNA is only assigned to one gene.

        Returns:
            Alignments: A new batch containing the filtered, non-redundant set of alignments.
        """
        if (n := len(self)) < 2:
            return self
        
        names = self.q_names if by_query else self.t_names
        starts = self.q_starts if by_query else self.t_starts
        ends = self.q_ends if by_query else self.t_ends
        scores = self.q_weighted_scores if by_query else self.t_weighted_scores
        
        # Map string names to integer IDs for C-level kernel processing
        _, name_ints = np.unique(names, return_inverse=True)

        if priority_mask is not None:
            scores[priority_mask] += 1e9
            
        # Deterministic tie-breaking for culling: Score -> Matches -> MAPQ
        order = np.lexsort((-self.qualities, -self.matches, -scores)).astype(np.int32)
        
        if group_by is None:
            group_by = np.zeros(n, dtype=np.int32)
        
        kept_mask = self.to_intervals(by_query=by_query).cull_overlaps(
            order=order, max_overlap_fraction=max_overlap_fraction,
            group_by=name_ints, secondary_group_by=group_by
        )
        return self[kept_mask]

    def swap_sides(self) -> 'Alignments':
        """Returns a new batch with the roles of query and target swapped.

        This operation effectively flips the alignment. The original query becomes the new target,
        and the original target becomes the new query. All corresponding data fields (names, lengths,
        coordinates) are swapped. CIGAR strings are also adjusted by swapping Insertion and Deletion
        operations via `Cigars.swap_sides`.

        Returns:
            Alignments: A new batch representing the swapped alignments.
        """
        return Alignments(
            q_names=self.t_names,
            q_lengths=self.t_lengths,
            q_starts=self.t_starts,
            q_ends=self.t_ends,
            t_names=self.q_names,
            t_lengths=self.q_lengths,
            t_starts=self.q_starts,
            t_ends=self.q_ends,
            strands=self.strands,
            lengths=self.lengths,
            matches=self.matches,
            mismatches=self.mismatches,
            scores=self.scores,
            qualities=self.qualities,
            cigars=self.cigars.swap_sides(),
            is_primary=self.is_primary,
            is_supplementary=self.is_supplementary,
            is_spliced=self.is_spliced,
            divergence=self.divergence,
            cs=self.cs,
            md=self.md
        )
        
    @classmethod
    def empty(cls) -> 'Alignments':
        return cls(
            q_names=np.empty(0, dtype=object), q_lengths=np.empty(0, dtype=np.int32),
            q_starts=np.empty(0, dtype=np.int32), q_ends=np.empty(0, dtype=np.int32),
            t_names=np.empty(0, dtype=object), t_lengths=np.empty(0, dtype=np.int32),
            t_starts=np.empty(0, dtype=np.int32), t_ends=np.empty(0, dtype=np.int32),
            strands=np.empty(0, dtype=np.int8), lengths=np.empty(0, dtype=np.int32),
            matches=np.empty(0, dtype=np.int32), mismatches=np.empty(0, dtype=np.int32),
            scores=np.empty(0, dtype=np.int32), qualities=np.empty(0, dtype=np.int8), cigars=Cigars.empty(),
            is_primary=np.empty(0, dtype=bool), is_supplementary=np.empty(0, dtype=bool),
            is_spliced=np.empty(0, dtype=bool), divergence=np.empty(0, dtype=np.float64),
            cs=np.empty(0, dtype=object), md=np.empty(0, dtype=object)
        )

    def to_intervals(self, by_query: bool = False) -> Intervals:
        """Converts the alignment coordinates into a high-performance `Intervals`.

        This is a key method for performing spatial queries on alignments. It extracts the start,
        end, and strand information into an `Intervals`, which has highly optimized methods
        for finding overlaps, nearest neighbors, etc. The `original_indices` of the `Intervals`
        are critically preserved as a direct mapping back to the indices of this `Alignments`.

        Args:
            by_query (bool, optional): If True, creates intervals from the query coordinates (`q_starts`, `q_ends`).
                If False, uses the target coordinates. Defaults to False.

        Returns:
            Intervals: An `Intervals` representing the genomic locations of the alignments.
        """
        starts = self.q_starts if by_query else self.t_starts
        ends = self.q_ends if by_query else self.t_ends

        return Intervals(
            starts=starts,
            ends=ends,
            strands=self.strands,
            # CRITICAL: Ensures we can map relational queries back to this alignment record!
            original_indices=np.arange(len(self), dtype=np.int32)
        )
    
    @property
    def is_partial_left(self) -> np.ndarray:
        return (self.t_starts == 0) & np.where(self.strands == 1, self.q_starts > 0, self.q_ends < self.q_lengths)
    
    @property
    def is_partial_right(self):
        return (self.t_ends == self.t_lengths) & np.where(
            self.strands == 1, self.q_ends < self.q_lengths, self.q_starts > 0
        )

    @property
    def is_partial(self) -> np.ndarray:
        """Creates a vectorized boolean mask of alignments where the query start/end hangs over the target start/end.

        Returns:
            np.ndarray: A 1D boolean array where `True` indicates a partial alignment.
        """
        return self.is_partial_left | self.is_partial_right

    @classmethod
    def from_records(cls, records: Iterable[Alignment]) -> 'Alignments':
        """Constructs an Alignments from an iterable of `AlignmentRecord` objects.

        This factory method is the inverse of `__getitem__` and `split`. It efficiently reconstructs
        the Structure-of-Arrays (SoA) layout from a list of record objects.

        Args:
            records (Iterable[Alignment]): An iterable of `AlignmentRecord` objects.

        Returns:
            Alignments: The newly constructed batch.
        """
        data = [
            (r.q_name, r.q_length, r.q_start, r.q_end, r.t_name, r.t_length,
             r.t_start, r.t_end, r.strand, r.length, r.match, r.mismatch,
             r.quality, r.is_primary, r.is_supplementary, r.is_spliced, r.divergence, r.cs, r.md)
            for r in records

        ]
        if not data:
            return cls.empty()
            
        qn, ql, qs, qe, tn, tl, ts, te, st, bl, ml, nm, mq, ip, isu, isp, div, cs_list, md_list = zip(*data, strict=True)
        return cls(
            q_names=np.array(qn, dtype=object),
            q_lengths=np.array(ql, dtype=np.int32),
            q_starts=np.array(qs, dtype=np.int32),
            q_ends=np.array(qe, dtype=np.int32),
            t_names=np.array(tn, dtype=object),
            t_lengths=np.array(tl, dtype=np.int32),
            t_starts=np.array(ts, dtype=np.int32),
            t_ends=np.array(te, dtype=np.int32),
            strands=np.array(st, dtype=np.int8),
            lengths=np.array(bl, dtype=np.int32),
            matches=np.array(ml, dtype=np.int32),
            mismatches=np.array(nm, dtype=np.int32),
            qualities=np.array(mq, dtype=np.uint8),
            cigars=Cigars.from_lists([r.cigar for r in records]),
            is_primary=np.array(ip, dtype=bool),
            is_supplementary=np.array(isu, dtype=bool),
            is_spliced=np.array(isp, dtype=bool),
            divergence=np.array(div, dtype=np.float64),
            cs=np.array(cs_list, dtype=object),
            md=np.array(md_list, dtype=object)
        )


# Kernels --------------------------------------------------------------------------------------------------------------
@njit(cache=True, nogil=True)
def parse_cigar_string(cigar_bytes: bytes) -> npt.NDArray[np.uint32]:
    """Fast Numba parser converting a CIGAR byte-string to a BAM-encoded uint32 array."""
    ops = 0
    # Iterate over the byte values directly
    for i in range(len(cigar_bytes)):
        char = cigar_bytes[i]
        # Check against ASCII values for valid CIGAR ops
        if char == 77 or char == 73 or char == 68 or char == 78 or char == 83 or char == 72 or char == 80 or char == 61 or char == 88 or char == 66:
            ops += 1

    out = np.empty(ops, dtype=np.uint32)
    idx = 0
    current_len = 0

    for i in range(len(cigar_bytes)):
        char = cigar_bytes[i]
        if 48 <= char <= 57:  # ASCII '0' to '9'
            current_len = current_len * 10 + (char - 48)
        else:
            op = 0
            if char == 77:  # 'M'
                op = 0
            elif char == 73:  # 'I'
                op = 1
            elif char == 68:  # 'D'
                op = 2
            elif char == 78:  # 'N'
                op = 3
            elif char == 83:  # 'S'
                op = 4
            elif char == 72:  # 'H'
                op = 5
            elif char == 80:  # 'P'
                op = 6
            elif char == 61:  # '='
                op = 7
            elif char == 88:  # 'X'
                op = 8
            elif char == 66:  # 'B'
                op = 9
            else:
                continue

            out[idx] = (np.uint32(current_len) << 4) | np.uint32(op)
            idx += 1
            current_len = 0

    return out


@njit(cache=True, nogil=True)
def _swap_cigar_kernel(cigar_data: npt.NDArray[np.uint32]) -> npt.NDArray[np.uint32]:
    """Swaps Insertions (1) with Deletions (2) at the bit-level."""
    out = np.empty_like(cigar_data)
    for i in range(len(cigar_data)):
        val = cigar_data[i]
        op = val & 0xF
        if op == 1:
            out[i] = (val & ~np.uint32(0xF)) | np.uint32(2)
        elif op == 2:
            out[i] = (val & ~np.uint32(0xF)) | np.uint32(1)
        else:
            out[i] = val
    return out
