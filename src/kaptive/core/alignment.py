r"""High-performance Structure-of-Arrays alignment and CIGAR data structures.

This module provides optimized representations for sequence alignments and CIGAR operation
encodings. Vectorized containers such as [`Alignments`][kaptive.core.alignment.Alignments] and
[`Cigars`][kaptive.core.alignment.Cigars] store alignment fields in flat NumPy arrays, enabling
fast bulk filtering, interval conversions via [`Intervals`][kaptive.core.interval.Intervals], and
strand manipulations using [`Strand`][kaptive.core.interval.Strand].
"""

from __future__ import annotations

from collections.abc import Iterable
from dataclasses import dataclass
from enum import IntEnum
from typing import Any, NamedTuple, Self

import numpy as np
import numpy.typing as npt
from numba import njit

from kaptive.core.collections import BatchedContainer, RaggedArrayContainer
from kaptive.core.interval import Intervals, Strand


# Classes --------------------------------------------------------------------------------------------------------------
class CigarOp(IntEnum):
    r"""BAM CIGAR operation encodings as a Python Enum.

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
    I = 1  # noqa: E741  # Insertion to the reference
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
        r"""Return the single-character string representation of the CIGAR operation.

        Returns:
            str: Single character such as `'M'`, `'I'`, or `'D'`.
        """
        return "MIDNSHP=XB"[self.value]


@dataclass(frozen=True, slots=True)
class Cigars(RaggedArrayContainer[npt.NDArray[np.uint32], "Cigars"]):
    r"""A high-performance, batched container for CIGAR data using a flat memory layout.

    Instead of storing CIGAR strings or lists of tuples for each alignment, this class concatenates all
    CIGAR operations into a single, large NumPy array (`data`). This "ragged array" is managed by `offsets`
    and `lengths` arrays, which define the slice of the `data` array corresponding to each individual alignment's
    CIGAR sequence.

    Each CIGAR operation is encoded into a single 32-bit unsigned integer, following the BAM specification:
    - The upper 28 bits store the length of the operation.
    - The lower 4 bits store the operation type (corresponding to [`CigarOp`][kaptive.core.alignment.CigarOp] values).

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
        r"""Return the number of CIGAR sequences in the batch.

        Returns:
            int: Number of CIGAR sequences.
        """
        return len(self.offsets)

    def __getitem__(self, item: int | slice | npt.NDArray | list) -> npt.NDArray[np.uint32] | Cigars:
        r"""Access CIGAR data by index, slice, or boolean mask.

        - If `item` is an integer, returns a NumPy array of encoded CIGAR operations for that alignment.
        - If `item` is a slice or mask, returns a new, smaller [`Cigars`][kaptive.core.alignment.Cigars]
          containing only the selected CIGARs.

        Args:
            item (int | slice | npt.NDArray | list): An integer index, slice, or mask/indices array.

        Returns:
            npt.NDArray[np.uint32] | Cigars: A single CIGAR array or a new
                [`Cigars`][kaptive.core.alignment.Cigars] batch.

        Raises:
            IndexError: If an integer index is out of range.
        """
        if isinstance(item, (int, np.integer)):
            if item < 0:
                item += len(self)
            if item < 0 or item >= len(self):
                raise IndexError("Batch index out of range")
            offset_val, length_val = self.offsets[item], self.lengths[item]
            return self.data[offset_val : offset_val + length_val]

        if isinstance(item, slice):
            indices = np.arange(len(self))[item]
        else:
            item_arr = np.asarray(item)
            indices = np.nonzero(item_arr)[0] if item_arr.dtype.kind == "b" else item_arr

        if len(indices) == 0:
            return self.empty()

        new_lengths = self.lengths[indices]
        new_offsets = np.zeros(len(new_lengths), dtype=np.int32)
        if len(new_lengths) > 1:
            np.cumsum(new_lengths[:-1], out=new_offsets[1:])

        extracted = np.concatenate([self.data[self.offsets[i] : self.offsets[i] + self.lengths[i]] for i in indices])
        return Cigars(extracted, new_offsets, new_lengths)

    @classmethod
    def empty(cls) -> Cigars:
        r"""Create an empty Cigars instance.

        Returns:
            Cigars: An empty [`Cigars`][kaptive.core.alignment.Cigars] batch.
        """
        return cls(
            np.empty(0, dtype=np.uint32),
            np.empty(0, dtype=np.int32),
            np.empty(0, dtype=np.int32),
        )

    @classmethod
    def concat(cls, batches: Iterable[Cigars]) -> Cigars:
        r"""Concatenate multiple Cigars objects into a single, larger batch.

        Args:
            batches (Iterable[Cigars]): An iterable of [`Cigars`][kaptive.core.alignment.Cigars] batches.

        Returns:
            Cigars: A new combined [`Cigars`][kaptive.core.alignment.Cigars] batch.
        """
        batches_list = list(batches)
        if not batches_list:
            return cls.empty()
        lengths = np.concatenate([b.lengths for b in batches_list])
        offsets = np.zeros(len(lengths), dtype=np.int32)
        if len(lengths) > 1:
            np.cumsum(lengths[:-1], out=offsets[1:])
        return cls(np.concatenate([b.data for b in batches_list]), offsets, lengths)

    def swap_sides(self) -> Cigars:
        r"""Return a new Cigars batch with Insertion (I) and Deletion (D) operations swapped.

        This is used when swapping the query and target roles of an alignment.

        Returns:
            Cigars: A new [`Cigars`][kaptive.core.alignment.Cigars] batch with I and D operations swapped.
        """
        return Cigars(_swap_cigar_kernel(self.data), self.offsets, self.lengths)

    @classmethod
    def from_lists(cls, cigar_lists: list[npt.NDArray[np.uint32]]) -> Cigars:
        r"""Construct a Cigars instance from a list of individual CIGAR NumPy arrays.

        Args:
            cigar_lists (list[npt.NDArray[np.uint32]]): List of 1D uint32 NumPy arrays of CIGAR operations.

        Returns:
            Cigars: The newly constructed [`Cigars`][kaptive.core.alignment.Cigars] batch.
        """
        if not cigar_lists:
            return cls.empty()
        lengths = np.array([len(c) for c in cigar_lists], dtype=np.int32)
        offsets = np.zeros(len(lengths), dtype=np.int32)
        if len(lengths) > 1:
            np.cumsum(lengths[:-1], out=offsets[1:])
        return cls(np.concatenate(cigar_lists), offsets, lengths)


class Alignment(NamedTuple):
    r"""A lightweight, read-only view of a single alignment record's data.

    This `NamedTuple` provides a convenient interface for accessing the data of a single alignment
    within an [`Alignments`][kaptive.core.alignment.Alignments] collection.

    Attributes:
        idx (int): Original index within the source [`Alignments`][kaptive.core.alignment.Alignments].
        q_name (str): Query sequence name.
        q_length (int): Total query sequence length.
        q_start (int): Start position on the query sequence (0-based, inclusive).
        q_end (int): End position on the query sequence (0-based, exclusive).
        t_name (str): Target sequence name.
        t_length (int): Total target sequence length.
        t_start (int): Start position on the target sequence (0-based, inclusive).
        t_end (int): End position on the target sequence (0-based, exclusive).
        strand (Strand): Alignment orientation ([`Strand.FORWARD`][kaptive.core.interval.Strand] or
            [`Strand.REVERSE`][kaptive.core.interval.Strand]).
        length (int): Alignment block length.
        match (int): Number of matching bases.
        mismatch (int): Number of mismatching bases.
        score (int): Alignment score.
        quality (int): Mapping quality score (MAPQ).
        cigar (npt.NDArray[np.uint32]): 1D NumPy array of BAM-encoded CIGAR operations.
        is_primary (bool): True if primary alignment.
        is_supplementary (bool): True if supplementary alignment.
        is_spliced (bool): True if spliced alignment.
        divergence (float): Estimated sequence divergence.
        cs (bytes | None): Optional cs tag byte string.
        md (bytes | None): Optional MD tag byte string.
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


@dataclass(frozen=True, slots=True)
class Alignments(BatchedContainer[Alignment, "Alignments"]):
    r"""A high-performance, vectorized batch of alignment records.

    This class stores all data for a collection of alignments in a Structure-of-Arrays (SoA) layout
    using NumPy arrays.

    Attributes:
        q_name_ids (npt.NDArray[np.int32]): Integer indices into `q_names_dict`.
        q_names_dict (tuple[str, ...]): Unique query sequence names dictionary.
        q_lengths (npt.NDArray[np.int32]): Total query sequence lengths.
        q_starts (npt.NDArray[np.int32]): Alignment start positions on query.
        q_ends (npt.NDArray[np.int32]): Alignment end positions on query.
        t_name_ids (npt.NDArray[np.int32]): Integer indices into `t_names_dict`.
        t_names_dict (tuple[str, ...]): Unique target sequence names dictionary.
        t_lengths (npt.NDArray[np.int32]): Total target sequence lengths.
        t_starts (npt.NDArray[np.int32]): Alignment start positions on target.
        t_ends (npt.NDArray[np.int32]): Alignment end positions on target.
        strands (npt.NDArray[np.int8]): Strand orientations (+1 or -1).
        lengths (npt.NDArray[np.int32]): Alignment block lengths.
        matches (npt.NDArray[np.int32]): Number of matching bases.
        mismatches (npt.NDArray[np.int32]): Number of mismatching bases.
        scores (npt.NDArray[np.int32]): Alignment scores.
        qualities (npt.NDArray[np.uint8]): Mapping quality scores (MAPQ).
        cigars (Cigars): Batched [`Cigars`][kaptive.core.alignment.Cigars] container.
        is_primary (npt.NDArray[np.bool_]): Boolean mask of primary alignments.
        is_supplementary (npt.NDArray[np.bool_]): Boolean mask of supplementary alignments.
        is_spliced (npt.NDArray[np.bool_]): Boolean mask of spliced alignments.
        divergence (npt.NDArray[np.float64]): Sequence divergence estimates.
        cs (npt.NDArray[np.object_]): Array of cs tags.
        md (npt.NDArray[np.object_]): Array of MD tags.
    """

    q_name_ids: npt.NDArray[np.int32]
    q_names_dict: tuple[str, ...]
    q_lengths: npt.NDArray[np.int32]
    q_starts: npt.NDArray[np.int32]
    q_ends: npt.NDArray[np.int32]
    t_name_ids: npt.NDArray[np.int32]
    t_names_dict: tuple[str, ...]
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
    def q_names(self) -> npt.NDArray[np.object_]:
        r"""Return array of query sequence names for each alignment.

        Returns:
            npt.NDArray[np.object_]: 1D object array of decoded query name strings.
        """
        return np.array([self.q_names_dict[i] for i in self.q_name_ids], dtype=object)

    @property
    def t_names(self) -> npt.NDArray[np.object_]:
        r"""Return array of target sequence names for each alignment.

        Returns:
            npt.NDArray[np.object_]: 1D object array of decoded target name strings.
        """
        return np.array([self.t_names_dict[i] for i in self.t_name_ids], dtype=object)

    @property
    def q_aln_lens(self) -> npt.NDArray[np.int32]:
        r"""Return query alignment span lengths.

        Returns:
            npt.NDArray[np.int32]: Alignment spans calculated as `q_ends - q_starts`.
        """
        return self.q_ends - self.q_starts

    @property
    def t_aln_lens(self) -> npt.NDArray[np.int32]:
        r"""Return target alignment span lengths.

        Returns:
            npt.NDArray[np.int32]: Alignment spans calculated as `t_ends - t_starts`.
        """
        return self.t_ends - self.t_starts

    @property
    def q_covs(self) -> npt.NDArray[np.float64]:
        r"""Return query alignment coverage fractions.

        Returns:
            npt.NDArray[np.float64]: Coverage ratios computed as `q_aln_lens / q_lengths`.
        """
        return np.divide(
            self.q_aln_lens,
            self.q_lengths,
            out=np.zeros_like(self.q_lengths, dtype=np.float64),
            where=self.q_lengths > 0,
        )

    @property
    def t_covs(self) -> npt.NDArray[np.float64]:
        r"""Return target alignment coverage fractions.

        Returns:
            npt.NDArray[np.float64]: Coverage ratios computed as `t_aln_lens / t_lengths`.
        """
        return np.divide(
            self.t_aln_lens,
            self.t_lengths,
            out=np.zeros_like(self.t_lengths, dtype=np.float64),
            where=self.t_lengths > 0,
        )

    def __len__(self) -> int:
        r"""Return the number of alignments in the batch.

        Returns:
            int: Number of alignment records.
        """
        return len(self.q_starts)

    @classmethod
    def from_mapping_iterators(
        cls, queries: list[tuple[str, int]], iterators: Iterable[Any]
    ) -> Self:
        r"""Construct an Alignments batch from mapping iterators.

        Args:
            queries (list[tuple[str, int]]): List of query tuples `(query_name, query_length)`.
            iterators (Iterable[Any]): Iterable of mapping hit iterators (`rammappy.align.MappingIterator`).

        Returns:
            Alignments: Vectorized [`Alignments`][kaptive.core.alignment.Alignments] batch.
        """
        ql, qs, qe, tl, ts, te, st, bl, ml, nm, sc, mq = [], [], [], [], [], [], [], [], [], [], [], []
        ip, isu, isp, div, cs_list, md_list = [], [], [], [], [], []
        cigar_lists = []
        qn_ids, tn_ids = [], []
        q_names_map, t_names_map = {}, {}
        q_names_list, t_names_list = [], []

        for (q_name, q_length), it in zip(queries, iterators):
            if q_name not in q_names_map:
                q_names_map[q_name] = len(q_names_list)
                q_names_list.append(q_name)
            q_id = q_names_map[q_name]
            for h in it:
                t_name = h.target_name.decode("ascii")
                if t_name not in t_names_map:
                    t_names_map[t_name] = len(t_names_list)
                    t_names_list.append(t_name)
                t_id = t_names_map[t_name]
                qn_ids.append(q_id)
                ql.append(q_length)
                qs.append(h.query_start)
                qe.append(h.query_end)
                tn_ids.append(t_id)
                tl.append(h.target_len)
                ts.append(h.target_start)
                te.append(h.target_end)
                st.append(1 if "Forward" in repr(h.strand) else -1)
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
        if not qn_ids:
            return cls.empty()

        return cls(
            q_name_ids=np.array(qn_ids, dtype=np.int32),
            q_names_dict=tuple(q_names_list),
            q_lengths=np.array(ql, dtype=np.int32),
            q_starts=np.array(qs, dtype=np.int32),
            q_ends=np.array(qe, dtype=np.int32),
            t_name_ids=np.array(tn_ids, dtype=np.int32),
            t_names_dict=tuple(t_names_list),
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
            md=np.array(md_list, dtype=object),
        )

    @classmethod
    def concat(cls, batches: Iterable[Alignments]) -> Self:
        r"""Concatenate multiple Alignments objects into a single larger batch.

        Args:
            batches (Iterable[Alignments]): Iterable of [`Alignments`][kaptive.core.alignment.Alignments] batches.

        Returns:
            Alignments: Combined [`Alignments`][kaptive.core.alignment.Alignments] batch.

        Raises:
            ValueError: If `batches` is empty or if batch fields cannot be concatenated.
        """
        batches_list = list(batches)
        if not batches_list:
            raise ValueError("Cannot concatenate an empty iterable of batches")

        kwargs = {}

        q_names_map, q_names_list = {}, []
        t_names_map, t_names_list = {}, []
        new_q_ids, new_t_ids = [], []

        for b in batches_list:
            q_remap = np.empty(len(b.q_names_dict), dtype=np.int32)
            for i, name in enumerate(b.q_names_dict):
                if name not in q_names_map:
                    q_names_map[name] = len(q_names_list)
                    q_names_list.append(name)
                q_remap[i] = q_names_map[name]
            new_q_ids.append(q_remap[b.q_name_ids])

            t_remap = np.empty(len(b.t_names_dict), dtype=np.int32)
            for i, name in enumerate(b.t_names_dict):
                if name not in t_names_map:
                    t_names_map[name] = len(t_names_list)
                    t_names_list.append(name)
                t_remap[i] = t_names_map[name]
            new_t_ids.append(t_remap[b.t_name_ids])

        kwargs["q_name_ids"] = np.concatenate(new_q_ids)
        kwargs["q_names_dict"] = tuple(q_names_list)
        kwargs["t_name_ids"] = np.concatenate(new_t_ids)
        kwargs["t_names_dict"] = tuple(t_names_list)

        for field_name in cls.__dataclass_fields__:
            if field_name in ("q_name_ids", "q_names_dict", "t_name_ids", "t_names_dict"):
                continue
            if field_name == "cigars":
                kwargs[field_name] = Cigars.concat([b.cigars for b in batches_list])
                continue
            first_val = getattr(batches_list[0], field_name)
            if isinstance(first_val, np.ndarray):
                kwargs[field_name] = np.concatenate([getattr(b, field_name) for b in batches_list])
            else:
                if any(getattr(b, field_name) != first_val for b in batches_list):
                    raise ValueError(f"Cannot concatenate batches with mismatched '{field_name}' values")
                kwargs[field_name] = first_val

        return cls(**kwargs)

    def __getitem__(self, item: int | slice | npt.NDArray | list) -> Alignment | Alignments:
        r"""Access alignment records by index, slice, or boolean array mask.

        Args:
            item (int | slice | npt.NDArray | list): Integer index, slice, or boolean NumPy array.

        Returns:
            Alignment | Alignments: A scalar [`Alignment`][kaptive.core.alignment.Alignment] record view
                if an integer is passed, or a new filtered [`Alignments`][kaptive.core.alignment.Alignments] batch.

        Raises:
            IndexError: If an integer index is out of range.
        """
        if isinstance(item, (int, np.integer)):
            if item < 0:
                item += len(self)
            if item < 0 or item >= len(self):
                raise IndexError("Batch index out of range")
            return Alignment(
                idx=item,
                q_name=self.q_names_dict[self.q_name_ids[item]],
                q_length=self.q_lengths[item],
                q_start=self.q_starts[item],
                q_end=self.q_ends[item],
                t_name=self.t_names_dict[self.t_name_ids[item]],
                t_length=self.t_lengths[item],
                t_start=self.t_starts[item],
                t_end=self.t_ends[item],
                strand=Strand(self.strands[item]),
                length=self.lengths[item],
                match=self.matches[item],
                mismatch=self.mismatches[item],
                score=self.scores[item],
                quality=self.qualities[item],
                cigar=self.cigars[item],
                is_primary=self.is_primary[item],
                is_supplementary=self.is_supplementary[item],
                is_spliced=self.is_spliced[item],
                divergence=self.divergence[item],
                cs=self.cs[item],
                md=self.md[item],
            )

        return Alignments(
            q_name_ids=self.q_name_ids[item],
            q_names_dict=self.q_names_dict,
            q_lengths=self.q_lengths[item],
            q_starts=self.q_starts[item],
            q_ends=self.q_ends[item],
            t_name_ids=self.t_name_ids[item],
            t_names_dict=self.t_names_dict,
            t_lengths=self.t_lengths[item],
            t_starts=self.t_starts[item],
            t_ends=self.t_ends[item],
            strands=self.strands[item],
            lengths=self.lengths[item],
            matches=self.matches[item],
            mismatches=self.mismatches[item],
            scores=self.scores[item],
            qualities=self.qualities[item],
            cigars=self.cigars[item],
            is_primary=self.is_primary[item],
            is_supplementary=self.is_supplementary[item],
            is_spliced=self.is_spliced[item],
            divergence=self.divergence[item],
            cs=self.cs[item],
            md=self.md[item],
        )

    def best(self, by_query: bool = True) -> Alignments:
        r"""Return an Alignments batch containing only the best alignment per query or target.

        Selection ranks by alignment score, matches, and mapping quality.

        Args:
            by_query (bool): If True, selects the best alignment per query sequence;
                if False, selects the best per target sequence. Defaults to True.

        Returns:
            Alignments: Filtered [`Alignments`][kaptive.core.alignment.Alignments] batch.
        """
        if (n := len(self)) == 0:
            return self

        name_ints = self.q_name_ids if by_query else self.t_name_ids

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

    def cull_overlaps(
        self,
        max_overlap_fraction: float = 0.1,
        group_by: np.ndarray | None = None,
        priority_mask: np.ndarray | None = None,
        by_query: bool = True,
    ) -> Alignments:
        r"""Greedily cull alignments that overlap significantly with higher-scoring alignments.

        Args:
            max_overlap_fraction (float): Maximum allowable overlap as a fraction of alignment length.
                Defaults to 0.1.
            group_by (np.ndarray | None): Optional integer array for grouping alignments.
                Overlaps are checked only within the same group. Defaults to None.
            priority_mask (np.ndarray | None): Optional boolean mask giving priority score boosts.
                Defaults to None.
            by_query (bool): If True, checks overlaps in query coordinates; if False, target coordinates.
                Defaults to True.

        Returns:
            Alignments: Filtered, non-overlapping [`Alignments`][kaptive.core.alignment.Alignments] batch.
        """
        if (n := len(self)) < 2:
            return self

        name_ints = self.q_name_ids if by_query else self.t_name_ids
        scores = self.scores.astype(np.float64)

        if priority_mask is not None:
            scores[priority_mask] += 1e9

        # Deterministic tie-breaking for culling: Score -> Matches -> MAPQ
        order = np.lexsort((-self.qualities, -self.matches, -scores)).astype(np.int32)

        if group_by is None:
            group_by = np.zeros(n, dtype=np.int32)

        kept_mask = self.to_intervals(by_query=by_query).cull_overlaps(
            order=order,
            max_overlap_fraction=max_overlap_fraction,
            group_by=name_ints,
            secondary_group_by=group_by,
        )
        return self[kept_mask]

    def swap_sides(self) -> Alignments:
        r"""Return a new Alignments batch with query and target roles swapped.

        Returns:
            Alignments: Swapped [`Alignments`][kaptive.core.alignment.Alignments] batch.
        """
        return Alignments(
            q_name_ids=self.t_name_ids,
            q_names_dict=self.t_names_dict,
            q_lengths=self.t_lengths,
            q_starts=self.t_starts,
            q_ends=self.t_ends,
            t_name_ids=self.q_name_ids,
            t_names_dict=self.q_names_dict,
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
            md=self.md,
        )

    @classmethod
    def empty(cls) -> Alignments:
        r"""Create an empty Alignments instance.

        Returns:
            Alignments: Empty [`Alignments`][kaptive.core.alignment.Alignments] batch.
        """
        return cls(
            q_name_ids=np.empty(0, dtype=np.int32),
            q_names_dict=(),
            q_lengths=np.empty(0, dtype=np.int32),
            q_starts=np.empty(0, dtype=np.int32),
            q_ends=np.empty(0, dtype=np.int32),
            t_name_ids=np.empty(0, dtype=np.int32),
            t_names_dict=(),
            t_lengths=np.empty(0, dtype=np.int32),
            t_starts=np.empty(0, dtype=np.int32),
            t_ends=np.empty(0, dtype=np.int32),
            strands=np.empty(0, dtype=np.int8),
            lengths=np.empty(0, dtype=np.int32),
            matches=np.empty(0, dtype=np.int32),
            mismatches=np.empty(0, dtype=np.int32),
            scores=np.empty(0, dtype=np.int32),
            qualities=np.empty(0, dtype=np.uint8),
            cigars=Cigars.empty(),
            is_primary=np.empty(0, dtype=bool),
            is_supplementary=np.empty(0, dtype=bool),
            is_spliced=np.empty(0, dtype=bool),
            divergence=np.empty(0, dtype=np.float64),
            cs=np.empty(0, dtype=object),
            md=np.empty(0, dtype=object),
        )

    def to_intervals(self, by_query: bool = False) -> Intervals:
        r"""Convert alignment coordinates into an Intervals collection.

        Args:
            by_query (bool): If True, uses query coordinates (`q_starts`, `q_ends`);
                if False, uses target coordinates (`t_starts`, `t_ends`). Defaults to False.

        Returns:
            Intervals: Spatial [`Intervals`][kaptive.core.interval.Intervals] collection.
        """
        starts = self.q_starts if by_query else self.t_starts
        ends = self.q_ends if by_query else self.t_ends

        return Intervals(
            starts=starts,
            ends=ends,
            strands=self.strands,
            # CRITICAL: Ensures we can map relational queries back to this alignment record!
            original_indices=np.arange(len(self), dtype=np.int32),
        )

    def is_partial_left(self, edge_tolerance: int = 0) -> npt.NDArray[np.bool_]:
        r"""Identify alignments that hang over the left edge of the target contig.

        Args:
            edge_tolerance (int): Allowed distance from edge in base pairs. Defaults to 0.

        Returns:
            npt.NDArray[np.bool_]: Boolean mask of partial left alignments.
        """
        return (self.t_starts <= edge_tolerance) & np.where(
            self.strands == 1, self.q_starts > 0, self.q_ends < self.q_lengths
        )

    def is_partial_right(self, edge_tolerance: int = 0) -> npt.NDArray[np.bool_]:
        r"""Identify alignments that hang over the right edge of the target contig.

        Args:
            edge_tolerance (int): Allowed distance from edge in base pairs. Defaults to 0.

        Returns:
            npt.NDArray[np.bool_]: Boolean mask of partial right alignments.
        """
        return (self.t_ends >= self.t_lengths - edge_tolerance) & np.where(
            self.strands == 1, self.q_ends < self.q_lengths, self.q_starts > 0
        )

    def is_partial(self, edge_tolerance: int = 0) -> npt.NDArray[np.bool_]:
        r"""Identify alignments that hang over either edge of the target contig.

        Args:
            edge_tolerance (int): Allowed distance from edge in base pairs. Defaults to 0.

        Returns:
            npt.NDArray[np.bool_]: Boolean mask of partial alignments.
        """
        return self.is_partial_left(edge_tolerance) | self.is_partial_right(edge_tolerance)

    @classmethod
    def from_records(cls, records: Iterable[Alignment]) -> Alignments:
        r"""Construct an Alignments batch from an iterable of Alignment record objects.

        Args:
            records (Iterable[Alignment]): Iterable of scalar [`Alignment`][kaptive.core.alignment.Alignment] objects.

        Returns:
            Alignments: Newly constructed [`Alignments`][kaptive.core.alignment.Alignments] batch.
        """
        records_list = list(records)
        if not records_list:
            return cls.empty()

        q_names_map: dict[str, int] = {}
        q_names_list: list[str] = []
        qn_ids: list[int] = []

        t_names_map: dict[str, int] = {}
        t_names_list: list[str] = []
        tn_ids: list[int] = []

        for r in records_list:
            if r.q_name not in q_names_map:
                q_names_map[r.q_name] = len(q_names_list)
                q_names_list.append(r.q_name)
            qn_ids.append(q_names_map[r.q_name])

            if r.t_name not in t_names_map:
                t_names_map[r.t_name] = len(t_names_list)
                t_names_list.append(r.t_name)
            tn_ids.append(t_names_map[r.t_name])

        return cls(
            q_name_ids=np.array(qn_ids, dtype=np.int32),
            q_names_dict=tuple(q_names_list),
            q_lengths=np.array([r.q_length for r in records_list], dtype=np.int32),
            q_starts=np.array([r.q_start for r in records_list], dtype=np.int32),
            q_ends=np.array([r.q_end for r in records_list], dtype=np.int32),
            t_name_ids=np.array(tn_ids, dtype=np.int32),
            t_names_dict=tuple(t_names_list),
            t_lengths=np.array([r.t_length for r in records_list], dtype=np.int32),
            t_starts=np.array([r.t_start for r in records_list], dtype=np.int32),
            t_ends=np.array([r.t_end for r in records_list], dtype=np.int32),
            strands=np.array([r.strand for r in records_list], dtype=np.int8),
            lengths=np.array([r.length for r in records_list], dtype=np.int32),
            matches=np.array([r.match for r in records_list], dtype=np.int32),
            mismatches=np.array([r.mismatch for r in records_list], dtype=np.int32),
            scores=np.array([r.score for r in records_list], dtype=np.int32),
            qualities=np.array([r.quality for r in records_list], dtype=np.uint8),
            cigars=Cigars.from_lists([r.cigar for r in records_list]),
            is_primary=np.array([r.is_primary for r in records_list], dtype=bool),
            is_supplementary=np.array([r.is_supplementary for r in records_list], dtype=bool),
            is_spliced=np.array([r.is_spliced for r in records_list], dtype=bool),
            divergence=np.array([r.divergence for r in records_list], dtype=np.float64),
            cs=np.array([r.cs for r in records_list], dtype=object),
            md=np.array([r.md for r in records_list], dtype=object),
        )


# Kernels --------------------------------------------------------------------------------------------------------------
@njit(cache=True, nogil=True)
def parse_cigar_string(cigar_bytes: bytes) -> npt.NDArray[np.uint32]:
    r"""Fast Numba parser converting a CIGAR byte-string to a BAM-encoded uint32 array.

    Args:
        cigar_bytes (bytes): CIGAR string encoded as ASCII bytes (e.g. b"100M5D20M").

    Returns:
        npt.NDArray[np.uint32]: 1D array of BAM-encoded 32-bit CIGAR operations.
    """
    ops = 0
    # Iterate over the byte values directly
    for i in range(len(cigar_bytes)):
        char = cigar_bytes[i]
        # Check against ASCII values for valid CIGAR ops
        if (
            char == 77
            or char == 73
            or char == 68
            or char == 78
            or char == 83
            or char == 72
            or char == 80
            or char == 61
            or char == 88
            or char == 66
        ):
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
    r"""Swap Insertions (1) and Deletions (2) at the bit level in a CIGAR uint32 array.

    Args:
        cigar_data (npt.NDArray[np.uint32]): 1D array of BAM-encoded CIGAR operations.

    Returns:
        npt.NDArray[np.uint32]: New uint32 array with I and D operations swapped.
    """
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
