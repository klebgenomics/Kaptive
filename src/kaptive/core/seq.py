r"""Biological sequence data structures, SoA containers, and translation utilities.

This module provides high-performance data structures for storing, manipulating, and
translating biological sequences (DNA, RNA, protein). Sequences are stored either as individual
immutable [`SeqRecord`][kaptive.core.seq.SeqRecord] instances or as contiguous, flat arrays in
[`Sequences`][kaptive.core.seq.Sequences] Structure-of-Arrays (SoA) containers.

Key Classes:
    - [`SeqRecord`][kaptive.core.seq.SeqRecord]: Immutable container for a single sequence record.
    - [`Sequences`][kaptive.core.seq.Sequences]: Memory-efficient, Numba-accelerated SoA container.
    - [`BacterialTranslationTable`][kaptive.core.seq.BacterialTranslationTable]: NCBI Translation Table 11 translator.
"""

from __future__ import annotations

from collections.abc import Generator, Iterable
from dataclasses import dataclass
from typing import Any, Self

import numpy as np
import numpy.typing as npt
from numba import njit, prange

from kaptive.core.collections import RaggedArrayContainer
from kaptive.core.interval import Interval, IntervalLike, Intervals, Strand


# Classes --------------------------------------------------------------------------------------------------------------
@dataclass(frozen=True, slots=True)
class SeqRecord:
    r"""A simple, immutable container for a single biological sequence.

    This class pairs a string identifier with sequence data stored as immutable bytes.
    It provides basic functionality for length checking, FASTA formatting, and sub-sequence extraction.

    Attributes:
        id (str): Unique identifier or name of the sequence.
        seq (bytes): Raw sequence data (e.g. DNA, RNA, or protein).
    """

    id: str
    seq: bytes

    def __len__(self) -> int:
        r"""Return the length of the sequence in bytes.

        Returns:
            int: Sequence length in bytes.
        """
        return len(self.seq)

    def to_fasta(self) -> bytes:
        r"""Format the sequence as a FASTA record byte string.

        Returns:
            bytes: A byte string containing the FASTA header (`>id\n`) followed by
                the sequence and a trailing newline.
        """
        return b">%b\n%b\n" % (self.id.encode(), self.seq)

    def extract(self, start: int | IntervalLike, end: int | None = None, strand: Strand = Strand.UNSTRANDED) -> bytes:
        r"""Extract a sub-sequence based on coordinates and orientation.

        If `strand` is negative (e.g. [`Strand`][kaptive.core.interval.Strand]),
        the extracted sub-sequence is automatically reverse-complemented.

        Args:
            start (int | IntervalLike): 0-based start coordinate, or an
                `IntervalLike` object (providing start, end, and strand).
            end (int | None): 0-based end coordinate (exclusive). If `start` is an interval,
                this must be None. Defaults to None.
            strand (Strand): Orientation of the extraction. Defaults to
                [`Strand`][kaptive.core.interval.Strand].

        Returns:
            bytes: Extracted (and potentially reverse-complemented) sub-sequence bytes.
        """
        if end is None:
            interval = Interval.from_item(start, strand=strand)
            start_val, end_val, strand_val = interval.start, interval.end, interval.strand
        else:
            start_val, end_val, strand_val = int(start), int(end), strand  # type: ignore

        new_seq = self.seq[start_val:end_val]
        if strand_val < 0:
            return bytes(new_seq.translate(BacterialTranslationTable._COMP)[::-1])
        return bytes(new_seq)


@dataclass(frozen=True, slots=True)
class Sequences(RaggedArrayContainer["SeqRecord", "Sequences"]):
    r"""A high-performance SoA container for biological sequences.

    Stores multiple sequences in a flat, contiguous memory layout using a 1D NumPy array
    of unsigned 8-bit integers (`np.uint8`). Individual sequences are accessed via parallel arrays
    of offsets and lengths.

    Attributes:
        ids (tuple[str, ...]): String identifiers for each sequence.
        seqs (npt.NDArray[np.uint8]): Single flat 1D array containing all sequence byte data concatenated.
        offsets (npt.NDArray[np.int32]): 1D array of start indices in `seqs`.
        lengths (npt.NDArray[np.int32]): 1D array of sequence lengths.
    """

    ids: tuple[str, ...]
    seqs: npt.NDArray[np.uint8]  # concatenated sequences
    offsets: npt.NDArray[np.int32]  # start index per sequence
    lengths: npt.NDArray[np.int32]  # length per sequence

    def __len__(self) -> int:
        r"""Return the total number of sequences in the batch.

        Returns:
            int: Number of sequence records in the container.
        """
        return len(self.ids)

    def to_dict(self) -> dict[str, Any]:
        r"""Convert sequence batch to a dictionary representation suitable for serialization.

        Returns:
            dict[str, Any]: Dictionary containing 'ids', ASCII string 'seqs', 'offsets', and 'lengths'.
        """
        return {
            "ids": self.ids,
            "seqs": self.seqs.tobytes().decode("ascii"),
            "offsets": self.offsets,
            "lengths": self.lengths,
        }

    @classmethod
    def from_dict(cls, data: dict[str, Any]) -> Sequences:
        r"""Deserialize a Sequences object from a dictionary representation.

        Args:
            data (dict[str, Any]): Dictionary containing 'ids', ASCII string 'seqs', 'offsets', and 'lengths'.

        Returns:
            Sequences: Deserialized sequence collection.
        """
        return cls(
            ids=tuple(data["ids"]),
            seqs=np.frombuffer(data["seqs"].encode("ascii"), dtype=np.uint8),
            offsets=np.array(data["offsets"], dtype=np.int32),
            lengths=np.array(data["lengths"], dtype=np.int32),
        )

    def to_fasta(self, use_indices: bool = False) -> bytes:
        r"""Format the entire batch as a single FASTA byte string.

        Args:
            use_indices (bool): If True, uses 0-based integer index as FASTA header (`>0`, `>1`, ...)
                instead of string `ids`. Defaults to False.

        Returns:
            bytes: Complete FASTA-formatted byte string containing all sequences.
        """
        if not self.ids and not use_indices:
            return b""

        seq_bytes = self.seqs.tobytes()
        if use_indices:
            return b"".join(
                [
                    b">%d\n%b\n" % (i, seq_bytes[o : o + length_val])
                    for i, (o, length_val) in enumerate(zip(self.offsets.tolist(), self.lengths.tolist()))
                ]
            )
        else:
            return b"".join(
                [
                    b">%b\n%b\n" % (i.encode(), seq_bytes[o : o + length_val])
                    for i, o, length_val in zip(self.ids, self.offsets.tolist(), self.lengths.tolist())
                ]
            )

    @classmethod
    def empty(cls) -> Sequences:
        r"""Create an empty Sequences object with zero-length arrays.

        Returns:
            Sequences: An empty sequences container.
        """
        return cls(
            (),
            np.empty(0, dtype=np.uint8),
            np.empty(0, dtype=np.int32),
            np.empty(0, dtype=np.int32),
        )

    @classmethod
    def concat(cls, batches: Iterable[Self]) -> Sequences:
        r"""Concatenate multiple Sequences containers into a single larger collection.

        Args:
            batches (Iterable[Sequences]): An iterable of sequence batches to combine.

        Returns:
            Sequences: A combined sequences container.
        """
        batches_list = list(batches)
        if not batches_list:
            return cls.empty()

        all_ids = sum((b.ids for b in batches_list), ())
        all_seqs = np.concatenate([b.seqs for b in batches_list])
        all_lengths = np.concatenate([b.lengths for b in batches_list])

        offsets = np.zeros(len(all_lengths), dtype=np.int32)
        if len(all_lengths) > 1:
            np.cumsum(all_lengths[:-1], out=offsets[1:])

        return cls(all_ids, all_seqs, offsets, all_lengths)

    @property
    def internal_stops(self) -> np.ndarray:
        r"""Vectorized check for internal stop codons in protein sequences.

        Scans each sequence for the stop codon character ('*') before the final position using a Numba kernel.

        Returns:
            np.ndarray: A 1D boolean array where True indicates presence of an internal stop.
        """
        return _internal_stops_kernel(self.seqs, self.offsets, self.lengths)

    def unique(self) -> Sequences:
        r"""Return a new Sequences container containing only unique sequences.

        Uses a Numba-accelerated 64-bit FNV-1a hash to identify unique sequences while preserving
        first-occurrence order.

        Returns:
            Sequences: A new container with duplicate sequences removed.
        """
        if len(self) <= 1:
            return self

        hashes = _hash_sequences_kernel(self.seqs, self.offsets, self.lengths)
        _, unique_indices = np.unique(hashes, return_index=True)
        unique_indices.sort()

        return self[unique_indices]  # type: ignore

    def __getitem__(self, item: int | slice | np.ndarray | list) -> SeqRecord | Sequences:
        r"""Access sequences by index, slice, or boolean mask.

        Args:
            item (int | slice | np.ndarray | list): An integer index, slice, or NumPy array mask/indices.

        Returns:
            SeqRecord | Sequences: A single [`SeqRecord`][kaptive.core.seq.SeqRecord] if integer index,
                or a sliced [`Sequences`][kaptive.core.seq.Sequences] collection.

        Raises:
            IndexError: If an integer index is out of bounds.
        """
        if isinstance(item, (int, np.integer)):
            item_idx = int(item)
            if item_idx < 0:
                item_idx += len(self)
            if item_idx < 0 or item_idx >= len(self):
                raise IndexError("Batch index out of range")
            offset_val = self.offsets[item_idx]
            length_val = self.lengths[item_idx]
            return SeqRecord(self.ids[item_idx], self.seqs[offset_val : offset_val + length_val].tobytes())

        if isinstance(item, slice):
            indices = np.arange(len(self))[item]
        else:
            indices = np.asarray(item)
            if indices.dtype == bool:
                indices = np.nonzero(indices)[0]

        starts = np.zeros(len(indices), dtype=np.int32)
        ends = self.lengths[indices].astype(np.int32)
        strands = np.ones(len(indices), dtype=np.int8)
        return self.extract(
            indices.astype(np.int32), starts, ends, strands, new_ids=tuple(self.ids[i] for i in indices)
        )

    def __iter__(self) -> Generator[SeqRecord, None, None]:
        r"""Iterate over the batch, yielding a SeqRecord for each sequence.

        Yields:
            SeqRecord: Scalar record for each sequence in the batch.
        """
        for i in range(len(self)):
            offset_val = self.offsets[i]
            length_val = self.lengths[i]
            yield SeqRecord(self.ids[i], self.seqs[offset_val : offset_val + length_val].tobytes())

    @classmethod
    def from_bytes(cls, seqs: list[bytes], ids: tuple[str, ...] | None = None) -> Sequences:
        r"""Construct Sequences from a list of byte strings.

        Args:
            seqs (list[bytes]): List of raw sequence byte strings.
            ids (tuple[str, ...] | None): Sequence identifiers. Defaults to string integer indices ("0", "1", ...).

        Returns:
            Sequences: Newly constructed sequences container.
        """
        ids = ids or tuple(str(i) for i in range(len(seqs)))
        return cls.from_records([SeqRecord(i, s) for i, s in zip(ids, seqs)])

    @classmethod
    def from_records(cls, records: list[SeqRecord]) -> Sequences:
        r"""Construct Sequences from a list of SeqRecord objects.

        Args:
            records (list[SeqRecord]): List of sequence records.

        Returns:
            Sequences: Newly constructed sequences container.
        """
        ids = tuple(r.id for r in records)
        seqs = [np.frombuffer(r.seq, dtype=np.uint8) for r in records]
        if not seqs:
            return cls.empty()
        out_seqs = np.concatenate(seqs, dtype=np.uint8)
        lengths = np.array([len(s) for s in seqs], dtype=np.int32)
        offsets = np.zeros(len(seqs), dtype=np.int32)
        if len(seqs) > 1:
            np.cumsum(lengths[:-1], out=offsets[1:])
        return cls(ids, out_seqs, offsets, lengths)

    def extract(
        self,
        indices: npt.NDArray[np.int32],
        starts: npt.NDArray[np.int32],
        ends: npt.NDArray[np.int32],
        strands: npt.NDArray[np.int8],
        new_ids: tuple[str, ...] | None = None,
    ) -> Sequences:
        r"""Vectorized sub-sequence extraction from the batch.

        Handles coordinates and reverse-complementation across multiple sequences in parallel using Numba.

        Args:
            indices (npt.NDArray[np.int32]): 1D array specifying parent sequence index for each extraction.
            starts (npt.NDArray[np.int32]): 0-based start coordinates relative to parent sequence.
            ends (npt.NDArray[np.int32]): 0-based end coordinates relative to parent sequence.
            strands (npt.NDArray[np.int8]): Strand orientations (1 for forward, -1 for reverse-complement).
            new_ids (tuple[str, ...] | None): Identifiers for extracted sequences. If None, auto-generates names.

        Returns:
            Sequences: New collection of extracted sub-sequences.
        """
        if len(indices) == 0:
            return self.empty()
        new_ids = new_ids or tuple(f"{self.ids[i]}_{x}_{y}_{z}" for i, x, y, z in zip(indices, starts, ends, strands))
        out_seqs, offsets, lengths = _extract_ragged_kernel(
            self.seqs, self.offsets, indices, starts, ends, strands, BacterialTranslationTable._COMP_MAP
        )
        return Sequences(new_ids, out_seqs, offsets, lengths)

    def extract_intervals(
        self,
        indices: npt.NDArray[np.integer],
        intervals: Intervals,
        new_ids: tuple[str, ...] | None = None,
    ) -> Sequences:
        r"""Wrapper around extract taking an Intervals collection.

        Args:
            indices (npt.NDArray[np.integer]): Target sequence indices for each interval.
            intervals (Intervals): Interval collection containing coordinates and strands.
            new_ids (tuple[str, ...] | None): New sequence identifiers. Defaults to None.

        Returns:
            Sequences: New collection of extracted sub-sequences.
        """
        return self.extract(
            indices.astype(np.int32),
            intervals.starts.astype(np.int32),
            intervals.ends.astype(np.int32),
            intervals.strands,
            new_ids=new_ids,
        )

    def translate(self, frames: npt.NDArray[np.int8] | None = None, to_stop: bool = False) -> Sequences:
        r"""Vectorized translation of nucleotide sequences into protein sequences.

        Translates sequences according to NCBI Translation Table 11 (Bacterial, Archaeal, and Plant Plastid Code)
        via [`BacterialTranslationTable`][kaptive.core.seq.BacterialTranslationTable].

        Args:
            frames (npt.NDArray[np.int8] | None): Reading frame offsets (0, 1, or 2) per sequence.
                Defaults to frame 0 for all.
            to_stop (bool): If True, truncates translation at the first stop codon (*). Defaults to False.

        Returns:
            Sequences: New collection containing translated protein sequences.
        """
        if len(self) == 0:
            return self.empty()
        if frames is None:
            frames = np.zeros(len(self), dtype=np.int8)
        out_seqs, offsets, lengths = _translate_ragged_kernel(
            self.seqs,
            self.offsets,
            self.lengths,
            frames,
            BacterialTranslationTable._CHAR_MAP,
            BacterialTranslationTable._CODON_MAP,
            to_stop,
        )
        return Sequences(self.ids, out_seqs, offsets, lengths)


class BacterialTranslationTable:
    r"""NCBI Translation Table 11 utilities for bacterial codon translation and complementation.

    Provides pre-computed character map and codon lookup tables for fast 3-base codon to amino acid conversion,
    reverse-complementation mapping, and start/stop codon validation.
    """

    _MAPPING = {
        b"TTT": b"F",
        b"TTC": b"F",
        b"TTA": b"L",
        b"TTG": b"L",
        b"TCT": b"S",
        b"TCC": b"S",
        b"TCA": b"S",
        b"TCG": b"S",
        b"TAT": b"Y",
        b"TAC": b"Y",
        b"TAA": b"*",
        b"TAG": b"*",
        b"TGT": b"C",
        b"TGC": b"C",
        b"TGA": b"*",
        b"TGG": b"W",
        b"CTT": b"L",
        b"CTC": b"L",
        b"CTA": b"L",
        b"CTG": b"L",
        b"CCT": b"P",
        b"CCC": b"P",
        b"CCA": b"P",
        b"CCG": b"P",
        b"CAT": b"H",
        b"CAC": b"H",
        b"CAA": b"Q",
        b"CAG": b"Q",
        b"CGT": b"R",
        b"CGC": b"R",
        b"CGA": b"R",
        b"CGG": b"R",
        b"ATT": b"I",
        b"ATC": b"I",
        b"ATA": b"I",
        b"ATG": b"M",
        b"ACT": b"T",
        b"ACC": b"T",
        b"ACA": b"T",
        b"ACG": b"T",
        b"AAT": b"N",
        b"AAC": b"N",
        b"AAA": b"K",
        b"AAG": b"K",
        b"AGT": b"S",
        b"AGC": b"S",
        b"AGA": b"R",
        b"AGG": b"R",
        b"GTT": b"V",
        b"GTC": b"V",
        b"GTA": b"V",
        b"GTG": b"V",
        b"GCT": b"A",
        b"GCC": b"A",
        b"GCA": b"A",
        b"GCG": b"A",
        b"GAT": b"D",
        b"GAC": b"D",
        b"GAA": b"E",
        b"GAG": b"E",
        b"GGT": b"G",
        b"GGC": b"G",
        b"GGA": b"G",
        b"GGG": b"G",
    }
    _START_CODONS = {b"TTG", b"CTG", b"ATT", b"ATC", b"ATA", b"ATG", b"GTG"}
    _STOP_CODONS = {b"TAA", b"TAG", b"TGA"}
    _COMP = bytes.maketrans(b"ACGTUacgtu", b"TGCAAtgcaa")
    _CHAR_MAP = np.full(256, 4, dtype=np.uint8)
    for _i, _c in enumerate(b"ACGT"):
        _CHAR_MAP[_c] = _i
        _CHAR_MAP[_c + 32] = _i  # Map lowercase automatically
    _CHAR_MAP[b"U"[0]] = 3
    _CHAR_MAP[b"u"[0]] = 3
    _CHAR_MAP.flags.writeable = False
    # 88 = ord('X')
    _CODON_MAP = np.full(125, 88, dtype=np.uint8)
    for _codon, _aa in _MAPPING.items():
        _idx = _CHAR_MAP[_codon[0]] * 25 + _CHAR_MAP[_codon[1]] * 5 + _CHAR_MAP[_codon[2]]
        _CODON_MAP[_idx] = _aa[0]
    _CODON_MAP.flags.writeable = False
    _COMP_MAP = np.arange(256, dtype=np.uint8)
    for _c, _comp in zip(b"ACGTUacgtu", b"TGCAAtgcaa"):
        _COMP_MAP[_c] = _comp
    _COMP_MAP.flags.writeable = False

    @classmethod
    def translate(
        cls, seq: bytes | bytearray | memoryview | npt.NDArray[np.uint8], to_stop: bool = False
    ) -> npt.NDArray[np.uint8]:
        r"""Translate a nucleotide sequence array or byte string to amino acid byte values.

        Args:
            seq (bytes | bytearray | memoryview | npt.NDArray[np.uint8]): Input nucleotide sequence.
            to_stop (bool): If True, truncates translation at the first stop codon ('*'). Defaults to False.

        Returns:
            npt.NDArray[np.uint8]: 1D array of ASCII uint8 values representing translated amino acids.
        """
        if len(seq) < 3:
            return np.array([], dtype=np.uint8)

        if not isinstance(seq, np.ndarray):
            seq = np.ascontiguousarray(np.frombuffer(seq, np.uint8))

        return _translate_kernel(seq, cls._CHAR_MAP, cls._CODON_MAP, to_stop)

    @classmethod
    def is_coding(cls, seq: bytes) -> bool:
        r"""Check if a byte sequence starts with a valid bacterial start codon and ends with a stop codon.

        Args:
            seq (bytes): Nucleotide sequence to check.

        Returns:
            bool: True if sequence length >= 3 and starts/ends with valid bacterial start/stop codons.
        """
        if len(seq) < 3:
            return False
        return seq[:3] in cls._START_CODONS and seq[-3:] in cls._STOP_CODONS


# Kernels --------------------------------------------------------------------------------------------------------------
@njit(parallel=True, cache=True, nogil=True)
def _hash_sequences_kernel(
    seqs: npt.NDArray[np.uint8], offsets: npt.NDArray[np.int32], lengths: npt.NDArray[np.int32]
) -> npt.NDArray[np.uint64]:
    r"""Compute a 64-bit FNV-1a hash for each sequence in parallel.

    Args:
        seqs (npt.NDArray[np.uint8]): Concatenated sequence bytes.
        offsets (npt.NDArray[np.int32]): Sequence start offsets.
        lengths (npt.NDArray[np.int32]): Sequence lengths.

    Returns:
        npt.NDArray[np.uint64]: 1D array of 64-bit hash values.
    """
    n = len(offsets)
    hashes = np.zeros(n, dtype=np.uint64)
    for i in prange(n):  # type: ignore
        s = offsets[i]
        length_val = lengths[i]
        h = np.uint64(14695981039346656037)  # FNV offset basis
        for j in range(length_val):
            h = (h ^ np.uint64(seqs[s + j])) * np.uint64(1099511628211)  # FNV prime
        # Mix the length to further reduce collision risk
        h = (h ^ np.uint64(length_val)) * np.uint64(1099511628211)
        hashes[i] = h
    return hashes


@njit(cache=True, nogil=True, parallel=True)
def _translate_kernel(
    seq_arr: npt.NDArray[np.uint8],
    char_map: npt.NDArray[np.uint8],
    codon_map: npt.NDArray[np.uint8],
    to_stop: bool,
) -> npt.NDArray[np.uint8]:
    r"""Parallel translation kernel for a single flat nucleotide sequence array.

    Args:
        seq_arr (npt.NDArray[np.uint8]): Nucleotide sequence byte array.
        char_map (npt.NDArray[np.uint8]): Character mapping table to 0..4 base indices.
        codon_map (npt.NDArray[np.uint8]): 125-element codon to amino acid ASCII byte map.
        to_stop (bool): Whether to stop translation at the first stop codon.

    Returns:
        npt.NDArray[np.uint8]: Array of translated amino acid ASCII byte values.
    """
    max_codons = len(seq_arr) // 3
    if to_stop:
        n_codons = 0
        for i in range(max_codons):
            c1 = char_map[seq_arr[i * 3]]
            c2 = char_map[seq_arr[i * 3 + 1]]
            c3 = char_map[seq_arr[i * 3 + 2]]
            idx = c1 * 25 + c2 * 5 + c3
            if codon_map[idx] == 42:
                break
            n_codons += 1
    else:
        n_codons = max_codons

    out = np.empty(n_codons, dtype=np.uint8)
    for i in prange(n_codons):  # type: ignore
        c1 = char_map[seq_arr[i * 3]]
        c2 = char_map[seq_arr[i * 3 + 1]]
        c3 = char_map[seq_arr[i * 3 + 2]]
        idx = c1 * 25 + c2 * 5 + c3
        out[i] = codon_map[idx]
    return out


@njit(cache=True, nogil=True, parallel=True)
def _extract_ragged_kernel(
    seqs: npt.NDArray[np.uint8],
    offsets: npt.NDArray[np.int32],
    indices: npt.NDArray[np.int32],
    starts: npt.NDArray[np.int32],
    ends: npt.NDArray[np.int32],
    strands: npt.NDArray[np.int8],
    comp_map: npt.NDArray[np.uint8],
) -> tuple[npt.NDArray[np.uint8], npt.NDArray[np.int32], npt.NDArray[np.int32]]:
    r"""Parallel Numba kernel for multi-sequence sub-region extraction and reverse-complementation.

    Args:
        seqs (npt.NDArray[np.uint8]): Concatenated parent sequence array.
        offsets (npt.NDArray[np.int32]): Parent sequence start offsets.
        indices (npt.NDArray[np.int32]): Index of parent sequence for each extraction region.
        starts (npt.NDArray[np.int32]): 0-based extraction start coordinates.
        ends (npt.NDArray[np.int32]): 0-based extraction end coordinates.
        strands (npt.NDArray[np.int8]): Strand orientations (1 forward, -1 reverse).
        comp_map (npt.NDArray[np.uint8]): Complement character lookup map.

    Returns:
        tuple[npt.NDArray[np.uint8], npt.NDArray[np.int32], npt.NDArray[np.int32]]:
            Tuple of (extracted concatenated sequence bytes, new offsets, new lengths).
    """
    n = len(indices)
    out_lengths = np.empty(n, dtype=np.int32)
    total_len = 0
    for i in range(n):
        out_lengths[i] = ends[i] - starts[i]
        total_len += out_lengths[i]

    out_offsets = np.empty(n, dtype=offsets.dtype)
    if n > 0:
        out_offsets[0] = 0
        for i in range(1, n):
            out_offsets[i] = out_offsets[i - 1] + out_lengths[i - 1]

    out_seqs = np.empty(total_len, dtype=seqs.dtype)
    for i in prange(n):  # type: ignore
        length_val = out_lengths[i]
        if length_val == 0:
            continue

        t_idx = indices[i]
        g_s = offsets[t_idx] + starts[i]
        g_e = offsets[t_idx] + ends[i]
        st = strands[i]
        out_s = out_offsets[i]

        if st >= 0:
            for c in range(length_val):
                out_seqs[out_s + c] = seqs[g_s + c]
        else:
            for c in range(length_val):
                out_seqs[out_s + c] = comp_map[seqs[g_e - 1 - c]]
    return out_seqs, out_offsets, out_lengths


@njit(cache=True, nogil=True, parallel=True)
def _translate_ragged_kernel(
    seqs: npt.NDArray[np.uint8],
    offsets: npt.NDArray[np.int32],
    lengths: npt.NDArray[np.int32],
    frames: npt.NDArray[np.int8],
    char_map: npt.NDArray[np.uint8],
    codon_map: npt.NDArray[np.uint8],
    to_stop: bool,
) -> tuple[npt.NDArray[np.uint8], npt.NDArray[np.int32], npt.NDArray[np.int32]]:
    r"""Parallel Numba kernel for translating a batch of ragged nucleotide sequences.

    Args:
        seqs (npt.NDArray[np.uint8]): Concatenated parent sequence array.
        offsets (npt.NDArray[np.int32]): Sequence start offsets.
        lengths (npt.NDArray[np.int32]): Sequence lengths.
        frames (npt.NDArray[np.int8]): Reading frame offsets per sequence (0, 1, 2).
        char_map (npt.NDArray[np.uint8]): Base-to-index mapping table.
        codon_map (npt.NDArray[np.uint8]): Codon-to-AA mapping table.
        to_stop (bool): Whether to truncate translation at first stop codon per sequence.

    Returns:
        tuple[npt.NDArray[np.uint8], npt.NDArray[np.int32], npt.NDArray[np.int32]]:
            Tuple of (translated concatenated amino acid bytes, new offsets, new lengths).
    """
    n = len(offsets)
    out_lengths = np.empty(n, dtype=np.int32)
    total_len = 0
    for i in range(n):
        length_val = lengths[i]
        f = frames[i]
        if length_val > f:
            adj_len = length_val - f
            max_codons = adj_len // 3 if adj_len >= 3 else 0
            if to_stop:
                n_codons = 0
                ptr = offsets[i] + f
                for c in range(max_codons):
                    c1, c2, c3 = char_map[seqs[ptr]], char_map[seqs[ptr + 1]], char_map[seqs[ptr + 2]]
                    if codon_map[c1 * 25 + c2 * 5 + c3] == 42:
                        break
                    n_codons += 1
                    ptr += 3
                out_lengths[i] = n_codons
            else:
                out_lengths[i] = max_codons
        else:
            out_lengths[i] = 0
        total_len += out_lengths[i]

    out_offsets = np.empty(n, dtype=np.int32)
    if n > 0:
        out_offsets[0] = 0
        for i in range(1, n):
            out_offsets[i] = out_offsets[i - 1] + out_lengths[i - 1]

    out_data = np.empty(total_len, dtype=np.uint8)
    for i in prange(n):  # type: ignore
        n_codons = out_lengths[i]
        if n_codons == 0:
            continue

        s = offsets[i] + frames[i]
        out_s = out_offsets[i]

        ptr = s
        for c in range(n_codons):
            c1, c2, c3 = char_map[seqs[ptr]], char_map[seqs[ptr + 1]], char_map[seqs[ptr + 2]]
            out_data[out_s + c] = codon_map[c1 * 25 + c2 * 5 + c3]
            ptr += 3
    return out_data, out_offsets, out_lengths


@njit(cache=True, nogil=True, parallel=True)
def _internal_stops_kernel(
    seqs: npt.NDArray[np.uint8], offsets: npt.NDArray[np.int32], lengths: npt.NDArray[np.int32]
) -> npt.NDArray[np.bool_]:
    r"""Parallel Numba kernel for detecting internal stop codons in protein sequences.

    Args:
        seqs (npt.NDArray[np.uint8]): Concatenated protein sequence bytes.
        offsets (npt.NDArray[np.int32]): Sequence start offsets.
        lengths (npt.NDArray[np.int32]): Sequence lengths.

    Returns:
        npt.NDArray[np.bool_]: 1D boolean array indicating presence of internal stop codons.
    """
    n = len(offsets)
    out = np.zeros(n, dtype=np.bool_)
    for i in prange(n):  # type: ignore
        s = offsets[i]
        length_val = lengths[i]
        if length_val > 0:  # Check up to the second-to-last character (s + length_val - 1)
            for j in range(s, s + length_val - 1):
                if seqs[j] == 42:  # 42 is ord('*')
                    out[i] = True
                    break
    return out
