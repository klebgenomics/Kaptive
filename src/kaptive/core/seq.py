from dataclasses import dataclass
from typing import Generator, Iterable

import numpy as np
import numpy.typing as npt
from numba import njit, prange

from kaptive.core.interval import Strand, Interval, IntervalLike, Intervals


# Classes --------------------------------------------------------------------------------------------------------------
@dataclass(frozen=True, slots=True)
class SeqRecord:
    """A simple, immutable container for a single biological sequence.

    This class pairs a string identifier with sequence data stored as immutable bytes.
    It provides basic functionality for length checking, formatting, and sub-sequence extraction.
    For high-performance operations on many sequences, `Sequences` should be preferred.

    Attributes:
        id (str): The unique identifier or name of the sequence.
        seq (bytes): The raw sequence data (e.g., DNA, RNA, or protein).
    """
    id: str
    seq: bytes

    def __len__(self) -> int:
        """Returns the length of the sequence in bytes."""
        return len(self.seq)

    def to_fasta(self) -> bytes:
        """Formats the sequence as a FASTA record.

        Returns:
            bytes: A byte string containing the FASTA header (`>id\\n`) followed by
                the sequence and a trailing newline.
        """
        return b">%b\n%b\n" % (self.id.encode(), self.seq)

    def extract(self, start: int | IntervalLike, end: int | None = None, strand: Strand = Strand.UNSTRANDED) -> bytes:
        """Extracts a sub-sequence based on coordinates and orientation.

        If the `strand` is set to `Strand.REVERSE` (-1), the extracted sub-sequence
        is automatically reverse-complemented using standard DNA/RNA base pairing.

        Args:
            start (int | IntervalLike): The 0-based start coordinate, or an `IntervalLike` object
                (which provides start, end, and strand).
            end (int | None, optional): The 0-based end coordinate (exclusive). If `start` is an
                interval, this should be None. Defaults to None.
            strand (Strand, optional): The orientation of the extraction. Defaults to UNSTRANDED.

        Returns:
            bytes: The extracted (and potentially reverse-complemented) sub-sequence.
        """
        # Ensure we are working with an immutable bytes object for consistent return types.
        # This avoids returning a list[bytearray] if a bytearray is passed in.
        if end is None:
            interval = Interval.from_item(start, strand=strand)
            start_val, end_val, strand_val = interval.start, interval.end, interval.strand
        else:
            start_val, end_val, strand_val = int(start), int(end), strand
            
        new_seq = self.seq[start_val:end_val]
        if strand_val < 0:
            return bytes(new_seq.translate(BacterialTranslationTable._COMP)[::-1])
        return bytes(new_seq)


@dataclass(frozen=True, slots=True)
class Sequences:
    """A high-performance SoA container for biological sequences.

    This class stores multiple sequences in a flat, contiguous memory layout using a 1D NumPy array
    of unsigned 8-bit integers (`np.uint8`). Individual sequences are accessed via parallel arrays
    of offsets and lengths. This Structure-of-Arrays (SoA) layout enables:

    1.  **Vectorized Operations:** Functions like sub-sequence extraction (`extract`) and translation
        (`translate`) are implemented using highly optimized, parallelized Numba kernels.
    2.  **Memory Efficiency:** Avoids the overhead of thousands of small string or byte objects.
    3.  **Fast I/O:** The flat array can be rapidly serialized/deserialized or written directly to formats
        like FASTA.

    Attributes:
        ids (tuple[str, ...]): A tuple containing the string identifiers for each sequence.
        seqs (npt.NDArray[np.uint8]): A single, flat 1D array containing all sequence data concatenated
            together. Characters are stored as ASCII byte values.
        offsets (npt.NDArray[np.int32]): A 1D array where `offsets[i]` gives the starting index
            in the `seqs` array for the i-th sequence.
        lengths (npt.NDArray[np.int32]): A 1D array where `lengths[i]` gives the length of the
            i-th sequence.
    """
    ids:     tuple[str, ...]
    seqs:    npt.NDArray[np.uint8]   # concatenated sequences
    offsets: npt.NDArray[np.int32]  # start index per sequence
    lengths: npt.NDArray[np.int32]  # length per sequence

    def __len__(self) -> int:
        """Returns the total number of sequences in the batch."""
        return len(self.ids)
    
    def __dict__(self) -> dict:
        """Serializes the batch into a dictionary, converting the sequence array to an ASCII string."""
        return {
            'ids': self.ids,
            'seqs': self.seqs.tobytes().decode('ascii'),
            'offsets': self.offsets.tolist(),
            'lengths': self.lengths.tolist()
        }
        
    @classmethod
    def from_dict(cls, data: dict) -> Sequences:
        """Deserializes a Sequences object from a dictionary representation."""
        return cls(
            ids=tuple(data['ids']),
            seqs=np.frombuffer(data['seqs'].encode('ascii'), dtype=np.uint8),
            offsets=np.array(data['offsets'], dtype=np.int32),
            lengths=np.array(data['lengths'], dtype=np.int32)
        )

    def to_fasta(self, use_indices: bool = False) -> bytes:
        """Formats the entire batch as a single FASTA byte string.

        This method efficiently iterates over the flat memory layout to construct the FASTA format.

        Args:
            use_indices (bool, optional): If True, uses the numeric index (0, 1, 2...) as the FASTA
                header instead of the string ID. Useful if IDs are empty or redundant. Defaults to False.

        Returns:
            bytes: The complete FASTA-formatted byte string for all sequences.
        """
        # Allow empty ids if we are explicitly using numeric indices
        if not self.ids and not use_indices:
            return b""
            
        seq_bytes = self.seqs.tobytes()
        if use_indices:
            return b"".join([
                b">%d\n%b\n" % (i, seq_bytes[o : o + l])
                for i, (o, l) in enumerate(zip(self.offsets.tolist(), self.lengths.tolist()))
            ])
        else:
            return b"".join([
                b">%b\n%b\n" % (i.encode(), seq_bytes[o : o + l])
                for i, o, l in zip(self.ids, self.offsets.tolist(), self.lengths.tolist())
            ])

    @classmethod
    def empty(cls) -> Sequences:
        """Creates an empty Sequences object with correctly typed, zero-length arrays."""
        return cls(
            (),
            np.empty(0, dtype=np.uint8),
            np.empty(0, dtype=np.int32),
            np.empty(0, dtype=np.int32),
        )

    @classmethod
    def concat(cls, batches: Iterable[Sequences]) -> Sequences:
        """Concatenates multiple Sequences objects into a single, larger collection.

        Args:
            batches (Iterable[Sequences]): A list or iterable of sequences to concatenate.

        Returns:
            Sequences: A new, combined sequences container.
        """
        batches = list(batches)
        if not batches:
            return cls.empty()

        all_ids = sum((b.ids for b in batches), ())
        all_seqs = np.concatenate([b.seqs for b in batches])
        all_lengths = np.concatenate([b.lengths for b in batches])

        offsets = np.zeros(len(all_lengths), dtype=np.int32)
        if len(all_lengths) > 1:
            np.cumsum(all_lengths[:-1], out=offsets[1:])
            
        return cls(all_ids, all_seqs, offsets, all_lengths)

    @property
    def internal_stops(self) -> np.ndarray:
        """Vectorized check for premature stop codons in protein sequences.

        This property scans each sequence for the stop codon character ('*') *before* the final
        position. It uses a fast Numba kernel for parallel execution.

        Returns:
            np.ndarray: A 1D boolean array where True indicates the presence of an internal stop.
        """
        return _internal_stops_kernel(self.seqs, self.offsets, self.lengths)

    def unique(self) -> Sequences:
        """Returns a new Sequences container containing only unique sequences.
        
        Uses a Numba-accelerated 64-bit FNV-1a hash to rapidly identify unique sequences.
        The first occurrence of each unique sequence is kept, preserving the original order.
        """
        if len(self) <= 1:
            return self
            
        hashes = _hash_sequences_kernel(self.seqs, self.offsets, self.lengths)
        _, unique_indices = np.unique(hashes, return_index=True)
        unique_indices.sort()
        
        return self[unique_indices]

    def __getitem__(self, item: int | slice | np.ndarray | list) -> 'SeqRecord | Sequences':
        """Accesses sequences by index, slice, or boolean mask.

        - If `item` is an integer, returns a single `SeqRecord` object.
        - If `item` is a slice or boolean mask, returns a new, smaller `Sequences` collection.

        Args:
            item: An integer index, slice, or NumPy array mask/indices.

        Returns:
            SeqRecord | Sequences: A single record or a new sequences collection.

        Raises:
            IndexError: If an integer index is out of range.
        """
        if isinstance(item, (int, np.integer)):
            if item < 0: item += len(self)
            if item < 0 or item >= len(self): raise IndexError("Index out of range")
            s = self.offsets[item]
            l = self.lengths[item]
            return SeqRecord(self.ids[item], self.seqs[s:s + l].tobytes())
            
        if isinstance(item, slice):
            indices = np.arange(len(self))[item]
        else:
            indices = np.asarray(item)
            if indices.dtype == bool:
                indices = np.nonzero(indices)[0]
                
        starts = np.zeros(len(indices), dtype=np.int32)
        ends = self.lengths[indices].astype(np.int32)
        strands = np.ones(len(indices), dtype=np.int8)
        return self.extract(indices.astype(np.int32), starts, ends, strands, new_ids=tuple(self.ids[i] for i in indices))

    def __iter__(self) -> Generator[SeqRecord, None, None]:
        """Iterates over the batch, yielding a single `SeqRecord` at a time.

        Note: Iterating over `Sequences` defeats its performance benefits. Use vectorized
        methods (`extract`, `translate`, etc.) whenever possible.
        """
        for i in range(len(self)):
            s = self.offsets[i]
            l = self.lengths[i]
            yield SeqRecord(self.ids[i], self.seqs[s:s + l].tobytes())

    @classmethod
    def from_bytes(cls, seqs: list[bytes], ids: tuple[str, ...] | None = None) -> Sequences:
        """Constructs Sequences from a list of byte strings.

        Args:
            seqs (list[bytes]): A list of biological sequences as bytes.
            ids (tuple[str, ...], optional): A tuple of identifiers. If not provided, numeric string
                indices ("0", "1", ...) will be used.

        Returns:
            Sequences: The newly constructed sequences container.
        """
        ids = ids or tuple(str(i) for i in range(len(seqs)))
        return cls.from_records([SeqRecord(i, s) for i, s in zip(ids, seqs)])

    @classmethod
    def from_records(cls, records: list[SeqRecord]) -> Sequences:
        """Constructs Sequences from a list of `SeqRecord` objects.

        This method handles flattening the individual byte sequences into the contiguous
        NumPy array and computing the necessary offsets.

        Args:
            records (list[SeqRecord]): A list of sequence records.

        Returns:
            Sequences: The newly constructed sequences container.
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

    def extract(self, indices: npt.NDArray[np.int32], starts: npt.NDArray[np.int32],
                ends: npt.NDArray[np.int32], strands: npt.NDArray[np.int8],
                new_ids: tuple[str, ...] | None = None) -> Sequences:
        """Vectorized sub-sequence extraction from the batch.

        This is a highly optimized method for extracting many regions from the sequences in the batch
        simultaneously. It handles reverse-complementation automatically if the strand is negative.
        The work is performed by a multi-threaded Numba kernel.

        Args:
            indices (npt.NDArray[np.int32]): An array of indices specifying *which* sequence in the batch
                each extraction refers to.
            starts (npt.NDArray[np.int32]): An array of 0-based start coordinates relative to each parent sequence.
            ends (npt.NDArray[np.int32]): An array of 0-based end coordinates relative to each parent sequence.
            strands (npt.NDArray[np.int8]): An array of strand orientations (1 for forward, -1 for reverse).
            new_ids (tuple[str, ...], optional): Identifiers for the new, extracted sequences.

        Returns:
            Sequences: A new collection containing the extracted sub-sequences.
        """
        if len(indices) == 0:
            return self.empty()
            
        new_ids = new_ids or tuple(str(i) for i in range(len(indices)))
        out_seqs, offsets, lengths = _extract_ragged_kernel(
            self.seqs, self.offsets, indices, starts, ends, strands, BacterialTranslationTable._COMP_MAP
        )
        return Sequences(new_ids, out_seqs, offsets, lengths)

    def extract_intervals(self, indices: npt.NDArray[np.integer], intervals: Intervals,
                          new_ids: tuple[str, ...] | None = None) -> Sequences:
        """Convenience wrapper around `extract` that takes an `Intervals` object directly.

        Args:
            indices (npt.NDArray[np.integer]): The sequence indices to extract from.
            intervals (Intervals): The coordinates to extract.
            new_ids (tuple[str, ...], optional): New identifiers. Defaults to None.

        Returns:
            Sequences: A new collection containing the extracted sub-sequences.
        """
        return self.extract(
            indices.astype(np.int32),
            intervals.starts.astype(np.int32),
            intervals.ends.astype(np.int32),
            intervals.strands,
            new_ids=new_ids
        )

    def translate(self, frames: npt.NDArray[np.int8] | None = None) -> Sequences:
        """Vectorized translation of nucleotide sequences into protein sequences.

        This method uses a highly optimized Numba kernel and lookup tables to translate
        all sequences in the batch simultaneously according to the standard bacterial genetic code
        (NCBI Translation Table 11).

        Args:
            frames (npt.NDArray[np.int8], optional): An array specifying the reading frame offset
                (0, 1, or 2) for each sequence. If None, frame 0 is used for all. Defaults to None.

        Returns:
            Sequences: A new collection containing the translated amino acid sequences.
        """
        if len(self) == 0:
            return self.empty()
        if frames is None:
            frames = np.zeros(len(self), dtype=np.int8)
        out_seqs, offsets, lengths = _translate_ragged_kernel(
            self.seqs, self.offsets, self.lengths, frames,
            BacterialTranslationTable._CHAR_MAP, BacterialTranslationTable._CODON_MAP
        )
        return Sequences(self.ids, out_seqs, offsets, lengths)


class BacterialTranslationTable:
    _MAPPING = {
        b'TTT': b'F', b'TTC': b'F', b'TTA': b'L', b'TTG': b'L',
        b'TCT': b'S', b'TCC': b'S', b'TCA': b'S', b'TCG': b'S',
        b'TAT': b'Y', b'TAC': b'Y', b'TAA': b'*', b'TAG': b'*',
        b'TGT': b'C', b'TGC': b'C', b'TGA': b'*', b'TGG': b'W',
        b'CTT': b'L', b'CTC': b'L', b'CTA': b'L', b'CTG': b'L',
        b'CCT': b'P', b'CCC': b'P', b'CCA': b'P', b'CCG': b'P',
        b'CAT': b'H', b'CAC': b'H', b'CAA': b'Q', b'CAG': b'Q',
        b'CGT': b'R', b'CGC': b'R', b'CGA': b'R', b'CGG': b'R',
        b'ATT': b'I', b'ATC': b'I', b'ATA': b'I', b'ATG': b'M',
        b'ACT': b'T', b'ACC': b'T', b'ACA': b'T', b'ACG': b'T',
        b'AAT': b'N', b'AAC': b'N', b'AAA': b'K', b'AAG': b'K',
        b'AGT': b'S', b'AGC': b'S', b'AGA': b'R', b'AGG': b'R',
        b'GTT': b'V', b'GTC': b'V', b'GTA': b'V', b'GTG': b'V',
        b'GCT': b'A', b'GCC': b'A', b'GCA': b'A', b'GCG': b'A',
        b'GAT': b'D', b'GAC': b'D', b'GAA': b'E', b'GAG': b'E',
        b'GGT': b'G', b'GGC': b'G', b'GGA': b'G', b'GGG': b'G',
    }
    _START_CODONS = {b'TTG', b'CTG', b'ATT', b'ATC', b'ATA', b'ATG', b'GTG'}
    _STOP_CODONS = {b'TAA', b'TAG', b'TGA'}
    _COMP = bytes.maketrans(b"ACGTUacgtu", b"TGCAAtgcaa")
    _CHAR_MAP = np.full(256, 4, dtype=np.uint8)
    for _i, _c in enumerate(b'ACGT'):
        _CHAR_MAP[_c] = _i
        _CHAR_MAP[_c + 32] = _i  # Map lowercase automatically
    _CHAR_MAP[b'U'[0]] = 3
    _CHAR_MAP[b'u'[0]] = 3
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
    def translate(cls, seq: bytes | bytearray | memoryview | npt.NDArray[np.uint8]) -> npt.NDArray[np.uint8]:
        if len(seq) < 3:
            return np.array([], dtype=np.uint8)

        if not isinstance(seq, np.ndarray):
            seq = np.ascontiguousarray(np.frombuffer(seq, np.uint8))

        return _translate_kernel(seq, cls._CHAR_MAP, cls._CODON_MAP)

    @classmethod
    def is_coding(cls, seq: bytes) -> bool:
        if len(seq) < 3:
            return False
        return seq[:3] in cls._START_CODONS and seq[-3:] in cls._STOP_CODONS


# class SeqEncoder:
#     __slots__ = ('_lut', 'alphabet_size')
#
#     def __init__(self, mapping: dict[bytes, int], unknown_int: int | None = None):
#         # Map unknowns to the next available integer to keep the alphabet contiguous
#         if unknown_int is None:
#             unknown_int = max(mapping.values()) + 1
#
#         self.alphabet_size = unknown_int + 1
#
#         # Pre-fill the 256-byte lookup table with unknown_int
#         lut = bytearray([unknown_int] * 256)
#
#         # Map both uppercase and lowercase source characters directly
#         for source_byte, target_int in mapping.items():
#             lut[source_byte[0]] = target_int
#             lut[source_byte.lower()[0]] = target_int
#
#         # Store as immutable bytes for optimal C-level .translate() performance
#         self._lut = bytes(lut)
#
#     def encode(self, seq: bytes) -> bytes:
#         return seq.translate(self._lut)
#
#     @classmethod
#     @cache
#     def dayhoff(cls) -> 'SeqEncoder':
#         return cls(
#             {
#                 b'C': 0,
#                 b'A': 1, b'G': 1, b'P': 1, b'S': 1, b'T': 1,
#                 b'D': 2, b'E': 2, b'N': 2, b'Q': 2,
#                 b'H': 3, b'K': 3, b'R': 3,
#                 b'I': 4, b'L': 4, b'M': 4, b'V': 4,
#                 b'F': 5, b'W': 5, b'Y': 5
#             }
#         )
#
#     @classmethod
#     @cache
#     def hydrophobic(cls) -> 'SeqEncoder':
#         return cls(
#             {
#                 b'A': 0, b'C': 0, b'F': 0, b'I': 0, b'L': 0, b'M': 0, b'V': 0, b'W': 0, b'Y': 0,
#                 b'R': 1, b'N': 1, b'D': 1, b'E': 1, b'Q': 1, b'G': 1, b'H': 1, b'K': 1, b'P': 1, b'S': 1, b'T': 1
#             }
#         )
#
#     @classmethod
#     @cache
#     def mmseqs12(cls) -> 'SeqEncoder':
#         return cls(
#             {
#                 b'A': 0, b'S': 0, b'T': 0,
#                 b'L': 1, b'M': 1,
#                 b'I': 2, b'V': 2,
#                 b'K': 3, b'R': 3,
#                 b'E': 4, b'Q': 4,
#                 b'N': 5, b'D': 5,
#                 b'F': 6, b'Y': 6,
#                 b'C': 7,
#                 b'G': 8,
#                 b'H': 9,
#                 b'P': 10,
#                 b'W': 11
#             }
#         )



# Kernels --------------------------------------------------------------------------------------------------------------
@njit(parallel=True, cache=True, nogil=True)
def _hash_sequences_kernel(seqs: npt.NDArray[np.uint8], offsets: npt.NDArray[np.int32],
                           lengths: npt.NDArray[np.int32]) -> npt.NDArray[np.uint64]:
    """Computes a 64-bit FNV-1a hash for each sequence in parallel."""
    n = len(offsets)
    hashes = np.zeros(n, dtype=np.uint64)
    for i in prange(n):
        s = offsets[i]
        l = lengths[i]
        h = np.uint64(14695981039346656037)  # FNV offset basis
        for j in range(l):
            h = (h ^ np.uint64(seqs[s + j])) * np.uint64(1099511628211)  # FNV prime
        # Mix the length to further reduce collision risk
        h = (h ^ np.uint64(l)) * np.uint64(1099511628211)
        hashes[i] = h
    return hashes


@njit(cache=True, nogil=True, parallel=True)
def _translate_kernel(seq_arr: npt.NDArray[np.uint8], char_map: npt.NDArray[np.uint8],
                      codon_map: npt.NDArray[np.uint8]) -> npt.NDArray[np.uint8]:
    n_codons = len(seq_arr) // 3
    out = np.empty(n_codons, dtype=np.uint8)
    for i in prange(n_codons):
        c1 = char_map[seq_arr[i * 3]]
        c2 = char_map[seq_arr[i * 3 + 1]]
        c3 = char_map[seq_arr[i * 3 + 2]]
        idx = c1 * 25 + c2 * 5 + c3
        out[i] = codon_map[idx]
    return out


@njit(cache=True, nogil=True, parallel=True)
def _extract_ragged_kernel(seqs: npt.NDArray[np.uint8], offsets: npt.NDArray[np.int32], indices: npt.NDArray[np.int32],
                           starts: npt.NDArray[np.int32], ends: npt.NDArray[np.int32], strands: npt.NDArray[np.int8],
                           comp_map: npt.NDArray[np.uint8]):
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
            out_offsets[i] = out_offsets[i-1] + out_lengths[i-1]
            
    out_seqs = np.empty(total_len, dtype=seqs.dtype)
    for i in prange(n):
        l = out_lengths[i]
        if l == 0: continue
        
        t_idx = indices[i]
        g_s = offsets[t_idx] + starts[i]
        g_e = offsets[t_idx] + ends[i]
        st = strands[i]
        out_s = out_offsets[i]
        
        if st >= 0:
            for c in range(l):
                out_seqs[out_s + c] = seqs[g_s + c]
        else:
            for c in range(l):
                out_seqs[out_s + c] = comp_map[seqs[g_e - 1 - c]]
    return out_seqs, out_offsets, out_lengths


@njit(cache=True, nogil=True, parallel=True)
def _translate_ragged_kernel(seqs: npt.NDArray[np.uint8], offsets: npt.NDArray[np.int32],
                             lengths: npt.NDArray[np.int32], frames: npt.NDArray[np.int8],
                             char_map: npt.NDArray[np.uint8], codon_map: npt.NDArray[np.uint8]):
    n = len(offsets)
    out_lengths = np.empty(n, dtype=np.int32)
    total_len = 0
    for i in range(n):
        l = lengths[i]
        f = frames[i]
        if l > f:
            adj_len = l - f
            out_lengths[i] = adj_len // 3 if adj_len >= 3 else 0
        else:
            out_lengths[i] = 0
        total_len += out_lengths[i]
        
    out_offsets = np.empty(n, dtype=np.int32)
    if n > 0:
        out_offsets[0] = 0
        for i in range(1, n):
            out_offsets[i] = out_offsets[i-1] + out_lengths[i-1]
            
    out_data = np.empty(total_len, dtype=np.uint8)
    for i in prange(n):
        n_codons = out_lengths[i]
        if n_codons == 0: continue
        
        s = offsets[i] + frames[i]
        out_s = out_offsets[i]
        
        ptr = s
        for c in range(n_codons):
            c1, c2, c3 = char_map[seqs[ptr]], char_map[seqs[ptr + 1]], char_map[seqs[ptr + 2]]
            out_data[out_s + c] = codon_map[c1 * 25 + c2 * 5 + c3]
            ptr += 3
    return out_data, out_offsets, out_lengths


@njit(cache=True, nogil=True, parallel=True)
def _internal_stops_kernel(seqs: npt.NDArray[np.uint8], offsets: npt.NDArray[np.int32],
                           lengths: npt.NDArray[np.int32]) -> npt.NDArray[np.bool_]:
    n = len(offsets)
    out = np.zeros(n, dtype=np.bool_)
    for i in prange(n):
        s = offsets[i]
        l = lengths[i]
        if l > 0:  # Check up to the second-to-last character (s + l - 1)
            for j in range(s, s + l - 1):
                if seqs[j] == 42:  # 42 is ord('*')
                    out[i] = True
                    break
    return out
