from dataclasses import dataclass
from typing import ClassVar, Generator

import numpy as np
import numpy.typing as npt
from numba import njit, prange

from kaptive.core.constants import Strand
from kaptive.core.interval import Interval, IntervalLike, IntervalBatch


# Classes --------------------------------------------------------------------------------------------------------------
@dataclass(frozen=True, slots=True)
class SeqRecord:
    id: str
    seq: bytes

    def __len__(self):
        return len(self.seq)

    def to_fasta(self) -> bytes:
        return b">%b\n%b\n" % (self.id.encode(), self.seq)

    def extract(self, start: int | IntervalLike, end: int | None = None, strand: Strand = Strand.UNSTRANDED) -> bytes:
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
class SeqBatch:
    ids:     tuple[str, ...]
    seqs:    npt.NDArray[np.uint8]   # concatenated sequences
    offsets: npt.NDArray[np.uint32]  # start index per sequence
    lengths: npt.NDArray[np.uint32]  # length per sequence
    
    _SEQ_DTYPE: ClassVar = np.uint8
    _OFFSET_DTYPE: ClassVar = np.uint32
    _LENGTH_DTYPE: ClassVar = np.uint32


    def __len__(self) -> int:
        return len(self.ids)

    def to_fasta(self) -> bytes:
        """Writes the record to the file handle in a flat fasta format"""
        return b''.join(i.to_fasta() for i in iter(self))

    @classmethod
    def empty(cls) -> 'SeqBatch':
        """Create an empty SeqBatch."""
        return cls(
            (),
            np.empty(0, dtype=np.uint8),
            np.empty(0, dtype=np.uint32),
            np.empty(0, dtype=np.uint32),
        )

    @property
    def internal_stops(self) -> np.ndarray:
        """Vectorized check for internal stop codons (*)."""
        return _internal_stops_kernel(self.seqs, self.offsets, self.lengths)

    def __getitem__(self, item: int | slice | np.ndarray | list) -> 'SeqRecord | SeqBatch':
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
                
        starts = np.zeros(len(indices), dtype=np.uint32)
        ends = self.lengths[indices].astype(np.uint32)
        strands = np.ones(len(indices), dtype=np.int8)
        return self.extract(indices.astype(np.uint32), starts, ends, strands, new_ids=tuple(self.ids[i] for i in indices))

    def __iter__(self) -> Generator[SeqRecord, None, None]:
        for i in range(len(self)):
            s = self.offsets[i]
            l = self.lengths[i]
            yield SeqRecord(self.ids[i], self.seqs[s:s + l].tobytes())

    @classmethod
    def from_bytes(cls, seqs: list[bytes], ids: tuple[str, ...] | None = None) -> 'SeqBatch':
        ids = ids or tuple(str(i) for i in range(len(seqs)))
        return cls.from_records([SeqRecord(i, s) for i, s in zip(ids, seqs)])

    @classmethod
    def from_records(cls, records: list[SeqRecord]) -> 'SeqBatch':
        ids = tuple(r.id for r in records)
        seqs = [np.frombuffer(r.seq, dtype=cls._SEQ_DTYPE) for r in records]
        if not seqs:
            return cls.empty()
        out_seqs = np.concatenate(seqs, dtype=cls._SEQ_DTYPE)
        lengths = np.array([len(s) for s in seqs], dtype=cls._LENGTH_DTYPE)
        offsets = np.zeros(len(seqs), dtype=cls._OFFSET_DTYPE)
        if len(seqs) > 1:
            np.cumsum(lengths[:-1], out=offsets[1:])
        return cls(ids, out_seqs, offsets, lengths)

    def extract(self, indices: npt.NDArray[np.uint32], starts: npt.NDArray[np.uint32],
                ends: npt.NDArray[np.uint32], strands: npt.NDArray[np.int8],
                new_ids: tuple[str, ...] | None = None) -> 'SeqBatch':
        if len(indices) == 0:
            return self.empty()
            
        new_ids = new_ids or tuple(str(i) for i in range(len(indices)))
        out_seqs, offsets, lengths = _extract_ragged_kernel(
            self.seqs, self.offsets, indices, starts, ends, strands, BacterialTranslationTable._COMP_MAP
        )
        return SeqBatch(new_ids, out_seqs, offsets, lengths)

    def extract_intervals(self, indices: npt.NDArray[np.integer], intervals: IntervalBatch,
                          new_ids: tuple[str, ...] | None = None) -> 'SeqBatch':
        """Convenience method to extract sequence blocks using an IntervalBatch directly."""
        return self.extract(
            indices.astype(np.uint32),
            intervals.starts.astype(np.uint32),
            intervals.ends.astype(np.uint32),
            intervals.strands,
            new_ids=new_ids
        )

    def translate(self, frames: npt.NDArray[np.uint8] | None = None, new_ids: tuple[str, ...] | None = None) -> 'SeqBatch':
        if len(self) == 0:
            return self.empty()
            
        new_ids = new_ids or self.ids
        if frames is None:
            frames = np.zeros(len(self), dtype=np.uint8)
        out_seqs, offsets, lengths = _translate_ragged_kernel(
            self.seqs, self.offsets, self.lengths, frames,
            BacterialTranslationTable._CHAR_MAP, BacterialTranslationTable._CODON_MAP
        )
        return SeqBatch(new_ids, out_seqs, offsets, lengths)


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
def _extract_ragged_kernel(seqs: npt.NDArray[np.uint8], offsets: npt.NDArray[np.uint32], indices: npt.NDArray[np.uint32],
                           starts: npt.NDArray[np.uint32], ends: npt.NDArray[np.uint32], strands: npt.NDArray[np.int8],
                           comp_map: npt.NDArray[np.uint8]):
    n = len(indices)
    out_lengths = np.empty(n, dtype=np.uint32)
    total_len = 0
    for i in range(n):
        out_lengths[i] = ends[i] - starts[i]
        total_len += out_lengths[i]
        
    out_offsets = np.empty(n, dtype=offsets.dtype)
    if n > 0:
        out_offsets[0] = 0
        for i in range(1, n):
            out_offsets[i] = out_offsets[i-1] + out_lengths[i-1]
            
    out_data = np.empty(total_len, dtype=seqs.dtype)
    for i in prange(n):
        l = out_lengths[i]
        if l == 0: continue
        
        t_idx = indices[i]
        g_s = offsets[t_idx] + starts[i]
        g_e = offsets[t_idx] + ends[i]
        st = strands[i]
        out_s = out_offsets[i]
        
        if st >= 0:
            for c in range(l): out_data[out_s + c] = seqs[g_s + c]
        else:
            for c in range(l): out_data[out_s + c] = comp_map[seqs[g_e - 1 - c]]
    return out_data, out_offsets, out_lengths


@njit(cache=True, nogil=True, parallel=True)
def _translate_ragged_kernel(seqs: npt.NDArray[np.uint8], offsets: npt.NDArray[np.uint32],
                             lengths: npt.NDArray[np.uint32], frames: npt.NDArray[np.uint8],
                             char_map: npt.NDArray[np.uint8], codon_map: npt.NDArray[np.uint8]):
    n = len(offsets)
    out_lengths = np.empty(n, dtype=lengths.dtype)
    total_len = 0
    for i in range(n):
        adj_len = lengths[i] - frames[i]
        out_lengths[i] = adj_len // 3 if adj_len >= 3 else 0
        total_len += out_lengths[i]
        
    out_offsets = np.empty(n, dtype=offsets.dtype)
    if n > 0:
        out_offsets[0] = 0
        for i in range(1, n):
            out_offsets[i] = out_offsets[i-1] + out_lengths[i-1]
            
    out_data = np.empty(total_len, dtype=seqs.dtype)
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
def _internal_stops_kernel(seqs: npt.NDArray[np.uint8], offsets: npt.NDArray[np.uint32],
                           lengths: npt.NDArray[np.uint32]) -> npt.NDArray[np.bool_]:
    n = len(offsets)
    out = np.zeros(n, dtype=np.bool_)
    for i in prange(n):
        s = offsets[i]
        l = lengths[i]
        if l > 0:
            # Check up to the second-to-last character (s + l - 1)
            for j in range(s, s + l - 1):
                if seqs[j] == 42:  # 42 is ord('*')
                    out[i] = True
                    break
    return out
