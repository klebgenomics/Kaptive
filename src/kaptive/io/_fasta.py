from collections.abc import Iterator
from typing import IO

from kaptive.core.seq import SeqRecord


# Readers --------------------------------------------------------------------------------------------------------------
class FastaReader(Iterator):
    """
    A high-performance FASTA file reader.

    Assumes the provided handle is opened in binary mode ('rb') to return the
    sequence as a byte string, which is highly efficient and ideal for wrapped
    FASTA lines.
    """
    def __init__(self, handle: IO[bytes]):
        self._handle = handle
        self._generator = self._parse_records()

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        self._handle.close()
        
    def __iter__(self) -> Iterator[SeqRecord]:
        return self

    def __next__(self) -> SeqRecord:
        return next(self._generator)
    
    def _parse_records(self) -> Iterator[SeqRecord]:
        header: str = ""
        seq_chunks: list[bytes] = []

        # Local variable caching speeds up lookups inside the tight loop
        join_bytes = b"".join

        for line in self._handle:
            # rstrip() is faster than strip() and removes \r\n or \n
            line = line.rstrip()
            if not line:
                continue

            # 62 is the ASCII integer value for '>'
            if line[0] == 62:
                if header:
                    # Yield the previous record
                    yield SeqRecord(seq=join_bytes(seq_chunks), id=header)

                seq_chunks.clear()
                
                # Decode the header and split into ID and description
                decoded_line = line[1:].decode("utf-8", errors="replace")
                header, _, description = decoded_line.partition(" ")

            else:
                seq_chunks.append(line)

        # Yield the final record once the loop ends
        if header:
            yield SeqRecord(seq=join_bytes(seq_chunks), id=header)