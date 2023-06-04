"""
Copyright 2023 Tom Stanton (tomdstanton@gmail.com)
https://github.com/tomdstanton/kaptive-mapper

This file is part of kaptive-mapper. kaptive-mapper is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by the Free Software Foundation,
either version 3 of the License, or (at your option) any later version. kaptive-mapper is distributed
in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
details. You should have received a copy of the GNU General Public License along with kaptive-mapper.
If not, see <https://www.gnu.org/licenses/>.
"""
from re import compile, MULTILINE, finditer
from io import DEFAULT_BUFFER_SIZE, BufferedReader
import json
import fcntl
from collections import OrderedDict
from pathlib import Path


# Record regexes ------------------------------------------------------------------------------------------------------
FASTA_REGEX = compile(rb'^>(?P<header>.+?)\n(?P<sequence>(?:\n|[^>])+)', MULTILINE)


# Classes -------------------------------------------------------------------------------------------------------------
class Assembly:
    def __init__(self, path: Path, contigs: dict):
        """Loads in an assembly and builds a BLAST database for it (if necessary)."""
        self.path = path
        self.name = path.name.strip('.gz').rsplit('.', 1)[0]
        self.contigs = contigs  # key = name, value = sequence

    def __repr__(self):
        return self.name

    def __str__(self):
        return str(self.path)

    def __len__(self):
        return self.path.stat().st_size

    # def type_gene_search(self, db: 'Database'):
    #     if db.type_gene_names and self.args.allelic_typing:
    #         all_gene_blast_hits = blast_hit_filter(
    #             get_blast_hits(self.fasta, str(self.args.allelic_typing), self.args, type_genes=True),
    #             query_cov=self.args.min_gene_cov, pident=self.args.min_gene_id
    #         )
    #         if not self.result:
    #             self.result = Result(db, self)
    #
    #         for gene_name in db.type_gene_names:
    #             blast_hits = sorted([x for x in all_gene_blast_hits if x.gene_name == gene_name],
    #                                 reverse=True, key=lambda z: z.bitscore)
    #             if blast_hits:
    #                 perfect_match = None
    #                 for hit in blast_hits:
    #                     if hit.pident == 100.0 and hit.query_cov == 100.0:
    #                         perfect_match = hit
    #                         break
    #                 if perfect_match:
    #                     hit = perfect_match
    #                     hit.result = str(perfect_match.allele_number)
    #                 else:
    #                     hit = blast_hits[0]
    #                     hit.result = str(blast_hits[0].allele_number) + '*'
    #                 self.result.type_gene_results[gene_name] = hit


class AssemblyPiece(object):
    """
    This class describes a piece of an assembly: which contig the piece is on and what the range is.
    """

    def __init__(self, assembly, contig_name, contig_start, contig_end, strand, blast_hits=None):
        self.assembly = assembly
        self.contig_name = contig_name
        self.start = contig_start
        self.end = contig_end
        self.strand = strand
        if not blast_hits:
            blast_hits = []
        self.blast_hits = blast_hits

    def __repr__(self):
        return f'{self.assembly.name}_{self.get_header()}'

    def __len__(self):
        """Returns the sequence length for this piece."""
        return self.end - self.start

    def get_header(self):
        """Returns a descriptive string for the FASTA header when saving this piece to file."""
        nice_contig_name = get_nice_contig_name(self.contig_name)
        return f'{nice_contig_name}_{self.start + 1}_to_{self.end}_{self.strand}_strand'

    def get_nice_header(self):
        """Like get_header, but uses spaces/parentheses instead of underscores for readability."""
        nice_contig_name = get_nice_contig_name(self.contig_name)
        return f'{nice_contig_name} ({self.start + 1}-{self.end}, {self.strand} strand)'

    def get_bandage_range(self):
        """Returns the assembly piece in a Bandage path format."""
        if is_contig_name_spades_format(self.contig_name):
            name = self.contig_name.split('_')[1]
        else:
            name = self.contig_name.split()[0]
        return f'({self.start + 1}) {name}+ ({self.end})'

    def get_sequence(self):
        """Returns the DNA sequence for this piece of the assembly."""
        seq = self.assembly.contigs[self.contig_name][self.start:self.end]
        return seq if self.strand == '+' else reverse_complement(seq)

    def combine(self, other):
        """
        If this assembly piece and the other can be combined, this function returns the combined
        piece. If they can't, it returns None.
        To be able to combine, pieces must be overlapping or directly adjacent and on the same
        strand.
        """
        if self.contig_name != other.contig_name or self.strand != other.strand:
            return None
        combined = IntRange([(self.start, self.end)])
        combined.add_range(other.start, other.end)
        if len(combined.ranges) == 1:
            new_start, new_end = combined.ranges[0]
            combined_hits = self.blast_hits + other.blast_hits
            return AssemblyPiece(self.assembly, self.contig_name, new_start, new_end, self.strand,
                                 combined_hits)
        else:
            return None

    def overlaps(self, contig_name, start, end):
        """Returns whether this assembly piece overlaps with the given parameters."""
        return self.contig_name == contig_name and self.start < end and start < self.end

    def earliest_hit_coordinate(self):
        """Returns the lowest query start coordinate in the BLAST hits."""
        if not self.blast_hits:
            return None
        return sorted([x.qstart for x in self.blast_hits])[0]

    def latest_hit_coordinate(self):
        """Returns the highest query end coordinate in the BLAST hits."""
        if not self.blast_hits:
            return None
        return sorted([x.qend for x in self.blast_hits])[-1]


def is_contig_name_spades_format(contig_name):
    """
    Returns whether or not the contig name appears to be in the SPAdes/Velvet format.
    Example: NODE_5_length_150905_cov_4.42519
    """
    contig_name_parts = contig_name.split('_')
    return len(contig_name_parts) > 5 and contig_name_parts[0] == 'NODE' and \
        contig_name_parts[2] == 'length' and contig_name_parts[4] == 'cov'


def get_nice_contig_name(contig_name):
    """
    For a contig with a SPAdes/Velvet format, this function returns a simplified string that is
    just NODE_XX where XX is the contig number.
    For any other format, this function trims off everything following the first whitespace.
    """
    if is_contig_name_spades_format(contig_name):
        return f'NODE_{contig_name.split("_")[1]}'
    else:
        return contig_name.split()[0]


def load_fasta(filename: Path) -> dict:
    """Returns the names and sequences for the given fasta file."""
    fasta_seqs = {}
    open_func = gzip.open if get_compression_type(filename) == 'gz' else open
    try:
        with open_func(filename, 'rt') as f:
            fasta_seqs = {i[0].split(' ')[0]: i[1] for i in SimpleFastaParser(f)}
    except Exception as e:
        LOGGER.error(e)
    return fasta_seqs


def create_table_file(output_prefix: Path, type_gene_names):
    """
    Creates the table file and writes a header line if necessary.
    If the file already exists and the header line is correct, then it does nothing (to allow
    multiple independent processes to append to the file).
    """
    path = Path(f'{output_prefix}_table.txt')

    # If the table already exists, we don't need to do anything.
    if path.exists():
        with open(path, 'rt') as existing_table:
            first_line = existing_table.readline().strip()
            if first_line.startswith('Assembly\tBest match locus'):
                return

    headers = ['Assembly',
               'Best match locus',
               'Best match type',
               'Match confidence',
               'Problems',
               'Coverage',
               'Identity',
               'Length discrepancy',
               'Expected genes in locus',
               'Expected genes in locus, details',
               'Missing expected genes',
               'Other genes in locus',
               'Other genes in locus, details',
               'Expected genes outside locus',
               'Expected genes outside locus, details',
               'Other genes outside locus',
               'Other genes outside locus, details']

    if type_gene_names:
        headers += type_gene_names

    with open(path, 'wt') as table:
        table.write('\t'.join(headers))
        table.write('\n')


def write_json_file(output_prefix: Path, json_record: OrderedDict):
    path = output_prefix.with_suffix('.json')
    json_records = [json_record]
    if path.exists():
        with open(path, 'rt') as json_out:
            fcntl.flock(json_out, fcntl.LOCK_EX)
            json_records += json.load(json_out, object_pairs_hook=OrderedDict)
            fcntl.flock(json_out, fcntl.LOCK_UN)

    with open(path, 'wt') as json_out:
        fcntl.flock(json_out, fcntl.LOCK_EX)
        json.dump(json_records, json_out, indent=4)
        fcntl.flock(json_out, fcntl.LOCK_UN)


def parse_fasta(filename, buffer_size=DEFAULT_BUFFER_SIZE):
    """
    Generator function that yields one record at a time from a FASTA file.
    """
    with open(filename, 'rb') as f:
        buffer = BufferedReader(f)
        while True:
            chunk = buffer.read(buffer_size)
            if not chunk:
                break
            matches = finditer(FASTA_REGEX, chunk)
            for match in matches:
                header = match.group('header').decode()
                sequence = match.group('sequence').replace(b'\n', b'').decode()
                yield header, sequence


def merge_assembly_pieces(pieces):
    """
    Takes a list of AssemblyPiece objects and returns another list of AssemblyPiece objects where
    the overlapping pieces have been merged.
    """
    while True:
        merged_pieces = []
        merge_count = 0
        while pieces:
            merged_piece = pieces[0]
            unmerged = []
            for other_piece in pieces[1:]:
                combined = merged_piece.combine(other_piece)
                if not combined:
                    unmerged.append(other_piece)
                else:
                    merged_piece = combined
                    merge_count += 1
            merged_pieces.append(merged_piece)
            pieces = unmerged
        if merge_count == 0:
            break
        else:
            pieces = merged_pieces
    return merged_pieces


def fill_assembly_piece_gaps(pieces, max_gap_fill_size, assembly):
    """
    This function takes a list of assembly pieces, and if any of them are close enough to each
    other, the gap will be merged in.
    It assumes that all given pieces are from the same assembly.
    """
    pieces_by_contig_and_strand = {}
    fixed_pieces = []
    for piece in pieces:
        contig = piece.contig_name
        strand = piece.strand
        if (contig, strand) not in pieces_by_contig_and_strand:
            pieces_by_contig_and_strand[(contig, strand)] = []
        pieces_by_contig_and_strand[(contig, strand)].append(piece)
    for (contig, strand), pieces_in_contig_and_strand in pieces_by_contig_and_strand.items():
        gap_filling_pieces = []
        sorted_pieces = sorted(pieces_in_contig_and_strand, key=lambda x: x.start)
        max_end = sorted_pieces[0].end
        gaps = []
        for piece in sorted_pieces[1:]:
            if piece.start > max_end and piece.start - max_end <= max_gap_fill_size:
                gaps.append((max_end, piece.start))
            max_end = max(max_end, piece.end)
        for gap in gaps:
            gap_filling_pieces.append(AssemblyPiece(assembly, contig, gap[0], gap[1], strand))
        before_merge = pieces_in_contig_and_strand + gap_filling_pieces
        filled_pieces = merge_assembly_pieces(before_merge)
        fixed_pieces += filled_pieces
    return fixed_pieces


def get_mean_identity(pieces):
    """Returns the mean identity (weighted by sequence length) for a list of assembly pieces."""
    identity_sum = 0.0
    length_sum = 0
    for piece in pieces:
        for hit in piece.blast_hits:
            length_sum += hit.length
            identity_sum += hit.length * hit.pident
    if identity_sum == 0.0:
        return 0.0
    else:
        return identity_sum / length_sum


