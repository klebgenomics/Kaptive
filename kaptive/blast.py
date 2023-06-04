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
from pathlib import Path
from collections import OrderedDict  # For Python 3.6 compatibility
# TODO: Remove this when Python 3.6 support is dropped

from kaptive.misc import quick_translate, run_command, find_files_with_suffixes
from kaptive.log import quit_with_error, log, warning
from kaptive.assembly import get_nice_contig_name, AssemblyPiece
from kaptive.database import Gene
from kaptive.intrange import IntRange


class BlastHit(object):
    """
    Stores the BLAST hit output mostly verbatim. However, it does convert the BLAST ranges
    (1-based, inclusive end) to Python ranges (0-based, exclusive end).
    Designed for:
    -outfmt 6 qseqid sseqid qstart qend sstart send evalue bitscore length pident qlen qseq sseq slen
    """

    def __init__(self, hit_string):
        parts = hit_string.split('\t')
        self.qseqid = parts[0]
        self.sseqid = parts[1]
        self.qstart = int(parts[2]) - 1
        self.qend = int(parts[3])
        self.sstart = int(parts[4])
        self.send = int(parts[5])
        if self.sstart <= self.send:
            self.strand = '+'
        else:
            self.sstart, self.send = self.send, self.sstart
            self.strand = '-'
        self.sstart -= 1
        self.evalue = float(parts[6])
        self.bitscore = float(parts[7])
        self.length = int(parts[8])
        self.pident = float(parts[9])
        self.qlen = int(parts[10])
        self.query_cov = 100.0 * len(parts[11]) / self.qlen
        self.qseq = parts[11]
        self.sseq = parts[12]
        self.slen = int(parts[13])

    def __repr__(self):
        return f'{self.qseqid}, {self.get_contig_details_string()}, {self.get_coverage_details_string()}, ' \
               f'{self.get_identity_details_string()}'

    def __len__(self):
        return self.length

    def get_contig_details_string(self):
        """Returns a string describing the hit's range and strand in the contig."""
        return f'Contig: {get_nice_contig_name(self.sseqid)} ({self.sstart}-{self.send}, {self.strand} strand)'

    def get_coverage_string(self):
        return '%.2f' % self.query_cov + '%'

    def get_coverage_details_string(self, extra_space=False):
        first_part = 'Cov: '
        second_part = self.get_coverage_string()
        if len(second_part) == 6 and extra_space:
            first_part += ' '
        return first_part + second_part

    def get_identity_string(self):
        return '%.2f' % self.pident + '%'

    def get_identity_details_string(self, extra_space=False):
        first_part = 'ID: '
        second_part = self.get_identity_string()
        if len(second_part) == 6 and extra_space:
            first_part += ' '
        return first_part + second_part

    def get_aligned_string(self, spacing_1, spacing_2, spacing_3, cov_space, id_space):
        """Returns a string describing the hit with spaced parts for alignment."""
        aligned_string = self.qseqid + '  '
        aligned_string = aligned_string.ljust(spacing_1)
        aligned_string += self.get_contig_details_string() + '  '
        aligned_string = aligned_string.ljust(spacing_2)
        aligned_string += self.get_coverage_details_string(cov_space) + '  '
        aligned_string = aligned_string.ljust(spacing_3)
        aligned_string += self.get_identity_details_string(id_space)
        return aligned_string

    def get_assembly_piece(self, assembly):
        """Returns the piece of the assembly which corresponds to this BLAST hit."""
        return AssemblyPiece(assembly, self.sseqid, self.sstart, self.send, self.strand, [self])

    def get_query_range(self):
        """Produces an IntRange object for the hit query."""
        return IntRange([(self.qstart, self.qend)])

    def in_assembly_pieces(self, assembly_pieces):
        """
        Returns True if the hit is in (or at least overlaps with) any of the given assembly pieces.
        """
        for piece in assembly_pieces:
            if piece.overlaps(self.sseqid, self.sstart, self.send):
                return True
        return False

    def get_blast_result_json_dict(self):
        blast_results = OrderedDict()
        blast_results['Coverage'] = self.get_coverage_string()
        blast_results['Identity'] = self.get_identity_string()
        blast_results['Contig name'] = self.sseqid
        blast_results['Contig start position'] = self.sstart
        blast_results['Contig end position'] = self.send
        blast_results['Contig strand'] = self.strand
        blast_results['Bit score'] = self.bitscore
        blast_results['E-value'] = self.evalue
        return blast_results


class GeneBlastHit(BlastHit):
    """This class adds a few gene-specific things to the BlastHit class."""

    def __init__(self, hit_string):
        BlastHit.__init__(self, hit_string)
        self.non_synonymous_substitutions = None
        self.prot_alignment = None
        self.gene = None
        self.locus = None
        self.neighbour_before = None
        self.neighbour_after = None
        self.over_identity_threshold = False
        self.prot_seq = None

    def add_reference(self, ref: Gene):
        self.prot_seq = quick_translate(self.sseq, to_stop=True)
        self.gene = ref
        self.locus = self.gene.locus

    def truncated_protein(self):
        # We can only evaluate in context with the reference
        return True if len(self.prot_seq) < len(self.gene.prot_seq) else False

    def fragmented(self, min_gene_cov=0.0):
        return True if self.query_cov < min_gene_cov else False

    def edge_of_contig(self):
        return True if self.sstart == 0 or self.send == self.slen else False

    def deletion(self):
        # Test for indels that might result in protein internal stops etc...
        return True if ("-" in self.sseq and not "-" in self.qseq and not self.fragmented) else False

    def insertion(self):
        # Test for indels that might result in protein internal stops etc...
        return True if ("-" in self.qseq and not "-" in self.sseq and not self.fragmented) else False

    def conflicts(self, other):
        """
        Returns whether this hit conflicts with the other hit.
        A conflict is defined as the hits overlapping by 50% or more of the shortest hit's length.
        A hit is not considered to conflict with itself.
        """
        if self is other:
            return False
        if self.sseqid != other.sseqid:
            return False
        max_start = max(self.sstart, other.sstart)
        min_end = min(self.send, other.send)
        if max_start < min_end:
            overlap = min_end - max_start
        else:
            overlap = 0
        min_length = min(self.send - self.sstart, other.send - other.sstart)
        frac_overlap = overlap / min_length
        return frac_overlap > 0.5

    def get_blast_result_json_dict(self):
        blast_results = super(GeneBlastHit, self).get_blast_result_json_dict()
        blast_results['Nucleotide length'] = len(self.sseq)
        blast_results['Protein length'] = len(self.prot_seq)
        blast_results['Nucleotide sequence'] = self.sseq
        blast_results['Protein sequence'] = self.prot_seq
        return blast_results

    def get_match_confidence(self):
        cov = self.query_cov
        ident = self.pident
        cov_thresholds = [100, 95.0, 97.0, 95.0]
        ident_thresholds = [99.0, 95.0, 95.0, 85.0]
        if cov == cov_thresholds[0] and ident >= ident_thresholds[0]:
            confidence = 'Very high'
        elif cov >= cov_thresholds[1] and ident >= ident_thresholds[1]:
            confidence = 'High'
        elif cov >= cov_thresholds[2] and ident >= ident_thresholds[2]:
            confidence = 'Good'
        elif cov >= cov_thresholds[3] and ident >= ident_thresholds[3]:
            confidence = 'Low'
        else:
            confidence = 'None'
        return confidence


class TypeGeneBlastHit(BlastHit):
    """This class adds a couple type gene-specific things to the BlastHit class."""

    def __init__(self, hit_string):
        BlastHit.__init__(self, hit_string)
        try:
            name_parts = self.qseqid.split('__')
            self.gene_name = name_parts[1]
            self.allele_number = int(name_parts[2])
        except (IndexError, ValueError):
            self.gene_name = ''
            self.allele_number = 0

    def get_blast_result_json_dict(self):
        blast_results = OrderedDict()
        blast_results['Coverage'] = self.get_coverage_string()
        blast_results['Identity'] = self.get_identity_string()
        blast_results['Assembly piece name'] = self.sseqid
        blast_results['Assembly piece start position'] = self.sstart
        blast_results['Assembly piece end position'] = self.send
        blast_results['Assembly piece strand'] = self.strand
        blast_results['Bit score'] = self.bitscore
        blast_results['E-value'] = self.evalue
        blast_results['Length'] = len(self.sseq)
        blast_results['Sequence'] = self.sseq
        return blast_results


def get_blast_hits(database: Path, query: Union[List[Path], str], threads: int, type_genes=False,
                   genes=False):
    """Returns a list BlastHit objects for a search of the given query in the given database."""
    command = f"blastn|-task|blastn|-db|{database}|-num_threads|{threads}|-dust|no|-evalue|1e-10|-culling_limit|1|" \
              f"-outfmt|6 qseqid sseqid qstart qend sstart send evalue bitscore length pident qlen qseq sseq slen"
    out = run_command(f"{command}|-query|-", cmd_split="|", pipe=query) if isinstance(
        query, str) else run_command(f"{command}|-query|{query}", cmd_split="|")
    if not out:
        return []
    if genes:
        return [GeneBlastHit(line) for line in out.splitlines()]
    elif type_genes:
        return [TypeGeneBlastHit(line) for line in out.splitlines()]
    else:
        return [BlastHit(line) for line in out.splitlines()]


def cull_conflicting_hits(hit_to_keep, blast_hits):
    """
    This function returns a (potentially) reduced set of BLAST hits which excludes BLAST hits that
    overlap too much (same part of assembly) with the hit to keep.
    """
    return [x for x in blast_hits if not x.conflicts(hit_to_keep)]


def cull_all_conflicting_hits(blast_hits):
    """
    This function returns a (potentially) reduced set of BLAST hits where none of the remaining
    hits conflict.
    """
    blast_hits.sort(key=lambda x: x.bitscore, reverse=True)
    kept_hits = []
    while blast_hits:
        kept_hits.append(blast_hits.pop(0))
        blast_hits = cull_conflicting_hits(kept_hits[-1], blast_hits)
    return kept_hits


def blast_hit_filter(hits, bitscore=0, pident=0, query_cov=0):
    """Simple BlastHit object filter"""
    return [i for i in hits if i.bitscore >= bitscore and i.pident >= pident and i.query_cov >= query_cov]


def makeblastdb(fasta: Union[Path, str], out_path: Path, blastdb_version=5) -> Path:
    command = f'makeblastdb -title kaptive_blastn -blastdb_version {blastdb_version} -dbtype nucl -out {out_path}'
    run_command(f"{command} -in -", pipe=fasta) if isinstance(fasta, str) else run_command(
        f"{command} -in {fasta}")
    if not blastn_db_exists:
        quit_with_error(f'No blastn database files found in {out_path}')
    else:
        return out_path


def blastn_db_exists(reference: Union[Path, str]) -> bool:
    """
    Check if minimap2 index exists for the given reference
    :param reference: Path to reference
    :return: bool
    """
    suffixes = ['.nhr', '.nin', '.nsq', '.not', '.ntf', '.nto', '.ndb']
    return len(find_files_with_suffixes(reference, suffixes)) == len(suffixes)


# def get_best_hit_per_query(hits: List[BlastHit], groups: dict, attr: str = 'bitscore',
#                            reverse: bool = True) -> List[BlastHit]:
#     return [sorted([
#         i for i in hits if i.qseqid in hit_names
#     ], key=lambda x: getattr(x, attr), reverse=reverse)[0] for group, hit_names in groups.items()]


