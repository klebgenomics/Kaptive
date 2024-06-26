"""
Copyright 2023 Tom Stanton (tomdstanton@gmail.com)
https://github.com/klebgenomics/Kaptive

This file is part of Kaptive. Kaptive is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by the Free Software Foundation,
either version 3 of the License, or (at your option) any later version. Kaptive is distributed
in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
details. You should have received a copy of the GNU General Public License along with Kaptive.
If not, see <https://www.gnu.org/licenses/>.
"""
from __future__ import annotations

from functools import cached_property
from pathlib import Path
from typing import Generator, TextIO
from itertools import chain
import re
from warnings import catch_warnings

from Bio import SeqIO
from Bio.SeqFeature import SeqFeature
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq

from kaptive.log import log, quit_with_error, warning
from kaptive.misc import check_file


# Constants -----------------------------------------------------------------------------------------------------------
_LOCUS_REGEX = re.compile(r'(?<=locus:)\w+|(?<=locus: ).*')
_TYPE_REGEX = re.compile(r'(?<=type:)\w+|(?<=type: ).*')
_DB_KEYWORDS = {
    'Klebsiella_k_locus_primary_reference': ['kpsc_k', 'kp_k', 'k_k'],
    'Klebsiella_k_locus_variant_reference': ['kpsc_k_variant', 'kp_k_variant', 'k_k_variant'],
    'Klebsiella_o_locus_primary_reference': ['kpsc_o', 'kp_o', 'k_o'],
    'Acinetobacter_baumannii_k_locus_primary_reference': ['ab_k'],
    'Acinetobacter_baumannii_OC_locus_primary_reference': ['ab_o']
}
_GENE_THRESHOLDS = {
    'Klebsiella_k_locus_primary_reference': 82.5,
    'Klebsiella_k_locus_variant_reference': 82.5,
    'Klebsiella_o_locus_primary_reference': 82.5,
    'Acinetobacter_baumannii_k_locus_primary_reference': 85,
    'Acinetobacter_baumannii_OC_locus_primary_reference': 85
}


# Classes -------------------------------------------------------------------------------------------------------------
class DatabaseError(Exception):
    pass


class Database(object):
    def __init__(self, path: Path | None = None, name: str | None = None, loci: dict[str, Locus] | None = None,
                 genes: dict[str, Gene] | None = None, is_elements: dict[str, Gene] | None = None,
                 extra_loci: dict[str, Locus] | None = None, extra_genes: dict[str, Gene] | None = None,
                 gene_threshold: float | None = None):
        self.path = path or Path()
        self.name = name or self.path.stem
        self.loci = loci or {}
        self.extra_loci = extra_loci or {}
        self.genes = genes or {}
        self.extra_genes = extra_genes or {}
        self.is_elements = is_elements or {}
        self.gene_threshold = gene_threshold or _GENE_THRESHOLDS.get(self.name, 0)

    def __repr__(self):
        return (f"{self.name} ({len(self.loci)} Loci) ({len(self.genes)} Genes) ({len(self.extra_loci)} Extra Loci) "
                f"({len(self.extra_genes)} Extra Genes)")

    def __str__(self) -> str:
        return self.name

    def __len__(self) -> int:
        return len(self.loci)

    def __iter__(self):
        return chain(self.loci.values(), self.extra_loci.values())

    def __getitem__(self, item: str | int) -> Locus | Gene:
        if isinstance(item, int):
            if not 0 <= item < len(self):
                raise DatabaseError(f'Index {item} out of range for database {self.name}')
            return list(self.loci.values())[item]
        if not (result := self.loci.get(item, self.extra_loci.get(item, self.genes.get(item, self.extra_genes.get(item, self.is_elements.get(item)))))):
            raise DatabaseError(f'Could not find {item} in database {self.name}')
        return result

    @cached_property
    def largest_locus(self) -> Locus:
        return max(self.loci.values(), key=len)

    def format(self, format_spec):
        # f"##gff-version 3\n{''.join([i.as_gff_string() for i in self.loci.values()])}"
        if format_spec in {'fna', 'ffn', 'faa'}:
            return ''.join([locus.format(format_spec) for locus in self])
        raise ValueError(f'Invalid format specifier: {format_spec}')

    def add_locus(self, locus: Locus):
        """
        Adds a locus and its genes to the database. Checks that the locus and genes don't already exist in the database.
        """
        locus_dict, gene_dict = (self.loci, self.genes) if not locus.extra() else (self.extra_loci, self.extra_genes)
        if locus.name in locus_dict:
            raise DatabaseError(f'Locus {locus.name} already exists in database {self.name}.')
        locus_dict[locus.name] = locus
        locus.db = self
        for gene in locus:
            if gene.name in gene_dict:
                raise DatabaseError(f'Gene {gene} already exists in database {self.name}.')
            gene_dict[gene.name] = gene

    def add_phenotype(self, loci: list[str], genes: dict[str, str], phenotype: str):
        extra_genes = {(g.name, 'present') for g in self.extra_genes.values() if g.gene_name in genes}
        for locus in (self.loci.keys() if loci == ["ALL"] else loci):
            if locus not in self.loci:
                raise PhenotypeError(f'Could not find {locus} in database {self.name}')
            if extra_genes:
                self.loci[locus].add_phenotype(None, extra_genes, phenotype)
            else:
                self.loci[locus].add_phenotype(genes, None, phenotype, strict=loci != ["ALL"])


class LocusError(Exception):
    pass


class PhenotypeError(Exception):
    pass


class Locus(object):
    def __init__(self, name: str | None = None, seq: Seq | None = Seq(''), genes: dict[str: Gene] | None = None,
                 type_label: str | None = None,
                 phenotypes: list[tuple[set[tuple[str, str]], str]] | None = None):
        self.name = name or ''
        self.seq = seq
        self._length = len(self.seq)
        self.genes = genes or {}
        self.type_label = type_label or ''
        self.phenotypes = phenotypes or []

    @classmethod
    def from_seqrecord(cls, record: SeqRecord, locus_name: str, type_name: str, load_seq: bool = True,
                       extract_translations: bool = False):
        if load_seq:
            self = cls(name=locus_name, seq=record.seq)
        else:
            self = cls(name=locus_name)
            self._length = len(record.seq)  # We are not loading the sequence, so we need to set the length manually
        n = 1
        for feature in record.features:  # type: SeqFeature
            if feature.type == 'CDS':
                gene = Gene.from_feature(record, feature, position_in_locus=n, locus=self)
                if gene.name in self.genes:
                    raise LocusError(f'Gene {gene} already exists in locus {self}')
                if gene.locus and gene.locus != self:
                    raise LocusError(f'Gene {gene} is from a different locus than locus {self}')
                if extract_translations:  # Force translation of the gene
                    gene.extract_translation()
                self.genes[gene.name] = gene
                n += 1
        self.type_label = type_name if not self.extra() else None  # Extra genes don't have a type
        return self

    def __hash__(self):
        return hash(self.name)  # The name of the locus is unique and immutable

    def __repr__(self):
        return self.name

    def __len__(self):
        return self._length or len(self.seq)

    def __getitem__(self, item) -> Gene:
        if not (result := self.genes.get(item)):
            raise LocusError(f'Could not find {item} in locus {self.name}')
        return result

    def __iter__(self):
        return iter(self.genes.values())

    def extra(self) -> bool:
        return self.name.startswith('Extra_genes')

    def add_phenotype(self, genes: dict[str, str] | None, extra_genes: set[tuple[str, str]] | None, phenotype: str, strict: bool = False):
        if extra_genes:
            self.phenotypes = sorted(self.phenotypes + [(extra_genes, phenotype)], key=lambda x: len(x[0]), reverse=True)
        else:  # Turn gene names into a set of tuples with the gene name and state
            genes = {(g.name, state) for g in self if (state := genes.get('ALL', genes.get(g.gene_name, genes.get(g.name, None))))}
            if genes:
                self.phenotypes = sorted(self.phenotypes + [(genes, phenotype)], key=lambda x: len(x[0]), reverse=True)
            elif strict:  # If strict, raise an error
                raise PhenotypeError(f"Phenotype ({phenotype}) based on ({genes}) does not apply to {self}")

    def format(self, format_spec):
        if format_spec == 'fna':
            if len(self.seq) == 0:
                warning(f'No DNA sequence for {self}')
                return ""
            return f'>{self.name}\n{self.seq}\n'
        if format_spec in {'ffn', 'faa'}:
            return ''.join([gene.format(format_spec) for gene in self])
        raise ValueError(f'Invalid format specifier: {format_spec}')

    def write(self, fna: Path | TextIO | None = None, ffn: Path | TextIO | None = None, faa: Path | TextIO | None = None):
        """Write the typing result to files or file handles."""
        [(fh / f'{self.name.replace("/", "_")}.{fmt}').write_text(self.format(fmt)) if isinstance(fh, Path) else
         fh.write(self.format(fmt)) for fh, fmt in [(fna, 'fna'), (ffn, 'ffn'), (faa, 'faa')] if fh]

    # def as_gff_record(self) -> GffRecord:
    #     return GffRecord(seqid=self.name, source='Kaptive', type_='region', start=1, end=len(self), score=0,
    #                      strand='+', phase=0, attributes={"ID": f"locus:{self.name}", "Name": self.name})
    #
    # def as_gff_string(self) -> str:
    #     return str(self.as_gff_record()) + ''.join([i.as_gff_string() for i in self.genes.values()])


class GeneError(Exception):
    pass


# TODO: consider renaming this to CDS
class Gene(object):
    """
    This class prepares and stores a CDS feature from a Kaptive reference genbank file.
    It is designed so that the Feature itself doesn't need to be stored, only the information required to
    extract it from the record.
    """
    def __init__(self, name: str | None = None, locus: Locus | None = None, position_in_locus: int | None = 0,
                 start: int | None = 0, end: int | None = 0, strand: str | None = None, protein_seq: Seq | None = None,
                 dna_seq: Seq | None = None, gene_name: str | None = None, product: str | None = None):
        self.name = name or ''
        self.locus = locus  # Keep reference to parent class for now
        self.position_in_locus = position_in_locus
        self.start = start  # 0-based
        self.end = end
        self.strand = strand  # Either + or -
        self.gene_name = gene_name or ''
        self.name = name or ''
        self.product = product or ''  # Can also be description
        self.dna_seq = dna_seq or Seq('')
        self.protein_seq = protein_seq or Seq('')

    @classmethod
    def from_feature(cls, record: SeqRecord, feature: SeqFeature, **kwargs):
        self = cls(
            start=feature.location.start, end=feature.location.end, strand='+' if feature.location.strand == 1 else '-',
            dna_seq=feature.extract(record.seq), product=feature.qualifiers.get('product', [''])[0], **kwargs)
        self.name = f"{self.locus.name}_{str(self.position_in_locus).zfill(2)}" + (f"_{x}" if (x := feature.qualifiers.get('gene', [''])[0]) else '')
        self.gene_name = x
        assert len(self.dna_seq) % 3 == 0, quit_with_error(f'DNA sequence of {self} is not a multiple of 3')
        return self

    @classmethod
    def from_seq(cls, record: tuple[str, str, str], **kwargs):
        """This function allows the creation of a gene from a fasta file. Currently used for storing IS elements"""
        return cls(name=record[0], product=record[1], dna_seq=Seq(record[2]), **kwargs)

    def __hash__(self):
        return hash(self.name)  # The name of the gene is unique and immutable

    def __repr__(self):
        return self.name

    def __len__(self):
        return len(self.dna_seq)

    def format(self, format_spec):
        if format_spec == 'ffn':
            if len(self.dna_seq) == 0:
                warning(f'No DNA sequence for {self}')
                return ""
            return f'>{self.name}\n{self.dna_seq}\n'
        if format_spec == 'faa':
            self.extract_translation()
            if len(self.protein_seq) == 0:
                warning(f'No protein sequence for {self.__repr__()}')
                return ""
            return f'>{self.name}\n{self.protein_seq}\n'
        raise ValueError(f'Invalid format specifier: {format_spec}')

    def extra(self) -> bool:
        return self.name.startswith('Extra_genes')

    def extract_translation(self, **kwargs):
        """
        Extracts the protein sequence from the DNA sequence of the gene. Implemented as a method so unnecessary
        translations are not performed.
        :param table: NCBI translation table number
        :param cds: if True, only translates the CDS
        :param to_stop: if True, stops translation at the first stop codon
        :param gap: gap character
        :param stop_symbol: stop codon character
        """
        if len(self.protein_seq) == 0:  # Only translate if the protein sequence is not already stored
            if len(self.dna_seq) == 0:
                raise GeneError(f'No DNA sequence for reference {self}')
            with catch_warnings(record=True) as w:
                self.protein_seq = self.dna_seq.translate(**kwargs)
                # for i in w:
                #     warning(f"{i.message}: {self.__repr__()}")
            if len(self.protein_seq) == 0:
                warning(f'No protein sequence for reference {self}')

    # def as_gff_record(self, ensembl_format: bool = False) -> list[GffRecord]:
    #     """
    #     Returns the gene, transcript and CDS GFF Records compatible with bcftools csq.
    #     See format specification here: https://samtools.github.io/bcftools/bcftools.html#csq
    #     """
    #     if ensembl_format:  # for bcftools csq
    #         return [
    #             GffRecord(
    #                 seqid=self.locus.name, source='Kaptive', type_='gene', start=self.start + 1,
    #                 end=self.end, score=0, strand=self.strand, phase=0,
    #                 attributes={"ID": f"gene:{self.name}", "biotype": "protein_coding", "Name": self.gene_name,
    #                             "description": self.product}),
    #             GffRecord(
    #                 seqid=self.locus.name, source='Kaptive', type_='transcript', start=self.start + 1,
    #                 end=self.end, score=0, strand=self.strand, phase=0,
    #                 attributes={"ID": f"transcript:{self.name}.t1", "Parent": f"gene:{self.name}",
    #                             "biotype": "protein_coding"}),
    #             GffRecord(seqid=self.locus.name, source='Kaptive', type_='CDS', start=self.start + 1,
    #                       end=self.end, score=0, strand=self.strand, phase=0,
    #                       attributes={"ID": f"CDS:{self.name}.cds1", "Parent": f"transcript:{self.name}.t1"})
    #         ]
    #     else:  # Usual NCBI format
    #         return [
    #             GffRecord(
    #                 seqid=self.locus.name, source='Kaptive', type_='gene', start=self.start + 1,
    #                 end=self.end, score=0, strand=self.strand, phase=0,
    #                 attributes={"ID": self.name, "Name": self.gene_name, "gene": self.gene_name, "locus_tag": self.name}
    #             ),
    #             GffRecord(
    #                 seqid=self.locus.name, source='Kaptive', type_='CDS', start=self.start + 1,
    #                 end=self.end, score=0, strand=self.strand, phase=0,
    #                 attributes={"ID": self.name, "Name": self.gene_name, "gene": self.gene_name, "Parent": self.name,
    #                             "product": self.product, "transl_table": '11', "locus_tag": self.name}),
    #         ]
    #
    # def as_gff_string(self) -> str:
    #     return ''.join([str(i) for i in self.as_gff_record()])


# class GffRecordError(Exception):
#     pass
#
#
# class GffRecord:
#     def __init__(self, seqid: str | None = None, source: str | None = None,
#                  type_: str | None = None, start: int | None = None, end: int | None = None, score: float | None = 0,
#                  strand: str | None = None, phase: float | None = 0, attributes: dict | None = None):
#         # https://github.com/The-Sequence-Ontology/Specifications/blob/master/gff3.md
#
#         self.seqid = seqid or ''
#         self.source = source or ''
#         self.type = type_ or ''
#         self.start = start or 0
#         self.end = end or 0
#         self.score = score or 0
#         self.strand = strand or ''
#         self.phase = phase or 0
#         self.attributes = attributes or {}
#
#     def __str__(self):
#         return f"{self.seqid}\t{self.source}\t{self.type}\t{self.start}\t{self.end}\t{self.score}\t{self.strand}\t" \
#                f"{self.phase}\t{';'.join(f'{key}={value}' for key, value in self.attributes.items() if value)}\n"


# Functions ------------------------------------------------------------------------------------------------------------
def name_from_record(record: SeqRecord, locus_regex: re.Pattern | None = _LOCUS_REGEX,
                     type_regex: re.Pattern | None = _TYPE_REGEX) -> tuple[str | None, str | None]:
    """
    This function extracts the locus and type names from a genbank record using regular expressions.
    If the locus_regex or type_regex are not provided, the default regexes are used.
    """

    locus_regex, type_regex = locus_regex or _LOCUS_REGEX, type_regex or _TYPE_REGEX  # If None, use the default regexes
    locus_name, type_name = set(), set()

    if not (source := next((f for f in record.features if f.type == 'source'), None)):
        quit_with_error(f'Could not find source feature in genbank record: {record.id}')
    if "note" not in source.qualifiers:
        quit_with_error(f'Could not find note qualifier in source feature of genbank record: {record.id}')

    for note in source.qualifiers['note']:
        if type_regex and (match := type_regex.search(note)):
            type_name.add(match.group())

        if note.startswith('Extra genes'):  # "Extra genes: gmlABD" -> "Extra_genes_gmlABD"
            locus_name.add(f"Extra_genes_{note.split(' ')[-1]}")

        elif locus_regex and (match := locus_regex.search(note)):
            locus_name.add(match.group())

    if len(locus_name) > 1:
        quit_with_error(f'Found multiple locus names in record: {record.id}\n\tNote: {source.qualifiers["note"]}')
    if len(type_name) > 1:
        quit_with_error(f'Found multiple type names in record: {record.id}\n\tNote: {source.qualifiers["note"]}')

    return locus_name.pop() if len(locus_name) == 1 else None, type_name.pop() if len(type_name) == 1 else None


def parse_logic(logic_file: Path, verbose: bool = False) -> Generator[tuple[list[str], dict[str, str], str], None, None]:
    log(f'Parsing logic {logic_file}', verbose=verbose)
    with logic_file.open() as f:
        if (line := f.readline()) != 'loci\tgenes\tphenotype\n':
            quit_with_error(f'Logic file {logic_file} has invalid header: {line}')
        for line in f:
            loci, genes, phenotype = line.strip().split('\t')
            yield loci.split(';'), dict(gene.split(",", 1) if "," in gene else (gene, 'present') for gene in genes.split(';')), phenotype


def get_database(argument: str) -> Path:
    """
    Returns the path to the database file.
    If an existing file is passed, it is returned, otherwise it will be treated as a keyword and used to
    find the respective database in the kaptive package.
    """
    if (db_path := Path(argument)).is_file():
        return check_file(db_path)

    expected_path = (Path(__file__).parent.parent / "reference_database")
    if not (dbs_in_package := list(expected_path.glob('*.gbk'))):
        quit_with_error(f'No databases found in expected path: {expected_path}')

    for db in dbs_in_package:
        if argument.lower() == db.name.lower():
            return db
        if argument.lower() == db.stem.lower():
            return db
        if argument.lower() in [i.lower() for i in _DB_KEYWORDS[db.stem]]:  # Check for keywords
            return db

    quit_with_error(f'No database found for {argument}\n'
                    f'Available databases: {", ".join([i.stem for i in dbs_in_package])}\n'
                    f'Valid keywords: {", ".join([i for x in _DB_KEYWORDS.values() for i in x])}')


def parse_database(genbank: Path, locus_filter: re.Pattern | None = None, load_locus_seqs: bool = True,
                   extract_translations: bool = False, verbose: bool = False, **kwargs) -> Generator[Locus, None, None]:
    """
    Wrapper around SeqIO.parse to parse a Kaptive database genbank file and return a generator of Locus objects
    """
    log(f'Parsing {genbank.name}', verbose=verbose)
    try:
        for record in SeqIO.parse(genbank, 'genbank'):
            locus_name, type_name = name_from_record(record, **kwargs)
            if not locus_name:
                quit_with_error(f'Could not parse locus name from {record.id}')
            if type_name == "unknown" or (not type_name and not locus_name.startswith('Extra')):
                type_name = f'unknown ({locus_name})'  # Add the locus name to the type name if it is unknown
            if locus_filter and not locus_filter.search(locus_name):
                continue
            yield Locus.from_seqrecord(record, locus_name, type_name, load_locus_seqs, extract_translations)
    except Exception as e:
        quit_with_error(f'Could not parse database {genbank.name}: {e}')


def load_database(argument: str | Path, gene_threshold: float | None = None, **kwargs):
    db = Database(path=get_database(argument), gene_threshold=gene_threshold)
    for locus in parse_database(db.path, **kwargs):
        db.add_locus(locus)
    if not db.loci:  # Check that loci were properly loaded
        quit_with_error(f'No loci found in database {db.name}')
    if (logic_file := db.path.with_suffix(".logic")).is_file():  # Load phenotype logic
        [db.add_phenotype(*i) for i in parse_logic(logic_file)]
    # if is_elements:  # Load IS elements
    #     self.is_elements |= {(x := Gene.from_seq(i, strand='+', gene_name='tnp')).name: x for i in parse_fasta(is_elements)}
    return db
