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

from pathlib import Path
from typing import Generator
import re
from itertools import chain

from Bio import SeqIO
from Bio.SeqFeature import SeqFeature
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq

from kaptive.log import quit_with_error
from kaptive.misc import check_file

# Constants -----------------------------------------------------------------------------------------------------------
LOCUS_REGEX = r'(?<=locus:)\w+|(?<=locus: ).*'
TYPE_REGEX = r'(?<=type:)\w+|(?<=type: ).*'
DB_KEYWORDS = {
    'Klebsiella_k_locus_primary_reference': ['kpsc_k', 'kp_k', 'k_k'],
    'Klebsiella_k_locus_variant_reference': ['kpsc_k_variant', 'kp_k_variant', 'k_k_variant'],
    'Klebsiella_o_locus_primary_reference': ['kpsc_o', 'kp_o', 'k_o'],
    'Acinetobacter_baumannii_k_locus_primary_reference': ['ab_k'],
    'Acinetobacter_baumannii_OC_locus_primary_reference': ['ab_o'],
    'VibrioPara_Kaptivedb_K': ['vp_k'],
    'VibrioPara_Kaptivedb_O': ['vp_o'],
}


# Classes -------------------------------------------------------------------------------------------------------------
class DatabaseError(Exception):
    pass


class Database(object):
    def __init__(self, path: Path | None = None, name: str | None = None, loci: dict[str, Locus] | None = None,
                 genes: dict[str, Gene] | None = None, is_elements: dict[str, Gene] | None = None,
                 # extra_loci: dict[str, Locus] | None = None, extra_genes: dict[str, Gene] | None = None
                 ):
        self.path = path or Path()
        self.name = name or self.path.stem.replace('_', ' ')
        self.loci = loci or {}
        # self.extra_loci = extra_loci or {}
        self.genes = genes or {}
        # self.extra_genes = extra_genes or {}
        self.is_elements = is_elements or {}

    @classmethod
    def from_genbank(cls, path: Path, is_elements: Path | None = None, *args, **kwargs):
        for locus in parse_database((self := cls(path)).path, *args, **kwargs):
            self.add_locus(locus)
        if not self.loci:  # Check that loci were properly loaded
            raise DatabaseError(f'No loci found in database {self.name}')
        for locus in self.loci.values():  # Load special logic
            if not locus.name.startswith('Extra') and locus.type == 'special logic':
                for locus_name, genes, type_name in load_special_logic(self.path):
                    if locus_name not in self.loci:
                        raise DatabaseError(f'Special logic locus {locus_name} not found in {self.name}')
                    else:
                        self.loci[locus_name].special_logic[type_name] = {} if genes == "None" else set(
                            genes.split(','))
                break
        if is_elements:
            try:
                self.is_elements |= {(x := Gene.from_seq(i, strand='+', gene_name='tnp')
                                      ).name: x for i in SeqIO.parse(is_elements, 'fasta')}
            except ValueError as e:
                raise DatabaseError(f'Could not parse IS element file {is_elements}: {e}')
        return self

    def __repr__(self):
        return self.name

    def __str__(self):
        return f"{self.name} - {len(self.loci)} Loci; {len(self.genes)} Genes; IS Elements {len(self.is_elements)}"

    def __len__(self):
        return len(self.loci)

    def __iter__(self):
        for locus_name, locus in self.loci.items():
            yield locus_name, locus

    def __getitem__(self, item):
        if item not in self.loci:
            raise KeyError(f'No locus named {item} in database {self.name}')
        return self.loci[item]

    def add_locus(self, locus: Locus):
        """
        Adds a locus and its genes to the database. Checks that the locus and genes don't already exist in the database.
        """
        # locus_dict, gene_dict = (self.loci, self.genes) if not locus.name.startswith('Extra_genes') else (
        #     self.extra_loci, self.extra_genes)
        if locus.name in self.loci:
            raise DatabaseError(f'Locus {locus.name} already exists in database {self.name}.')
        self.loci[locus.name] = locus
        locus.db = self
        for gene_name, gene in locus:
            if gene_name in self.genes:
                raise DatabaseError(f'Gene {gene_name} already exists in database {self.name}.')
            self.genes[gene_name] = gene

    def as_fasta(self) -> str:
        if not self.loci:
            raise DatabaseError(f'No loci in database {self.name}')
        return ''.join(locus.as_fasta() for locus in self.loci.values())

    def as_gene_fasta(self) -> str:
        if not self.genes:
            raise DatabaseError(f'No genes in database {self.name}')
        return ''.join(gene.as_fasta() for gene in chain(self.genes.values(), self.is_elements.values()))

    def as_protein_fasta(self) -> str:
        if not self.genes:
            raise DatabaseError(f'No genes in database {self.name}')
        return ''.join(gene.as_protein_fasta() for gene in chain(self.genes.values(), self.is_elements.values()))

    def as_gff_string(self) -> str:
        if not self.loci:
            raise DatabaseError(f'No loci in database {self.name}')
        return f"##gff-version 3\n{''.join([i.as_gff_string() for i in self.loci.values()])}"


class LocusError(Exception):
    pass


class Locus(object):
    def __init__(self, name: str | None = None, type_: str | None = None, seq: Seq | None = None,
                 genes: dict[str: Gene] | None = None, db: Database | None = None,
                 special_logic: dict | None = None, length: int | None = None):
        self.name = name or ''
        self.db = db
        self.type = type_ or ''
        self.seq = seq or Seq('')
        self.length = length or len(self.seq)
        self.genes = genes or {}
        self.special_logic = special_logic or {}

    @classmethod
    def from_seqrecord(cls, record: SeqRecord, locus_name: str, type_name: str, load_seq: bool = True):
        if load_seq:
            self = cls(name=locus_name, type_=type_name, seq=record.seq)
        else:
            self = cls(name=locus_name, type_=type_name, length=len(record.seq))
        [self.add_gene(Gene.from_feature(record, feature)) for feature in record.features if feature.type == 'CDS']
        return self

    def __hash__(self):
        return hash(self.name)  # The name of the locus is unique and immutable

    def __repr__(self):
        return self.name

    def __len__(self):
        return self.length

    def __iter__(self):
        for gene_name, gene in self.genes.items():
            yield gene_name, gene

    def add_gene(self, gene: Gene):
        if gene.name in self.genes:
            raise LocusError(f'Gene {gene.name} already exists in locus {self.name}')
        if gene.locus and gene.locus != self:
            raise LocusError(f'Gene {gene.name} is from a different locus than locus {self.name}')
        else:
            gene.locus = self
        self.genes[gene.name] = gene

    def as_fasta(self):
        """Returns the locus sequence as a FASTA string."""
        if len(self.seq) == 0:
            raise LocusError(f'Locus {self.name} has no sequence.')
        return f'>{self.name}\n{self.seq}\n'

    def as_gene_fasta(self) -> str:
        if not self.genes:
            raise LocusError(f'Locus {self.name} has no genes.')
        return ''.join(gene.as_fasta() for gene in self.genes.values())

    def as_protein_fasta(self) -> str:
        if not self.genes:
            raise LocusError(f'Locus {self.name} has no genes.')
        return ''.join(gene.as_protein_fasta() for gene in self.genes.values())

    def as_gff_record(self) -> GffRecord:
        return GffRecord(seqid=self.name, source='Kaptive', type_='region', start=1, end=self.length, score=0,
                         strand='+', phase=0, attributes={"ID": f"locus:{self.name}", "Name": self.name})

    def as_gff_string(self) -> str:
        return str(self.as_gff_record()) + ''.join([i.as_gff_string() for i in self.genes.values()])


class GeneError(Exception):
    pass


class Gene(object):
    """
    This class prepares and stores a CDS feature from a Kaptive reference genbank file.
    It is designed so that the Feature itself doesn't need to be stored, only the information required to
    extract it from the record.
    """

    # TODO: consider renaming this to CDS

    def __init__(self, name: str | None = None, locus: Locus | None = None, position_in_locus: int | None = 0,
                 start: int | None = 0, end: int | None = 0, strand: str | None = None, protein_seq: Seq | None = None,
                 dna_seq: Seq | None = None, gene_name: str | None = None, product: str | None = None):
        self.name = name or ''
        self.locus = locus
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
        self = cls(start=feature.location.start, end=feature.location.end, strand='+' if feature.strand == 1 else '-',
                   **kwargs)
        if 'locus_tag' not in feature.qualifiers:
            raise GeneError(f'{feature} does not have a locus tag.')
        self.name = feature.qualifiers['locus_tag'][0]
        self.gene_name = feature.qualifiers['gene'][0] if 'gene' in feature.qualifiers else ''
        self.product = feature.qualifiers['product'][0] if 'product' in feature.qualifiers else ''
        self.dna_seq = feature.extract(record.seq)
        if 'translation' in feature.qualifiers:
            self.protein_seq = Seq(feature.qualifiers['translation'][0])
        return self

    @classmethod
    def from_seq(cls, record: SeqRecord, **kwargs):
        """This function allows the creation of a gene from a fasta file. Currently used for storing IS elements"""
        return cls(name=record.id, product=record.description, dna_seq=record.seq, **kwargs)

    def __hash__(self):
        return hash(self.name)

    def __repr__(self):
        return self.name

    def __len__(self):
        return len(self.dna_seq)

    def extract_translation(self, table: int = 11, cds: bool = False, to_stop: bool = True, gap: str = '-',
                            stop_symbol: str = '*'):
        """
        Extracts the protein sequence from the DNA sequence of the gene. Implemented as a method so unnecessary
        translations are not performed.
        :param table: NCBI translation table number
        :param cds: if True, only translates the CDS
        :param to_stop: if True, stops translation at the first stop codon
        :param gap: gap character
        :param stop_symbol: stop codon character
        """
        # First try extracting the translation from the feature
        if len(self.protein_seq) == 0:  # Only translate if the protein sequence is not already stored
            if len(self.dna_seq) == 0:
                raise GeneError(f'No DNA sequence for reference {self}')
            self.protein_seq = self.dna_seq.translate(table=table, cds=cds, to_stop=to_stop, gap=gap,
                                                      stop_symbol=stop_symbol)
            if len(self.protein_seq) == 0:
                raise GeneError(f'No protein sequence for reference {self}')

    def as_fasta(self) -> str:
        if len(self.dna_seq) == 0:
            raise GeneError(f'No DNA sequence for reference {self}')
        return f'>{self.name}\n{self.dna_seq}\n'

    def as_protein_fasta(self) -> str:
        self.extract_translation()
        return f'>{self.name}\n{self.protein_seq}\n'

    def as_gff_record(self, ensembl_format: bool = False) -> list[GffRecord]:
        """
        Returns the gene, transcript and CDS GFF Records compatible with bcftools csq.
        See format specification here: https://samtools.github.io/bcftools/bcftools.html#csq
        """
        if ensembl_format:  # for bcftools csq
            return [
                GffRecord(
                    seqid=self.locus.name, source='Kaptive', type_='gene', start=self.start + 1,
                    end=self.end, score=0, strand=self.strand, phase=0,
                    attributes={"ID": f"gene:{self.name}", "biotype": "protein_coding", "Name": self.gene_name,
                                "description": self.product}),
                GffRecord(
                    seqid=self.locus.name, source='Kaptive', type_='transcript', start=self.start + 1,
                    end=self.end, score=0, strand=self.strand, phase=0,
                    attributes={"ID": f"transcript:{self.name}.t1", "Parent": f"gene:{self.name}",
                                "biotype": "protein_coding"}),
                GffRecord(seqid=self.locus.name, source='Kaptive', type_='CDS', start=self.start + 1,
                          end=self.end, score=0, strand=self.strand, phase=0,
                          attributes={"ID": f"CDS:{self.name}.cds1", "Parent": f"transcript:{self.name}.t1"})
            ]
        else:  # Usual NCBI format
            return [
                GffRecord(
                    seqid=self.locus.name, source='Kaptive', type_='gene', start=self.start + 1,
                    end=self.end, score=0, strand=self.strand, phase=0,
                    attributes={"ID": self.name, "Name": self.gene_name, "gene": self.gene_name, "locus_tag": self.name}
                ),
                GffRecord(
                    seqid=self.locus.name, source='Kaptive', type_='CDS', start=self.start + 1,
                    end=self.end, score=0, strand=self.strand, phase=0,
                    attributes={"ID": self.name, "Name": self.gene_name, "gene": self.gene_name, "Parent": self.name,
                                "product": self.product, "transl_table": '11', "locus_tag": self.name}),
            ]

    def as_gff_string(self) -> str:
        return ''.join([str(i) for i in self.as_gff_record()])

    def get_reference_info_json_dict(self):
        reference_dict = {}
        if self.gene_name:
            reference_dict['Gene'] = self.gene_name
        if self.product:
            reference_dict['Product'] = self.product
        # if self.ec_number:
        #     reference_dict['EC number'] = self.ec_number
        reference_dict['Nucleotide length'] = len(self.dna_seq)
        reference_dict['Protein length'] = len(self.protein_seq)
        reference_dict['Nucleotide sequence'] = self.dna_seq
        reference_dict['Protein sequence'] = self.protein_seq
        return reference_dict


class GffRecordError(Exception):
    pass


class GffRecord:
    def __init__(self, seqid: str | None = None, source: str | None = None,
                 type_: str | None = None, start: int | None = None, end: int | None = None, score: float | None = 0,
                 strand: str | None = None, phase: float | None = 0, attributes: dict | None = None):
        # https://github.com/The-Sequence-Ontology/Specifications/blob/master/gff3.md

        self.seqid = seqid or ''
        self.source = source or ''
        self.type = type_ or ''
        self.start = start or 0
        self.end = end or 0
        self.score = score or 0
        self.strand = strand or ''
        self.phase = phase or 0
        self.attributes = attributes or {}

    # @classmethod
    # def from_tsv(cls, line: str | bytes, **kwargs):
    #     if not line:
    #         raise GffRecordError("Empty line passed to GffRecord")
    #
    #     line = line.decode().strip() if isinstance(line, bytes) else line.strip()
    #     if line.startswith('#'):
    #         raise GffRecordError("Comment line passed to GffRecord")
    #
    #     self = cls(**kwargs)
    #     if len(line := line.split('\t')) < 9:
    #         raise GffRecordError(f"Gff line has less than 9 columns: {line}")
    #
    #     self.seqid = line[0]
    #     self.source = line[1]
    #     self.type = line[2]
    #     self.start = int(line[3])
    #     self.end = int(line[4])
    #     self.score = line[5]
    #     self.strand = line[6]
    #     self.phase = line[7]
    #     self.attributes = parse_attributes(line[8])
    #     return self

    def __str__(self):
        return f"{self.seqid}\t{self.source}\t{self.type}\t{self.start}\t{self.end}\t{self.score}\t{self.strand}\t" \
               f"{self.phase}\t{';'.join(f'{key}={value}' for key, value in self.attributes.items() if value)}\n"


# Functions ------------------------------------------------------------------------------------------------------------
def name_from_record(record: SeqRecord, locus_regex: re.Pattern, type_regex: re.Pattern
                     ) -> tuple[str | None, str | None]:
    locus_name, type_name = set(), set()
    source = next((f for f in record.features if f.type == 'source'), None)
    if not source:
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


def load_special_logic(genbank: Path) -> Generator[tuple[str, str, str], None, None]:
    """
    This function loads special logic file if  any of the reference loci have a type of 'special logic',
    that implies that a corresponding file exists to describe that logic.
    Logic files are tab-delimited files with 3 columns: Locus name, Extra genes, Type.
    Returns a dictionary with the locus name as the key and the extra genes and type as the value.
    """
    special_logic_filename = check_file(genbank.with_suffix(".logic"))
    lines = [i.split('\t') for i in special_logic_filename.read_text().splitlines()]
    if lines[0] != ['locus', 'extra_genes', 'type']:
        quit_with_error(f'Invalid special logic file: {special_logic_filename}')
    for line in lines[1:]:
        yield line[0], line[1], line[2]


# def parse_attributes(attr_string: str) -> dict:
#     # attr_dict = {match[0]: match[1] for match in re.findall(GFF_ATTRIBUTE_REGEX, attr_string)}
#     attr_dict = {}
#     for attr in attr_string.split(';'):
#         if len(key_value := attr.split('=')) == 2:
#             attr_dict[key_value[0]] = key_value[1]
#     return attr_dict


def join_attributes(attr_dict: dict) -> str:
    return ';'.join([f'{key}={value}' for key, value in attr_dict.items()])


def extract(args):
    if args.format in ["gbk", "ids"]:
        for record in SeqIO.parse(args.db, 'genbank'):
            locus_name, type_name = name_from_record(record, args.locus_regex, args.type_regex)
            record.name = locus_name  # Set the name of the record to the locus name
            if args.filter and not args.filter.search(locus_name):
                continue
            if args.format == 'ids':
                if type_name == "unknown" or (not type_name and not locus_name.startswith('Extra')):
                    type_name = f'unknown ({locus_name})'
                args.out.write(f"{locus_name}\t{type_name}\n")
            else:
                args.out.write(record.format('genbank'))  # Do we need to add a newline?
    else:
        for locus in parse_database(args.db, args.locus_regex, args.type_regex, args.filter,
                                    load_seq=args.format == 'loci'):
            if args.format == 'loci':
                args.out.write(locus.as_fasta())
            elif args.format == 'genes':
                args.out.write(locus.as_gene_fasta())
            elif args.format == 'proteins':
                args.out.write(locus.as_protein_fasta())
            elif args.format == 'gff':
                args.out.write(locus.as_gff_string())


def get_database(argument: str) -> Path:
    """
    Returns the path to the database file.
    If an existing file is passed, it is returned, otherwise it will be treated as a keyword and used to
    find the respective database in the kaptive package.
    """
    if (db_path := Path(argument)).is_file():
        return db_path
    if not (dbs_in_package := list((Path(__file__).parent.parent / "reference_database").glob('*.gbk'))):
        quit_with_error('No databases found in kaptive package')

    for db in dbs_in_package:
        if argument.lower() == db.name.lower():
            return db
        if argument.lower() == db.stem.lower():
            return db
        if argument.lower() in [i.lower() for i in DB_KEYWORDS[db.stem]]:
            return db

    quit_with_error(f'No database found for {argument}')


def parse_database(
        genbank: Path, locus_regex: re.Pattern = re.compile(LOCUS_REGEX),
        type_regex: re.Pattern = re.compile(TYPE_REGEX),
        locus_filter: re.Pattern | None = None, load_seq: bool = True) -> Generator[Locus, None, None]:
    """
    Wrapper around SeqIO.parse to parse a Kaptive database genbank file and return a generator of Locus objects
    """
    try:
        for record in SeqIO.parse(genbank, 'genbank'):
            locus_name, type_name = name_from_record(record, locus_regex, type_regex)
            if not locus_name:
                raise DatabaseError(f'Could not parse locus name from {record.id}')
            if type_name == "unknown" or (not type_name and not locus_name.startswith('Extra')):
                type_name = f'unknown ({locus_name})'
            if locus_filter and not locus_filter.search(locus_name):
                continue
            yield Locus.from_seqrecord(record, locus_name, type_name, load_seq)
    except ValueError as e:
        raise DatabaseError(f'Could not parse database {genbank.name}: {e}')
