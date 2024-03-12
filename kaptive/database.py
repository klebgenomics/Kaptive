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
from typing import Generator
import re
import warnings

from Bio import SeqIO
from Bio.SeqFeature import SeqFeature
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq

from kaptive.log import log, quit_with_error, warning
from kaptive.misc import parse_fasta

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

    @classmethod
    def from_genbank(cls, path: Path, locus_filter: re.Pattern | None = None, load_seq: bool = True,
                     gene_threshold: float | None = None, verbose: bool = False, **kwargs):
        for locus in parse_database((self := cls(path=path, gene_threshold=gene_threshold)).path, locus_filter, load_seq, verbose, **kwargs):
            self.add_locus(locus)
        if not self.loci:  # Check that loci were properly loaded
            raise DatabaseError(f'No loci found in database {self.name}')
        if (logic_file := self.path.with_suffix(".logic")).is_file():  # Load phenotype logic
            [self.add_phenotype(*i) for i in parse_logic(logic_file)]
        # if is_elements:  # Load IS elements
        #     self.is_elements |= {(x := Gene.from_seq(i, strand='+', gene_name='tnp')).name: x for i in parse_fasta(is_elements)}
        return self

    def add_phenotype(self, loci: list[str], genes: dict[str, str], phenotype: str):
        if 0 < len(extra_genes := {g.name: 'present' for g in self.extra_genes.values() if g.gene_name in genes}):
            if len(extra_genes) < len(genes):  # If there are extra genes, they must be the only genes
                raise PhenotypeError(f'Extra genes {extra_genes} combined with specific genes {genes} in phenotype {phenotype}')
            genes = {'ALL': 'present'}  # If there are only extra genes, set genes to ALL
        for locus in self.loci.keys() if loci == ["ALL"] else loci:
            if locus not in self.loci:
                raise PhenotypeError(f'Could not find {locus} in database {self.name}')
            self.loci[locus].add_phenotype(genes, extra_genes, phenotype, strict=loci != ["ALL"])

    def __repr__(self):
        return (f"{self.name} ({len(self.loci)} Loci) ({len(self.genes)} Genes) ({len(self.extra_loci)} Extra Loci) "
                f"({len(self.extra_genes)} Extra Genes) ({len(self.is_elements)} IS Elements)")

    def __str__(self) -> str:
        return self.name

    def __len__(self) -> int:
        return len(self.loci)

    def __iter__(self) -> Generator[Locus, None, None]:
        return self.get_loci()

    def __getitem__(self, item) -> Locus | Gene:
        if not (result := self.loci.get(item, self.extra_loci.get(item, self.genes.get(item, self.extra_genes.get(item, self.is_elements.get(item)))))):
            raise DatabaseError(f'Could not find {item} in database {self.name}')
        return result

    @cached_property
    def largest_locus(self) -> Locus:
        return max(self.loci.values(), key=len)

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

    def get_loci(self) -> Generator[Locus, None, None]:
        yield from self.loci.values()
        if self.extra_loci:
            yield from self.extra_loci.values()

    def get_genes(self) -> Generator[Gene, None, None]:
        yield from self.genes.values()
        if self.extra_genes:
            yield from self.extra_genes.values()
        if self.is_elements:
            yield from self.is_elements.values()

    def as_fasta(self) -> str:
        if not self.loci:
            raise DatabaseError(f'No loci in database {self.name}')
        return ''.join(locus.as_fasta() for locus in self.get_loci())

    def as_gene_fasta(self) -> str:
        return ''.join(gene.as_fasta() for gene in self.get_genes())

    def as_protein_fasta(self, **kwargs) -> str:
        return ''.join(gene.as_protein_fasta(**kwargs) for gene in self.get_genes())

    def as_gff_string(self) -> str:
        return f"##gff-version 3\n{''.join([i.as_gff_string() for i in self.loci.values()])}"


class LocusError(Exception):
    pass


class PhenotypeError(Exception):
    pass


class Locus(object):
    def __init__(self, name: str | None = None, seq: Seq | None = None, genes: dict[str: Gene] | None = None,
                 db: Database | None = None, phenotypes: list[tuple[set[tuple[str, str]], str]] | None = None,
                 length: int | None = None):
        self.name = name or ''
        self.db = db
        self.seq = seq or Seq('')
        self.length = length or len(self.seq)
        self.genes = genes or {}
        self.phenotypes = phenotypes or []

    @classmethod
    def from_seqrecord(cls, record: SeqRecord, locus_name: str, type_name: str, load_seq: bool = True):
        if load_seq:
            self = cls(name=locus_name, seq=record.seq)
        else:
            self = cls(name=locus_name, length=len(record.seq))
        [self.add_gene(Gene.from_feature(record, feature)) for feature in record.features if feature.type == 'CDS']
        if not self.extra():
            self.add_phenotype({'ALL': 'present'}, {}, type_name)  # Add the default phenotype
        return self

    def __hash__(self):
        return hash(self.name)  # The name of the locus is unique and immutable

    def __repr__(self):
        return self.name

    def __len__(self):
        return self.length

    def __getitem__(self, item) -> Gene:
        if not (result := self.genes.get(item)):
            raise LocusError(f'Could not find {item} in locus {self.name}')
        return result

    def __iter__(self) -> Generator[Gene, None, None]:
        return self.get_genes()

    def extra(self) -> bool:
        return self.name.startswith('Extra_genes')

    def get_genes(self) -> Generator[Gene, None, None]:
        yield from self.genes.values()

    def add_gene(self, gene: Gene):
        if gene.name in self.genes:
            raise LocusError(f'Gene {gene.name} already exists in locus {self.name}')
        if gene.locus and gene.locus != self:
            raise LocusError(f'Gene {gene.name} is from a different locus than locus {self.name}')
        else:
            gene.locus = self
        self.genes[gene.name] = gene

    def add_phenotype(self, genes: dict[str, str], extra_genes: dict[str, str], phenotype: str, strict: bool = False):
        """
        Genes is a dictionary of gene names/gene_names and states, e.g. {'wcaJ': 'present', 'wbaP': 'absent'}.
        Locus genes will be checked to see if their gene name/gene_name is present in the genes dictionary, and
        their full gene name will be used if it is (e.g. KL1_01_galF will be used if galF is present in the genes dictionary).
        If extra_genes is not None, it is a dictionary of gene names and states, and genes is set to {'ALL': 'present'}
        This method is aware of the ALL gene and will raise an error if it is used in conjunction with specific genes.
        """
        if not (genes := {(gene.name, state) for gene in self if (state := genes.get('ALL', genes.get(
            gene.gene_name, genes.get(gene.name, None))))}):
            if strict:
                raise PhenotypeError(f"Phenotype ({phenotype}) based on ({genes}) does not apply to {self}")
        else:
            genes |= set(extra_genes.items())   # Add extra genes to the genes set
            self.phenotypes = sorted(self.phenotypes + [(genes, phenotype)], key=lambda x: len(x[0]), reverse=True)

    def as_fasta(self):
        """Returns the locus sequence as a FASTA string."""
        if len(self.seq) == 0:
            raise LocusError(f'Locus {self.name} has no sequence.')
        return f'>{self.name}\n{self.seq}\n'

    def as_gene_fasta(self) -> str:
        if not self.genes:
            raise LocusError(f'Locus {self.name} has no genes.')
        return ''.join(gene.as_fasta() for gene in self.genes.values())

    def as_protein_fasta(self, **kwargs) -> str:
        if not self.genes:
            raise LocusError(f'Locus {self.name} has no genes.')
        return ''.join(gene.as_protein_fasta(**kwargs) for gene in self.genes.values())

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
        self = cls(start=feature.location.start, end=feature.location.end, strand='+' if feature.location.strand == 1 else '-',
                   **kwargs)
        if not (locus_tag := feature.qualifiers.get('locus_tag', [None])[0]):
            raise GeneError(f'{feature} does not have a locus tag.')
        self.name = locus_tag
        self.gene_name = feature.qualifiers.get('gene', [self.name])[0]  # Use locus tag if gene name is not present
        self.product = feature.qualifiers.get('product', [''])[0]  # Use empty string if product is not present
        self.dna_seq = feature.extract(record.seq)
        # TODO: force on-the-fly translation or use existing translation?
        #  Existing translation is faster but may be incorrect
        # if 'translation' in feature.qualifiers:
        #     self.protein_seq = Seq(feature.qualifiers['translation'][0])
        return self

    @classmethod
    def from_seq(cls, record: tuple[str, str, str], **kwargs):
        """This function allows the creation of a gene from a fasta file. Currently used for storing IS elements"""
        return cls(name=record[0], product=record[1], dna_seq=Seq(record[2]), **kwargs)

    def __hash__(self):
        return hash(self.name)

    def __repr__(self):
        return self.name

    def __len__(self):
        return len(self.dna_seq)

    def extra(self) -> bool:
        return self.name.startswith('Extra_genes')

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
            with warnings.catch_warnings(record=True) as w:
                self.protein_seq = self.dna_seq.translate(table=table, cds=cds, to_stop=to_stop, gap=gap,
                                                          stop_symbol=stop_symbol)
                # for i in w:
                #     warning(f"{i.message}: {self.__repr__()}")
            if len(self.protein_seq) == 0:
                # raise GeneError(f'No protein sequence for reference {self}\nDNA sequence: {self.dna_seq}')
                warning(f'No protein sequence for reference {self}')

    def as_fasta(self) -> str:
        if len(self.dna_seq) > 0:
            return f'>{self.name}\n{self.dna_seq}\n'
        else:
            raise GeneError(f'No DNA sequence for reference {self}')

    def as_protein_fasta(self, min_proportion: int = 0.1) -> str:
        """
        Returns the protein sequence as a FASTA string. If the proportion of the translated sequence to the DNA sequence
        is less than min_proportion, an empty string is returned.
        """
        self.extract_translation()
        if len(self.protein_seq) / (len(self.dna_seq) / 3) >= min_proportion:
            return f'>{self.name}\n{self.protein_seq}\n'
        else:
            return ''

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

    def __str__(self):
        return f"{self.seqid}\t{self.source}\t{self.type}\t{self.start}\t{self.end}\t{self.score}\t{self.strand}\t" \
               f"{self.phase}\t{';'.join(f'{key}={value}' for key, value in self.attributes.items() if value)}\n"


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
    else:  # load seqs == True if format == 'loci' else False
        for locus in parse_database(args.db, args.filter, args.format == 'loci', locus_regex=args.locus_regex,
                                    type_regex=args.type_regex):
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
    if not (dbs_in_package := list(expected_path := (Path(__file__).parent.parent / "reference_database").glob('*.gbk'))):
        quit_with_error(f'No databases found in expected path: {expected_path}')

    for db in dbs_in_package:
        if argument.lower() == db.name.lower():
            return db
        if argument.lower() == db.stem.lower():
            return db
        if argument.lower() in [i.lower() for i in _DB_KEYWORDS[db.stem]]:  # Check for keywords
            return db

    quit_with_error(f'No database found for {argument}')


def parse_database(genbank: Path, locus_filter: re.Pattern | None = None, load_seq: bool = True, verbose: bool = False,
                   **kwargs) -> Generator[Locus, None, None]:
    """
    Wrapper around SeqIO.parse to parse a Kaptive database genbank file and return a generator of Locus objects
    """
    log(f'Parsing database {genbank}', verbose=verbose)
    try:
        for record in SeqIO.parse(genbank, 'genbank'):
            locus_name, type_name = name_from_record(record, **kwargs)
            if not locus_name:
                raise DatabaseError(f'Could not parse locus name from {record.id}')
            if type_name == "unknown" or (not type_name and not locus_name.startswith('Extra')):
                type_name = f'unknown ({locus_name})'  # Add the locus name to the type name if it is unknown
            if locus_filter and not locus_filter.search(locus_name):
                continue
            yield Locus.from_seqrecord(record, locus_name, type_name, load_seq)
    except ValueError as e:
        raise DatabaseError(f'Could not parse database {genbank.name}: {e}')
