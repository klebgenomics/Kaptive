from itertools import chain
from functools import cached_property
from pathlib import Path
from typing import Generator, IO, Union, Optional
import re
from warnings import catch_warnings, warn
from io import TextIOBase
from json import load as json_load

import numpy as np

from Bio import SeqIO
from Bio.SeqFeature import SeqFeature
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq

from kaptive.utils import Resources, LiteralFile, GitRepo

# Constants -----------------------------------------------------------------------------------------------------------
_LOCUS_REGEX = re.compile(r'(?<=locus:)\w+|(?<=locus: ).*')
_TYPE_REGEX = re.compile(r'(?<=type:)\w+|(?<=type: ).*')


# Exceptions and warnings ----------------------------------------------------------------------------------------------
class LocusError(Exception): pass
class PhenotypeError(Exception): pass
class DatabaseError(Exception): pass
class DatabaseManagerError(Exception): pass


# Classes -------------------------------------------------------------------------------------------------------------
class DatabaseManager:
    __slots__ = ('_databases', '_metadata', '_keywords', '_database_files', '_unit_files', '_metadata_files', '_logic_files')
    def __init__(self):
        self._databases: list[str] = []
        self._metadata: dict[str, dict] = {}
        self._keywords: dict[str, str] = {}
        self._database_files: dict[str, Path] = {}
        self._unit_files: dict[str, Path] = {}
        self._metadata_files: dict[str, Path] = {}
        self._logic_files: dict[str, Path] = {}
        for args in self._scan(self._PATH): self._load(*args)

    def __len__(self): return len(self._databases)
    def __iter__(self): return iter(self._databases)
    def __contains__(self, name: str): return name in self._databases or name in self._keywords
    def __getitem__(self, name_or_kwd: str) -> 'Database':
        """Load an installed database using its name or keyword"""
        name = self._keywords.get(name_or_kwd, name_or_kwd)
        if name not in self._database_files:
            raise DatabaseError(f'Invalid database name or keyword: {name_or_kwd}')
        db = Database.from_file(self._database_files[name], self._logic_files.get(name), self._unit_files.get(name))
        db._metadata = self._metadata[name]
        return db

    @staticmethod
    def _scan(path: Path) -> Generator[tuple[Path, Path, Optional[Path], Optional[Path]], None, None]:
        if not path or not path.is_dir(): raise DatabaseManagerError('Path must be a directory')
        for database_file in path.glob('*.gbk'):
            if (metadata_file := database_file.with_suffix('.json')) is None:
                raise FileNotFoundError(f'No metadata file found for {database_file}')
            if not (logic_file := database_file.with_suffix('.logic')).is_file(): logic_file = None
            if not (unit_file := database_file.with_suffix('.units')).is_file(): unit_file = None
            yield database_file, metadata_file, logic_file, unit_file

    def _load(self, database_file: Path, metadata_file: Path, logic_file: Path = None,
              unit_file: Path = None, add_to_package: bool = False):
            self._databases.append(database := database_file.stem)
            if add_to_package:
                database_file = self._add_to_package(database_file)
                metadata_file = self._add_to_package(metadata_file)
                logic_file = self._add_to_package(logic_file)
                unit_file = self._add_to_package(unit_file)
            self._database_files[database] = database_file
            self._metadata_files[database] = metadata_file
            if logic_file: self._logic_files[database] = logic_file
            if unit_file: self._unit_files[database] = unit_file
            self._metadata[database] = (metadata := json_load(metadata_file.open('r', encoding='utf-8')))
            for keyword in metadata['keywords']: self._keywords[keyword] = database

    def _add_to_package(self, path: Path = None) -> Optional[Path]:
        if path is None: return None
        new_path = self._PATH / path.name
        new_path.write_text(path.read_text(encoding='utf-8'), encoding='utf-8')
        return new_path

    def fetch(self, repo: GitRepo, add_to_package: bool = False):
        """Fetches databases from a GitHub repository, optionally adding them to the package"""
        for args in self._scan(repo.local_path): self._load(*args, add_to_package=add_to_package)

class Database:
    def __init__(self, name: str, loci: dict[str, Locus] = None, genes: dict[str, Gene] = None,
                 extra_genes: dict[str, Gene] = None, gene_threshold: float = 0):
        self.name = name
        self.loci = loci or {}
        self.genes = genes or {}
        self.extra_genes = extra_genes or {}
        self.gene_threshold = gene_threshold
        self._expected_gene_counts = None

    def __repr__(self):
        return f"{self.name} ({len(self.loci)} Loci) ({len(self.genes)} Genes) ({len(self.extra_genes)} Extra Genes)"

    def __str__(self) -> str:
        return self.name

    def __len__(self) -> int:
        return len(self.loci)

    def __iter__(self):
        return iter(self.loci.values())

    def __getitem__(self, item: str | int) -> Locus | Gene:
        if isinstance(item, int):
            if not 0 <= item < len(self):
                raise DatabaseError(f'Index {item} out of range for database {self.name}')
            return list(self.loci.values())[item]
        if not (
        result := self.loci.get(item, self.genes.get(item, self.extra_genes.get(item)))):
            raise DatabaseError(f'Could not find {item} in database {self.name}')
        return result

    @cached_property
    def largest_locus(self) -> Locus:
        return max(self.loci.values(), key=len)

    @property
    def expected_gene_counts(self) -> np.ndarray:
        if self._expected_gene_counts is None:
            self._expected_gene_counts = np.array([len(l.genes) for l in self.loci.values()])
        return self._expected_gene_counts

    def format(self, format_spec):
        if format_spec in {'fna'}:
            return ''.join([locus.format(format_spec) for locus in self.loci.values()])
        elif format_spec in {'ffn', 'faa'}:
            return ''.join([gene.format(format_spec) for gene in chain(self.genes.values(), self.extra_genes.values())])
        else:
            raise ValueError(f'Invalid format specifier: {format_spec}')

    def add_locus(self, locus: Locus):
        """
        Adds a locus and its genes to the database. Checks that the locus and genes don't already exist in the database.
        """
        if not locus.name.startswith('Extra_genes'):
            if locus.name in self.loci:
                raise DatabaseError(f'Locus {locus.name} already exists in database {self.name}.')
            self.loci[locus.name] = locus
            for gene in locus:
                if gene.name in self.genes:
                    raise DatabaseError(f'Gene {gene} already exists in database {self.name}.')
                self.genes[gene.name] = gene
        else:
            for gene in locus:
                if gene.name in self.extra_genes:
                    raise DatabaseError(f'Gene {gene} already exists in database {self.name}.')
                self.extra_genes[gene.name] = gene

    def add_phenotype(self, loci: list[str], genes: dict[str, str], phenotype: str):
        extra_genes = {(g.name, 'present') for g in self.extra_genes.values() if g.gene_name in genes}
        for locus in (self.loci.keys() if loci == ["ALL"] else loci):
            if locus in self.loci:
                if extra_genes:
                    self.loci[locus].add_phenotype(None, extra_genes, phenotype)
                else:
                    self.loci[locus].add_phenotype(genes, None, phenotype, strict=loci != ["ALL"])
            # else:
            #     raise PhenotypeError(f'Could not find {locus} in database {self.name}')


class Locus:
    def __init__(self, name: str = None, seq: Seq | None = Seq(''), genes: dict[str, Gene] = None,
                 type_label: str = None, phenotypes: list[tuple[set[tuple[str, str]], str]] = None,
                 index: int | None = 0):
        self.name = name or ''
        self.seq = seq
        self._length = len(self.seq)
        self.genes = genes or {}
        self.type_label = type_label or ''
        self.phenotypes = phenotypes or []
        self.index = index

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

                name = f"{locus_name}_{str(n).zfill(2)}" + (
                    f"_{gene_name}" if (gene_name := feature.qualifiers.get('gene', [''])[0]) else '')

                gene = Gene(
                    name=name, gene_name=gene_name, dna_seq=feature.extract(record.seq), start=feature.location.start,
                    end=feature.location.end, strand='+' if feature.location.strand == 1 else '-',
                    product=feature.qualifiers.get('product', [''])[0]
                )

                if not len(gene.dna_seq) % 3 == 0:  # Check the gene is a multiple of 3 (complete CDS)
                    # TODO: this is quite strict, but enforces the inclusion of complete CDS in Kaptive databases
                    raise GeneError(f'DNA sequence of {name} is not a multiple of 3')

                if name in self.genes:
                    raise LocusError(f'Gene {name} already exists in locus {self}')

                if extract_translations:  # Force translation of the gene
                    gene.extract_translation()

                self.genes[name] = gene
                n += 1
        self.type_label = type_name if not self.name.startswith('Extra_genes') else None  # Extra genes don't have a type
        return self

    def __hash__(self):  # TODO: Check if this is used at all
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

    def add_phenotype(self, genes: dict[str, str], extra_genes: set[tuple[str, str]], phenotype: str,
                      strict: bool = False):
        if extra_genes:
            self.phenotypes = sorted(self.phenotypes + [(extra_genes, phenotype)], key=lambda x: len(x[0]),
                                     reverse=True)
        else:  # Turn gene names into a set of tuples with the gene name and state
            genes = {(g.name, state) for g in self if
                     (state := genes.get('ALL', genes.get(g.gene_name, genes.get(g.name, None))))}
            if genes:
                self.phenotypes = sorted(self.phenotypes + [(genes, phenotype)], key=lambda x: len(x[0]), reverse=True)
            elif strict:  # If strict, raise an error
                raise PhenotypeError(f"Phenotype ({phenotype}) based on ({genes}) does not apply to {self}")

    def format(self, format_spec):
        if format_spec == 'fna':
            if len(self.seq) == 0:
                warn(f'No DNA sequence for {self}', UserWarning)
                return ""
            return f'>{self.name}\n{self.seq}\n'
        if format_spec in {'ffn', 'faa'}:
            return ''.join([gene.format(format_spec) for gene in self])
        raise ValueError(f'Invalid format specifier: {format_spec}')

    def write(self, fna: Union[str, Path, IO[str]] = None, ffn: Union[str, Path, IO[str]] = None,
              faa: Union[str, Path, IO[str]] = None):
        """Write the typing result to files or file handles."""
        for f, fmt in [(fna, 'fna'), (ffn, 'ffn'), (faa, 'faa')]:
            if f:
                if isinstance(f, TextIOBase):
                    f.write(self.format(fmt))
                elif isinstance(f, (Path, str)):
                    Path(f / f'{self.name.replace("/", "_")}.{fmt}').write_text(self.format(fmt))


class GeneError(Exception):
    pass


class Gene:
    """
    This class prepares and stores a CDS feature from a Kaptive reference genbank file.
    It is designed so that the Feature itself doesn't need to be stored, only the information required to
    extract it from the record.
    """

    def __init__(self, name: str, start: int = 0, end: int = 0, strand: str = "+",
                 protein_seq: Seq = None, dna_seq: Seq = None, gene_name: str = None, product: str = None):
        self.name = name or ''
        self.start = start  # 0-based
        self.end = end
        self.strand = strand  # Either + or -
        self.gene_name = gene_name or ''
        self.product = product or ''  # Can also be description
        self.dna_seq = dna_seq or Seq('')
        self.protein_seq = protein_seq or Seq('')

    def __hash__(self):
        return hash(self.name)  # The name of the gene is unique and immutable

    def __repr__(self):
        return self.name

    def __len__(self):
        return len(self.dna_seq)

    def format(self, format_spec):
        if format_spec == 'ffn':
            if len(self.dna_seq) == 0:
                warn(f'No DNA sequence for {self}', UserWarning)
                return ""
            return f'>{self.name}\n{self.dna_seq}\n'
        if format_spec == 'faa':
            self.extract_translation()
            if len(self.protein_seq) == 0:
                warn(f'No protein sequence for {self.__repr__()}', UserWarning)
                return ""
            return f'>{self.name}\n{self.protein_seq}\n'
        raise ValueError(f'Invalid format specifier: {format_spec}')

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
                warn(f'No protein sequence for reference {self}', UserWarning)


# Functions ------------------------------------------------------------------------------------------------------------
def name_from_record(record: SeqRecord, locus_regex: re.Pattern = _LOCUS_REGEX,
                     type_regex: re.Pattern = _TYPE_REGEX) -> tuple[Optional[str], Optional[str]]:
    """
    This function extracts the locus and type names from a genbank record using regular expressions.
    If the locus_regex or type_regex are not provided, the default regexes are used.
    """

    locus_regex, type_regex = locus_regex or _LOCUS_REGEX, type_regex or _TYPE_REGEX  # If None, use the default regexes
    locus_name, type_name = set(), set()

    if not (source := next((f for f in record.features if f.type == 'source'), None)):
        raise RuntimeError(f'Could not find source feature in genbank record: {record.id}')
    if "note" not in source.qualifiers:
        raise RuntimeError(f'Could not find note qualifier in source feature of genbank record: {record.id}')

    for note in source.qualifiers['note']:
        if type_regex and (match := type_regex.search(note)):
            type_name.add(match.group())

        if note.startswith('Extra genes'):  # "Extra genes: gmlABD" -> "Extra_genes_gmlABD"
            locus_name.add(f"Extra_genes_{note.split(' ')[-1]}")

        elif locus_regex and (match := locus_regex.search(note)):
            locus_name.add(match.group())

    if len(locus_name) > 1:
        raise RuntimeError(f'Found multiple locus names in record: {record.id}\n\tNote: {source.qualifiers["note"]}')
    if len(type_name) > 1:
        raise RuntimeError(f'Found multiple type names in record: {record.id}\n\tNote: {source.qualifiers["note"]}')

    return locus_name.pop() if len(locus_name) == 1 else None, type_name.pop() if len(type_name) == 1 else None


def parse_logic(logic_file: Union[str, Path]) -> Generator[tuple[list[str], dict[str, str], str], None, None]:
    with open(logic_file, 'rt') as f:
        if (line := f.readline()) != 'loci\tgenes\tphenotype\n':
            raise RuntimeError(f'Logic file {logic_file} has invalid header: {line}')
        for line in f:
            loci, genes, phenotype = line.strip().split('\t')
            yield loci.split(';'), dict(
                gene.split(",", 1) if "," in gene else (gene, 'present') for gene in genes.split(';')), phenotype


def get_database(argument: Union[str, Path]) -> tuple[str, Path]:
    """
    Returns the path to the database file.
    If an existing file is passed, it is returned, otherwise it will be treated as a keyword and used to
    find the respective database in the kaptive package.
    """
    
    if path := LiteralFile(argument):
        return path.stem, path

    if not (dbs_in_package := list(Resources.data.glob('*.gbk'))):
        raise RuntimeError(f'No databases found in expected path: {Resources.data}')

    # Check keywords
    for db in dbs_in_package:
        db_stem, _ = path.splitext(db)
        if argument == db_stem or argument in _DB_KEYWORDS[db_stem]:
            return db_stem, LiteralFile(path.join(Resources.data, db))

    raise RuntimeError(f'No database found for {argument}\n'
                    f'Available databases: {", ".join(dbs_in_package)}\n'
                    f'Valid keywords: {", ".join([i for x in _DB_KEYWORDS.values() for i in x])}')


def parse_database(db: Union[str, Path], locus_filter: re.Pattern = None, load_locus_seqs: bool = True,
                   extract_translations: bool = False, verbose: bool = False, **kwargs) -> Generator[Locus, None, None]:
    """
    Wrapper around SeqIO.parse to parse a Kaptive database genbank file and return a generator of Locus objects
    """
    db_name, db_path = get_database(db)
    try:
        for record in SeqIO.parse(db_path, 'genbank'):
            locus_name, type_name = name_from_record(record, **kwargs)
            if not locus_name:
                raise RuntimeError(f'Could not parse locus name from {record.id}')
            if type_name == "unknown" or (not type_name and not locus_name.startswith('Extra_genes')):
                type_name = f'unknown ({locus_name})'  # Add the locus name to the type name if it is unknown
            if locus_filter and not locus_filter.search(locus_name):
                continue
            yield Locus.from_seqrecord(record, locus_name, type_name, load_locus_seqs, extract_translations)
    except Exception as e:
        raise RuntimeError(f'Could not parse database {db_name}: {e}')


def load_database(argument: Union[str, Path], gene_threshold: float = None, **kwargs) -> Database:
    db_name, db_path = get_database(argument)
    db = Database(db_name, gene_threshold=gene_threshold)
    for locus in parse_database(db_path, **kwargs):
        db.add_locus(locus)
    if not db.loci:  # Check that loci were properly loaded
        raise RuntimeError(f'No loci found in database {db.name}')
    if path.isfile(logic_file := f'{path.splitext(db_path)[0]}.logic'):  # Load phenotype logic
        logic = list(parse_logic(logic_file))

        [db.add_phenotype(*i) for i in parse_logic(logic_file)]
    for n, locus in enumerate(db.loci.values()):
        locus.index = n
    return db
