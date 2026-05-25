"""
Updated database module with dataclasses, enums and SoA batch layouts for vectorized math and machine-learning.
"""
from re import compile as re_compile
from fnmatch import filter as fnmatch_filter
from dataclasses import dataclass, field
from typing import Generator
import pickle
import tomllib
from pathlib import Path
from tempfile import TemporaryDirectory

import numpy as np
import numpy.typing as npt

from kaptive.core.seq import SeqRecord, SeqBatch, BacterialTranslationTable
from kaptive.core.interval import IntervalBatch
import sys
from kaptive import db
sys.modules['kaptive.database'] = db

from kaptive.utils import GitRepo, check_file


# Exceptions and warnings ----------------------------------------------------------------------------------------------
class DatabaseError(Exception):
    pass


# Data Models ----------------------------------------------------------------------------------------------------------
@dataclass(frozen=True, slots=True)
class DatabaseMetadata:
    """
    Strict schema for Database metadata.
    Provides dependency-free validation and ergonomic attribute access.
    """
    name: str = field(
        metadata={"description": "The name of the database, e.g. 'Klebsiella pneumoniae Species Complex K'"})
    keyword: str = field(metadata={"description": "The database keyword, e.g. 'kpsc_k'"})
    genbank: str = field(metadata={
        "description": "The name of the main database file, e.g. 'Klebsiella_pneumoniae_Species_Complex_K.gbk'"})
    organism: str = field(
        metadata={"description": "The name of the database organism, e.g. 'Klebsiella pneumoniae Species Complex'"})
    taxon: int = field(metadata={"description": "The NCBI Taxonomy ID of the database organism, e.g. 3390273"})
    antigen: str = field(metadata={"description": "The name of the database antigen, e.g. 'Capsular polysaccharide'"})
    pathway: str = field(
        metadata={"description": "The name of the database antigen synthesis pathway, e.g. 'Wzx/Wzy-dependent'"})
    version: str = field(metadata={"description": "The version of the database, e.g. '3.2.1'"})
    id_threshold: float = field(metadata={"description": "The identity threshold of the database, e.g. 82.5"})
    doi: list[str] = field(metadata={"description": "A list of DOIs associated with the database, e.g. ['TBD']"})
    owner: str = field(metadata={"description": "The owner of the database Github repo, e.g. 'klebgenomics'"})
    repo: str = field(
        metadata={"description": "The name of the database Github repo, e.g. 'KpSC_surface_antigen_loci'"})
    branch: str = field(metadata={"description": "The branch of the database Github repo, e.g. 'main'"})
    contact: dict = field(metadata={
        "description": "The details of the database curators, e.g. {'Kelly Wyres': 'kaptive.typing@gmail.com'}"})
    phenotype_logic: dict = field(default_factory=dict,
                                  metadata={"description": "Phenotype logic rules defining required loci and genes"})
    antigenic_units: dict = field(default_factory=dict, metadata={"description": "Antigenic unit mappings"})

    @property
    def parsed_version(self) -> tuple[int, ...]:
        """Parses a Semantic Version string into a tuple for safe numeric comparison."""
        import re
        return tuple(int(x) for x in re.findall(r'\d+', str(self.version)))

    @classmethod
    def from_dict(cls, data: dict) -> 'DatabaseMetadata':
        if not isinstance(data, dict):
            raise DatabaseError("Metadata must be a dictionary.")

        try:
            meta = cls(
                name=data['name'],
                keyword=data['keyword'],
                genbank=data['genbank'],
                organism=data['organism'],
                taxon=int(data['taxon']),
                antigen=data['antigen'],
                pathway=data['pathway'],
                version=data['version'],
                id_threshold=float(data['id_threshold']),
                doi=data['doi'],
                owner=data['owner'],
                repo=data['repo'],
                branch=data['branch'],
                contact=data['contact'],
                phenotype_logic=data.get('phenotype_logic', data.get('logic', {})),
                antigenic_units=data.get('antigenic_units', data.get('units', {}))
            )
        except KeyError as e:
            raise DatabaseError(f"Metadata is missing required field: {e.args[0]!r}")
        except ValueError as e:
            raise DatabaseError(f"Metadata has an invalid value type: {e}")

        return meta


@dataclass(slots=True, frozen=True)
class Phenotype:
    id: str
    loci: set[str]
    extra_genes: set[str]
    inactive_genes: set[str]
    priority: int = 50


@dataclass(frozen=True, slots=True)
class Database:
    """
    Optimised, flat representation of a Kaptive antigen database in memory utilising a SoA-layout for math operations.
    """
    metadata: DatabaseMetadata = field(
        metadata={"description": "The strict schema metadata associated with the database."})

    # --- LOCUS-LEVEL DATA ---
    loci: SeqBatch  # Locus nucleotide sequences and IDs mapping to locus indices.
    serotypes: tuple[str, ...]  # Tuple of serotypes mapping to locus indices.
    locus_gene_slices: tuple[slice, ...]  # Slices replacing index arrays for fast locus-to-gene mapping.
    locus_intervals: tuple[IntervalBatch, ...]  # Interval batches mapping to locus indices.

    # --- GENE-LEVEL DATA (Flat SoA) ---
    genes: SeqBatch  # Gene nucleotide sequences mapping to global gene indices.
    translations: SeqBatch  # Gene protein sequences mapping to global gene indices.
    extra_genes: npt.NDArray[np.bool_]  # Boolean mask indicating if a gene belongs to an extra locus (dtype=bool).
    gene_locus_indices: npt.NDArray[np.int32]  # Maps every global gene index back to its locus index (dtype=np.int32).

    # Integer-Encoded Math Arrays
    cluster_keys: tuple[str, ...]  # Tuple mapping integer IDs back to cluster names for O(1) lookups.
    gene_cluster_ids: npt.NDArray[np.int32]  # Integer cluster IDs for each global gene index (dtype=np.int32).
    description_keys: tuple[str, ...]  # Tuple mapping integer IDs back to descriptions for O(1) lookups.
    gene_description_ids: npt.NDArray[np.int32]  # Integer description IDs for each global gene index (dtype=np.int32).
    gene_positions: npt.NDArray[np.int32]  # Expected biological positions for each global gene index (dtype=np.int32).

    phenotypes: tuple[Phenotype, ...]

    @property
    def max_locus_length(self) -> int:
        """Returns the length of the longest locus sequence in the database."""
        return int(np.max(self.loci.lengths)) if len(self.loci) > 0 else 0

    @property
    def cluster_vocab(self) -> dict[str, int]:
        """Vocabulary mapping gene cluster names to integer IDs."""
        return {k: i for i, k in enumerate(self.cluster_keys)}

    @property
    def description_vocab(self) -> dict[str, int]:
        """Vocabulary mapping gene descriptions to integer IDs."""
        return {k: i for i, k in enumerate(self.description_keys)}

    @staticmethod
    def _parse_phenotype(id_: str, data: dict, locus_iterable, cluster_iterable) -> Phenotype:
        loci, inactive, extra = [], [], []
        for token, result, iterable in (('loci', loci, locus_iterable),
                                        ('extra_genes', extra, cluster_iterable),
                                        ('inactive_genes', inactive, cluster_iterable)):
            for t in data.get(token, []):
                if '*' in t:
                    result += fnmatch_filter(iterable, t)
                else:
                    if t in iterable:
                        result.append(t)

        return Phenotype(id_, set(loci), set(extra), set(inactive), data.get('priority', 50))

    @classmethod
    def load(cls, file_or_keyword: str | Path) -> 'Database':
        try:
            file = check_file(file_or_keyword)
            if file.suffix == '.gbk':
                return cls.from_genbank(file)
            elif file.suffix == '.pkl':
                return cls.from_pickle(file)
            raise DatabaseError(f'File {file_or_keyword} not supported')
        except FileNotFoundError:
            try:
                file = DatabaseManager._get_existing_db_path(str(file_or_keyword))
                return pickle.loads(file.read_bytes())
            except DatabaseError:
                return DatabaseManager.install(str(file_or_keyword))

    @classmethod
    def from_pickle(cls, file: str | Path) -> 'Database':
        return pickle.loads(check_file(file).read_bytes())

    @classmethod
    def from_genbank(cls, file: str | Path) -> 'Database':
        """
        Parses a Kaptive database from a legacy Genbank-formatted file and a metadata TOML.
        """
        file = check_file(file)
        from gb_io import iter as GenbankIterator
        _LOCUS_REGEX = re_compile(r'locus:\s?(.*)$')
        _SEROTYPE_REGEX = re_compile(r'type:\s?(.*)$')
        _EXTRA_REGEX = re_compile(r'Extra genes:\s?(.*)$')

        global_gene_idx = 0

        # Locus trackers
        locus_records, serotype_names, locus_gene_slices, locus_intervals = [], [], [], []

        # Gene trackers
        gene_ids, gene_seqs, protein_seqs, extra_genes = [], [], [], []
        gene_cluster_ids, gene_description_ids, gene_expected_positions = [], [], []

        # Global Vocabulary tracker
        cluster_vocab, description_vocab = {}, {}

        with file.open('rb') as fh:
            for rec in GenbankIterator(fh):
                locus_name, serotype, extra = None, None, False

                if not (notes := [i.value for i in rec.features[0].qualifiers if i.key == 'note']):
                    raise DatabaseError(f'Locus has no "note" qualifiers: {rec.name}')

                # Iterate over notes to extract locus and type names, or whether the locus is an "Extra genes" locus
                for note in notes:  # type: str
                    if match := _EXTRA_REGEX.search(note):
                        extra = True
                        locus_name = match.group(1)
                        break

                    if not locus_name and (match := _LOCUS_REGEX.search(note)):
                        locus_name = match.group(1)

                    if not serotype and (match := _SEROTYPE_REGEX.search(note)):
                        serotype = match.group(1)

                if not locus_name:
                    raise DatabaseError(f'Locus has no valid "locus" qualifiers: {rec.name}')

                locus_record = SeqRecord(locus_name, rec.sequence.upper())

                # Local trackers for the current locus
                starts, ends, strands = [], [], []

                local_gene_idx = 0
                locus_start_idx = global_gene_idx

                for feat in rec.features[1:]:
                    if feat.kind != 'CDS':
                        continue

                    cluster, description = '', ''
                    for i in feat.qualifiers:
                        if not cluster and i.key == 'gene':
                            cluster = i.value
                        if not description and i.key == 'product':
                            description = i.value

                    if not extra:
                        gene_id = f'{locus_name}_{local_gene_idx + 1:02}_{cluster}'
                    else:
                        gene_id = cluster

                    if cluster not in cluster_vocab:
                        cluster_vocab[cluster] = len(cluster_vocab)
                    cluster_id = cluster_vocab[cluster]

                    if description not in description_vocab:
                        description_vocab[description] = len(description_vocab)
                    description_id = description_vocab[description]

                    expected_pos = local_gene_idx + 1  # Standardizing biological position

                    # Coordinate parsing
                    loc = feat.location
                    start, end = sorted((loc.start, loc.end))
                    strand_val = -1 if loc.strand in (-1, '-') else 1

                    # Append to flat gene arrays
                    starts.append(start)
                    ends.append(end)
                    strands.append(strand_val)
                    gene_ids.append(gene_id)
                    gene_cluster_ids.append(cluster_id)
                    gene_description_ids.append(description_id)
                    gene_seqs.append(dna_seq := locus_record.extract(start, end, strand_val))
                    protein_seqs.append(BacterialTranslationTable.translate(dna_seq).tobytes())


                    # Extra genes do not have expected biological positions or synteny matrices
                    gene_expected_positions.append(0 if extra else expected_pos)

                    local_gene_idx += 1
                    global_gene_idx += 1

                if local_gene_idx == 0:
                    continue

                # 1. Store the slice instead of an array of indices
                locus_gene_slices.append(slice(locus_start_idx, global_gene_idx))

                # 2. Handle Sequence Extraction via IntervalBatch
                intervals = IntervalBatch(
                    np.array(starts, dtype=np.int32),
                    np.array(ends, dtype=np.int32),
                    np.array(strands, dtype=np.int8)
                )

                # 4. Append final locus metadata
                locus_records.append(locus_record)
                serotype_names.append(serotype or "")
                locus_intervals.append(intervals)
                extra_genes.extend([extra] * local_gene_idx)

            # Pre-compute an array mapping every global gene index back to its locus index
            gene_locus_indices = np.zeros(global_gene_idx, dtype=np.int32)
            for i, s in enumerate(locus_gene_slices):
                gene_locus_indices[s] = i

        db_gene_ids = tuple(gene_ids)
        loci = SeqBatch.from_records(locus_records)
        cluster_keys = tuple(cluster_vocab.keys())
        phenotypes = []
        if (metadata_file := file.with_suffix('.toml')).is_file():
            with metadata_file.open('rb') as fp:
                metadata = DatabaseMetadata.from_dict(tomllib.load(fp))
                for k, v in metadata.phenotype_logic.items():
                    phenotypes.append(cls._parse_phenotype(k, v, loci.ids, cluster_keys))
        else:
            raise DatabaseError("Missing required TOML metadata file alongside Genbank file.")

        return cls(
            metadata=metadata,
            loci=loci,
            serotypes=tuple(serotype_names),
            locus_gene_slices=tuple(locus_gene_slices),
            locus_intervals=tuple(locus_intervals),
            genes=SeqBatch.from_bytes(gene_seqs, ids=db_gene_ids),
            translations=SeqBatch.from_bytes(protein_seqs, ids=db_gene_ids),
            extra_genes=np.array(extra_genes, dtype=bool),
            gene_locus_indices=gene_locus_indices,
            cluster_keys=cluster_keys,
            gene_cluster_ids=np.array(gene_cluster_ids, dtype=np.int32),
            description_keys=tuple(description_vocab.keys()),
            gene_description_ids=np.array(gene_description_ids, dtype=np.int32),
            gene_positions=np.array(gene_expected_positions, dtype=np.int32),
            phenotypes=tuple(phenotypes),
        )


class DatabaseManager:
    _KNOWN = {  # Tuple of username, repo and the name of the database - key is the database keyword from metadata
        "kpsc_k": ("klebgenomics", "KpSC_surface_antigen_loci", "Klebsiella_pneumoniae_Species_Complex_K"),
        "kpsc_o": ("klebgenomics", "KpSC_surface_antigen_loci", "Klebsiella_pneumoniae_Species_Complex_O")
    }
    _DB_DIR = Path.home() / '.kaptive' / 'databases' / 'data'
    _DB_DIR.mkdir(parents=True, exist_ok=True)

    @classmethod
    def _get_db_path(cls, kwd: str) -> Path:
        """Returns the expected file path for a given database keyword."""
        return cls._DB_DIR / f"{kwd}.pkl"

    @classmethod
    def _get_existing_db_path(cls, kwd: str) -> Path:
        """Returns the file path for an installed database, raising an error if it doesn't exist."""
        db_path = cls._get_db_path(kwd)
        if not db_path.is_file():
            raise DatabaseError(f'Database "{kwd}" has not been installed.')
        return db_path

    @classmethod
    def reset(cls):
        """Removes all installed databases by deleting their compiled files."""
        if cls._DB_DIR.exists():
            for file_path in cls._DB_DIR.glob('*.pkl'):
                file_path.unlink()

    @classmethod
    def uninstall(cls, kwd: str):
        """Uninstalls a specific database."""
        cls._get_existing_db_path(kwd).unlink()

    @classmethod
    def installed(cls) -> list[str]:
        """Returns a list of keywords for all currently installed databases."""
        if not cls._DB_DIR.exists():
            return []
        return [p.stem for p in cls._DB_DIR.glob('*.pkl')]

    @classmethod
    def known(cls) -> list[str]:
        """Returns a list of keywords for all currently known databases."""
        return list(cls._KNOWN.keys())

    @classmethod
    def update(cls) -> Generator['Database', None, None]:
        """Updates all installed databases"""
        for db_path in cls._DB_DIR.glob('*.pkl'):
            meta = pickle.loads(db_path.read_bytes()).metadata  # type: DatabaseMetadata
            yield cls.add(meta.owner, meta.repo, meta.name, branch=meta.branch, local_meta=meta)

    @classmethod
    def install(cls, kwd: str) -> 'Database':
        """Installs a known database."""
        if (known_info := cls._KNOWN.get(kwd, None)) is None:
            raise DatabaseError(f'"{kwd}" is not a known database, choose from {list(cls._KNOWN.keys())}')
        return cls.add(*known_info)

    @classmethod
    def add(cls, owner: str, repo_name: str, db_name: str, branch: str = 'main',
            local_meta: DatabaseMetadata | None = None) -> 'Database':
        """Add a database from a known Git repository."""

        repo = GitRepo(owner, repo_name, branch=branch)
        toml_filename = f"{db_name}.toml"
        gbk_filename = f"{db_name}.gbk"

        # Fetch TOML first to quickly check version
        toml_bytes = repo.fetch(toml_filename)
        remote_meta_dict = tomllib.loads(toml_bytes.decode('utf-8'))
        remote_meta = DatabaseMetadata.from_dict(remote_meta_dict)

        db_keyword = remote_meta.keyword
        db_path = cls._get_db_path(db_keyword)

        if local_meta is None and db_path.is_file():
            local_db = pickle.loads(db_path.read_bytes())
            local_meta = getattr(local_db, 'metadata', None)

        if local_meta and local_meta.parsed_version >= remote_meta.parsed_version:
            return  # Up to date

        # Fetch GBK if update is needed
        gbk_bytes = repo.fetch(gbk_filename)

        # Use temp directory for parsing since Database.from_genbank expects files
        with TemporaryDirectory() as tmpdir:
            tmp_path = Path(tmpdir)
            gbk_path = tmp_path / gbk_filename
            toml_path = tmp_path / toml_filename
            gbk_path.write_bytes(gbk_bytes)
            toml_path.write_bytes(toml_bytes)
            db_obj = Database.from_genbank(gbk_path)

        db_path.write_bytes(pickle.dumps(db_obj, protocol=pickle.HIGHEST_PROTOCOL))
        return db_obj

    @classmethod
    def load(cls, kwd: str) -> Database:
        """Load an installed database using its keyword"""
        return pickle.loads(cls._get_existing_db_path(kwd).read_bytes())

    @classmethod
    def save(cls, db: Database):
        """Dump an installed database using its keyword"""
        return cls._get_db_path(db.metadata.keyword).write_bytes(pickle.dumps(db, protocol=pickle.HIGHEST_PROTOCOL))
