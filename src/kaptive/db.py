"""
Updated database module with dataclasses, enums and SoA batch layouts for vectorized math and machine-learning.
"""
import os
import pickle
import tomllib
import json
import urllib.error
from collections.abc import Generator
from dataclasses import dataclass, asdict
from fnmatch import filter as fnmatch_filter
from pathlib import Path
from re import compile as re_compile
from tempfile import TemporaryDirectory
from urllib.request import urlopen

import numpy as np
import numpy.typing as npt

from kaptive.core.interval import Intervals
from kaptive.core.seq import SeqRecord, Sequences


# Exceptions and warnings ----------------------------------------------------------------------------------------------
class DatabaseError(Exception):
    pass

# Data Models ----------------------------------------------------------------------------------------------------------
@dataclass(frozen=True, slots=True)
class DatabaseMetadata:
    """Strict schema for Database metadata with dependency-free validation and ergonomic attribute access.
    
    Attributes:
        name (str): The name of the database, e.g. 'Klebsiella pneumoniae Species Complex K'
        keyword (str): The database keyword, e.g. 'kpsc_k'
        genbank (str): The name of the main database file, e.g. 'Klebsiella_pneumoniae_Species_Complex_K.gbk'
        organism (str): The name of the database organism, e.g. 'Klebsiella pneumoniae Species Complex'
        taxon (int): The NCBI Taxonomy ID of the database organism, e.g. 3390273
        antigen (str): The name of the database antigen, e.g. 'Capsular polysaccharide'
        pathway (str): The name of the database antigen synthesis pathway, e.g. 'Wzx/Wzy-dependent'
        version (str): The version of the database, e.g. '3.2.1'
        id_threshold (float): The identity threshold of the database, e.g. 82.5
        doi (list[str]): A list of DOIs associated with the database, e.g. ['TBD']
        owner (str): The owner of the database Github repo, e.g. 'klebgenomics'
        repo (str): The name of the database Github repo, e.g. 'KpSC_surface_antigen_loci'
        branch (str): The branch of the database Github repo, e.g. 'main'
        contact (dict): The details of the database curators, e.g. {'Kelly Wyres': 'kaptive.typing@gmail.com'}
        phenotype_logic (dict): Phenotype logic rules defining required loci and genes
        antigenic_units (dict): Antigenic unit mappings
    """
    name: str
    keyword: str
    genbank: str
    organism: str
    taxon: int
    antigen: str
    pathway: str
    version: str
    id_threshold: float
    doi: list[str]
    owner: str
    repo: str
    branch: str
    contact: dict
    phenotype_logic: dict
    antigenic_units: dict

    @property
    def parsed_version(self) -> tuple[int, ...]:
        """Parses a Semantic Version string into a tuple for safe numeric comparison."""
        pat = re_compile(r'\d+')
        return tuple(int(x) for x in pat.findall(str(self.version)))

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
    """
    Class to represent a Kaptive phenotype which in most cases is either an antigen or serotype.

    Attributes:
        id (str): The id of the phenotype which is used as the serotyping phenotype.
        loci (set[str]): The loci in the database this phenotype applies to.
        extra_genes (set[str]): Cluster names of extra genes that must ALL be present in this phenotype.
        inactive_genes (set[str]): Cluster names of ANY genes that must be inactive in this phenotype.
        priority (int, optional): Used for sorting multiple matching phenotypes. Defaults to 50.
        as_suffix (bool, optional): Append the phenotype as a suffix to all matching phenotypes. Defaults to False.
    """
    id: str
    loci: set[str]
    extra_genes: set[str]
    inactive_genes: set[str]
    priority: int = 50
    as_suffix: bool = False


@dataclass(frozen=True, slots=True)
class Phenotypes:
    """A Structure-of-Arrays (SoA) container for vectorized phenotype evaluation."""
    ids: tuple[str, ...]
    locus_masks: npt.NDArray[np.bool_]
    extra_masks: npt.NDArray[np.bool_]
    inactive_masks: npt.NDArray[np.bool_]
    priorities: npt.NDArray[np.int32]
    as_suffix: npt.NDArray[np.bool_]

    def __len__(self) -> int:
        return len(self.ids)

    def to_dict(self) -> dict:
        return {
            'ids': self.ids,
            'locus_masks': self.locus_masks,
            'extra_masks': self.extra_masks,
            'inactive_masks': self.inactive_masks,
            'priorities': self.priorities,
            'as_suffix': self.as_suffix
        }

    @classmethod
    def from_dict(cls, data: dict) -> 'Phenotypes':
        return cls(
            ids=tuple(data['ids']),
            locus_masks=np.array(data['locus_masks'], dtype=bool),
            extra_masks=np.array(data['extra_masks'], dtype=bool),
            inactive_masks=np.array(data['inactive_masks'], dtype=bool),
            priorities=np.array(data['priorities'], dtype=np.int32),
            as_suffix=np.array(data['as_suffix'], dtype=bool)
        )


@dataclass(frozen=True, slots=True)
class Database:
    """Optimized, flat representation of a Kaptive antigen database in memory utilizing a Structure-of-Arrays (SoA) layout.

    This class eschews traditional object-oriented hierarchies (e.g., Locus -> Gene -> Sequence) in favor of flat,
    parallel arrays (Structure of Arrays). This layout provides significant performance benefits:

    1.  **Memory Locality:** Data of the same type (e.g., all gene cluster IDs) is stored contiguously in memory,
        drastically reducing cache misses during iterative processing.
    2.  **Vectorization:** Numpy arrays allow operations to be performed on thousands of genes simultaneously
        (e.g., masking extra genes, looking up descriptions) without native Python loops, enabling the use of
        optimized C/Fortran libraries under the hood.
    3.  **Fast Lookups:**  Vocabularies (like `cluster_keys`) convert strings into integer IDs (`gene_cluster_ids`).
        Comparing integers is orders of magnitude faster than comparing strings.

    The database manages two main entities: Loci (full locus sequences) and Genes (individual coding sequences within loci).
    Mappings between these entities are maintained via indices and slices, not object references.

    Attributes:
        metadata (DatabaseMetadata): The strict, validated metadata schema associated with the database.
        loci (Sequences): A vectorized batch of all locus nucleotide sequences and their IDs.
            The length of this batch defines the total number of loci (`N_loci`).
        serotypes (tuple[str, ...]): A tuple of length `N_loci` mapping each locus index to its associated serotype.
        locus_gene_offsets (npt.NDArray[np.uint32]): A 1D array of length `N_loci`. Each element is the start index
            for the genes belonging to that locus in the global gene arrays.
        locus_gene_lengths (npt.NDArray[np.uint32]): A 1D array of length `N_loci`. Each element is the number of
            genes belonging to that locus.
        gene_intervals (Intervals): A vectorized batch of all gene physical coordinates relative to their parent loci.
            The length of this batch is `N_genes`.
        genes (Sequences): A vectorized batch of all gene nucleotide sequences across all loci.
            The length of this batch defines the total number of genes in the database (`N_genes`).
        translations (Sequences): A vectorized batch of all gene protein sequences (translated from `genes`).
            Has length `N_genes`.
        extra_genes (npt.NDArray[np.bool_]): A 1D boolean mask of length `N_genes`. True if the gene belongs to an
            "Extra genes" locus (which do not have expected positions or synteny rules).
        gene_locus_indices (npt.NDArray[np.uint16]): A 1D array of length `N_genes`. Contains the integer index of the
            parent locus for every gene, allowing fast reverse-lookups from gene to locus.
        cluster_keys (tuple[str, ...]): The vocabulary of unique gene cluster names (e.g., 'wzi'). The index in
            this tuple corresponds to the integer ID used in `gene_cluster_ids`.
        gene_cluster_ids (npt.NDArray[np.uint16]): A 1D array of length `N_genes`. Each element is the integer ID
            representing the cluster name of that gene.
        description_keys (tuple[str, ...]): The vocabulary of unique gene product descriptions. The index in this
            tuple corresponds to the integer ID used in `gene_description_ids`.
        gene_description_ids (npt.NDArray[np.uint16]): A 1D array of length `N_genes`. Each element is the integer ID
            representing the product description of that gene.
        gene_positions (npt.NDArray[np.uint16]): A 1D array of length `N_genes`. The expected biological position
            (1-indexed) of the gene within its locus. Set to 0 for extra genes.
        phenotypes (Phenotypes): A vectorized SoA batch dictating how serotypes are assigned based on
            locus and gene presence/absence logic.
    """
    metadata: DatabaseMetadata
    loci: Sequences
    serotypes: tuple[str, ...]
    locus_gene_offsets: npt.NDArray[np.uint32]
    locus_gene_lengths: npt.NDArray[np.uint32]
    gene_intervals: Intervals
    genes: Sequences
    translations: Sequences
    extra_genes: npt.NDArray[np.bool_]
    gene_locus_indices: npt.NDArray[np.uint16]
    cluster_keys: tuple[str, ...]
    gene_cluster_ids: npt.NDArray[np.uint16]
    description_keys: tuple[str, ...]
    gene_description_ids: npt.NDArray[np.uint16]
    gene_positions: npt.NDArray[np.uint16]
    phenotypes: Phenotypes

    @property
    def max_locus_length(self) -> int:
        """Returns the length of the longest locus sequence in the database.

        Useful for pre-allocating memory buffers of sufficient size when aligning against the database.

        Returns:
            int: The maximum sequence length among all loci, or 0 if the database contains no loci.
        """
        return int(np.max(self.loci.lengths)) if len(self.loci) > 0 else 0

    @property
    def cluster_vocab(self) -> dict[str, int]:
        """A dictionary mapping gene cluster names (strings) to their corresponding integer IDs.

        This is a convenience property for O(1) string-to-ID lookups. The reverse lookup
        (ID-to-string) should use the `cluster_keys` tuple directly (e.g., `cluster_keys[id]`).

        Returns:
            dict[str, int]: The cluster vocabulary mapping.
        """
        return {k: i for i, k in enumerate(self.cluster_keys)}

    @property
    def description_vocab(self) -> dict[str, int]:
        """A dictionary mapping gene product descriptions (strings) to their corresponding integer IDs.

        This is a convenience property for O(1) string-to-ID lookups. The reverse lookup
        (ID-to-string) should use the `description_keys` tuple directly (e.g., `description_keys[id]`).

        Returns:
            dict[str, int]: The description vocabulary mapping.
        """
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

    @staticmethod
    def _check_file(file: str | Path, min_size: int = 1) -> Path:
        """Validates that a file exists and meets a minimum size requirement."""
        if isinstance(file, str):
            file = Path(file)

        if file.is_file() and file.stat().st_size >= min_size:
            return file
        raise FileNotFoundError(file)

    @classmethod
    def load(cls, file_or_keyword: str | Path) -> 'Database':
        """Loads a Database from a file path (Genbank or Pickle) or a known database keyword.

        This is the primary entry point for instantiating a Database object. It acts as a factory,
        delegating to the appropriate loader based on the input:

        1. If a valid file path is provided ending in '.gbk', it calls `from_genbank`.
        2. If a valid file path is provided ending in '.pkl', it calls `from_pickle`.
        3. If a file is not found, it assumes the input is a keyword for an installed database
           and attempts to load it via `DatabaseManager`. If the database isn't installed, it attempts
           to install it first.

        Args:
            file_or_keyword (str | Path): The path to a .gbk or .pkl file, or a registered database keyword
                (e.g., 'kpsc_k').

        Returns:
            Database: The loaded and initialized Database object.

        Raises:
            DatabaseError: If the file extension is unsupported, or if the keyword is not recognized
                by the DatabaseManager.
        """
        try:
            file = cls._check_file(file_or_keyword)
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
        """Loads a pre-compiled Database object from a serialized pickle file.

        This method is significantly faster than parsing from Genbank, as all the arrays
        and vocabularies are already constructed and optimized.

        Args:
            file (str | Path): The path to the `.pkl` file containing the serialized Database.

        Returns:
            Database: The unpickled Database object.
        """
        return pickle.loads(cls._check_file(file).read_bytes())

    @classmethod
    def from_genbank(cls, file: str | Path) -> 'Database':
        """Compiles a Database object by parsing a legacy GenBank format file and its associated TOML metadata.

        This is an expensive compilation step that converts the nested, string-heavy GenBank structure
        into the flat, integer-based SoA layout defined by this class. 
        
        The method performs the following tasks:

        1.  Iterates through loci in the `.gbk` file, extracting sequence data.
        2.  Parses CDS (Coding Sequence) features to extract gene coordinates, clusters, and descriptions.
        3.  Builds the string-to-integer vocabularies (`cluster_vocab`, `description_vocab`).
        4.  Constructs the flat numpy arrays (e.g., `gene_cluster_ids`, `gene_positions`).
        5.  Creates `Sequences` and `Intervals` objects for vectorized sequence operations.
        6.  Translates all extracted gene nucleotide sequences into protein sequences.
        7.  Loads and validates metadata from a companion `.toml` file (which must have the same name as the `.gbk` file).
        8.  Parses complex phenotype logic from the metadata.

        Args:
            file (str | Path): The path to the `.gbk` database file. A companion `.toml` file
                must exist in the same directory.

        Returns:
            Database: The newly compiled, optimized Database object.

        Raises:
            DatabaseError: If loci lack required qualifiers ('note', 'locus'), or if the associated
                `.toml` metadata file is missing.
        """
        file = cls._check_file(file)
        from gb_io import iter as GenbankIterator
        _LOCUS_REGEX = re_compile(r'locus:\s?(.*)$')
        _SEROTYPE_REGEX = re_compile(r'type:\s?(.*)$')
        _EXTRA_REGEX = re_compile(r'Extra genes:\s?(.*)$')

        global_gene_idx = 0

        # Locus trackers
        locus_records, serotype_names, locus_gene_offsets, locus_gene_lengths, locus_intervals = [], [], [], [], []

        # Gene trackers
        gene_ids, extra_genes = [], []
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

                    # Extra genes do not have expected biological positions or synteny matrices
                    gene_expected_positions.append(0 if extra else expected_pos)

                    local_gene_idx += 1
                    global_gene_idx += 1

                if local_gene_idx == 0:
                    continue

                locus_gene_offsets.append(locus_start_idx)
                locus_gene_lengths.append(local_gene_idx)

                # 2. Handle Sequence Extraction via Intervals
                intervals = Intervals(
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
            gene_locus_indices = np.zeros(global_gene_idx, dtype=np.uint16)
            for i, (o, l) in enumerate(zip(locus_gene_offsets, locus_gene_lengths)):
                gene_locus_indices[o : o + l] = i

        db_gene_ids = tuple(gene_ids)
        loci = Sequences.from_records(locus_records)
        cluster_keys = tuple(cluster_vocab.keys())
        phenotype_objs = []
        if (metadata_file := file.with_suffix('.toml')).is_file():
            with metadata_file.open('rb') as fp:
                metadata = DatabaseMetadata.from_dict(tomllib.load(fp))
                for k, v in metadata.phenotype_logic.items():
                    phenotype_objs.append(cls._parse_phenotype(k, v, loci.ids, cluster_keys))
        else:
            raise DatabaseError("Missing required TOML metadata file alongside Genbank file.")

        # Initialize Phenotype SoA arrays
        n_pheno, n_loci, n_clusters = len(phenotype_objs), len(loci), len(cluster_keys)
        pheno_ids = []
        locus_masks = np.zeros((n_pheno, n_loci), dtype=bool)
        extra_masks = np.zeros((n_pheno, n_clusters), dtype=bool)
        inactive_masks = np.zeros((n_pheno, n_clusters), dtype=bool)
        priorities = np.zeros(n_pheno, dtype=np.int32)
        as_suffix = np.zeros(n_pheno, dtype=bool)
        
        locus_vocab = {name: i for i, name in enumerate(loci.ids)}
        for i, p in enumerate(phenotype_objs):
            pheno_ids.append(p.id)
            for loc in p.loci:
                locus_masks[i, locus_vocab[loc]] = True
            for ext in p.extra_genes:
                extra_masks[i, cluster_vocab[ext]] = True
            for ina in p.inactive_genes:
                inactive_masks[i, cluster_vocab[ina]] = True
            priorities[i] = p.priority
            as_suffix[i] = p.as_suffix
            
        global_intervals = Intervals.concat(locus_intervals)
        genes = loci.extract_intervals(gene_locus_indices, global_intervals, new_ids=db_gene_ids)
        translations = genes.translate()

        return cls(
            metadata=metadata,
            loci=loci,
            serotypes=tuple(serotype_names),
            locus_gene_offsets=np.array(locus_gene_offsets, dtype=np.uint32),
            locus_gene_lengths=np.array(locus_gene_lengths, dtype=np.uint32),
            gene_intervals=global_intervals,
            genes=genes,
            translations=translations,
            extra_genes=np.array(extra_genes, dtype=bool),
            gene_locus_indices=gene_locus_indices,
            cluster_keys=cluster_keys,
            gene_cluster_ids=np.array(gene_cluster_ids, dtype=np.uint16),
            description_keys=tuple(description_vocab.keys()),
            gene_description_ids=np.array(gene_description_ids, dtype=np.uint16),
            gene_positions=np.array(gene_expected_positions, dtype=np.uint16),
            phenotypes=Phenotypes(
                ids=tuple(pheno_ids),
                locus_masks=locus_masks,
                extra_masks=extra_masks,
                inactive_masks=inactive_masks,
                priorities=priorities,
                as_suffix=as_suffix
            ),
        )


# Classes --------------------------------------------------------------------------------------------------------------
class DatabaseManager:
    """Class for managing Kaptive databases both on the user's disk and in curator Github repositories.

    This class provides a comprehensive mechanism for downloading, compiling, and managing Kaptive databases.
    Databases are maintained as source files (GenBank and TOML) in Git repositories. The `DatabaseManager`
    fetches these files, compiles them into optimized, flat `Database` objects (using a Structure of Arrays
    layout for vectorized operations), and stores them locally as serialized pickle files (`.pkl`) in the user's
    `~/.kaptive` directory.

    The manager handles:
        *   **Installation:** Fetching a known database from its remote repository (`install`), or a custom database
            from any Git repository (`add`), compiling it, and caching the result locally.
        *   **Updates:** Checking the local compiled database against the remote repository's version (specified in
            the TOML metadata) and downloading/recompiling if a newer version exists (`update`).
        *   **Storage & Retrieval:** Saving and loading these compiled `.pkl` files efficiently (`save`, `load`).
        *   **Lifecycle Management:** Uninstalling specific databases (`uninstall`) or completely resetting the
            local cache (`reset`).
    """
    _KNOWN = {  # Tuple of username, repo and the name of the database - key is the database keyword from metadata
        "kpsc_k": ("klebgenomics", "KpSC_surface_antigen_loci", "Klebsiella_pneumoniae_Species_Complex_K"),
        "kpsc_o": ("klebgenomics", "KpSC_surface_antigen_loci", "Klebsiella_pneumoniae_Species_Complex_O"),
        "kosc_k": ("klebgenomics", "KoSC-surface-antigen-loci", "Klebsiella_oxytoca_Species_Complex_K_locus_database"),
        "kosc_o": ("klebgenomics", "KoSC-surface-antigen-loci", "Klebsiella_oxytoca_Species_Complex_O_locus_database"),
        "kaer_k": ("tomdstanton", "KAC_cps_reference", "KAC_k_locus_primary_reference"),
        "ab_k": ("tomdstanton", "A_baumannii_Surface_Antigen_Loci", "Acinetobacter_baumannii_K"),
        "ab_o": ("tomdstanton", "A_baumannii_Surface_Antigen_Loci", "Acinetobacter_baumannii_OC"),
        "ecoli_kps": ("tomdstanton", "EC-K-typing", "EC-K-typing_group2and3"),
        "vpar_k": ("tomdstanton", "vibrio_parahaemolyticus_genomoserotyping", "VibrioPara_Kaptivedb_K"),
        "vpar_o": ("tomdstanton", "vibrio_parahaemolyticus_genomoserotyping", "VibrioPara_Kaptivedb_O"),
    }
    _DB_DIR = Path(os.environ.get("KAPTIVE_DB_DIR", Path.home() / '.kaptive'))
    _DB_DIR.mkdir(parents=True, exist_ok=True)

    @classmethod
    def _get_db_path(cls, kwd: str) -> 'Path':
        """Returns the expected file path for a given database keyword."""
        return cls._DB_DIR / f"{kwd}.pkl"

    @classmethod
    def _get_existing_db_path(cls, kwd: str) -> 'Path':
        """Returns the file path for an installed database, raising an error if it doesn't exist."""
        db_path = cls._get_db_path(kwd)
        if not db_path.is_file():
            raise DatabaseError(f'Database "{kwd}" has not been installed.')
        return db_path

    @classmethod
    def reset(cls) -> 'None':
        """Removes all installed databases by deleting their compiled files from the local directory.
        
        This clears the user's `~/.kaptive` directory of any `.pkl` files, effectively uninstalling
        all downloaded and compiled databases.
        """
        if cls._DB_DIR.exists():
            for file_path in cls._DB_DIR.glob('*.pkl'):
                file_path.unlink()
            for file_path in cls._DB_DIR.glob('*.json'):
                file_path.unlink()

    @classmethod
    def uninstall(cls, kwd: str) -> 'None':
        """Uninstalls a specific database by removing its compiled local file.

        Args:
            kwd (str): The keyword of the database to uninstall (e.g., 'kpsc_k').

        Raises:
            DatabaseError: If the specified database is not currently installed.
        """
        db_path = cls._get_existing_db_path(kwd)
        db_path.unlink()
        if db_path.with_suffix('.json').exists():
            db_path.with_suffix('.json').unlink()

    @classmethod
    def installed(cls) -> 'list[str]':
        """Returns a list of keywords for all currently installed databases.

        Returns:
            list[str]: A list of database keywords corresponding to the `.pkl` files present
                in the local cache directory. Returns an empty list if no databases are installed.
        """
        if not cls._DB_DIR.exists():
            return []
        return [p.stem for p in cls._DB_DIR.glob('*.pkl')]

    @classmethod
    def known(cls) -> 'list[str]':
        """Returns a list of keywords for all currently known, officially supported databases.

        These databases can be installed simply by providing their keyword to the `install` method.

        Returns:
            list[str]: A list of known database keywords.
        """
        return list(cls._KNOWN.keys())

    @classmethod
    def update(cls, kwd: str | list[str] = "all") -> 'Generator[Database, None, None]':
        """Updates installed databases by checking against their remote repositories.

        This method extracts metadata to determine the source repository and current version.
        It then calls `add` to check the remote repository for a newer version. If a newer version
        is found, it is downloaded, compiled, and the local file is replaced. If 'all' or a list
        of keywords is provided, updates are fetched concurrently.

        Args:
            kwd (str | list[str]): The keyword(s) of the database to update (e.g., 'kpsc_k' or 'all').
                Defaults to "all", which updates all currently installed databases.

        Yields:
            Database: The updated compiled `Database` object for each database that required an update.
                If a database was already up-to-date, nothing is yielded for that database.
        """
        import json
        
        if kwd == "all":
            kwd = cls.installed()
            if not kwd:
                return
                
        if isinstance(kwd, list):
            import concurrent.futures
            
            def _update_one(k: str) -> 'Database | None':
                db_path = cls._get_existing_db_path(k)
                json_path = db_path.with_suffix('.json')
                if json_path.is_file():
                    meta = DatabaseMetadata.from_dict(json.loads(json_path.read_text()))
                else:
                    meta = pickle.loads(db_path.read_bytes()).metadata
                db_name = Path(meta.genbank).with_suffix('').name
                return cls.add(meta.owner, meta.repo, db_name, branch=meta.branch, local_meta=meta)
                
            with concurrent.futures.ThreadPoolExecutor() as executor:
                future_to_db = {executor.submit(_update_one, k): k for k in kwd}
                for future in concurrent.futures.as_completed(future_to_db):
                    if (res := future.result()) is not None:
                        yield res
        else:
            db_path = cls._get_existing_db_path(kwd)
            json_path = db_path.with_suffix('.json')
            if json_path.is_file():
                meta = DatabaseMetadata.from_dict(json.loads(json_path.read_text()))
            else:
                meta = pickle.loads(db_path.read_bytes()).metadata
            db_name = Path(meta.genbank).with_suffix('').name
            if (res := cls.add(meta.owner, meta.repo, db_name, branch=meta.branch, local_meta=meta)) is not None:
                yield res

    @classmethod
    def install(cls, kwd: str | list[str]) -> 'Database | list[Database | None]':
        """Installs a known, officially supported database by its keyword.

        This looks up the repository details (owner, repo, database name) associated with the provided
        keyword in the internal registry of known databases and delegates the installation to `add`.
        If 'all' or a list of keywords is provided, they are installed concurrently.

        Args:
            kwd (str | list[str]): The keyword(s) of the known database to install (e.g., 'kpsc_k' or 'all').

        Returns:
            Database | list[Database | None]: The newly compiled and installed `Database` object(s).

        Raises:
            DatabaseError: If the provided keyword is not found in the list of known databases.
        """
        if kwd == "all":
            kwd = list(cls._KNOWN.keys())

        if isinstance(kwd, list):
            import concurrent.futures
            
            def _install_one(k: str) -> 'Database | None':
                if (known_info := cls._KNOWN.get(k, None)) is None:
                    raise DatabaseError(f'"{k}" is not a known database, choose from {list(cls._KNOWN.keys())}')
                return cls.add(*known_info)

            with concurrent.futures.ThreadPoolExecutor() as executor:
                return list(executor.map(_install_one, kwd))

        if (known_info := cls._KNOWN.get(kwd, None)) is None:
            raise DatabaseError(f'"{kwd}" is not a known database, choose from {list(cls._KNOWN.keys())}')
        return cls.add(*known_info)

    @classmethod
    def add(cls, owner: str, repo_name: str, db_name: str, branch: str = 'main',
            local_meta: DatabaseMetadata | None = None) -> 'Database | None':
        """Adds or updates a database from a specific Git repository.

        This is the core method for retrieving and compiling databases. The installation mechanism
        follows these steps:
        1.  Initializes a connection to the specified Git repository.
        2.  Fetches the database's TOML metadata file from the repository.
        3.  Compares the remote version (from the fetched TOML) against the local version (either provided
            via `local_meta` or by reading the existing local database file). If the local version is
            equal to or newer than the remote version, the installation/update is skipped.
        4.  If an update or new installation is required, it fetches the potentially large GenBank (`.gbk`) file.
        5.  Saves the fetched `.toml` and `.gbk` files to a temporary directory.
        6.  Compiles the database by parsing the GenBank and TOML files using `Database.from_genbank`.
        7.  Serializes (pickles) the compiled `Database` object and saves it to the local `~/.kaptive` cache directory.

        Args:
            owner (str): The owner or organization of the Git repository (e.g., 'klebgenomics').
            repo_name (str): The name of the Git repository (e.g., 'KpSC_surface_antigen_loci').
            db_name (str): The base name of the database files within the repository (e.g.,
                'Klebsiella_pneumoniae_Species_Complex_K' without extensions).
            branch (str, optional): The Git branch to fetch from. Defaults to 'main'.
            local_meta (DatabaseMetadata, optional): Pre-loaded metadata of the local installation,
                used to skip unnecessary reading of the local file during updates. Defaults to None.

        Returns:
            Database | None: The compiled `Database` object if an installation or update occurred.
                Returns `None` if the local version was already up-to-date with the remote repository.
        """

        base_url = f"https://raw.githubusercontent.com/{owner}/{repo_name}/{branch}"
        toml_url = f"{base_url}/{db_name}.toml"
        gbk_url = f"{base_url}/{db_name}.gbk"

        # Fetch TOML first to quickly check version
        try:
            with urlopen(toml_url) as response:
                toml_bytes = response.read()
        except urllib.error.HTTPError as e:
            if e.code == 404:
                raise DatabaseError(f"Remote file not found: {toml_url}\n"
                                    f"Ensure the repository branch, name, and owner are correct.") from e
            raise DatabaseError(f"HTTP Error {e.code} fetching {toml_url}: {e.reason}") from e
        except urllib.error.URLError as e:
            raise DatabaseError(f"Network error: Failed to fetch {toml_url}."
                                f"Ensure you have an active internet connection. ({e.reason})") from e
            
        remote_meta_dict = tomllib.loads(toml_bytes.decode('utf-8'))
        remote_meta = DatabaseMetadata.from_dict(remote_meta_dict)

        db_keyword = remote_meta.keyword
        db_path = cls._get_db_path(db_keyword)
        json_path = db_path.with_suffix('.json')

        if local_meta is None and db_path.is_file():
            import json
            if json_path.is_file():
                local_meta = DatabaseMetadata.from_dict(json.loads(json_path.read_text()))
            else:
                local_db = pickle.loads(db_path.read_bytes())
                local_meta = getattr(local_db, 'metadata', None)

        if local_meta and local_meta.parsed_version >= remote_meta.parsed_version:
            return  # Up to date

        # Fetch GBK if update is needed
        try:
            with urlopen(gbk_url) as response:
                gbk_bytes = response.read()
        except urllib.error.HTTPError as e:
            if e.code == 404:
                raise DatabaseError(f"Remote file not found: {gbk_url}\n"
                                    f"Ensure the repository branch, name, and owner are correct.") from e
            raise DatabaseError(f"HTTP Error {e.code} fetching {gbk_url}: {e.reason}") from e
        except urllib.error.URLError as e:
            raise DatabaseError(f"Network error: Failed to fetch {gbk_url}. "
                                f"Ensure you have an active internet connection. ({e.reason})") from e

        # Use temp directory for parsing since Database.from_genbank expects files
        with TemporaryDirectory() as tmpdir:
            tmp_path = Path(tmpdir)
            gbk_path = tmp_path / f"{db_name}.gbk"
            toml_path = tmp_path / f"{db_name}.toml"
            gbk_path.write_bytes(gbk_bytes)
            toml_path.write_bytes(toml_bytes)
            db_obj = Database.from_genbank(gbk_path)

        cls.save(db_obj)
        return db_obj

    @classmethod
    def load(cls, kwd: str) -> 'Database':
        """Loads a locally installed, compiled database using its keyword.

        Args:
            kwd (str): The keyword of the database to load.

        Returns:
            Database: The unpickled, compiled `Database` object.

        Raises:
            DatabaseError: If the database is not installed locally.
        """
        return pickle.loads(cls._get_existing_db_path(kwd).read_bytes())

    @classmethod
    def save(cls, db: Database) -> 'int':
        """Serializes and saves a compiled Database object to the local cache directory.

        The database is saved as a `.pkl` file named after its keyword in the `~/.kaptive` directory.
        A companion `.json` file containing the metadata is also saved for fast version checking.

        Args:
            db (Database): The compiled database object to save.

        Returns:
            int: The number of bytes written to the file.
        """
        db_path = cls._get_db_path(db.metadata.keyword)
        db_path.with_suffix('.json').write_text(json.dumps(asdict(db.metadata)))
        return db_path.write_bytes(pickle.dumps(db, protocol=pickle.HIGHEST_PROTOCOL))
