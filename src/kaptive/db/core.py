r"""Core in-memory representation of Kaptive antigen reference databases.

This module defines the [`Database`][kaptive.db.core.Database] class, which maintains
all locus sequences, gene coordinates, protein translations, phenotypic logic rules,
and search indexes in a memory-efficient Structure-of-Arrays (SoA) layout.

Classes:
    Database: Optimized SoA in-memory reference database ([`Database`][kaptive.db.core.Database]).
"""

import pickle
import tomllib
from collections.abc import Iterable
from dataclasses import dataclass
from fnmatch import filter as fnmatch_filter
from pathlib import Path
from re import compile as re_compile
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from kaptive.compare import LocusData

import numpy as np
import numpy.typing as npt

from kaptive.core.interval import Intervals
from kaptive.core.kmers import FracMinHashIndex
from kaptive.core.seq import SeqRecord, Sequences
from kaptive.db.models import DatabaseError, DatabaseMetadata, Phenotype, Phenotypes


@dataclass(frozen=True, slots=True)
class Database:
    r"""Optimized representation of a Kaptive antigen database in memory using Structure-of-Arrays (SoA).

    This class eschews traditional object hierarchies (e.g., Locus -> Gene -> Sequence) in favor of flat,
    parallel arrays (Structure of Arrays). This layout provides significant performance benefits:

    1.  **Memory Locality:** Data of the same type (e.g., all gene cluster IDs) is stored contiguously,
        drastically reducing cache misses during iterative alignment.
    2.  **Vectorization:** NumPy arrays allow operations to be performed on thousands of genes simultaneously
        without Python loop overhead.
    3.  **Fast Lookups:** Vocabularies convert strings into integer IDs (`gene_cluster_ids`). Comparing
        integers is orders of magnitude faster than string comparison.

    The database manages two main entities: Loci (full locus sequences) and Genes (individual coding sequences).
    Mappings between these entities are maintained via indices and slices rather than object references.

    Attributes:
        metadata (DatabaseMetadata): Strict, validated metadata schema associated with the database
            ([`DatabaseMetadata`][kaptive.db.models.DatabaseMetadata]).
        loci (Sequences): Vectorized batch of locus nucleotide sequences and IDs
            ([`Sequences`][kaptive.core.seq.Sequences]). Length is `N_loci`.
        serotypes (tuple[str, ...]): Tuple of length `N_loci` mapping each locus index to its serotype name.
        locus_gene_offsets (npt.NDArray[np.uint32]): 1D array of length `N_loci` containing global gene start indices.
        locus_gene_lengths (npt.NDArray[np.uint32]): 1D array of length `N_loci` containing gene counts per locus.
        gene_intervals (Intervals): Vectorized batch of gene coordinates relative to parent loci
            ([`Intervals`][kaptive.core.interval.Intervals]). Length is `N_genes`.
        genes (Sequences): Vectorized batch of gene nucleotide sequences
            ([`Sequences`][kaptive.core.seq.Sequences]). Length is `N_genes`.
        translations (Sequences): Vectorized batch of translated protein sequences
            ([`Sequences`][kaptive.core.seq.Sequences]). Length is `N_genes`.
        extra_genes (npt.NDArray[np.bool_]): Boolean mask of length `N_genes` indicating extra genes without synteny.
        gene_locus_indices (npt.NDArray[np.uint16]): 1D array of length `N_genes` mapping gene back to locus index.
        cluster_keys (tuple[str, ...]): Vocabulary of unique gene cluster names.
        gene_cluster_ids (npt.NDArray[np.uint16]): 1D array of length `N_genes` storing cluster integer IDs.
        description_keys (tuple[str, ...]): Vocabulary of unique gene product descriptions.
        gene_description_ids (npt.NDArray[np.uint16]): 1D array of length `N_genes` storing product description IDs.
        gene_positions (npt.NDArray[np.uint16]): 1D array of length `N_genes` storing expected gene position
            (1-indexed).
        phenotypes (Phenotypes): Vectorized batch dictating serotype assignment logic
            ([`Phenotypes`][kaptive.db.models.Phenotypes]).
        loci_sketches (FracMinHashIndex): Precomputed FracMinHash sketches for locus containment testing
            ([`FracMinHashIndex`][kaptive.core.kmers.FracMinHashIndex]).

    See Also:
        [`DatabaseManager`][kaptive.db.manager.DatabaseManager],
        [`DatabaseMetadata`][kaptive.db.models.DatabaseMetadata],
        [`Phenotypes`][kaptive.db.models.Phenotypes]
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
    loci_sketches: "FracMinHashIndex"

    def get_locus_data(self, locus_name: str) -> "LocusData":
        r"""Extracts locus data including protein sequences and coordinate backbone for a specific locus.

        Args:
            locus_name (str): Identifier of the locus sequence (e.g., 'K1').

        Returns:
            LocusData: Container object populated with locus proteins and coordinate backbone
                ([`LocusData`][kaptive.compare.LocusData]).

        Raises:
            ValueError: If `locus_name` is not present in `self.loci.ids`.

        See Also:
            [`LocusData`][kaptive.compare.LocusData]
        """
        from kaptive.compare import LocusData

        locus_idx = self.loci.ids.index(locus_name)
        start = self.locus_gene_offsets[locus_idx]
        length = self.locus_gene_lengths[locus_idx]

        return LocusData(
            proteins=self.translations[start : start + length],
            name=locus_name,
            backbone=self.gene_intervals[start : start + length],
            pieces=None,
            gene_ctg_indices=None,
        )

    @property
    def max_locus_length(self) -> int:
        r"""Returns the length of the longest locus sequence in the database.

        Useful for pre-allocating memory buffers of sufficient size when aligning against the database.

        Returns:
            int: The maximum sequence length among all loci, or 0 if the database contains no loci.
        """
        return int(np.max(self.loci.lengths)) if len(self.loci) > 0 else 0

    @property
    def cluster_vocab(self) -> dict[str, int]:
        r"""Dictionary mapping gene cluster names (strings) to integer IDs.

        Convenience property for O(1) string-to-ID lookups.

        Returns:
            dict[str, int]: Mapping from cluster string name to integer ID.

        See Also:
            [`Database`][kaptive.db.core.Database]
        """
        return {k: i for i, k in enumerate(self.cluster_keys)}

    @property
    def description_vocab(self) -> dict[str, int]:
        r"""Dictionary mapping gene product descriptions (strings) to integer IDs.

        Convenience property for O(1) string-to-ID lookups.

        Returns:
            dict[str, int]: Mapping from product description string to integer ID.

        See Also:
            [`Database`][kaptive.db.core.Database]
        """
        return {k: i for i, k in enumerate(self.description_keys)}

    @staticmethod
    def _parse_phenotype(
        id_: str, data: dict, locus_iterable: Iterable[str], cluster_iterable: Iterable[str]
    ) -> Phenotype:
        r"""Parses phenotype definition rule dictionary into a structured Phenotype dataclass.

        Supports wildcards (`*`) in locus and gene names using pattern matching.

        Args:
            id_ (str): Identifier for the phenotype (e.g. serotype name).
            data (dict): Dictionary containing phenotype rule specifications (loci, extra_genes,
                inactive_genes, priority).
            locus_iterable (Iterable[str]): Valid locus names in database.
            cluster_iterable (Iterable[str]): Valid gene cluster names in database.

        Returns:
            Phenotype: Validated Phenotype object ([`Phenotype`][kaptive.db.models.Phenotype]).

        See Also:
            [`Phenotype`][kaptive.db.models.Phenotype]
        """
        loci, inactive, extra = [], [], []
        for token, result, iterable in (
            ("loci", loci, locus_iterable),
            ("extra_genes", extra, cluster_iterable),
            ("inactive_genes", inactive, cluster_iterable),
        ):
            for t in data.get(token, []):
                if "*" in t:
                    result += fnmatch_filter(iterable, t)
                else:
                    if t in iterable:
                        result.append(t)

        return Phenotype(id_, set(loci), set(extra), set(inactive), data.get("priority", 50))

    @staticmethod
    def _check_file(file: str | Path, min_size: int = 1) -> Path:
        r"""Validates that a file exists and meets a minimum size requirement.

        Args:
            file (str | Path): File path to validate.
            min_size (int): Minimum file size in bytes. Defaults to 1.

        Returns:
            Path: Validated Path object.

        Raises:
            FileNotFoundError: If the file does not exist or size is less than `min_size`.
        """
        if isinstance(file, str):
            file = Path(file)

        if file.is_file() and file.stat().st_size >= min_size:
            return file
        raise FileNotFoundError(file)

    @classmethod
    def load(cls, file: str | Path) -> "Database":
        r"""Loads a Database from a file path (GenBank or Pickle).

        Factory entry point that delegates loading based on file extension:
        1. If `.gbk` file path is provided, invokes [`from_genbank`][kaptive.db.core.Database.from_genbank].
        2. If `.pkl` file path is provided, invokes [`from_pickle`][kaptive.db.core.Database.from_pickle].

        Args:
            file (str | Path): Path to `.gbk`/`.pkl` file.

        Returns:
            Database: Loaded and initialized [`Database`][kaptive.db.core.Database] instance.

        Raises:
            DatabaseError: If file extension is unsupported.
            FileNotFoundError: If the file does not exist.

        See Also:
            [`from_genbank`][kaptive.db.core.Database.from_genbank],
            [`from_pickle`][kaptive.db.core.Database.from_pickle]
        """
        file_path = cls._check_file(file)
        if file_path.suffix == ".gbk":
            return cls.from_genbank(file_path)
        elif file_path.suffix == ".pkl":
            return cls.from_pickle(file_path)
        raise DatabaseError(f"File {file} not supported")

    @classmethod
    def from_pickle(cls, file: str | Path) -> "Database":
        r"""Loads a pre-compiled Database object from a serialized pickle file.

        Args:
            file (str | Path): Path to `.pkl` database file.

        Returns:
            Database: Deserialized [`Database`][kaptive.db.core.Database] instance.

        Raises:
            FileNotFoundError: If the file does not exist.

        See Also:
            [`load`][kaptive.db.core.Database.load]
        """
        return pickle.loads(cls._check_file(file).read_bytes())

    @classmethod
    def from_genbank(cls, file: str | Path) -> "Database":
        r"""Compiles a Database object by parsing a legacy GenBank format file and its associated TOML metadata.

        This is an expensive compilation step that converts the nested, string-heavy GenBank structure
        into the flat, integer-based SoA layout defined by this class.

        The method performs the following tasks:

        1.  Iterates through loci in the `.gbk` file, extracting sequence data.
        2.  Parses CDS (Coding Sequence) features to extract gene coordinates, clusters, and descriptions.
        3.  Builds the string-to-integer vocabularies (`cluster_vocab`, `description_vocab`).
        4.  Constructs the flat numpy arrays (e.g., `gene_cluster_ids`, `gene_positions`).
        5.  Creates `Sequences` and `Intervals` objects for vectorized sequence operations.
        6.  Translates all extracted gene nucleotide sequences into protein sequences.
        7.  Loads and validates metadata from a companion `.toml` file (which must have the same name as
            the `.gbk` file).
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

        _LOCUS_REGEX = re_compile(r"locus:\s?(.*)$")
        _SEROTYPE_REGEX = re_compile(r"type:\s?(.*)$")
        _EXTRA_REGEX = re_compile(r"Extra genes:\s?(.*)$")

        global_gene_idx = 0

        # Locus trackers
        locus_records, serotype_names, locus_gene_offsets, locus_gene_lengths, locus_intervals = (
            [],
            [],
            [],
            [],
            [],
        )

        # Gene trackers
        gene_ids, extra_genes = [], []
        gene_cluster_ids, gene_description_ids, gene_expected_positions = [], [], []

        # Global Vocabulary tracker
        cluster_vocab, description_vocab = {}, {}

        with file.open("rb") as fh:
            for rec in GenbankIterator(fh):
                locus_name, serotype, extra = None, None, False

                if not (notes := [i.value for i in rec.features[0].qualifiers if i.key == "note"]):
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
                    if feat.kind != "CDS":
                        continue

                    cluster, description = "", ""
                    for i in feat.qualifiers:
                        if not cluster and i.key == "gene":
                            cluster = i.value
                        if not description and i.key == "product":
                            description = i.value

                    if not extra:
                        gene_id = f"{locus_name}_{local_gene_idx + 1:02}_{cluster}"
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
                    strand_val = -1 if loc.strand in (-1, "-") else 1

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
                    np.array(strands, dtype=np.int8),
                )

                # 4. Append final locus metadata
                locus_records.append(locus_record)
                serotype_names.append(serotype or "")
                locus_intervals.append(intervals)
                extra_genes.extend([extra] * local_gene_idx)

            # Pre-compute an array mapping every global gene index back to its locus index
            gene_locus_indices = np.zeros(global_gene_idx, dtype=np.uint16)
            for i, (o, length) in enumerate(zip(locus_gene_offsets, locus_gene_lengths)):
                gene_locus_indices[o : o + length] = i

        db_gene_ids = tuple(gene_ids)
        loci = Sequences.from_records(locus_records)
        cluster_keys = tuple(cluster_vocab.keys())
        phenotype_objs = []
        if (metadata_file := file.with_suffix(".toml")).is_file():
            with metadata_file.open("rb") as fp:
                metadata = DatabaseMetadata.from_dict(tomllib.load(fp))
                for k, v in metadata.phenotype_logic.items():
                    phenotype_objs.append(cls._parse_phenotype(k, v, loci.ids, cluster_keys))
        else:
            raise DatabaseError("Missing required TOML metadata file alongside Genbank file.")

        # Initialize Phenotype SoA arrays
        n_pheno, n_loci, n_clusters = len(phenotype_objs), len(loci), len(cluster_keys)
        pheno_ids = []
        locus_vocab = {name: i for i, name in enumerate(loci.ids)}
        locus_masks = np.zeros((n_pheno, n_loci), dtype=bool)
        extra_masks = np.zeros((n_pheno, n_clusters), dtype=np.int8)
        inactive_masks = np.zeros((n_pheno, n_clusters), dtype=np.int8)
        priorities = np.zeros(n_pheno, dtype=np.int8)
        as_suffix = np.zeros(n_pheno, dtype=bool)

        for i, p in enumerate(phenotype_objs):
            pheno_ids.append(p.id)
            for loc in p.loci:
                locus_masks[i, locus_vocab[loc]] = True
            for ext in p.extra_genes:
                extra_masks[i, cluster_vocab[ext]] = 1
            for ina in p.inactive_genes:
                inactive_masks[i, cluster_vocab[ina]] = 1
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
                ids=np.array(pheno_ids, dtype="S32"),
                locus_masks=locus_masks,
                extra_masks=extra_masks,
                inactive_masks=inactive_masks,
                extra_counts=extra_masks.sum(axis=1, dtype=np.int8),
                priorities=priorities,
                as_suffix=as_suffix,
            ),
            loci_sketches=FracMinHashIndex.build(loci, sort_by_hash=False),
        )
