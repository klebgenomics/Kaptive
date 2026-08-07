r"""Data models and custom exception classes for Kaptive database representation.

This module provides data models, structured metadata schemas, and exception classes used by
[`Database`][kaptive.db.core.Database] and [`DatabaseManager`][kaptive.db.manager.DatabaseManager]
to encapsulate surface antigen reference database specifications, phenotype definitions,
and vectorized phenotype matching containers.

Classes:
    DatabaseError: Custom exception for database loading, validation, and format errors.
    DatabaseMetadata: Strict schema for reference database metadata with attributes and validation.
    Phenotype: Single phenotype rule mapping loci and gene requirements to a serotype identifier.
    Phenotypes: Structure-of-Arrays (SoA) batch container for vectorized phenotype evaluation.
"""

from collections.abc import Iterable
from dataclasses import dataclass
from re import compile as re_compile
from typing import Any, Self

import numpy as np
import numpy.typing as npt

from kaptive.core.collections import BatchedContainer


# Exceptions and warnings ----------------------------------------------------------------------------------------------
class DatabaseError(Exception):
    r"""Exception raised for database loading, metadata validation, or format errors.

    This exception is raised when database metadata is invalid, required files are missing,
    or reference database files fail validation in [`DatabaseMetadata`][kaptive.db.models.DatabaseMetadata],
    [`Database`][kaptive.db.core.Database], or [`DatabaseManager`][kaptive.db.manager.DatabaseManager].
    """

    pass


# Data Models ----------------------------------------------------------------------------------------------------------
@dataclass(frozen=True, slots=True)
class DatabaseMetadata:
    r"""Strict schema for Database metadata with dependency-free validation and ergonomic attribute access.

    Represents the metadata associated with a Kaptive reference database, including organism details,
    locus pathway classifications, repository location, curator contact details, and phenotype logic rules.
    Used by [`Database`][kaptive.db.core.Database] and [`DatabaseManager`][kaptive.db.manager.DatabaseManager].

    Attributes:
        name (str): The name of the database, e.g. 'Klebsiella pneumoniae Species Complex K'.
        keyword (str): The database keyword, e.g. 'kpsc_k'.
        genbank (str): The name of the main database file, e.g. 'Klebsiella_pneumoniae_Species_Complex_K.gbk'.
        organism (str): The name of the database organism, e.g. 'Klebsiella pneumoniae Species Complex'.
        taxon (int): The NCBI Taxonomy ID of the database organism, e.g. 3390273.
        antigen (str): The name of the database antigen, e.g. 'Capsular polysaccharide'.
        pathway (str): The name of the database antigen synthesis pathway, e.g. 'Wzx/Wzy-dependent'.
        version (str): The version of the database, e.g. '3.2.1'.
        id_threshold (float): The identity threshold of the database, e.g. 82.5.
        doi (list[str]): A list of DOIs associated with the database, e.g. ['TBD'].
        owner (str): The owner of the database Github repo, e.g. 'klebgenomics'.
        repo (str): The name of the database Github repo, e.g. 'KpSC_surface_antigen_loci'.
        branch (str): The branch of the database Github repo, e.g. 'main'.
        contact (dict): The details of the database curators, e.g. {'Kelly Wyres': 'kaptive.typing@gmail.com'}.
        phenotype_logic (dict): Phenotype logic rules defining required loci and genes.
        antigenic_units (dict): Antigenic unit mappings.
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
    contact: dict  # type: ignore
    phenotype_logic: dict  # type: ignore
    antigenic_units: dict  # type: ignore

    @property
    def parsed_version(self) -> tuple[int, ...]:
        r"""Parses the semantic version string into a tuple of integers for numeric comparison.

        Extracts numeric digit sequences from the [`DatabaseMetadata`][kaptive.db.models.DatabaseMetadata]
        attribute and converts them into a tuple of integers (e.g., `'3.2.1'` becomes `(3, 2, 1)`).

        Returns:
            tuple[int, ...]: A tuple of extracted integer components representing the database version.
        """
        pat = re_compile(r"\d+")
        return tuple(int(x) for x in pat.findall(str(self.version)))

    @classmethod
    def from_dict(cls, data: dict) -> "DatabaseMetadata":  # type: ignore
        r"""Instantiates a `DatabaseMetadata` object from a dictionary.

        Validates required fields, casts numeric types, and sets default fallback dictionaries for
        phenotype logic and antigenic units.

        Args:
            data (dict): Dictionary containing metadata fields (e.g., parsed from JSON or TOML).

        Returns:
            DatabaseMetadata: Validated [`DatabaseMetadata`][kaptive.db.models.DatabaseMetadata] instance.

        Raises:
            DatabaseError: If `data` is not a dict, missing required keys, or contains invalid attribute types.
        """
        if not isinstance(data, dict):
            raise DatabaseError("Metadata must be a dictionary.")

        try:
            meta = cls(
                name=data["name"],
                keyword=data["keyword"],
                genbank=data["genbank"],
                organism=data["organism"],
                taxon=int(data["taxon"]),
                antigen=data["antigen"],
                pathway=data["pathway"],
                version=data["version"],
                id_threshold=float(data["id_threshold"]),
                doi=data["doi"],
                owner=data["owner"],
                repo=data["repo"],
                branch=data["branch"],
                contact=data["contact"],
                phenotype_logic=data.get("phenotype_logic", data.get("logic", {})),
                antigenic_units=data.get("antigenic_units", data.get("units", {})),
            )
        except KeyError as e:
            raise DatabaseError(f"Metadata is missing required field: {e.args[0]!r}")
        except ValueError as e:
            raise DatabaseError(f"Metadata has an invalid value type: {e}")

        return meta


@dataclass(slots=True, frozen=True)
class Phenotype:
    r"""Single locus phenotype rule mapping loci and gene requirements to a serotype identifier.

    Defines criteria for assigning a specific phenotype (e.g., K-type or O-type serotype)
    based on identified reference loci, required extra genes, and forbidden inactive genes.
    Processed by [`Database`][kaptive.db.core.Database] into vectorized [`Phenotypes`][kaptive.db.models.Phenotypes]
    batches.

    Attributes:
        id (str): Unique phenotype or serotype identifier string.
        loci (set[str]): Locus names in the database to which this phenotype applies.
        extra_genes (set[str]): Set of gene cluster names that must all be present for this phenotype match.
        inactive_genes (set[str]): Set of gene cluster names that must not be inactivated for this phenotype match.
        priority (int): Sorting priority when resolving multiple matching phenotypes. Defaults to `50`.
        as_suffix (bool): Whether to append this phenotype identifier as a suffix to matching phenotypes.
            Defaults to `False`.
    """

    id: str
    loci: set[str]
    extra_genes: set[str]
    inactive_genes: set[str]
    priority: int = 50
    as_suffix: bool = False


@dataclass(frozen=True, slots=True)
class Phenotypes(BatchedContainer[Any, "Phenotypes"]):
    r"""Structure-of-Arrays (SoA) container for vectorized phenotype evaluation.

    Encapsulates boolean matrix masks and priority arrays across a batch of [`Phenotype`][kaptive.db.models.Phenotype]
    definitions for high-performance vectorized evaluation during serotyping. Inherits from
    [`BatchedContainer`][kaptive.core.collections.BatchedContainer].

    Attributes:
        ids (npt.NDArray[np.bytes_]): 1D byte string array (e.g. `S32`) of phenotype identifier strings.
        locus_masks (npt.NDArray[np.bool_]): 2D boolean array of shape `(N, num_loci)` indicating locus requirements.
        extra_masks (npt.NDArray[np.int8]): 2D integer array of shape `(N, num_extra_genes)` for required extra genes.
        inactive_masks (npt.NDArray[np.int8]): 2D integer array of shape `(N, num_inactive_genes)` for forbidden
            inactive genes.
        extra_counts (npt.NDArray[np.int8]): 1D integer array storing the sum of extra required genes per phenotype.
        priorities (npt.NDArray[np.int8]): 1D integer array of shape `(N,)` indicating resolution priority values.
        as_suffix (npt.NDArray[np.bool_]): 1D boolean array of shape `(N,)` indicating if phenotype is used as a suffix.
    """

    ids: npt.NDArray[np.bytes_]
    locus_masks: npt.NDArray[np.bool_]
    extra_masks: npt.NDArray[np.int8]
    inactive_masks: npt.NDArray[np.int8]
    extra_counts: npt.NDArray[np.int8]
    priorities: npt.NDArray[np.int8]
    as_suffix: npt.NDArray[np.bool_]

    def __len__(self) -> int:
        r"""Returns the number of phenotype records in the container batch.

        Returns:
            int: Number of phenotype records.
        """
        return len(self.ids)

    def __getitem__(self, item: int | slice | npt.NDArray[Any] | list[int]) -> "Any | Phenotypes":
        r"""Slices or masks the `Phenotypes` container batch along the primary dimension.

        Args:
            item (int | slice | npt.NDArray | list): Slice object, boolean mask array, or list of indices to select.

        Returns:
            Any | Phenotypes: A new [`Phenotypes`][kaptive.db.models.Phenotypes] container containing
                the selected subset of records.

        Raises:
            NotImplementedError: If a single integer index is provided, as scalar indexing is not supported
                on SoA containers.
        """
        if isinstance(item, (int, np.integer)):
            raise NotImplementedError("Single item access not implemented for Phenotypes")
        return Phenotypes(
            ids=self.ids[item],
            locus_masks=self.locus_masks[item],
            extra_masks=self.extra_masks[item],
            inactive_masks=self.inactive_masks[item],
            extra_counts=self.extra_counts[item],
            priorities=self.priorities[item],
            as_suffix=self.as_suffix[item],
        )

    @classmethod
    def empty(cls) -> "Phenotypes":
        r"""Constructs an empty `Phenotypes` batch container.

        Returns:
            Phenotypes: An empty [`Phenotypes`][kaptive.db.models.Phenotypes] instance with 0 elements
                and 2D empty arrays.
        """
        return cls(
            ids=np.empty(0, dtype="S32"),
            locus_masks=np.empty((0, 0), dtype=bool),
            extra_masks=np.empty((0, 0), dtype=np.int8),
            inactive_masks=np.empty((0, 0), dtype=np.int8),
            extra_counts=np.empty(0, dtype=np.int8),
            priorities=np.empty(0, dtype=np.int8),
            as_suffix=np.empty(0, dtype=bool),
        )

    @classmethod
    def concat(cls, batches: Iterable[Self]) -> Self:  # type: ignore
        r"""Concatenates multiple `Phenotypes` batch containers into a single `Phenotypes` instance.

        Args:
            batches (Iterable[Phenotypes]): An iterable of [`Phenotypes`][kaptive.db.models.Phenotypes]
                container batches to concatenate.

        Returns:
            Phenotypes: A single combined [`Phenotypes`][kaptive.db.models.Phenotypes] batch container.
                Returns an empty container if `batches` is empty.
        """
        batches = list(batches)
        if not batches:
            return cls.empty()  # type: ignore
        return cls(
            ids=np.concatenate([b.ids for b in batches]),
            locus_masks=np.concatenate([b.locus_masks for b in batches]),
            extra_masks=np.concatenate([b.extra_masks for b in batches]),
            inactive_masks=np.concatenate([b.inactive_masks for b in batches]),
            extra_counts=np.concatenate([b.extra_counts for b in batches]),
            priorities=np.concatenate([b.priorities for b in batches]),
            as_suffix=np.concatenate([b.as_suffix for b in batches]),
        )

    def to_dict(self) -> dict:  # type: ignore
        r"""Converts the `Phenotypes` container attributes into a dictionary representation.

        Returns:
            dict: Dictionary mapping attribute names (`ids`, `locus_masks`, `extra_masks`,
                `inactive_masks`, `priorities`, `as_suffix`) to their stored values/arrays.
        """
        return {
            "ids": np.char.decode(self.ids, "utf-8").tolist(),
            "locus_masks": self.locus_masks,
            "extra_masks": self.extra_masks,
            "inactive_masks": self.inactive_masks,
            "extra_counts": self.extra_counts,
            "priorities": self.priorities,
            "as_suffix": self.as_suffix,
        }

    @classmethod
    def from_dict(cls, data: dict) -> "Phenotypes":  # type: ignore
        r"""Reconstructs a `Phenotypes` batch container from a dictionary of array data.

        Args:
            data (dict): Dictionary containing array and tuple entries corresponding to `Phenotypes` attributes.

        Returns:
            Phenotypes: Reconstructed [`Phenotypes`][kaptive.db.models.Phenotypes] container instance.
        """
        return cls(  # type: ignore
            ids=np.array([p.encode("utf-8") for p in data["ids"]], dtype="S32"),
            locus_masks=np.array(data["locus_masks"], dtype=bool),
            extra_masks=np.array(data["extra_masks"], dtype=bool),
            inactive_masks=np.array(data["inactive_masks"], dtype=bool),
            priorities=np.array(data["priorities"], dtype=np.int8),
            as_suffix=np.array(data["as_suffix"], dtype=bool),
        )
