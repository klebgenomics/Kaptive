r"""Database management module for downloading, compiling, and managing Kaptive databases.

This module defines the [`DatabaseManager`][kaptive.db.manager.DatabaseManager] class, which acts as the core controller
for managing Kaptive surface antigen database assets both locally and remotely.

Databases are curated as source GenBank (`.gbk`) and metadata (`.toml`) files in remote GitHub repositories.
[`DatabaseManager`][kaptive.db.manager.DatabaseManager] handles remote retrieval over HTTP, version evaluation via
[`DatabaseMetadata`][kaptive.db.models.DatabaseMetadata], on-the-fly parsing and compilation into optimized
[`Database`][kaptive.db.core.Database] instances, local caching in `~/.kaptive`, fast disk serialization (pickling),
and lifecycle commands (`install`, `update`, `uninstall`, `reset`, `add`, `load`, `save`).

Classes:
    DatabaseManager: Class managing local database storage and remote GitHub database retrieval.
"""

import concurrent.futures
import json
import os
import pickle
import tomllib
import urllib.error
from collections.abc import Generator
from dataclasses import asdict
from pathlib import Path
from tempfile import TemporaryDirectory
from urllib.request import urlopen

from kaptive.db.core import Database
from kaptive.db.models import DatabaseError, DatabaseMetadata


class DatabaseManager:
    r"""Class for managing Kaptive databases both on the user's disk and in curator GitHub repositories.

    This class provides a comprehensive mechanism for downloading, compiling, and managing Kaptive databases.
    Databases are maintained as source files (GenBank and TOML) in Git repositories. The
    [`DatabaseManager`][kaptive.db.manager.DatabaseManager] fetches these files, compiles them into optimized, flat
    [`Database`][kaptive.db.core.Database] objects (using a Structure-of-Arrays layout for vectorized operations),
    and stores them locally as serialized pickle files (`.pkl`) alongside `.json` metadata sidecars in the user's
    local directory (defaults to `~/.kaptive` or `$KAPTIVE_DB_DIR`).

    The manager handles:

    - **Installation:** Fetching a known database from its remote repository
      ([`install`][kaptive.db.manager.DatabaseManager.install]), or a custom database from any GitHub repository
      ([`add`][kaptive.db.manager.DatabaseManager.add]), compiling it, and caching the result locally.
    - **Updates:** Checking the local compiled database against the remote repository's version (specified in
      the TOML metadata) and downloading/recompiling if a newer version exists
      ([`update`][kaptive.db.manager.DatabaseManager.update]).
    - **Storage & Retrieval:** Saving ([`save`][kaptive.db.manager.DatabaseManager.save]) and loading
      ([`load`][kaptive.db.manager.DatabaseManager.load]) these compiled `.pkl` files efficiently.
    - **Lifecycle Management:** Uninstalling specific databases
      ([`uninstall`][kaptive.db.manager.DatabaseManager.uninstall]) or completely resetting the local cache
      ([`reset`][kaptive.db.manager.DatabaseManager.reset]).

    Attributes:
        _KNOWN (dict[str, tuple[str, str, str]]): Internal lookup mapping of officially supported database keywords
            to tuples of `(repository_owner, repository_name, database_base_name)`.
        _DB_DIR (pathlib.Path): Local cache directory path where `.pkl` database files and `.json` metadata
            sidecars are stored.
    """

    _KNOWN = {  # Tuple of username, repo and the name of the database - key is the database keyword from metadata
        "kpsc_k": ("klebgenomics", "KpSC_surface_antigen_loci", "Klebsiella_pneumoniae_Species_Complex_K"),
        "kpsc_o": ("klebgenomics", "KpSC_surface_antigen_loci", "Klebsiella_pneumoniae_Species_Complex_O"),
        "kosc_k": ("klebgenomics", "KoSC-surface-antigen-loci", "Klebsiella_oxytoca_Species_Complex_K_locus_database"),
        "kosc_o": ("klebgenomics", "KoSC-surface-antigen-loci", "Klebsiella_oxytoca_Species_Complex_O_locus_database"),
        "ab_k": ("johannajkenyon", "Abaumannii_surface_polysaccharide_loci", "Acinetobacter_baumannii_K"),
        "ab_o": ("johannajkenyon", "Abaumannii_surface_polysaccharide_loci", "Acinetobacter_baumannii_OC"),
        "ecoli_kps": ("rgladstone", "EC-K-typing", "EC-K-typing_group2and3"),
    }
    _DB_DIR = Path(os.environ.get("KAPTIVE_DB_DIR", Path.home() / ".kaptive"))
    _DB_DIR.mkdir(parents=True, exist_ok=True)

    @classmethod
    def _get_db_path(cls, kwd: str) -> Path:
        r"""Return the expected local `.pkl` file path for a database keyword.

        Args:
            kwd (str): The database keyword identifier (e.g., `'kpsc_k'`).

        Returns:
            Path: The expected `pathlib.Path` pointing to `_DB_DIR / "{kwd}.pkl"`.

        See Also:
            [`DatabaseManager`][kaptive.db.manager.DatabaseManager]
        """
        return cls._DB_DIR / f"{kwd}.pkl"

    @classmethod
    def _get_existing_db_path(cls, kwd: str) -> Path:
        r"""Return the file path for an installed database, validating its existence.

        Args:
            kwd (str): The database keyword identifier (e.g., `'kpsc_k'`).

        Returns:
            Path: The `pathlib.Path` to the existing `.pkl` database file.

        Raises:
            DatabaseError: If the specified database `.pkl` file does not exist in the local cache directory.

        See Also:
            [`DatabaseManager`][kaptive.db.manager.DatabaseManager],
            [`DatabaseError`][kaptive.db.models.DatabaseError]
        """
        db_path = cls._get_db_path(kwd)
        if not db_path.is_file():
            raise DatabaseError(f'Database "{kwd}" has not been installed.')
        return db_path

    @classmethod
    def reset(cls) -> None:
        r"""Remove all installed databases by deleting their compiled files from the local directory.

        This clears the user's `~/.kaptive` cache directory of any `.pkl` database files and `.json` metadata
        sidecar files, effectively uninstalling all downloaded and compiled databases.

        Returns:
            None

        See Also:
            [`DatabaseManager`][kaptive.db.manager.DatabaseManager],
            [`uninstall`][kaptive.db.manager.DatabaseManager.uninstall]
        """
        if cls._DB_DIR.exists():
            for file_path in cls._DB_DIR.glob("*.pkl"):
                file_path.unlink()
            for file_path in cls._DB_DIR.glob("*.json"):
                file_path.unlink()

    @classmethod
    def uninstall(cls, kwd: str) -> None:
        r"""Uninstall a specific database by removing its compiled local `.pkl` and `.json` files.

        Args:
            kwd (str): The keyword of the database to uninstall (e.g., `'kpsc_k'`).

        Returns:
            None

        Raises:
            DatabaseError: If the specified database is not currently installed locally.

        See Also:
            [`DatabaseManager`][kaptive.db.manager.DatabaseManager],
            [`DatabaseError`][kaptive.db.models.DatabaseError]
        """
        db_path = cls._get_existing_db_path(kwd)
        db_path.unlink()
        if db_path.with_suffix(".json").exists():
            db_path.with_suffix(".json").unlink()

    @classmethod
    def installed(cls) -> list[str]:
        r"""Return a list of keywords for all currently installed databases.

        Scans the local storage directory for `.pkl` files and extracts their keywords from the file stems.

        Returns:
            list[str]: A list of database keywords corresponding to installed `.pkl` database files.
                Returns an empty list if no databases are installed or if the storage directory does not exist.

        See Also:
            [`DatabaseManager`][kaptive.db.manager.DatabaseManager],
            [`known`][kaptive.db.manager.DatabaseManager.known]
        """
        if not cls._DB_DIR.exists():
            return []
        return [p.stem for p in cls._DB_DIR.glob("*.pkl")]

    @classmethod
    def known(cls) -> list[str]:
        r"""Return a list of keywords for all currently known, officially supported databases.

        These databases can be installed directly by providing their keyword to
        [`install`][kaptive.db.manager.DatabaseManager.install].

        Returns:
            list[str]: A list of known database keywords (e.g., `['kpsc_k', 'kpsc_o', 'kosc_k', ...]`).

        See Also:
            [`DatabaseManager`][kaptive.db.manager.DatabaseManager],
            [`install`][kaptive.db.manager.DatabaseManager.install]
        """
        return list(cls._KNOWN.keys())

    @classmethod
    def update(cls, kwd: str | list[str] = "all") -> Generator[Database, None, None]:
        r"""Update installed databases by checking against their remote GitHub repositories.

        Extracts local metadata to determine the source repository and version, then checks the remote GitHub
        repository for a newer version. If a newer version is available, the source files (`.gbk` and `.toml`)
        are fetched, compiled into a new [`Database`][kaptive.db.core.Database] object, and saved to disk. When
        updating multiple databases or `"all"`, remote fetches are executed concurrently using a thread pool executor.

        Args:
            kwd (str | list[str]): The keyword(s) of the database to update (e.g., `'kpsc_k'`,
                `['kpsc_k', 'ab_k']`, or `'all'`). Defaults to `"all"`, which updates all currently installed databases.

        Yields:
            Database: The newly compiled [`Database`][kaptive.db.core.Database] object for each database that
                required an update. Databases that are already up-to-date yield nothing.

        Raises:
            DatabaseError: If a requested database is not installed locally, or if network/parsing failures occur
                during update.

        See Also:
            [`installed`][kaptive.db.manager.DatabaseManager.installed],
            [`add`][kaptive.db.manager.DatabaseManager.add],
            [`DatabaseManager`][kaptive.db.manager.DatabaseManager],
            [`Database`][kaptive.db.core.Database],
            [`DatabaseMetadata`][kaptive.db.models.DatabaseMetadata],
            [`DatabaseError`][kaptive.db.models.DatabaseError]
        """
        if kwd == "all":
            kwd = cls.installed()
            if not kwd:
                return

        if isinstance(kwd, list):

            def _fetch_update_one(k: str):
                db_path = cls._get_existing_db_path(k)
                json_path = db_path.with_suffix(".json")
                if json_path.is_file():
                    meta = DatabaseMetadata.from_dict(json.loads(json_path.read_text()))
                else:
                    meta = pickle.loads(db_path.read_bytes()).metadata
                db_name = Path(meta.genbank).with_suffix("").name
                return cls._fetch_files(meta.owner, meta.repo, db_name, branch=meta.branch, local_meta=meta)

            with concurrent.futures.ThreadPoolExecutor() as executor:
                fetched_list = list(executor.map(_fetch_update_one, kwd))

            for fetched in fetched_list:
                if fetched is not None:
                    yield cls._compile_and_save(*fetched)
        else:
            db_path = cls._get_existing_db_path(kwd)
            json_path = db_path.with_suffix(".json")
            if json_path.is_file():
                meta = DatabaseMetadata.from_dict(json.loads(json_path.read_text()))
            else:
                meta = pickle.loads(db_path.read_bytes()).metadata
            db_name = Path(meta.genbank).with_suffix("").name
            if (res := cls.add(meta.owner, meta.repo, db_name, branch=meta.branch, local_meta=meta)) is not None:
                yield res

    @classmethod
    def install(cls, kwd: str | list[str]) -> Database | list[Database | None]:
        r"""Install known, officially supported databases by keyword.

        Looks up the repository details (owner, repo, database name) associated with the provided keyword(s)
        in the internal registry (`_KNOWN`) and delegates file retrieval
        and compilation to [`add`][kaptive.db.manager.DatabaseManager.add]. If `'all'` or a list of keywords is
        supplied, fetching is performed concurrently via a thread pool executor.

        Args:
            kwd (str | list[str]): The keyword(s) of the known database(s) to install (e.g., `'kpsc_k'`,
                `['kpsc_k', 'ab_k']`, or `'all'`).

        Returns:
            Database | list[Database | None]: For a single keyword, returns the compiled
                [`Database`][kaptive.db.core.Database] object (or `None` if already up-to-date). For a list of
                keywords or `'all'`, returns a list of compiled [`Database`][kaptive.db.core.Database] objects
                (or `None` for entries that were up-to-date).

        Raises:
            DatabaseError: If any keyword is not recognized in the list of known databases, or if network/parsing
                errors occur.

        See Also:
            [`known`][kaptive.db.manager.DatabaseManager.known],
            [`DatabaseManager`][kaptive.db.manager.DatabaseManager],
            [`add`][kaptive.db.manager.DatabaseManager.add],
            [`Database`][kaptive.db.core.Database],
            [`DatabaseError`][kaptive.db.models.DatabaseError]
        """
        if kwd == "all":
            kwd = list(cls._KNOWN.keys())

        if isinstance(kwd, list):

            def _fetch_one(k: str):
                if (known_info := cls._KNOWN.get(k, None)) is None:
                    raise DatabaseError(f'"{k}" is not a known database, choose from {list(cls._KNOWN.keys())}')
                return cls._fetch_files(*known_info)

            with concurrent.futures.ThreadPoolExecutor() as executor:
                fetched_list = list(executor.map(_fetch_one, kwd))

            results = []
            for fetched in fetched_list:
                if fetched is None:
                    results.append(None)
                else:
                    results.append(cls._compile_and_save(*fetched))
            return results

        if (known_info := cls._KNOWN.get(kwd, None)) is None:
            raise DatabaseError(f'"{kwd}" is not a known database, choose from {list(cls._KNOWN.keys())}')
        return cls.add(*known_info)

    @classmethod
    def _fetch_files(
        cls,
        owner: str,
        repo_name: str,
        db_name: str,
        branch: str = "main",
        local_meta: DatabaseMetadata | None = None,
    ) -> tuple[str, bytes, bytes] | None:
        r"""Fetch database source files (`.toml` and `.gbk`) from a remote GitHub repository over HTTP.

        Fetches the database metadata (`.toml`) first to determine the remote database version. Compares the remote
        version against `local_meta` (or existing local database metadata on disk). If the local version is already
        up-to-date ([`parsed_version`][kaptive.db.models.DatabaseMetadata.parsed_version] `>=` remote version), file
        downloading is skipped and `None` is returned. If the remote version is newer or no local copy exists, the
        GenBank (`.gbk`) sequence file is downloaded and returned along with the metadata.

        Args:
            owner (str): The owner or organization of the GitHub repository (e.g., `'klebgenomics'`).
            repo_name (str): The name of the GitHub repository (e.g., `'KpSC_surface_antigen_loci'`).
            db_name (str): The base name of the database files within the repository
                (e.g., `'Klebsiella_pneumoniae_Species_Complex_K'`).
            branch (str): The Git branch to fetch from. Defaults to `'main'`.
            local_meta (DatabaseMetadata | None): Pre-loaded metadata of the local installation for version comparison.
                If `None`, attempts to load local metadata from existing cache files. Defaults to `None`.

        Returns:
            tuple[str, bytes, bytes] | None: A tuple containing `(db_name, gbk_bytes, toml_bytes)` if a remote
                download is required, or `None` if the local version is already up-to-date.

        Raises:
            DatabaseError: If remote files are not found (HTTP 404), an HTTP error occurs, or network connectivity
                fails.

        See Also:
            [`DatabaseMetadata`][kaptive.db.models.DatabaseMetadata],
            [`parsed_version`][kaptive.db.models.DatabaseMetadata.parsed_version],
            [`DatabaseError`][kaptive.db.models.DatabaseError],
            [`add`][kaptive.db.manager.DatabaseManager.add]
        """
        base_url = f"https://raw.githubusercontent.com/{owner}/{repo_name}/{branch}"
        toml_url = f"{base_url}/{db_name}.toml"
        gbk_url = f"{base_url}/{db_name}.gbk"

        try:
            with urlopen(toml_url) as response:
                toml_bytes = response.read()
        except urllib.error.HTTPError as e:
            if e.code == 404:
                raise DatabaseError(
                    f"Remote file not found: {toml_url}\nEnsure the repository branch, name, and owner are correct."
                ) from e
            raise DatabaseError(f"HTTP Error {e.code} fetching {toml_url}: {e.reason}") from e
        except urllib.error.URLError as e:
            raise DatabaseError(
                f"Network error: Failed to fetch {toml_url}.Ensure you have an active internet connection. ({e.reason})"
            ) from e

        remote_meta_dict = tomllib.loads(toml_bytes.decode("utf-8"))
        remote_meta = DatabaseMetadata.from_dict(remote_meta_dict)

        db_keyword = remote_meta.keyword
        db_path = cls._get_db_path(db_keyword)
        json_path = db_path.with_suffix(".json")

        if local_meta is None and db_path.is_file():
            if json_path.is_file():
                local_meta = DatabaseMetadata.from_dict(json.loads(json_path.read_text()))
            else:
                local_db = pickle.loads(db_path.read_bytes())
                local_meta = getattr(local_db, "metadata", None)

        if local_meta and local_meta.parsed_version >= remote_meta.parsed_version:
            return None

        try:
            with urlopen(gbk_url) as response:
                gbk_bytes = response.read()
        except urllib.error.HTTPError as e:
            if e.code == 404:
                raise DatabaseError(
                    f"Remote file not found: {gbk_url}\nEnsure the repository branch, name, and owner are correct."
                ) from e
            raise DatabaseError(f"HTTP Error {e.code} fetching {gbk_url}: {e.reason}") from e
        except urllib.error.URLError as e:
            raise DatabaseError(
                f"Network error: Failed to fetch {gbk_url}. Ensure you have an active internet connection. ({e.reason})"
            ) from e

        return db_name, gbk_bytes, toml_bytes

    @classmethod
    def _compile_and_save(cls, db_name: str, gbk_bytes: bytes, toml_bytes: bytes) -> Database:
        r"""Compile raw GenBank and TOML file bytes into a Database instance and save to local storage.

        Writes the provided GenBank and TOML byte buffers into temporary files, invokes
        [`from_genbank`][kaptive.db.core.Database.from_genbank] to construct the vectorized Structure-of-Arrays
        [`Database`][kaptive.db.core.Database] object, saves it to the local cache via
        [`save`][kaptive.db.manager.DatabaseManager.save], and returns the compiled instance.

        Args:
            db_name (str): Base name of the database used for temporary file naming.
            gbk_bytes (bytes): Raw byte contents of the GenBank (`.gbk`) file.
            toml_bytes (bytes): Raw byte contents of the TOML (`.toml`) metadata file.

        Returns:
            Database: The compiled, initialized [`Database`][kaptive.db.core.Database] instance.

        Raises:
            DatabaseError: If GenBank or TOML parsing fails during compilation.

        See Also:
            [`from_genbank`][kaptive.db.core.Database.from_genbank],
            [`save`][kaptive.db.manager.DatabaseManager.save],
            [`Database`][kaptive.db.core.Database]
        """
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
    def add(
        cls,
        owner: str,
        repo_name: str,
        db_name: str,
        branch: str = "main",
        local_meta: DatabaseMetadata | None = None,
    ) -> Database | None:
        r"""Add or update a database directly from a specified remote Git repository.

        This is the primary method for adding custom or official databases from GitHub. The procedure:

        1. Constructs raw GitHub URL endpoints for the repository's `.toml` metadata and `.gbk` GenBank files.
        2. Downloads and parses the remote TOML metadata to extract version information.
        3. Compares the remote version against local metadata (if available). If up-to-date, skips remaining steps and
           returns `None`.
        4. Downloads the raw GenBank file content over HTTP.
        5. Writes source files into a temporary directory and compiles them using
           [`from_genbank`][kaptive.db.core.Database.from_genbank].
        6. Serializes and caches the compiled [`Database`][kaptive.db.core.Database] object into the local storage
           directory.

        Args:
            owner (str): Owner or organization of the GitHub repository (e.g., `'klebgenomics'`).
            repo_name (str): Name of the GitHub repository (e.g., `'KpSC_surface_antigen_loci'`).
            db_name (str): Base name of the database files in the repository
                (e.g., `'Klebsiella_pneumoniae_Species_Complex_K'`).
            branch (str): Git branch name to fetch from. Defaults to `'main'`.
            local_meta (DatabaseMetadata | None): Pre-loaded metadata of local database installation.
                Defaults to `None`.

        Returns:
            Database | None: The newly compiled [`Database`][kaptive.db.core.Database] object if installed or updated,
                or `None` if the local version was already up-to-date.

        Raises:
            DatabaseError: If repository files are not found, network issues occur, or file compilation fails.

        See Also:
            [`DatabaseManager`][kaptive.db.manager.DatabaseManager],
            [`Database`][kaptive.db.core.Database],
            [`DatabaseMetadata`][kaptive.db.models.DatabaseMetadata],
            [`DatabaseError`][kaptive.db.models.DatabaseError]
        """
        fetched = cls._fetch_files(owner, repo_name, db_name, branch=branch, local_meta=local_meta)
        if fetched is None:
            return None
        return cls._compile_and_save(*fetched)

    @classmethod
    def load(cls, kwd: str) -> Database:
        r"""Load a locally installed, compiled database using its keyword.

        Reads and unpickles the serialized `.pkl` database file from the local cache directory.

        Args:
            kwd (str): The keyword identifier of the database to load (e.g., `'kpsc_k'`).

        Returns:
            Database: The deserialized [`Database`][kaptive.db.core.Database] instance.

        Raises:
            DatabaseError: If the specified database is not installed locally.

        See Also:
            [`DatabaseManager`][kaptive.db.manager.DatabaseManager],
            [`Database`][kaptive.db.core.Database],
            [`DatabaseError`][kaptive.db.models.DatabaseError]
        """
        return pickle.loads(cls._get_existing_db_path(kwd).read_bytes())

    @classmethod
    def get(cls, file_or_keyword: str | Path) -> Database:
        r"""Load a Database from a file path or resolve and load it by keyword.

        If `file_or_keyword` points to an existing file, it is loaded directly.
        Otherwise, it is treated as a keyword. If the keyword is not installed locally,
        it will be automatically downloaded and installed.

        Args:
            file_or_keyword (str | Path): File path or recognized database keyword.

        Returns:
            Database: The loaded [`Database`][kaptive.db.core.Database] instance.
        """
        from kaptive.db.core import Database

        try:
            file_path = Path(file_or_keyword)
            if file_path.is_file():
                return Database.load(file_path)
        except (TypeError, ValueError, OSError):
            pass

        try:
            return cls.load(str(file_or_keyword))
        except DatabaseError:
            result = cls.install(str(file_or_keyword))
            if isinstance(result, list):
                result = result[0]
            if result is None:
                return cls.load(str(file_or_keyword))
            return result

    @classmethod
    def save(cls, db: Database) -> int:
        r"""Serialize and save a compiled Database object and its metadata to local storage.

        Saves the database as a `.pkl` file named `{keyword}.pkl` in the local cache directory (`_DB_DIR`).
        Also writes a companion `{keyword}.json` file containing the serialized metadata for fast version checking.

        Args:
            db (Database): The compiled [`Database`][kaptive.db.core.Database] object to save.

        Returns:
            int: The total number of bytes written to the `.pkl` database file.

        See Also:
            [`DatabaseManager`][kaptive.db.manager.DatabaseManager],
            [`Database`][kaptive.db.core.Database],
            [`DatabaseMetadata`][kaptive.db.models.DatabaseMetadata]
        """
        db_path = cls._get_db_path(db.metadata.keyword)
        db_path.with_suffix(".json").write_text(json.dumps(asdict(db.metadata)))
        return db_path.write_bytes(pickle.dumps(db, protocol=pickle.HIGHEST_PROTOCOL))
