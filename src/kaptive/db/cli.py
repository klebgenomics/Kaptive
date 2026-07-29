r"""Command-line interface commands for managing and querying Kaptive reference databases.

This module provides CLI command implementations for listing installed databases, downloading
or adding new databases, updating installed databases, extracting FASTA sequences (loci,
genes, proteins), printing database metadata, and resetting local database caches.

Classes:
    List: List installed local databases ([`List`][kaptive.db.cli.List]).
    Database: Parent command managing database subcommands ([`Database`][kaptive.db.cli.Database]).
    Install: Install remote database by keyword ([`Install`][kaptive.db.cli.Install]).
    Update: Update installed databases ([`Update`][kaptive.db.cli.Update]).
    Reset: Uninstall local databases and clear cache ([`Reset`][kaptive.db.cli.Reset]).
    Add: Add database from GitHub repository ([`Add`][kaptive.db.cli.Add]).
    Metadata: Display database metadata ([`Metadata`][kaptive.db.cli.Metadata]).
    Extract: Parent command for sequence extraction ([`Extract`][kaptive.db.cli.Extract]).
    Loci: Extract locus nucleotide sequences ([`Loci`][kaptive.db.cli.Loci]).
    Genes: Extract gene nucleotide sequences ([`Genes`][kaptive.db.cli.Genes]).
    Proteins: Extract translated protein sequences ([`Proteins`][kaptive.db.cli.Proteins]).
"""

import argparse
import sys

from kaptive.cli import Colors, Command


class List(Command):
    r"""📋 Command for listing all currently installed local databases.

    Displays the keywords of all compiled `.pkl` databases found in the user's
    local Kaptive directory (`~/.kaptive`).

    Aliases:
        ls

    See Also:
        [`DatabaseManager.installed`][kaptive.db.manager.DatabaseManager.installed]
    """

    aliases = ["ls"]

    def __call__(self, args: argparse.Namespace) -> None:
        r"""Executes the `list` command to print installed database keywords.

        Args:
            args (argparse.Namespace): Parsed command-line arguments.

        See Also:
            [`DatabaseManager.installed`][kaptive.db.manager.DatabaseManager.installed]
        """
        from kaptive.db import DatabaseManager

        if installed := DatabaseManager.installed():
            print("\n".join(installed))
        else:
            self.cli.msg("❌ No databases installed")


# Database command -----------------------------------------------------------------------------------------------------
class Database(Command):
    r"""📦 Parent command for managing local and remote Kaptive databases.

    Aggregates subcommands for listing, installing, updating, resetting, adding,
    extracting, and displaying metadata for Kaptive databases.

    Aliases:
        db
    """

    aliases = ["db"]

    def register_subcommands(self) -> None:
        r"""Registers subcommands for database management."""
        self.subcommands = [
            List(),
            Add(),
            Install(),
            Update(),
            Reset(),
            Extract(),
            Metadata(),
        ]


class Install(Command):
    r"""📦 Command for installing known reference databases via keyword.

    Downloads GenBank and TOML definition files from official repositories, compiles
    them into vectorized [`Database`][kaptive.db.core.Database] objects, and caches
    them locally.

    See Also:
        [`DatabaseManager.install`][kaptive.db.manager.DatabaseManager.install]
    """

    def setup_arguments(self) -> None:
        r"""Configures argument parser options for the `install` subcommand."""
        opts = self.parser.add_argument_group("📥 Inputs")
        opts.add_argument("database", help="Database keyword (see: `kaptive db list`) or 'all'")

    def __call__(self, args: argparse.Namespace) -> None:
        r"""Executes database installation for a specific keyword or 'all'.

        Args:
            args (argparse.Namespace): Parsed CLI arguments containing `database`.

        Raises:
            DatabaseError: If the specified keyword is not a known database.

        See Also:
            [`DatabaseManager.install`][kaptive.db.manager.DatabaseManager.install],
            [`DatabaseError`][kaptive.db.models.DatabaseError]
        """
        if args.database == "all":
            self.cli.msg("📥 Installing all known databases concurrently...")
        else:
            self.cli.msg(f"📥 Installing database '{args.database}'...")

        from kaptive.db import DatabaseManager

        DatabaseManager.install(args.database)

        if args.database == "all":
            self.cli.msg("✅ Successfully installed all known databases.")
        else:
            self.cli.msg(f"✅ Successfully installed '{args.database}'.")


class Update(Command):
    r"""🔄 Command for updating installed local databases from remote repositories.

    Checks installed databases against their source GitHub repositories for newer
    versions defined in TOML metadata and re-compiles modified databases.

    See Also:
        [`DatabaseManager.update`][kaptive.db.manager.DatabaseManager.update]
    """

    def setup_arguments(self) -> None:
        r"""Configures argument parser options for the `update` subcommand."""
        opts = self.parser.add_argument_group("📥 Inputs")
        opts.add_argument(
            "database",
            nargs="?",
            default="all",
            help="Database keyword (see: `kaptive db list`) or 'all' (default: all)",
        )

    def __call__(self, args: argparse.Namespace) -> None:
        r"""Executes database updates for a specific keyword or 'all' installed databases.

        Args:
            args (argparse.Namespace): Parsed CLI arguments containing `database`.

        See Also:
            [`DatabaseManager.update`][kaptive.db.manager.DatabaseManager.update]
        """
        if args.database == "all":
            self.cli.msg("🔄 Checking all installed databases for updates concurrently...")
        else:
            self.cli.msg(f"🔄 Checking '{args.database}' for updates...")

        from kaptive.db import DatabaseManager

        updated = False
        for db in DatabaseManager.update(args.database):
            self.cli.msg(f"✅ Updated {db.metadata.name} to version {db.metadata.version}")
            updated = True

        if not updated:
            self.cli.msg("🎉 All databases are already up to date.")


class Reset(Command):
    r"""🧹 Command for uninstalling all local databases and resetting local cache.

    Deletes all compiled `.pkl` and metadata `.json` files from the local Kaptive cache directory.

    See Also:
        [`DatabaseManager.reset`][kaptive.db.manager.DatabaseManager.reset]
    """

    def __call__(self, args: argparse.Namespace) -> None:
        r"""Executes local cache reset and database removal.

        Args:
            args (argparse.Namespace): Parsed command-line arguments.

        See Also:
            [`DatabaseManager.reset`][kaptive.db.manager.DatabaseManager.reset]
        """
        self.cli.msg("🧹 Uninstalling all local databases...")
        from kaptive.db import DatabaseManager

        DatabaseManager.reset()
        self.cli.msg("✅ All local databases have been uninstalled and reset.")


class Add(Command):
    r"""🔗 Command for adding a custom reference database from a GitHub repository.

    Fetches GenBank and TOML metadata files from any specified GitHub owner/repo/branch,
    compiles the database, and registers it in the local cache.

    See Also:
        [`DatabaseManager.add`][kaptive.db.manager.DatabaseManager.add]
    """

    def setup_arguments(self) -> None:
        r"""Configures argument parser options for the `add` subcommand."""
        opts = self.parser.add_argument_group("📥 Inputs")
        opts.add_argument("database", help="Name for the new database")

        opts = self.parser.add_argument_group(Colors.wrap("🌐 GitHub Details", Colors.BOLD))
        opts.add_argument("owner", help="GitHub repository owner")
        opts.add_argument("repo_name", help="GitHub repository name")
        opts.add_argument(
            "-b",
            "--branch",
            help="GitHub repository branch (default: main)",
            default="main",
            nargs="?",
        )

    def __call__(self, args: argparse.Namespace) -> None:
        r"""Executes retrieval, compilation, and registration of a custom GitHub database.

        Args:
            args (argparse.Namespace): Parsed CLI arguments containing `database`, `owner`,
                `repo_name`, and optional `branch`.

        See Also:
            [`DatabaseManager.add`][kaptive.db.manager.DatabaseManager.add]
        """
        from kaptive.db import DatabaseManager

        self.cli.msg(f"⤵️ Adding {args.database} from {args.owner}/{args.repo_name}/{args.branch}")
        if db := DatabaseManager.add(args.owner, args.repo_name, args.database, args.branch):
            self.cli.msg(f"✅ Added {db.metadata.name} v{db.metadata.version} successfully!")
        else:
            self.cli.msg("❌ Failed to add database! Is it already installed?")


class Metadata(Command):
    r"""📊 Command for printing detailed metadata of a Kaptive database.

    Displays summary information including organism, taxon ID, antigen type, synthesis pathway,
    version, identity threshold, GenBank filename, DOIs, repository URL, and curator contacts.

    Aliases:
        info

    See Also:
        [`Database.load`][kaptive.db.core.Database.load],
        [`DatabaseMetadata`][kaptive.db.models.DatabaseMetadata]
    """

    aliases = ["info"]

    def setup_arguments(self) -> None:
        r"""Configures argument parser options for the `metadata` subcommand."""
        opts = self.parser.add_argument_group("📥 Inputs")
        opts.add_argument("database", help="Database path or keyword (see: `kaptive db list`)")

    def __call__(self, args: argparse.Namespace) -> None:
        r"""Loads database and prints formatted metadata table to standard output.

        Args:
            args (argparse.Namespace): Parsed CLI arguments containing `database`.

        Raises:
            DatabaseError: If the specified database cannot be loaded.

        See Also:
            [`Database.load`][kaptive.db.core.Database.load]
        """
        from kaptive.db import Database

        db = Database.load(args.database)
        meta = db.metadata
        fields = [
            ("Organism", meta.organism),
            ("Taxon", str(meta.taxon)),
            ("Antigen", meta.antigen),
            ("Pathway", meta.pathway),
            ("Version", meta.version),
            ("Keyword", meta.keyword),
            ("Threshold", f"{meta.id_threshold}%"),
            ("GenBank", meta.genbank),
            ("DOIs", ", ".join(meta.doi) if meta.doi else "None"),
            ("Repository", f"https://github.com/{meta.owner}/{meta.repo}/tree/{meta.branch}"),
            ("Contact", ", ".join(f"{k} <{v}>" for k, v in meta.contact.items())),
        ]

        max_len = max(len(k) for k, v in fields)
        print(
            Colors.wrap(f"\n📊 Metadata for {meta.name}\n", Colors.BOLD_CYAN)
            + "\n".join(f"  {Colors.wrap(k.ljust(max_len), Colors.BOLD)}  {v}" for k, v in fields)
            + "\n"
        )


class Extract(Command):
    r"""📤 Parent command for extracting database records in FASTA format.

    Aggregates subcommands for extracting locus nucleotide sequences (`loci`), gene
    coding sequences (`genes`), and translated amino acid sequences (`proteins`).
    """

    def register_subcommands(self) -> None:
        r"""Registers extraction subcommands (`Loci`, `Genes`, `Proteins`)."""
        self.subcommands = [Loci(), Genes(), Proteins()]

    def get_shared_parser(self) -> argparse.ArgumentParser:
        r"""Creates shared parent parser containing common output arguments.

        Returns:
            argparse.ArgumentParser: Parser configured with `--out` and `--use-indices` flags.
        """
        parser = argparse.ArgumentParser(add_help=False)

        opts = parser.add_argument_group("📥 Inputs")
        opts.add_argument("database", help="Database path or keyword (see: `kaptive db list`)")

        opts = parser.add_argument_group("📤 Outputs")
        opts.add_argument(
            "-o",
            "--out",
            default=sys.stdout.buffer,
            metavar="FILE",
            help="Output file to write fasta to (default: stdout)",
        )
        opts.add_argument(
            "--use-indices",
            action="store_true",
            help="Use numeric indices instead of string IDs for fasta headers",
        )
        return parser


class Loci(Command):
    r"""🧬 Command for extracting locus nucleotide sequences in FASTA format.

    See Also:
        [`Database.load`][kaptive.db.core.Database.load],
        [`Sequences.to_fasta`][kaptive.core.seq.Sequences.to_fasta]
    """

    def __call__(self, args: argparse.Namespace) -> None:
        r"""Extracts locus FASTA sequences to specified output file or stdout.

        Args:
            args (argparse.Namespace): Parsed CLI arguments containing `database`, `out`,
                and `use_indices`.

        See Also:
            [`Database.load`][kaptive.db.core.Database.load]
        """
        self.cli.msg(f"💽 Loading database {args.database}...")
        from kaptive.db import Database

        db = Database.load(args.database)
        out_handle = self.cli.open_file(args.out, "wb")
        self.cli.msg("📤 Extracting loci...")
        out_handle.write(db.loci.to_fasta(args.use_indices))
        self.cli.msg(f"✅ Written locus sequences to {args.out}.")


class Genes(Command):
    r"""🧩 Command for extracting gene coding sequences in FASTA format.

    See Also:
        [`Database.load`][kaptive.db.core.Database.load],
        [`Sequences.to_fasta`][kaptive.core.seq.Sequences.to_fasta]
    """

    def __call__(self, args: argparse.Namespace) -> None:
        r"""Extracts gene nucleotide FASTA sequences to output stream.

        Args:
            args (argparse.Namespace): Parsed CLI arguments containing `database`, `out`,
                and `use_indices`.

        See Also:
            [`Database.load`][kaptive.db.core.Database.load]
        """
        self.cli.msg(f"💽 Loading database {args.database}...")
        from kaptive.db import Database

        db = Database.load(args.database)
        out_handle = self.cli.open_file(args.out, "wb")
        self.cli.msg("📤 Extracting genes...")
        out_handle.write(db.genes.to_fasta(args.use_indices))
        self.cli.msg(f"✅ Written gene sequences to {args.out}.")


class Proteins(Command):
    r"""🧶 Command for extracting translated protein sequences in FASTA format.

    See Also:
        [`Database.load`][kaptive.db.core.Database.load],
        [`Sequences.to_fasta`][kaptive.core.seq.Sequences.to_fasta]
    """

    def __call__(self, args: argparse.Namespace) -> None:
        r"""Extracts translated protein FASTA sequences to output stream.

        Args:
            args (argparse.Namespace): Parsed CLI arguments containing `database`, `out`,
                and `use_indices`.

        See Also:
            [`Database.load`][kaptive.db.core.Database.load]
        """
        self.cli.msg(f"💽 Loading database {args.database}...")
        from kaptive.db import Database

        db = Database.load(args.database)
        out_handle = self.cli.open_file(args.out, "wb")
        self.cli.msg("📤 Extracting proteins...")
        out_handle.write(db.translations.to_fasta(args.use_indices))
        self.cli.msg(f"✅ Written protein sequences to {args.out}.")
