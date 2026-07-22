import argparse
import os
import re
import sys
from abc import ABC
from collections.abc import Iterable, Sized
from importlib import metadata
from pathlib import Path
from typing import IO


# Classes --------------------------------------------------------------------------------------------------------------
class Colors:
    """A non-instantiable namespace for ANSI escape sequences.

    This class provides a collection of string constants for terminal text styling.
    It is designed as a pure namespace and cannot be instantiated.
    """

    ENABLED = sys.stdout.isatty() and not os.environ.get("NO_COLOR")

    def __init__(self):
        """Prevents instantiation of this namespace class."""
        raise TypeError("The Colors class is a namespace and cannot be instantiated.")

    # --- Static Attributes ---
    RESET = "\033[0m"
    BOLD = "\033[1m"
    BOLD_RED = "\033[1;31m"
    BOLD_CYAN = "\033[1;36m"

    @classmethod
    def wrap(cls, text: str, *styles: str) -> "str":
        """Wraps text in the specified color(s)/style(s) and appends the reset sequence."""
        if not cls.ENABLED:
            return text
        return f"{''.join(styles)}{text}{cls.RESET}"


class KaptiveHelpFormatter(argparse.RawTextHelpFormatter):
    """Custom formatter to add ANSI colors to argparse output!"""

    def _format_usage(self, usage, actions, groups, prefix):
        # Filter out optional actions to cleanly generate the positional usage string
        positionals = [a for a in actions if not a.option_strings]
        result = super()._format_usage(usage, positionals, groups, prefix)

        # Replace the subcommand set {add,install,...} with [subcommand]
        result = re.sub(r'\{[a-zA-Z0-9_,\.-]+\}', Colors.wrap('[subcommand]', Colors.BOLD_CYAN), result)
        
        actual_prefix = prefix if prefix is not None else "usage: "
        target = f"{actual_prefix}{self._prog}"
        
        if result.startswith(target):
            # Inject [options] cleanly if any optional actions exist
            if any(a.option_strings for a in actions):
                colored_options = Colors.wrap('[options]', Colors.BOLD_CYAN)
                result = result.replace(target, f"{target} {colored_options}", 1)
            # Colorize prefix
            result = result.replace(actual_prefix, Colors.wrap(actual_prefix, Colors.BOLD_CYAN), 1)
            
        return result

    def start_section(self, heading):
        if heading:
            heading = Colors.wrap(heading, Colors.BOLD_CYAN)
        super().start_section(heading)

    def _format_action(self, action):
        result = super()._format_action(action)
        # Remove the set string (e.g. {add,install,...}) from the subcommands header
        if type(action).__name__ == '_SubParsersAction':
            lines = result.split('\n', 1)
            if len(lines) > 1:
                result = lines[1]
        return result


class HelpOnErrorParser(argparse.ArgumentParser):
    """An ArgumentParser that prints full help on error."""
    
    def error(self, message):
        if match := re.search(r"invalid choice: '?([^']+)'? \(choose from (.*)\)", message):
            invalid = match.group(1)
            choices = [c.strip("'").strip() for c in match.group(2).split(", ")]
            from difflib import get_close_matches
            if matches := get_close_matches(invalid, choices):
                message += f"\n    💡 Did you mean '{Colors.wrap(matches[0], Colors.BOLD_CYAN)}'?"
                
        self.print_help(sys.stderr)
        self.exit(2,
                  f"\n{Colors.wrap('❌ Error:', Colors.BOLD_RED)} {message}\n")


class Cli:
    """
    Class defining the root Kaptive CLI
    """

    def __init__(self, description: str | None = None, epilog: str | None = None):
        self.verbose = False
        self.global_parser = HelpOnErrorParser(add_help=False)
        self.global_parser.add_argument(
            "-V", "--verbose", action="store_true", help="Enable verbose output/progress"
        )

        self.parser = HelpOnErrorParser(
            description=Colors.wrap(description, Colors.BOLD)
            if description
            else description,
            epilog=Colors.wrap(epilog, Colors.BOLD) if epilog else epilog,
            parents=[self.global_parser],
            formatter_class=KaptiveHelpFormatter,
        )

        try:
            __version__ = metadata.version("kaptive")
        except metadata.PackageNotFoundError:
            __version__ = "unknown"

        self.parser.add_argument(
            "-v",
            "--version",
            action="version",
            version=f"%(prog)s {__version__}",
            help="Show program's version number and exit",
        )

        # Rename the default options group for the root parser
        if hasattr(self.parser, "_optionals"):
            self.parser._optionals.title = Colors.wrap("🌎 Global options", Colors.BOLD)

        self.subparsers = self.parser.add_subparsers(
            title=Colors.wrap("💬 Commands", Colors.BOLD), dest="command", required=True
        )
        self._open_handles = []

    def add_command(self, command: "Command"):
        """Builds and attaches a top-level Command to the root CLI."""
        command.cli = self
        command.build(self.subparsers, parent_parsers=[self.global_parser])

    def run(self, args: list[str] | None = None):
        """Parses arguments and automatically executes the bound command."""
        parsed_args = self.parser.parse_args(args)
        self.verbose = getattr(parsed_args, "verbose", False)

        if hasattr(parsed_args, "func"):
            from kaptive.client import KaptiveWebClientError
            from kaptive.db import DatabaseError
            try:
                parsed_args.func(parsed_args)
            except (DatabaseError, KaptiveWebClientError) as e:
                self.msg(f"❌ {e}")
                sys.exit(1)
        else:
            self.parser.print_help()

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.cleanup()
        if exc_type is KeyboardInterrupt:
            self.msg("\n🛑 Cancelled by user.")
            sys.exit(1)
        elif exc_type is BrokenPipeError:
            devnull = os.open(os.devnull, os.O_WRONLY)
            os.dup2(devnull, sys.stdout.fileno())
            sys.exit(130)
        elif exc_type is PermissionError:
            self.msg(f"🔒 Permission denied: {exc_val}")
            sys.exit(1)
        elif exc_type is FileNotFoundError:
            self.msg(f"📄 File not found: {exc_val}")
            sys.exit(1)

    def exit(self, msg: str, code: int = 1):
        """Prints a message to stderr and cleanly exits."""
        self.msg(f"❌ {msg}")
        sys.exit(code)

    def __del__(self):
        self.cleanup()

    def cleanup(self):
        """Safely close any open file handles managed by this CLI."""
        for handle in self._open_handles:
            if handle not in (sys.stdout, sys.stdin, sys.stderr):
                handle.close()
        self._open_handles.clear()

    def msg(self, msg: str | None, **kwargs) -> None:
        """
        Prints an informational message to stderr.

        Informational messages are routed to stderr by default to keep stdout clean
        for data piping and redirection. The message is only printed if the --verbose flag was provided.
        """
        if self.verbose:
            print(msg, file=sys.stderr, **kwargs)

    def progress(self, iterable: Iterable, msg: str) -> Iterable:
        """
        Wraps an iterable to show a progress message when verbose is True.
        """
        try:
            total = len(iterable)
        except TypeError:
            total = "?"

        for i, item in enumerate(iterable, start=1):
            if self.verbose:
                print(f"\r{msg} {i}/{total}", end="", file=sys.stderr, flush=True)
            yield item

        if self.verbose:
            print(file=sys.stderr)

    def open_file(self, file: str, mode: str = "rb") -> IO:
        if file == "-" or file == "stdout":  # If the path is '-', return stdout
            return sys.stdout.buffer if "b" in mode else sys.stdout
            # We don't add stdout/stdin to _open_handles

        handle = open(file, mode)
        self._open_handles.append(handle)
        return handle


class Command(ABC):
    """
    Abstract base class for Kaptive CLI Commands.
    Provides a Click-like interface for nesting commands and sharing arguments using argparse.
    """

    name: str = ""
    aliases: list[str] = []
    description: str = ""
    help_text: str = ""

    def __init__(self):
        self.parser: argparse.ArgumentParser | None = None
        self.subcommands: list["Command"] = []
        self.cli: Cli | None = None

        # Auto-populate name from the class name
        if not self.name:
            self.name = type(self).__name__.lower()

        # Auto-populate description from class docstring
        if not self.description:
            if type(self).__doc__ and type(self).__doc__ != Command.__doc__:
                self.description = type(self).__doc__  # cleandoc(type(self).__doc__)

        # Auto-populate short help text from the first line of the description
        if not self.help_text and self.description:
            self.help_text = self.description.strip().split("\n")[0]

        self.register_subcommands()

    def register_subcommands(self):
        """Override to populate self.subcommands with child Command instances."""
        pass

    def setup_arguments(self):
        """Override to add specific arguments to self.parser."""
        pass

    def get_shared_parser(self) -> argparse.ArgumentParser | None:
        """Override to return an ArgumentParser(add_help=False) containing arguments shared with children."""
        return None

    def __call__(self, args: argparse.Namespace):
        """Override to __call__ the command's logic. If not overridden, command acts as a group/folder."""
        pass

    def build(
        self,
        subparsers: argparse._SubParsersAction,
        parent_parsers: list[argparse.ArgumentParser] | None = None,
    ):
        """Wires the command and its children into the argparse structure."""
        parents = parent_parsers or []

        self.parser = subparsers.add_parser(
            name=self.name,
            aliases=self.aliases,
            description=Colors.wrap(self.description, Colors.BOLD),
            help=self.help_text or self.description,
            parents=parents,
            formatter_class=KaptiveHelpFormatter,
        )

        # 1. Add specific arguments for this command
        self.setup_arguments()

        # Rename the default options group and move it to the bottom of the help menu
        if hasattr(self.parser, "_optionals"):
            self.parser._optionals.title = Colors.wrap("🌎 Global options", Colors.BOLD)
            # Pop the group out of the internal list and append it to the end
            groups = self.parser._action_groups
            if self.parser._optionals in groups:
                groups.append(groups.pop(groups.index(self.parser._optionals)))

        # 2. Bind the execution function (only if __call__ was actually overridden)
        if type(self).__call__ != Command.__call__:
            self.parser.set_defaults(func=self.__call__)

        # 3. Process subcommands (if any)
        if self.subcommands:
            # If this command doesn't do anything itself, it MUST require a subcommand
            is_required = type(self).__call__ == Command.__call__
            sub_action = self.parser.add_subparsers(
                title=Colors.wrap(f"'{self.name}' subcommands", Colors.BOLD),
                dest=f"{self.name}_subcommand",
                required=is_required,
            )

            # Collect shared arguments to pass down
            child_parents = parents.copy()
            if shared := self.get_shared_parser():
                child_parents.append(shared)

            for cmd in self.subcommands:
                cmd.cli = self.cli
                cmd.build(sub_action, parent_parsers=child_parents)


# List command ---------------------------------------------------------------------------------------------------------
class List(Command):
    """📋 Lists installed databases"""

    aliases = ["ls"]

    def __call__(self, args: argparse.Namespace):
        from kaptive.db import DatabaseManager

        if installed := DatabaseManager.installed():
            print("\n".join(installed))
        else:
            self.cli.msg("❌ No databases installed")


# Database command -----------------------------------------------------------------------------------------------------
class Database(Command):
    """📦 Manage Kaptive databases"""

    aliases = ["db"]

    def register_subcommands(self):
        self.subcommands = [List(), Add(), Install(), Update(), Reset(), Extract(), Metadata()]


class Install(Command):
    """📦 Install a known database via keyword (or 'all')."""

    def setup_arguments(self):
        opts = self.parser.add_argument_group("📥 Inputs")
        opts.add_argument("database", help="Database keyword (see: `kaptive db list`) or 'all'")

    def __call__(self, args: argparse.Namespace):
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
    """🔄 Update installed databases from their remote repositories."""

    def setup_arguments(self):
        opts = self.parser.add_argument_group("📥 Inputs")
        opts.add_argument(
            "database",
            nargs="?",
            default="all",
            help="Database keyword (see: `kaptive db list`) or 'all' (default: all)"
        )

    def __call__(self, args: argparse.Namespace):
        if args.database == "all":
            self.cli.msg("🔄 Checking all installed databases for updates concurrently...")
        else:
            self.cli.msg(f"🔄 Checking '{args.database}' for updates...")
            
        from kaptive.db import DatabaseManager
        updated = False
        for db in DatabaseManager.update(args.database):
            self.cli.msg(
                f"✅ Updated {db.metadata.name} to version {db.metadata.version}"
            )
            updated = True
            
        if not updated:
            self.cli.msg("🎉 All databases are already up to date.")


class Reset(Command):
    """🧹 Uninstall all local databases and clear the cache."""

    def __call__(self, args: argparse.Namespace):
        self.cli.msg("🧹 Uninstalling all local databases...")
        from kaptive.db import DatabaseManager

        DatabaseManager.reset()
        self.cli.msg("✅ All local databases have been uninstalled and reset.")


class Add(Command):
    """🔗 Add a new database from a GitHub repository."""

    def setup_arguments(self):
        opts = self.parser.add_argument_group("📥 Inputs")
        opts.add_argument("database", help="Name for the new database")

        opts = self.parser.add_argument_group(
            Colors.wrap("🌐 GitHub Details", Colors.BOLD)
        )
        opts.add_argument("owner", help="GitHub repository owner")
        opts.add_argument("repo_name", help="GitHub repository name")
        opts.add_argument(
            "-b",
            "--branch",
            help="GitHub repository branch (default: main)",
            default="main",
            nargs="?",
        )

    def __call__(self, args: argparse.Namespace):
        from kaptive.db import DatabaseManager

        self.cli.msg( f"⤵️ Adding {args.database} from {args.owner}/{args.repo_name}/{args.branch}")
        if db := DatabaseManager.add(
            args.owner, args.repo_name, args.database, args.branch
        ):
            self.cli.msg(f"✅ Added {db.metadata.name} v{db.metadata.version} successfully!")
        else:
            self.cli.msg("❌ Failed to add database! Is it already installed?")


class Metadata(Command):
    """📊 Print metadata of a database"""

    aliases = ['info']

    def setup_arguments(self):
        opts = self.parser.add_argument_group("📥 Inputs")
        opts.add_argument("database", help="Database path or keyword (see: `kaptive db list`)")

    def __call__(self, args: argparse.Namespace):
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
            Colors.wrap(f"\n📊 Metadata for {meta.name}\n", Colors.BOLD_CYAN) +
            '\n'.join(f"  {Colors.wrap(k.ljust(max_len), Colors.BOLD)}  {v}" for k, v in fields) +
            '\n'
        )


class Extract(Command):
    """📤 Extract records from a known database."""

    def register_subcommands(self):
        self.subcommands = [Loci(), Genes(), Proteins()]

    def get_shared_parser(self) -> argparse.ArgumentParser:
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
    """🧬 Extracts locus sequences from the database in fasta format."""

    def __call__(self, args: argparse.Namespace):
        self.cli.msg(f"💽 Loading database {args.database}...")
        from kaptive.db import Database

        db = Database.load(args.database)
        out_handle = self.cli.open_file(args.out, "wb")
        self.cli.msg("📤 Extracting loci...")
        out_handle.write(db.loci.to_fasta(args.use_indices))
        self.cli.msg(f"✅ Written locus sequences to {args.out}.")


class Genes(Command):
    """🧩 Extracts gene sequences from the database in fasta format."""

    def __call__(self, args: argparse.Namespace):
        self.cli.msg(f"💽 Loading database {args.database}...")
        from kaptive.db import Database

        db = Database.load(args.database)
        out_handle = self.cli.open_file(args.out, "wb")
        self.cli.msg("📤 Extracting genes...")
        out_handle.write(db.genes.to_fasta(args.use_indices))
        self.cli.msg(f"✅ Written gene sequences to {args.out}.")


class Proteins(Command):
    """🧶 Extracts protein sequences from the database in fasta format."""

    def __call__(self, args: argparse.Namespace):
        self.cli.msg(f"💽 Loading database {args.database}...")
        from kaptive.db import Database

        db = Database.load(args.database)
        out_handle = self.cli.open_file(args.out, "wb")
        self.cli.msg("📤 Extracting proteins...")
        out_handle.write(db.translations.to_fasta(args.use_indices))
        self.cli.msg(f"✅ Written protein sequences to {args.out}.")


# Exporter -------------------------------------------------------------------------------------------------------------
class ResultExporter:
    """
    Evaluates CLI arguments once and sets up a pipeline of output writers.
    This eliminates conditional branching inside the processing loop and allows
    reusability between the 'Type' and 'Convert' commands.
    """
    file_suffix = 'kaptive_results'
    def __init__(self, cli: Cli, args: argparse.Namespace):
        self.writers = []

        if tsv_file := getattr(args, "out", getattr(args, "tsv", None)):
            from kaptive.serotyping import KaptiveRow
            tsv_handle = cli.open_file(tsv_file, mode="wb")
            tsv_handle.write(KaptiveRow.header())
            self.writers.append(lambda r: tsv_handle.write(bytes(KaptiveRow.from_result(r))))

        if pha4ge_file := getattr(args, "pha4ge", None):
            from kaptive.serotyping import Pha4geRow
            pha4ge_handle = cli.open_file(pha4ge_file, mode="wb")
            pha4ge_handle.write(Pha4geRow.header())
            self.writers.append(lambda r: pha4ge_handle.write(bytes(Pha4geRow.from_result(r))))

        if json_file := getattr(args, "json", None):
            try:
                from orjson import dumps, OPT_SERIALIZE_NUMPY, OPT_APPEND_NEWLINE
            except ImportError:
                cli.exit("orjson not installed. Please run: pip install kaptive[json]")
            json_handle = cli.open_file(json_file, mode="wb")
            self.writers.append(
                lambda r: json_handle.write(dumps(r.to_dict(), option=OPT_SERIALIZE_NUMPY | OPT_APPEND_NEWLINE))
            )

        if loci_dir := getattr(args, "loci", None):  # type: Path
            self.writers.append(lambda r: (loci_dir / f"{r.genome}_{self.file_suffix}.fna"
                                           ).write_bytes(r.locus_seqs.to_fasta()))

        if genes_dir := getattr(args, "genes", None):  # type: Path
            self.writers.append(lambda r: (genes_dir / f"{r.genome}_{self.file_suffix}.ffn"
                                           ).write_bytes(r.gene_seqs.to_fasta()))

        if proteins_dir := getattr(args, "proteins", None):  # type: Path
            self.writers.append(lambda r: (proteins_dir / f"{r.genome}_{self.file_suffix}.faa"
                                           ).write_bytes(r.translations.to_fasta()))

        if plot_dir := getattr(args, "plots", None):  # type: Path
            try:
                from kaptive.plotting import SerotypingResultPlotter
            except ImportError:
                cli.exit("plotly not installed. Please run: pip install kaptive[plot]")
            self.writers.append(lambda r: SerotypingResultPlotter(r).write_html(
                plot_dir / f"{r.genome}_{self.file_suffix}.html", include_plotlyjs=False, full_html=True))

    def __call__(self, result):
        """Pass the result to all registered output writers."""
        for write in self.writers:
            write(result)


# Type command ---------------------------------------------------------------------------------------------------------
class Type(Command):
    """💉 In silico serotyping of assemblies."""

    aliases = ["assembly"]  # Backwards compatibility with < v3.3

    def setup_arguments(self):
        opts = self.parser.add_argument_group(Colors.wrap("📥 Inputs", Colors.BOLD))
        opts.add_argument("database", help="Database path or keyword (see: `kaptive db list`)")
        opts.add_argument(
            "genomes",
            nargs="+",
            help="Genome assemblies in fasta format; can be compressed",
        )

        opts = self.parser.add_argument_group(Colors.wrap("📤 Outputs", Colors.BOLD))
        opts.add_argument(
            "-o",
            "--out",
            metavar="FILE",
            default="stdout",
            help="Write serotyping results as a TSV report to a file (default: %(default)s)",
        )
        opts.add_argument(
            "-l",
            "--loci",
            metavar="DIR",
            nargs="?",
            const="./",
            type=Path,
            help="Write locus nucleotide fasta files to a directory (default: %(const)s)",
        )
        opts.add_argument(
            "-g",
            "--genes",
            metavar="DIR",
            nargs="?",
            const="./",
            type=Path,
            help="Write gene nucleotide fasta files to a directory (default: %(const)s)",
        )
        opts.add_argument(
            "-p",
            "--proteins",
            metavar="DIR",
            nargs="?",
            const="./",
            type=Path,
            help="Write translation amino-acid fasta files to a directory (default: %(const)s)",
        )
        opts.add_argument(
            "-j",
            "--json",
            metavar="FILE",
            nargs="?",
            const="kaptive_results.jsonl",
            help="Write serialised results to a newline-delimited JSON (default: %(const)s)",
        )
        opts.add_argument(
            "--pha4ge",
            metavar="FILE",
            nargs="?",
            const="kaptive_results.pha4ge",
            type=Path,
            help="Write PHA4GE-compliant serotyping report to a TSV file (default: %(const)s)",
        )
        opts.add_argument(
            "--plots",
            metavar="DIR",
            nargs="?",
            const="./",
            type=Path,
            help="Generate interactive locus plots to a directory (default: %(const)s)",
        )

        opts = self.parser.add_argument_group(
            Colors.wrap("🔬 Confidence options", Colors.BOLD)
        )
        opts.add_argument(
            "--max-other-genes",
            type=int,
            metavar="",
            default=1,
            help="Typeable if <= other genes (default: %(default)s)",
        )
        opts.add_argument(
            "--min-completeness",
            type=float,
            metavar="",
            default=0.5,
            help="Typeable if >= completeness (default: %(default)s)",
        )
        opts.add_argument(
            "--below-threshold",
            action="store_true",
            help="Typeable if any genes in locus are below threshold (default: False)",
        )

        opts = self.parser.add_argument_group(
            Colors.wrap("🔧 Other options", Colors.BOLD)
        )
        opts.add_argument(
            "-t",
            "--threads",
            type=int,
            default=0,
            metavar="",
            help="Number threads or 0 for all available (default: 0)",
        )

    def __call__(self, args: argparse.Namespace):
        self.cli.msg(f"💽 Loading database {args.database}...")
        from kaptive.db import Database
        from kaptive.serotyping import Serotyper

        db = Database.load(args.database)
        exporter = ResultExporter(self.cli, args)

        serotyper = Serotyper(
            db=db,
            max_other_genes=args.max_other_genes,
            min_completeness=args.min_completeness,
            allow_below_threshold=args.below_threshold,
        )
        for genome in self.cli.progress(args.genomes, "💉 Serotyping genomes..."):
            if result := serotyper(genome):
                exporter(result)

        self.cli.msg(f"✅ Serotyping complete. Results written to '{args.out}'.")


# Client command -------------------------------------------------------------------------------------------------------
class Convert(Command):
    """💱 Convert serialized Kaptive results into different formats."""

    def setup_arguments(self):
        opts = self.parser.add_argument_group(Colors.wrap("📥 Inputs", Colors.BOLD))
        opts.add_argument(
            "jsonl",
            default='stdin',
            help="Serialised results in JSON-lines format (default: stdin)",
        )

        opts = self.parser.add_argument_group(Colors.wrap("📤 Outputs", Colors.BOLD))
        opts.add_argument(
            "-t",
            "--tsv",
            metavar="FILE",
            nargs="?",
            const="stdout",
            help="Write serotyping results as a TSV report to a file (default: %(const)s)",
        )
        opts.add_argument(
            "-l",
            "--loci",
            metavar="DIR",
            nargs="?",
            const="./",
            type=Path,
            help="Write locus nucleotide fasta files to a directory (default: %(const)s)",
        )
        opts.add_argument(
            "-g",
            "--genes",
            metavar="DIR",
            nargs="?",
            const="./",
            type=Path,
            help="Write gene nucleotide fasta files to a directory (default: %(const)s)",
        )
        opts.add_argument(
            "-p",
            "--proteins",
            metavar="DIR",
            nargs="?",
            const="./",
            type=Path,
            help="Write translation amino-acid fasta files to a directory (default: %(const)s)",
        )
        opts.add_argument(
            "--pha4ge",
            metavar="FILE",
            nargs="?",
            const="kaptive_results.pha4ge",
            type=Path,
            help="Write PHA4GE-compliant serotyping report to a TSV file (default: %(const)s)",
        )
        opts.add_argument(
            "--plots",
            metavar="DIR",
            nargs="?",
            const="./",
            type=Path,
            help="Generate interactive locus plots to a directory (default: %(const)s)",
        )

    def __call__(self, args: argparse.Namespace):
        try:
            from orjson import loads
        except ImportError:
            self.cli.exit("orjson not installed. Please run: pip install kaptive[json]")

        from kaptive.serotyping import SerotypingResult

        exporter = ResultExporter(self.cli, args)

        handle = self.cli.open_file(args.jsonl, mode="rb")
        for line in self.cli.progress(handle, "💱 Converting results..."):
            line = line.strip()
            if not line:
                continue
            
            result = SerotypingResult.from_dict(loads(line))
            exporter(result)
            
        self.cli.msg("✅ Conversion complete.")


def main():
    """Main entry point for the Kaptive CLI."""
    description = (
        "🦠 Kaptive: The tool for in silico serotyping of surface antigen loci."
    )
    epilog = "📚 For more information and documentation, visit https://klebgenomics.github.io/Kaptive/"

    with Cli(description=description, epilog=epilog) as app:
        app.add_command(Database())
        app.add_command(Type())
        app.add_command(Convert())
        app.run()


if __name__ == "__main__":
    main()
