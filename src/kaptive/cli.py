r"""Command line interface framework and root runner for Kaptive.

This module provides ANSI color utilities ([`Colors`][kaptive.cli.Colors]), customized
argparse formatters ([`KaptiveHelpFormatter`][kaptive.cli.KaptiveHelpFormatter]), error-handling
parsers ([`HelpOnErrorParser`][kaptive.cli.HelpOnErrorParser]), the core application host
([`Cli`][kaptive.cli.Cli]), and the abstract command specification ([`Command`][kaptive.cli.Command]).
"""

import argparse
import os
import re
import sys
from abc import ABC
from collections.abc import Iterable
from pathlib import Path
from typing import IO, Any, Self

from kaptive import __version__


# Classes --------------------------------------------------------------------------------------------------------------
class Colors:
    r"""Namespace for ANSI terminal color escape sequences.

    Provides constants and styling utility methods for terminal text formatting.
    This class is designed as a pure namespace and cannot be instantiated.

    Attributes:
        ENABLED (bool): Flag indicating whether ANSI color output is active.
        RESET (str): ANSI reset escape sequence.
        BOLD (str): ANSI bold text formatting sequence.
        BOLD_RED (str): ANSI bold red text sequence.
        BOLD_CYAN (str): ANSI bold cyan text sequence.
    """

    ENABLED = sys.stdout.isatty() and not os.environ.get("NO_COLOR")

    def __init__(self) -> None:
        r"""Prevent instantiation of namespace class.

        Raises:
            TypeError: Always raised when instantiated.
        """
        raise TypeError("The Colors class is a namespace and cannot be instantiated.")

    # --- Static Attributes ---
    RESET = "\033[0m"
    BOLD = "\033[1m"
    BOLD_RED = "\033[1;31m"
    BOLD_CYAN = "\033[1;36m"

    @classmethod
    def wrap(cls, text: str, *styles: str) -> str:
        r"""Wrap text in specified ANSI styling escape sequences.

        Appends reset sequence at the end if color output is enabled.

        Args:
            text (str): Input text string to format.
            *styles (str): ANSI style escape sequences to apply.

        Returns:
            str: Styled text string, or original text if ANSI color is disabled.
        """
        if not cls.ENABLED:
            return text
        return f"{''.join(styles)}{text}{cls.RESET}"


class KaptiveHelpFormatter(argparse.RawTextHelpFormatter):
    r"""Custom argparse help formatter adding ANSI colors and layout enhancements.

    Inherits from [`argparse.RawTextHelpFormatter`][argparse.RawTextHelpFormatter] to preserve
    raw multiline descriptions while adding styled section headings and usage placeholders.
    """

    def _format_usage(
        self,
        usage: str | None,
        actions: Iterable[argparse.Action],
        groups: Iterable[argparse._ArgumentGroup],
        prefix: str | None,
    ) -> str:
        r"""Format usage syntax string with ANSI color highlights.

        Args:
            usage (str | None): Custom usage message pattern.
            actions (Iterable[argparse.Action]): Sequence of argument actions.
            groups (Iterable[argparse._ArgumentGroup]): Sequence of argument groups.
            prefix (str | None): Usage line prefix string (defaults to "usage: ").

        Returns:
            str: Colorized usage string.
        """
        # Filter out optional actions to cleanly generate the positional usage string
        positionals = [a for a in actions if not a.option_strings]
        result = super()._format_usage(usage, positionals, groups, prefix)  # type: ignore

        # Replace the subcommand set {add,install,...} with [subcommand]
        result = re.sub(r"\{[a-zA-Z0-9_,\.-]+\}", Colors.wrap("[subcommand]", Colors.BOLD_CYAN), result)

        actual_prefix = prefix if prefix is not None else "usage: "
        target = f"{actual_prefix}{self._prog}"

        if result.startswith(target):
            # Inject [options] cleanly if any optional actions exist
            if any(a.option_strings for a in actions):
                colored_options = Colors.wrap("[options]", Colors.BOLD_CYAN)
                result = result.replace(target, f"{target} {colored_options}", 1)
            # Colorize prefix
            result = result.replace(actual_prefix, Colors.wrap(actual_prefix, Colors.BOLD_CYAN), 1)

        return result

    def start_section(self, heading: str | None) -> None:
        r"""Start a new help output section with a colorized heading.

        Args:
            heading (str | None): Title heading string for the section.
        """
        if heading:
            heading = Colors.wrap(heading, Colors.BOLD_CYAN)
        super().start_section(heading)

    def _format_action(self, action: argparse.Action) -> str:
        r"""Format help output text for a single argument action.

        Args:
            action (argparse.Action): Argument action instance being formatted.

        Returns:
            str: Formatted action description string.
        """
        result = super()._format_action(action)
        # Remove the set string (e.g. {add,install,...}) from the subcommands header
        if type(action).__name__ == "_SubParsersAction":
            lines = result.split("\n", 1)
            if len(lines) > 1:
                result = lines[1]
        return result


class HelpOnErrorParser(argparse.ArgumentParser):
    r"""Argument parser that outputs full help text and suggestions on error.

    Extends [`argparse.ArgumentParser`][argparse.ArgumentParser] to print help menu to
    `sys.stderr` when invalid arguments or choices are encountered, suggesting close matches.
    """

    def error(self, message: str) -> None:  # type: ignore
        r"""Print help text and error message before terminating CLI session.

        Args:
            message (str): Error message string reported by parser.

        Raises:
            SystemExit: Exits process with status code 2.
        """
        if match := re.search(r"invalid choice: '?([^']+)'? \(choose from (.*)\)", message):
            invalid = match.group(1)
            choices = [c.strip("'").strip() for c in match.group(2).split(", ")]
            from difflib import get_close_matches

            if matches := get_close_matches(invalid, choices):
                message += f"\n    💡 Did you mean '{Colors.wrap(matches[0], Colors.BOLD_CYAN)}'?"

        self.print_help(sys.stderr)
        self.exit(2, f"\n{Colors.wrap('❌ Error:', Colors.BOLD_RED)} {message}\n")


class Cli:
    r"""Root Command Line Interface runner and context manager for Kaptive.

    Coordinates global option parsing, subcommand registration, safe file handle
    tracking, error trapping, and progress updates.

    Attributes:
        verbose (bool): Flag indicating whether verbose logging is enabled.
        global_parser (HelpOnErrorParser): Parser handling shared global options.
        parser (HelpOnErrorParser): Main argument parser for application commands.
        subparsers (argparse._SubParsersAction): Subparser registry for subcommands.
    """

    def __init__(self, description: str | None = None, epilog: str | None = None) -> None:
        r"""Initialize root CLI application instance.

        Args:
            description (str | None): Header description text displayed in top-level help.
            epilog (str | None): Footer epilog text displayed in top-level help.
        """
        self.verbose = False
        self.global_parser = HelpOnErrorParser(add_help=False)
        self.global_parser.add_argument("-V", "--verbose", action="store_true", help="Enable verbose output/progress")

        self.parser = HelpOnErrorParser(
            description=Colors.wrap(description, Colors.BOLD) if description else description,
            epilog=Colors.wrap(epilog, Colors.BOLD) if epilog else epilog,
            parents=[self.global_parser],
            formatter_class=KaptiveHelpFormatter,
        )

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
        self._open_handles: list[IO] = []  # type: ignore

    def add_command(self, command: "Command") -> None:
        r"""Build and attach top-level command to root CLI subparsers.

        Args:
            command (Command): Subcommand instance implementing [`Command`][kaptive.cli.Command].
        """
        command.cli = self
        command.build(self.subparsers, parent_parsers=[self.global_parser])

    def run(self, args: list[str] | None = None) -> None:
        r"""Parse arguments and execute bound subcommand handler.

        Args:
            args (list[str] | None): Command line arguments to parse. Defaults to `sys.argv[1:]`.

        Raises:
            SystemExit: If database error or web client error occurs during command execution.
        """
        parsed_args = self.parser.parse_args(args)
        self.verbose = getattr(parsed_args, "verbose", False)

        if hasattr(parsed_args, "func"):
            from kaptive.client import KaptiveWebClientError
            from kaptive.db import DatabaseError

            try:
                parsed_args.func(parsed_args)
            except (DatabaseError, KaptiveWebClientError) as e:
                print(f"❌ {e}", file=sys.stderr)
                sys.exit(1)
        else:
            self.parser.print_help()

    def __enter__(self) -> Self:
        r"""Enter context manager scope returning CLI instance.

        Returns:
            Cli: The active CLI host instance.
        """
        return self

    def __exit__(
        self,
        exc_type: type[BaseException] | None,
        exc_val: BaseException | None,
        exc_tb: Any | None,
    ) -> None:
        r"""Exit context manager scope and clean up allocated resources.

        Args:
            exc_type (type[BaseException] | None): Exception type raised within context, if any.
            exc_val (BaseException | None): Exception instance raised within context, if any.
            exc_tb (Any | None): Traceback object associated with exception, if any.

        Raises:
            SystemExit: Terminates process with standard exit codes for handled signals.
        """
        self.cleanup()
        if exc_type is KeyboardInterrupt:
            print("\n🛑 Cancelled by user.", file=sys.stderr)
            sys.exit(1)
        elif exc_type is BrokenPipeError:
            devnull = os.open(os.devnull, os.O_WRONLY)
            os.dup2(devnull, sys.stdout.fileno())
            sys.exit(130)
        elif exc_type is PermissionError:
            print(f"🔒 Permission denied: {exc_val}", file=sys.stderr)
            sys.exit(1)
        elif exc_type is FileNotFoundError:
            print(f"📄 File not found: {exc_val}", file=sys.stderr)
            sys.exit(1)

    def exit(self, msg: str, code: int = 1) -> None:
        r"""Print error message to stderr and terminate process.

        Args:
            msg (str): Error message string to report.
            code (int): Process exit code to pass to `sys.exit`. Defaults to 1.

        Raises:
            SystemExit: Exits process with specified status code.
        """
        print(f"❌ {msg}", file=sys.stderr)
        sys.exit(code)

    def __del__(self) -> None:
        r"""Clean up resource handles upon garbage collection."""
        self.cleanup()

    def cleanup(self) -> None:
        r"""Close any open file handles tracked by this CLI instance."""
        for handle in self._open_handles:
            if handle not in (sys.stdout, sys.stdin, sys.stderr):
                handle.close()
        self._open_handles.clear()

    def msg(self, msg: str | None, **kwargs: Any) -> None:
        r"""Print informational message to stderr if verbose mode is enabled.

        Args:
            msg (str | None): Message string to print.
            **kwargs (Any): Additional keyword arguments passed to `print`.
        """
        if self.verbose:
            print(msg, file=sys.stderr, **kwargs)

    def progress(self, iterable: Iterable, msg: str) -> Iterable:  # type: ignore
        r"""Wrap iterable to display live progress counter on stderr when verbose mode is active.

        Args:
            iterable (Iterable): Source collection or generator to iterate over.
            msg (str): Prefix status message displayed during iteration.

        Returns:
            Iterable: Generator yielding items from original iterable while rendering progress.
        """
        try:
            total = len(iterable)  # type: ignore
        except TypeError:
            total = "?"

        for i, item in enumerate(iterable, start=1):
            if self.verbose:
                print(f"\r{msg} {i}/{total}", end="", file=sys.stderr, flush=True)
            yield item

        if self.verbose:
            print(file=sys.stderr)

    def open_file(self, file: str, mode: str = "rb") -> IO:  # type: ignore
        r"""Open file path or return standard stdout stream handle with tracking.

        Args:
            file (str): File path string or `"-"` / `"stdout"` for standard output stream.
            mode (str): File open mode flag. Defaults to `"rb"`.

        Returns:
            IO: Open file handle or system stream object.
        """
        if file == "-" or file == "stdout":  # If the path is '-', return stdout
            return sys.stdout.buffer if "b" in mode else sys.stdout

        handle = open(file, mode)
        self._open_handles.append(handle)
        return handle


class Command(ABC):
    r"""Abstract base class for Kaptive CLI subcommands.

    Provides a declarative structure for defining CLI argument options, nested subcommands,
    and command execution handlers.

    Attributes:
        name (str): Identifier name used to trigger command on CLI.
        aliases (list[str]): Alternative trigger names or abbreviations for command.
        description (str): Detailed text explanation shown in command help.
        help_text (str): Brief single-line summary shown in parent subparser tables.
        parser (argparse.ArgumentParser | None): Bound argument parser instance.
        subcommands (list[Command]): Child subcommand instances.
        cli (Cli | None): Parent root CLI host application reference.
    """

    name: str = ""
    aliases: list[str] = []
    description: str = ""
    help_text: str = ""

    def __init__(self) -> None:
        r"""Initialize command instance and populate metadata attributes."""
        self.parser: argparse.ArgumentParser | None = None
        self.subcommands: list[Command] = []
        self.cli: Cli | None = None

        # Auto-populate name from the class name
        if not self.name:
            self.name = type(self).__name__.lower()

        # Auto-populate description from class docstring
        if not self.description:
            if type(self).__doc__ and type(self).__doc__ != Command.__doc__:
                self.description = type(self).__doc__  # type: ignore

        # Auto-populate short help text from the first line of the description
        if not self.help_text and self.description:
            self.help_text = self.description.strip().split("\n")[0]

        self.register_subcommands()

    def register_subcommands(self) -> None:
        r"""Register child subcommand instances into `subcommands` list."""
        pass

    def setup_arguments(self) -> None:
        r"""Configure argument flags and positional parameters on `self.parser`."""
        pass

    def get_shared_parser(self) -> argparse.ArgumentParser | None:
        r"""Return shared parent parser containing options passed to subcommands.

        Returns:
            argparse.ArgumentParser | None: Shared non-help argument parser or `None`.
        """
        return None

    def add_output_arguments(
        self,
        opts: argparse._ArgumentGroup,
        tsv_flags: tuple[str, str] = ("-o", "--out"),
        include_json: bool = True,
    ) -> None:
        r"""Add standard report and FASTA file output options to argument group.

        Args:
            opts (argparse._ArgumentGroup): Target argument group to populate.
            tsv_flags (tuple[str, str]): Short and long flag options for TSV report output.
                Defaults to `("-o", "--out")`.
            include_json (bool): Flag indicating whether to include `--json` argument option.
                Defaults to `True`.
        """
        help_msg = (
            "Write serotyping results as a TSV report to a file (default: %(default)s)"
            if tsv_flags[0] == "-o"
            else "Write serotyping results as a TSV report to a file (default: %(const)s)"
        )
        opts.add_argument(
            tsv_flags[0],
            tsv_flags[1],
            metavar="FILE",
            nargs="?" if tsv_flags[0] == "-t" else None,
            default="stdout" if tsv_flags[0] == "-o" else None,
            const="stdout" if tsv_flags[0] == "-t" else None,
            help=help_msg,
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
        if include_json:
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

    def __call__(self, args: argparse.Namespace) -> None:
        r"""Execute command business logic for parsed arguments.

        Args:
            args (argparse.Namespace): Parsed command line argument namespace.
        """
        pass

    def build(
        self,
        subparsers: argparse._SubParsersAction,  # type: ignore
        parent_parsers: list[argparse.ArgumentParser] | None = None,
    ) -> None:
        r"""Wire command parser and subcommands into parent argparse hierarchy.

        Args:
            subparsers (argparse._SubParsersAction): Target subparser action registry.
            parent_parsers (list[argparse.ArgumentParser] | None): Parent shared parsers to inherit.
        """
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


def main() -> None:
    r"""Execute main entry point for the Kaptive CLI application."""
    from kaptive.db.cli import Database
    from kaptive.serotyping.cli import Convert, Type

    description = "🦠 Kaptive: The tool for in silico serotyping of surface antigen loci."
    epilog = "📚 For more information and documentation, visit https://klebgenomics.github.io/Kaptive/"

    with Cli(description=description, epilog=epilog) as app:
        app.add_command(Database())
        app.add_command(Type())
        app.add_command(Convert())
        app.run()


if __name__ == "__main__":
    main()
