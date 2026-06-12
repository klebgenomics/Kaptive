import sys
import argparse
from abc import ABC, abstractmethod
from pathlib import Path
from typing import IO
from io import IOBase


# Classes --------------------------------------------------------------------------------------------------------------
class Colors:
    """A non-instantiable namespace for ANSI escape sequences.

    This class provides a collection of string constants for terminal text styling.
    It is designed as a pure namespace and cannot be instantiated. The color codes
    are generated programmatically to ensure consistency and prevent typos.
    """
    def __init__(self):
        """Prevents instantiation of this namespace class."""
        raise TypeError("The Colors class is a namespace and cannot be instantiated.")

    # --- Static Attributes ---
    RESET = "\033[0m"
    BOLD = "\033[1m"
    DIM = "\033[2m"
    ITALIC = "\033[3m"
    UNDERLINE = "\033[4m"

    # --- Code Generation ---
    _STYLES = {'': '0', 'BOLD': '1'}
    _COLORS = {
        'BLACK': '30', 'RED': '31', 'GREEN': '32', 'YELLOW': '33',
        'BLUE': '34', 'PURPLE': '35', 'CYAN': '36', 'WHITE': '37'
    }
    _BG_COLORS = {
        'BLACK': '40', 'RED': '41', 'GREEN': '42', 'YELLOW': '43',
        'BLUE': '44', 'PURPLE': '45', 'CYAN': '46', 'WHITE': '47'
    }

    # Dynamically generate attributes for foreground colors
    for _style_name, _style_code in _STYLES.items():
        for _color_name, _color_code in _COLORS.items():
            _attr_name = f"{_style_name}_{_color_name}" if _style_name else _color_name
            vars()[_attr_name] = f"\033[{_style_code};{_color_code}m"

    # Dynamically generate attributes for background colors
    for _color_name, _color_code in _BG_COLORS.items():
        _attr_name = f"BG_{_color_name}"
        vars()[_attr_name] = f"\033[{_color_code}m"

    # Cleanup namespace to avoid polluting the class with temporary variables
    del _STYLES, _COLORS, _BG_COLORS, _style_name, _style_code, _color_name, _color_code, _attr_name

    @classmethod
    def wrap(cls, text: str, *styles: str) -> str:
        """Wraps text in the specified color(s)/style(s) and appends the reset sequence."""
        return f"{''.join(styles)}{text}{cls.RESET}"
    

class Cli:
    """
    Class defining the root Kaptive CLI 
    """
    def __init__(self, description: str = "", epilog: str = ""):
        self.parser = argparse.ArgumentParser(
            description=description,
            epilog=epilog,
            formatter_class=argparse.RawTextHelpFormatter
        )
        self.subparsers = self.parser.add_subparsers(title="Commands", dest="command", required=True)
        self._open_handles = []

    def add_command(self, command: 'Command'):
        """Builds and attaches a top-level Command to the root CLI."""
        command.cli = self
        command.build(self.subparsers)

    def run(self, args: list[str] | None = None):
        """Parses arguments and automatically executes the bound command."""
        parsed_args = self.parser.parse_args(args)
        if hasattr(parsed_args, 'func'):
            parsed_args.func(parsed_args)
        else:
            self.parser.print_help()

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.cleanup()

    def __del__(self):
        self.cleanup()

    def cleanup(self):
        """Safely close any open file handles managed by this CLI."""
        for handle in self._open_handles:
            if handle not in (sys.stdout, sys.stdin, sys.stderr):
                handle.close()
        self._open_handles.clear()

    def write_to_file_or_directory(self, path: str | Path | IO, mode: str = 'at') -> Path | IO:
        """
        Writes to a file or creates a directory based on the provided path.

        If the path is '-' or 'stdout', it returns stdout.
        If the path has a suffix, it's treated as a file and opened for appending.
        If the path has no suffix, it's treated as a directory and created if it doesn't exist.

        :param path: The path to the file or directory.
        :param mode: The mode to open the file in if it's a file.
        :return: A file handle (IO) if a file is specified, or a Path object if a directory is specified.
        """
        if isinstance(path, IOBase):
            return path
        if path == '-' or path == 'stdout':  # If the path is '-', return stdout
            return sys.stdout
        if not isinstance(path, Path):  # Coerce to Path object
            path = Path(path)
        if path.suffix:  # If the path has an extension, it's probably a file
            # NB: We can't use is_file or is_dir because it may not exist yet, open() will create or append
            handle = open(path, mode)
            self._open_handles.append(handle)
            return handle
        else:
            path.mkdir(exist_ok=True, parents=True)  # Create the directory if it doesn't exist
        return path


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
        self.subcommands: list['Command'] = []
        self.cli: Cli | None = None
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

    def execute(self, args: argparse.Namespace):
        """Override to execute the command's logic. If not overridden, command acts as a group/folder."""
        pass

    def build(self, subparsers: argparse._SubParsersAction, parent_parsers: list[argparse.ArgumentParser] | None = None):
        """Wires the command and its children into the argparse structure."""
        parents = parent_parsers or []

        self.parser = subparsers.add_parser(
            name=self.name,
            aliases=self.aliases,
            description=self.description,
            help=self.help_text or self.description,
            parents=parents,
            formatter_class=argparse.RawTextHelpFormatter
        )

        # 1. Add specific arguments for this command
        self.setup_arguments()

        # 2. Bind the execution function (only if execute was actually overridden)
        if type(self).execute != Command.execute:
            self.parser.set_defaults(func=self.execute)

        # 3. Process subcommands (if any)
        if self.subcommands:
            # If this command doesn't do anything itself, it MUST require a subcommand
            is_required = (type(self).execute == Command.execute)
            sub_action = self.parser.add_subparsers(
                title=f"'{self.name}' subcommands",
                dest=f"{self.name}_subcommand",
                required=is_required
            )

            # Collect shared arguments to pass down
            child_parents = parents.copy()
            if shared := self.get_shared_parser():
                child_parents.append(shared)

            for cmd in self.subcommands:
                cmd.cli = self.cli
                cmd.build(sub_action, parent_parsers=child_parents)


class DatabaseCommand(Command):
    """Group command for database operations."""
    name = 'database'
    aliases = ['db']
    description = Colors.wrap('Manage and extract records from databases', Colors.BOLD)

    def register_subcommands(self):
        self.subcommands = [
            DatabaseInstallCommand(),
            DatabaseUpdateCommand(),
            DatabaseResetCommand()
        ]

    def get_shared_parser(self) -> argparse.ArgumentParser:
        # Creates a parent parser for shared arguments. `add_help=False` is critical!
        parser = argparse.ArgumentParser(add_help=False)
        parser.add_argument('db_path_or_keyword', help='Database path or keyword (e.g., kpsc_k)')
        return parser


class DatabaseInstallCommand(Command):
    name = 'install'
    description = 'Install a known database via keyword'

    def execute(self, args: argparse.Namespace):
        from kaptive.db import DatabaseManager
        print(f"Installing database '{args.db_path_or_keyword}'...")
        DatabaseManager.install(args.db_path_or_keyword)
        print(f"Successfully installed '{args.db_path_or_keyword}'.")


class DatabaseUpdateCommand(Command):
    name = 'update'
    description = 'Update installed databases from remote'

    def execute(self, args: argparse.Namespace):
        from kaptive.db import DatabaseManager
        print("Checking for updates...")
        updated = False
        for db in DatabaseManager.update():
            print(f"Updated {db.metadata.name} to version {db.metadata.version}")
            updated = True
        if not updated:
            print("All databases are already up to date.")


class DatabaseResetCommand(Command):
    name = 'reset'
    description = 'Uninstall all local databases'

    def execute(self, args: argparse.Namespace):
        from kaptive.db import DatabaseManager
        DatabaseManager.reset()
        print("All local databases have been uninstalled and reset.")
