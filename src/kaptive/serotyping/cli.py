r"""Command line interface commands and exporter for serotyping.

This module provides CLI command implementations for performing *in silico* serotyping
on genome assemblies and converting serialized serotyping results to various tabular,
JSON, or graphical output formats.

Classes:
    [`ResultExporter`][kaptive.serotyping.cli.ResultExporter]: Evaluates output options and dispatches
        serialization tasks.
    [`Type`][kaptive.serotyping.cli.Type]: Subcommand for *in silico* serotyping of genome assemblies.
    [`Convert`][kaptive.serotyping.cli.Convert]: Subcommand for converting serialized JSON-lines results.
"""

import argparse
from typing import Any

from kaptive.cli import Cli, Colors, Command


class ResultExporter:
    r"""Evaluates CLI arguments once and sets up a pipeline of output writers.

    This eliminates conditional branching inside the processing loop and allows
    reusability between the [`Type`][kaptive.serotyping.cli.Type] and
    [`Convert`][kaptive.serotyping.cli.Convert] commands.

    Attributes:
        file_suffix (str): Default filename suffix used when writing output files (default: `'kaptive_results'`).
        writers (list[Callable[[SerotypingResult], None]]): List of registered writer callback functions.
    """

    file_suffix = "kaptive_results"

    def __init__(self, cli: Cli, args: argparse.Namespace) -> None:
        r"""Initialize output writers based on parsed command-line arguments.

        Inspects output flags in `args` and registers appropriate serialization callbacks
        for TSV, PHA4GE TSV, JSON, locus nucleotide FASTA, gene nucleotide FASTA,
        translated protein FASTA, and interactive HTML plots.

        The PHA4GE TSV output adheres to the Public Health Alliance for Genomic Epidemiology
        genotyping specification (https://github.com/pha4ge/genotyping-specification).

        Args:
            cli (Cli): Parent `Cli` execution context.
            args (argparse.Namespace): Parsed command-line arguments containing output flags.

        Raises:
            SystemExit: If `--json` is set but `orjson` is not installed, or `--plots`
                is set but `plotly` is not installed.
        """
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
                from orjson import OPT_APPEND_NEWLINE, OPT_SERIALIZE_NUMPY, dumps
            except ImportError:
                cli.exit("orjson not installed. Please run: pip install kaptive[json]")
            json_handle = cli.open_file(json_file, mode="wb")
            self.writers.append(
                lambda r: json_handle.write(dumps(r.to_dict(), option=OPT_SERIALIZE_NUMPY | OPT_APPEND_NEWLINE))
            )

        if loci_dir := getattr(args, "loci", None):
            self.writers.append(
                lambda r: (loci_dir / f"{r.genome}_{self.file_suffix}.fna").write_bytes(r.locus_seqs.to_fasta())
            )

        if genes_dir := getattr(args, "genes", None):
            self.writers.append(
                lambda r: (genes_dir / f"{r.genome}_{self.file_suffix}.ffn").write_bytes(r.gene_seqs.to_fasta())
            )

        if proteins_dir := getattr(args, "proteins", None):
            self.writers.append(
                lambda r: (proteins_dir / f"{r.genome}_{self.file_suffix}.faa").write_bytes(r.translations.to_fasta())
            )

        if plot_dir := getattr(args, "plots", None):
            try:
                from kaptive.plotting import SerotypingResultPlotter
            except ImportError:
                cli.exit("plotly not installed. Please run: pip install kaptive[plot]")
            self.writers.append(
                lambda r: SerotypingResultPlotter()(r).write_html(
                    plot_dir / f"{r.genome}_{self.file_suffix}.html",
                    include_plotlyjs="cdn",
                    full_html=True,
                )
            )

    def __call__(self, result: Any) -> None:
        r"""Pass the serotyping result to all registered output writers.

        Args:
            result (SerotypingResult): The [`SerotypingResult`][kaptive.serotyping.models.SerotypingResult]
                instance to serialize and write out.
        """
        for write in self.writers:
            write(result)


# Type command ---------------------------------------------------------------------------------------------------------
class Type(Command):
    r"""💉 In silico serotyping of genome assemblies

    Aliases:
        assembly
    """

    aliases = ["assembly"]  # Backwards compatibility with < v3.3

    def setup_arguments(self) -> None:
        r"""Configure argument parser options for the type subcommand.

        Defines input database/genome arguments, output formatting flags via
        `add_output_arguments`, confidence options,
        and thread count parameters.
        """
        opts = self.parser.add_argument_group(Colors.wrap("📥 Inputs", Colors.BOLD))
        opts.add_argument("database", help="Database path or keyword (see: `kaptive db list`)")
        opts.add_argument(
            "genomes",
            nargs="+",
            help="Genome assemblies in fasta format; can be compressed",
        )

        opts = self.parser.add_argument_group(Colors.wrap("📤 Outputs", Colors.BOLD))
        self.add_output_arguments(opts, tsv_flags=("-o", "--out"), include_json=True)

        opts = self.parser.add_argument_group(Colors.wrap("🔬 Confidence options", Colors.BOLD))
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

        opts = self.parser.add_argument_group(Colors.wrap("🔧 Other options", Colors.BOLD))
        opts.add_argument(
            "-t",
            "--threads",
            type=int,
            default=0,
            metavar="",
            help="Number threads or 0 for all available (default: 0)",
        )
        opts.add_argument(
            "--partial-edge-tolerance",
            type=int,
            default=5,
            metavar="",
            help="Tolerance in bases from contig edge to call a partial gene (default: %(default)s)",
        )

    def __call__(self, args: argparse.Namespace) -> None:
        r"""Execute the serotyping workflow on input genome assemblies.

        Loads the requested locus database, initializes the serotyping engine, iterates
        through genome assemblies to call serotypes, and streams results to configured output targets.

        Args:
            args (argparse.Namespace): Parsed command-line arguments.
        """
        self.cli.msg(f"💽 Loading database {args.database}...")
        from kaptive.db import DatabaseManager
        from kaptive.serotyping import Serotyper

        db = DatabaseManager.get(args.database)
        exporter = ResultExporter(self.cli, args)

        serotyper = Serotyper(
            db=db,
            max_other_genes=args.max_other_genes,
            min_completeness=args.min_completeness,
            allow_below_threshold=args.below_threshold,
            partial_edge_tolerance=args.partial_edge_tolerance,
        )
        for genome in self.cli.progress(args.genomes, "💉 Serotyping genomes..."):
            if result := serotyper(genome):
                exporter(result)

        self.cli.msg(f"✅ Serotyping complete. Results written to '{args.out}'.")


# Client command -------------------------------------------------------------------------------------------------------
class Convert(Command):
    r"""🔄 Convert serialized Kaptive results into different formats

    Reads serialized JSON-lines serotyping output records and converts them into tabular
    TSV, PHA4GE TSV, or sequence FASTA files without re-running the serotyping pipeline.
    """

    def setup_arguments(self) -> None:
        r"""Configure argument parser options for the convert subcommand.

        Defines the input JSON-lines source parameter and output target flags via
        `add_output_arguments`.
        """
        opts = self.parser.add_argument_group(Colors.wrap("📥 Inputs", Colors.BOLD))
        opts.add_argument(
            "jsonl",
            default="stdin",
            help="Serialised results in JSON-lines format (default: stdin)",
        )

        opts = self.parser.add_argument_group(Colors.wrap("📤 Outputs", Colors.BOLD))
        self.add_output_arguments(opts, tsv_flags=("-t", "--tsv"), include_json=False)

    def __call__(self, args: argparse.Namespace) -> None:
        r"""Execute result format conversion from serialized JSON-lines input.

        Deserializes line-delimited JSON records into [`SerotypingResult`][kaptive.serotyping.models.SerotypingResult]
        objects and dispatches them to registered output writers.

        Args:
            args (argparse.Namespace): Parsed command-line arguments.

        Raises:
            SystemExit: If `orjson` is not installed in the current Python environment.
        """
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
