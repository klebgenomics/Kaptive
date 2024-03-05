"""
This module contains classes for interacting with bacterial genome assemblies and contigs and a pipeline
to type them.

Copyright 2023 Tom Stanton (tomdstanton@gmail.com)
https://github.com/klebgenomics/Kaptive

This file is part of Kaptive. Kaptive is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by the Free Software Foundation,
either version 3 of the License, or (at your option) any later version. Kaptive is distributed
in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
details. You should have received a copy of the GNU General Public License along with Kaptive.
If not, see <https://www.gnu.org/licenses/>.
"""
from __future__ import annotations

from pathlib import Path
from statistics import stdev, mean
from subprocess import Popen, PIPE, DEVNULL

from Bio.Seq import Seq
from matplotlib import pyplot as plt

from kaptive.log import warning, log, quit_with_error
from kaptive.database import Locus
from kaptive.alignment import Alignment
from kaptive.typing import AssemblyTypingResult
from kaptive.misc import parse_fasta


# Classes -------------------------------------------------------------------------------------------------------------
class AssemblyError(Exception):
    pass


class Assembly:
    def __init__(self, path: Path | None = None, name: str | None = None, contigs: dict[str: Contig] | None = None):
        self.path = path or Path()
        self.name = name or ''
        self.contigs = contigs or {}

    @classmethod
    def from_path(cls, path: Path, allow_plasmids: bool = True) -> Assembly:
        self = cls(path=path, name=path.name.strip('.gz').rsplit('.', 1)[0])
        for n, d, s in parse_fasta(path):
            self.contigs[c.name] = (c := Contig(self, n, d, Seq(s)))
        if not self.contigs:
            raise AssemblyError(f"No contigs found in {path}")
        return self

    def __repr__(self):
        return self.name

    def __hash__(self):
        return hash(self.name)

    def __len__(self):
        return self.path.stat().st_size

    def __iter__(self) -> tuple[str, Contig]:
        for contig_name, contig in self.contigs.items():
            yield contig_name, contig

    def __getitem__(self, contig_name):
        return self.contigs[contig_name]

    def __contains__(self, contig_name):
        return contig_name in self.contigs


class ContigError(Exception):
    pass


class Contig(object):
    """
    This class describes a contig in an assembly: the name, length, and sequence.
    """

    def __init__(self, assembly: Assembly | None = None, name: str | None = None, description: str | None = None,
                 sequence: Seq | None = Seq(''), replicon: str | None = None, topology: str | None = None):
        self.assembly = assembly
        self.name = name or ''
        self.description = description or ''
        self.sequence = sequence
        self.replicon = replicon or ''
        self.topology = topology or ''
        if self.description:
            self._get_replicon()
            self._get_topology()

    def __repr__(self):
        return self.name

    def __len__(self):
        return len(self.sequence)

    def __hash__(self):
        return hash(self.name)

    def _get_replicon(self):
        if any(i in self.name + self.description for i in {'plasmid', '__pl'}):
            self.replicon = "plasmid"
        elif any(i in self.name + self.description for i in {'chromosome', '__ch'}):
            self.replicon = "chromosome"
        else:
            self.replicon = "unknown"

    def _get_topology(self):
        if any(i in self.name + self.description for i in {'circular=true', 'complete', 'circular=Y'}):
            self.topology = "circular"
        elif any(i in self.name + self.description for i in {'circular=false', 'circular=N', 'linear'}):
            self.topology = "linear"
        else:
            self.topology = "unknown"


# Functions -----------------------------------------------------------------------------------------------------------
def get_scores(results: list[tuple[Locus: list[Alignment]]], score: str, weight: str,
               reverse_sort: bool = True) -> list[tuple[Locus, float]]:
    if weight == 'none':
        results = [(l, sum(getattr(x, score) for x in r)) for l, r in results]
    elif weight == "locus_length":
        results = [(l, sum(getattr(x, score) for x in r) / len(l)) for l, r in results]
    elif weight == "genes_expected":
        results = [(l, sum(getattr(x, score) for x in r) / len(l.genes)) for l, r in results]
    elif weight == "genes_found":
        results = [(l, sum(getattr(x, score) for x in r) / len(r)) for l, r in results]
    elif weight == "prop_genes_found":
        results = [(l, sum(getattr(x, score) for x in r) * (len(r) / len(l.genes))) for l, r in results]
    else:
        quit_with_error(f"Invalid weight {weight}")
    return sorted(results, key=lambda x: x[1], reverse=reverse_sort)


def score_stats(scores: list[tuple[Locus, float]]) -> list[tuple[Locus, float, float]]:
    """Returns locus, score and zscore"""
    _scores = {i[1] for i in scores}  # Get unique scores
    if len(_scores) == 1:
        return [(scores[0][0], scores[0][1], 0)]  # If only one score, return it with a zscore of 0
    m, sd = mean(_scores), stdev(_scores)  # Calculate mean and standard deviation
    return [(scores[0][0], scores[0][1], 0)] if sd == 0 else [  # If standard deviation is 0, return zscore of 0
        (scores[i][0], scores[i][1], (scores[i][1] - m) / sd) for i in range(len(scores))]


def typing_pipeline(assembly: Path, args):
    """
    Typing function for assemblies. Aligns genes to the assembly and scores the results.
    """
    try:
        assembly = Assembly.from_path(assembly)
        # assembly = Assembly.from_path(args.input[0])
    except AssemblyError as e:
        warning(f'Error parsing {assembly}: {e}')
        return
    log(f'Typing {assembly}', args.verbose)
    minimap_cmd = f"minimap2 -c -t {args.threads} {assembly.path} -"
    minimap_cmd += f" -x {args.preset}" if args.preset else ""
    minimap_cmd += f" {args.args}" if args.args else ""
    gene_alignments = {}
    extra_gene_alignments = {}  # Init dict of extra gene alignments {gene: alignment}
    is_element_alignments = []

    with Popen(minimap_cmd.split(), stdin=PIPE, stdout=PIPE, stderr=DEVNULL) as proc:
        for a in proc.communicate(args.db.as_gene_fasta().encode())[0].splitlines():  # Write gene fasta to proc.stdin
            if not (a := Alignment.from_paf_line(a)).target_name in assembly.contigs:
                continue  # Skip alignment if it's not to a contig in the assembly, happens if plasmids aren't allowed
            if a.query_name in args.db.is_elements:  # Iterate over lines in proc.stdout
                is_element_alignments.append(a)  # If the alignment is to an IS element, add to list
            elif a.query_name in args.db.genes:  # If the alignment is to a gene

                if (x := gene_alignments.get(a.query_name, None)) and a.AS < x.AS:
                    continue  # Skip this alignment if it's worse than the one we already have
                gene_alignments[a.query_name] = a
            elif a.query_name in args.db.extra_genes:  # If the alignment is to an extra gene
                if (x := extra_gene_alignments.get(a.query_name, None)) and a.AS < x.AS:
                    continue
                extra_gene_alignments[a.query_name] = a
            else:
                warning(f'Alignment to unknown gene {a.query_name} in {assembly}')

    if not gene_alignments:  # If no alignments, skip this assembly
        warning(f'No gene alignments found for {assembly}')
        args.out.write(f"{assembly}: No gene alignments found\n")
        return

    for_scoring = [(l, y) for l in args.db.loci.values() if (y := [x for g in l.genes.values() if (
        x := gene_alignments.get(g.name)) and x.percent_query_coverage >= args.min_cov])]

    scores = {  # Get scores for each metric and weight
        metric: {
            weight: score_stats(get_scores(for_scoring, metric, weight)) for weight in {
                'none', 'locus_length', 'genes_expected', 'genes_found', 'prop_genes_found', args.weight}
        } for metric in {'percent_identity', 'matching_bases', 'AS', 'percent_query_coverage', args.score}
    }

    result = AssemblyTypingResult.from_alignments(assembly, gene_alignments, extra_gene_alignments,
                                                  is_element_alignments, scores, args)

    args.out.write(result.as_table(args.debug))  # type: 'TextIOWrapper'
    if args.fasta:  # type: 'Path'
        (args.fasta / f'{result.assembly.name}_kaptive_results.fna').write_text(result.as_fasta())
    if args.json:  # type: 'TextIOWrapper'
        args.json.write(result.as_json())
    if args.figures:  # type: 'Path'
        ax = result.as_GraphicRecord().plot()[0]
        # Add title
        ax.set_title(f"{result.assembly.name} {result.best_match} ({result.phenotype}) - {result.confidence}")
        # Save figure
        ax.figure.savefig(args.figures / f'{result.assembly.name}_kaptive_results.png')
        plt.close()
