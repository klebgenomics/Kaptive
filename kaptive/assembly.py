"""
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
from gzip import open as gzopen
from statistics import stdev, mean
from functools import cached_property
from subprocess import Popen, PIPE
# import fcntl

from Bio import SeqIO
from Bio.Seq import Seq

from kaptive.log import warning, log
from kaptive.database import Locus
from kaptive.alignment import Alignment
from kaptive.misc import check_file


# Classes -------------------------------------------------------------------------------------------------------------
class AssemblyError(Exception):
    pass


class Assembly:
    def __init__(self, path: Path | None = None, name: str | None = None, contigs: dict[str: Contig] | None = None):
        self.path = path or Path()
        self.name = name or ''
        self.contigs = contigs or {}

    @classmethod
    def from_path(cls, path: str):
        path = check_file(path)
        return cls(path=path, name=path.name.strip('.gz').rsplit('.', 1)[0])

    def __repr__(self):
        return self.name

    def __hash__(self):
        return hash(self.name)

    def __len__(self):
        return self.path.stat().st_size

    def __iter__(self):
        for contig_name, contig in self.contigs.items():
            yield contig_name, contig

    def __getitem__(self, contig_name):
        if contig_name not in self.contigs:
            raise KeyError(f'Contig {contig_name} not found in assembly {self.name}')
        return self.contigs[contig_name]

    def __contains__(self, contig_name):
        return contig_name in self.contigs

    def open(self, mode: str = 'rt'):
        if not self.path.suffix == '.gz':
            try:
                return open(self.path, mode)
            except Exception as e:
                raise AssemblyError(f'Could not open assembly {self}: {e}')
        else:
            try:
                return gzopen(self.path, mode)
            except Exception as e:
                raise AssemblyError(f'Could not open assembly {self}: {e}')

    def parse(self):
        with self.open() as handle:
            try:
                for name, seq in SeqIO.FastaIO.SimpleFastaParser(handle):
                    name, description = name.split(' ', 1) if ' ' in name else (name, '')
                    self.add_contig(Contig(name, description, Seq(seq), self))
            except Exception as e:
                raise AssemblyError(f'Could not parse FASTA file {self}: {e}')

    def add_contig(self, contig: Contig):
        """Adds a contig to the assembly."""
        if contig.name in self.contigs:
            raise AssemblyError(f'Contig {contig.name} already exists in assembly {self.name}')
        if contig.assembly is not None and contig.assembly != self:
            raise AssemblyError(f'Contig {contig.name} already belongs to assembly {contig.assembly.name}')
        self.contigs[contig.name] = contig
        contig.assembly = self


class ContigError(Exception):
    pass


class Contig(object):
    """
    This class describes a contig in an assembly: the name, length, and sequence.
    """

    def __init__(self, name: str | None = None, description: str | None = None, sequence: Seq | None = Seq(''),
                 assembly: Assembly | None = None):
        self.name = name or ''
        self.description = description or ''
        self.sequence = sequence
        self.assembly = assembly

    def __repr__(self):
        return self.name

    def __len__(self):
        return len(self.sequence)

    def __hash__(self):
        return hash(self.name)

    @cached_property
    def replicon(self):
        if 'plasmid' in self.name + self.description or '__pl' in self.name + self.description:
            return "plasmid"
        elif 'chromosome' in self.name + self.description or '__ch' in self.name + self.description:
            return "chromosome"
        else:
            return "unknown"

    @cached_property
    def topology(self):
        if 'circular=true' in self.description or 'circular=Y' in self.description or 'complete' in self.description:
            return "circular"
        elif 'circular=false' in self.description or 'circular=N' in self.description or 'linear' in self.description:
            return "linear"
        else:
            return "unknown"


# Functions -----------------------------------------------------------------------------------------------------------
def get_scores(results: dict[Locus: list[Alignment]], score: str, weight: str | None, reverse_sort: bool = True) -> \
        list[tuple[Locus, float]]:
    if not weight:
        results = [(l, sum(getattr(x, score) for x in r)) for l, r in results.items()]  # Sum scores for each locus
    elif weight == "locus_length":
        results = [(l, sum(getattr(x, score) for x in r) / len(l)) for l, r in results.items()]
    elif weight == "genes_expected":
        results = [(l, sum(getattr(x, score) for x in r) / len(l.genes)) for l, r in results.items()]
    elif weight == "genes_found":
        results = [(l, sum(getattr(x, score) for x in r) / len(r)) for l, r in results.items()]
    else:
        results = [(l, sum(getattr(x, score) / getattr(x, weight) for x in r)) for l, r in results.items()]
    return sorted(results, key=lambda x: x[1], reverse=reverse_sort)


def score_stats(scores: list[tuple[Locus, float]]) -> list[tuple[Locus, float, float], None, None]:
    """Returns locus, score and zscore"""
    if len(scores) == 1:
        return [(scores[0][0], scores[0][1], 0)]  # If only one score, return it with a zscore of 0
    m, sd = mean(s := [i[1] for i in scores]), stdev(s)  # Calculate mean and standard deviation
    return [(scores[0][0], scores[0][1], 0)] if sd == 0 else [  # If standard deviation is 0, return zscore of 0
        (scores[i][0], scores[i][1], (scores[i][1] - m) / sd) for i in range(len(scores))]


def type_assembly(assembly, args):
    try:
        assembly.parse()
    except AssemblyError as e:
        warning(f'Error parsing {assembly}: {e}')
        return
    log(f'Typing:\t{assembly}', args.verbose)

    minimap_cmd = f"minimap2 -c -t {args.threads} {assembly.path} -"
    minimap_cmd += f" -x {args.preset}" if args.preset else ""
    minimap_cmd += f" {args.args}" if args.args else ""

    alignments = {}  # Init dict of best gene alignments {gene: alignment}
    is_element_alignments = []  # Init list of IS element alignments
    with Popen(minimap_cmd.split(), stdin=PIPE, stdout=PIPE, stderr=PIPE) as proc:
        for a in proc.communicate(args.db.as_gene_fasta().encode())[0].splitlines():  # Get alignments from minimap2
            if (a := Alignment.from_paf_line(a)).query_name in args.db.is_elements:
                is_element_alignments.append(a)   # If the alignment is to an IS element, add to list
            elif a.query_name in args.db.genes:  # If the alignment is to a gene
                if a.query_name in alignments and a.alignment_score < alignments[a.query_name].alignment_score:
                    continue  # Skip this alignment if it's worse than the one we already have
                alignments[a.query_name] = a
            else:
                warning(f'Alignment to unknown gene {a.query_name} in {assembly}')

    alignments_per_locus = {}  # Init per-locus gene alignment dict {locus: {gene: alignment}}
    not_used_for_scoring = {}  # Used for extra genes and anything else not used for scoring
    for a in alignments.values():
        locus = args.db.genes[a.query_name].locus
        alignment_dict = alignments_per_locus if not locus.name.startswith(
            'Extra') and a.percent_query_coverage >= args.min_cov else not_used_for_scoring
        if locus not in alignment_dict:
            alignment_dict[locus] = {}
        alignment_dict[locus][a.query_name] = a

    if not alignments_per_locus:  # If no alignments, skip this assembly
        warning(f'No gene alignments found for {assembly}')
        return

    # Convert to dict of lists for get_locus_scores, will add back in the alignments not used for scoring later
    locus_results = {k: list(v.values()) for k, v in alignments_per_locus.items()}
    scores = {  # Get scores for each metric and weight
        metric: {
            weight: score_stats(get_scores(locus_results, metric, weight)) for weight in
            {None, 'locus_length', 'genes_expected', 'genes_found', args.weight}
        } for metric in {'percent_identity', 'matching_bases', 'alignment_score', 'percent_query_coverage', args.score}
    }

    # Add back in the alignments not used for scoring
    for locus, gene_alignments in not_used_for_scoring.items():
        if locus not in locus_results:
            locus_results[locus] = list(gene_alignments.values())
        else:
            locus_results[locus].extend(list(gene_alignments.values()))

    # Get the best locus and its alignments
    best_locus, best_score, best_zscore = scores[args.score][args.weight][0]  # Best locus is the first in the list
    from kaptive.result import AssemblyResult  # Import here to avoid circular import
    result = AssemblyResult.from_best_locus(  # Create an AssemblyResult object
        args.db, best_locus, assembly,
        locus_results[best_locus],  # Best locus alignments,
        [x for locus in locus_results if locus != best_locus for x in locus_results[locus]],  # Other locus alignments
        is_element_alignments, args.gene_threshold, args.gene_distance
    )

    args.out.write(
        '\t'.join(result.as_list() + [
            ';'.join([i.gene.name for i in result.truncated_genes]),  # Truncated genes
            ';'.join([str(i) for i in result.extra_genes]),           # Extra genes
            str(len({i.contig.name for i in result.contig_pieces})),  # Number of contigs
            str(len(result.contig_pieces)),                           # Number of contig pieces
            ';'.join(str(i) for i in result.contig_pieces),           # Contig pieces
            f"{best_score:.1f}",                                      # Best score
            f"{best_zscore:.1f}",                                     # Best zscore
            args.score,                                               # Score metric
            args.weight if args.weight else 'None',                   # Weight metric
            args.preset if args.preset else 'None',                   # Minimap2 preset
            args.args if args.args else 'None',                       # Minimap2 args
            str(args.gene_threshold),                                 # Gene threshold
            str(args.gene_distance),                                  # Gene distance
            str(args.min_cov),                                        # Min coverage
            ';'.join(f"{i.query_name},{i.percent_identity:.2f}%,{i.percent_query_coverage:.2f}%" for i in
                     not_used_for_scoring.values() for i in i.values()),  # Not used for scoring
            ';'.join(
                f"{s}_{w}:{'|'.join(f'{l},{x:.4f},{y:.4f}' for l, x, y, in scores[s][w])}" for s in scores for w in scores[s]
            )  # Scores
        ] if args.debug else result.as_list()) + '\n'
    )

    if args.fasta:
        (args.fasta / f'{assembly.name}_kaptive_results.fna').write_text(result.as_fasta(), encoding='utf-8', mode='wt')

    if args.json:  # Use a JSONLines file to append and read from in chunks, may still need fcntl.flock for cluster
        args.json.write(result.as_json())
