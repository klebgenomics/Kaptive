r"""Core engine for in silico serotyping of bacterial genome assemblies.

This module implements the primary serotyping workflow in the [`Serotyper`][kaptive.serotyping.core.Serotyper]
class. It performs sequence alignment of reference locus genes against target genome assemblies using `rammappy`,
scores candidate loci, reconstructs spatial locus pieces ([`LocusPieces`][kaptive.serotyping.models.LocusPieces]),
classifies gene hits ([`GeneHits`][kaptive.serotyping.models.GeneHits] and
[`GeneState`][kaptive.serotyping.models.GeneState]), translates nucleotide alignments to amino acid sequences,
evaluates phenotype modifications, and produces comprehensive
[`SerotypingResult`][kaptive.serotyping.models.SerotypingResult] objects.
"""

from pathlib import Path

import numpy as np
import rammappy
from rammappy.align import Aligner

from kaptive import __version__
from kaptive.core.alignment import Alignments
from kaptive.core.genome import GenomeAssembly
from kaptive.core.pairwise import PairwiseAligner
from kaptive.core.seq import Sequences
from kaptive.db import Database
from kaptive.serotyping.models import (
    GeneHits,
    GeneState,
    LocusPieces,
    SerotypingResult,
)


class Serotyper:
    r"""High-performance *in silico* serotyping engine for bacterial genome assemblies.

    The [`Serotyper`][kaptive.serotyping.core.Serotyper] utilizes a reference database
    ([`Database`][kaptive.db.Database]) containing surface antigen locus definitions, reference gene sequences,
    and phenotypic rules to evaluate input assemblies ([`GenomeAssembly`][kaptive.core.genome.GenomeAssembly]).

    It executes a four-phase serotyping pipeline:

    1. **Mapping & Scoring**: Maps reference genes to assembly contigs using `rammappy`, culls overlapping hits,
       and ranks candidate loci based on gene coverage and locus completeness.
    2. **Locus Reconstruction**: Clusters gene hits spatially to bound locus regions into
       [`LocusPieces`][kaptive.serotyping.models.LocusPieces] and identifies missing or unexpected genes inside/outside
       locus boundaries.
    3. **Gene State & Identity Evaluation**: Translates gene alignments, performs protein-level pairwise alignment
       with [`PairwiseAligner`][kaptive.core.pairwise.PairwiseAligner], assesses frame shifts and truncations, and
       assigns [`GeneState`][kaptive.serotyping.models.GeneState] (NORMAL, PARTIAL, TRUNCATED, NOVEL).
    4. **Phenotype & Confidence Scoring**: Applies phenotypic rules (e.g. active/inactive gene clusters) and
       determines overall serotype typeability.

    Attributes:
        max_other_genes (int): Maximum allowed unexpected genes inside locus before classifying sample as untypeable.
        min_completeness (float): Minimum locus completeness fraction required for typeability call.
        allow_below_threshold (bool): Whether to permit genes falling below identity threshold while remaining typeable.
        preset (rammappy.Preset | None): Custom `rammappy` alignment preset, if specified.
        scoring_metric (str): Scoring metric used for locus scoring (default `"scores"`).
        min_gene_coverage (float): Minimum query coverage fraction required for gene alignments to be considered valid.
        partial_edge_tolerance (int): Base pair distance tolerance from contig boundaries for identifying partial genes.
    """

    def __init__(
        self,
        db: Database,
        max_other_genes: int = 1,
        min_completeness: float = 0.5,
        allow_below_threshold: bool = False,
        preset: rammappy.Preset | None = None,
        scoring_metric: str = "scores",
        min_gene_coverage: float = 0.20,
        partial_edge_tolerance: int = 5,
    ) -> None:
        r"""Initialize the serotyping engine with a reference database and search parameters.

        Args:
            db (Database): The reference surface antigen database containing loci, genes, and phenotype definitions.
                See [`Database`][kaptive.db.Database].
            max_other_genes (int): Maximum allowed unexpected genes inside the locus boundary before flagging as
                untypeable. Defaults to `1`.
            min_completeness (float): Minimum proportion of expected locus genes required to consider the call
                typeable. Defaults to `0.5`.
            allow_below_threshold (bool): If `False`, any gene inside the locus falling below the identity threshold
                makes the result untypeable. Defaults to `False`.
            preset (rammappy.Preset | None): Optional `rammappy` mapping preset. Defaults to `None`.
            scoring_metric (str): Scoring metric used for candidate locus ranking. Defaults to `"scores"`.
            min_gene_coverage (float): Minimum gene alignment query coverage fraction (0.0 to 1.0) for valid scoring.
                Defaults to `0.20`.
            partial_edge_tolerance (int): Distance tolerance in base pairs from contig edges to classify a hit
                as partial. Defaults to `5`.
        """
        self._db: Database = db
        self.max_other_genes: int = max_other_genes
        self.min_completeness: float = min_completeness
        self.allow_below_threshold: bool = allow_below_threshold
        self.preset: rammappy.Preset | None = preset
        self.scoring_metric: str = scoring_metric
        self.min_gene_coverage: float = min_gene_coverage
        self.partial_edge_tolerance: int = partial_edge_tolerance
        self._protein_aligner = PairwiseAligner()

        # Count expected genes per locus for weighting
        self._expected_genes_per_locus = np.zeros(len(self._db.loci), dtype=np.float32)
        np.add.at(
            self._expected_genes_per_locus,
            self._db.gene_locus_indices[~self._db.extra_genes],
            1.0,
        )
        self._expected_genes_per_locus = np.maximum(self._expected_genes_per_locus, 1.0)

        # Prepare metadata for alignment iterator parsing
        self._gene_seqs = [
            (
                str(i).encode(),
                bytes(
                    self._db.genes.seqs[
                        self._db.genes.offsets[i] : self._db.genes.offsets[i] + self._db.genes.lengths[i]
                    ]
                ),
            )
            for i in range(len(self._db.genes))
        ]
        self._gene_meta = [(str(i), self._db.genes.lengths[i]) for i in range(len(self._db.genes))]

    def __call__(self, genome: GenomeAssembly | str | Path) -> SerotypingResult | None:
        r"""Perform *in silico* serotyping on a target bacterial genome assembly.

        Maps reference locus genes against the provided genome assembly, ranks candidate loci,
        reconstructs locus boundaries, evaluates gene integrity and amino acid identity, and resolves the predicted
        serotype phenotype.

        Args:
            genome (GenomeAssembly | str | Path): Target genome assembly as a
                [`GenomeAssembly`][kaptive.core.genome.GenomeAssembly] instance or filesystem path (`str` or `Path`)
                to a FASTA file.

        Returns:
            SerotypingResult | None: Complete serotyping analysis result containing best matching locus, predicted
                phenotype, gene hit classifications, spatial locus pieces, and confidence metrics.
                See [`SerotypingResult`][kaptive.serotyping.models.SerotypingResult].

        Raises:
            FileNotFoundError: If `genome` is passed as a file path that does not exist on disk.
            ValueError: If the genome assembly contains no valid contigs or sequence data cannot be parsed.
        """
        genome = GenomeAssembly.ensure(genome)

        contig_index = genome.get_rammappy_index()
        aligner = Aligner(index=contig_index, preset=None, do_cigar=True, do_cs=False, do_md=False)
        opts = aligner.options
        opts.filtering.best_n = 50000
        opts.filtering.pri_ratio = 0.0
        aligner.options = opts

        gene_aln_batch = aligner.map_batch(self._gene_seqs)
        gene_alns = Alignments.from_mapping_iterators(self._gene_meta, gene_aln_batch)

        # Calculate total coverage per gene across all alignments for reporting
        q_indices = gene_alns.q_names.astype(np.int32)
        q_lengths = gene_alns.q_aln_lens
        total_q_covs = np.zeros(len(self._db.genes), dtype=np.float32)
        np.add.at(total_q_covs, q_indices, q_lengths)
        total_q_covs /= self._db.genes.lengths

        # Scoring phase ------------------------------------------------------------------------------------------------
        # Calculate alignment query coverage per alignment
        q_covs = gene_alns.q_covs
        valid_cov_mask = q_covs >= self.min_gene_coverage

        valid_alns = gene_alns[valid_cov_mask]
        valid_q_covs = q_covs[valid_cov_mask]
        valid_gene_indices = valid_alns.q_names.astype(np.int32)  # type: ignore

        # Sort to find the best alignment per gene (highest q_cov, tie-breaker: highest score)
        order = np.lexsort((-valid_alns.scores, -valid_q_covs, valid_gene_indices))  # type: ignore
        valid_alns = valid_alns[order]
        valid_gene_indices = valid_gene_indices[order]
        valid_q_covs = valid_q_covs[order]

        # Take the first (best) alignment for each gene
        _, unique_indices = np.unique(valid_gene_indices, return_index=True)
        best_gene_indices = valid_gene_indices[unique_indices]
        best_q_covs = valid_q_covs[unique_indices]

        valid_locus_indices = self._db.gene_locus_indices[best_gene_indices]
        valid_not_extra = ~self._db.extra_genes[best_gene_indices]

        # Calculate locus scores by summing best q_covs of valid expected genes
        locus_scores = np.zeros(len(self._db.loci), dtype=np.float64)
        np.add.at(
            locus_scores,
            valid_locus_indices[valid_not_extra],
            best_q_covs[valid_not_extra],
        )

        # Calculate completeness count per locus (unique expected genes matched per locus)
        locus_counts = np.zeros(len(self._db.loci), dtype=np.float32)
        matched_expected_genes = best_gene_indices[valid_not_extra]
        np.add.at(locus_counts, self._db.gene_locus_indices[matched_expected_genes], 1.0)

        locus_completeness = locus_counts / self._expected_genes_per_locus
        final_locus_scores = locus_scores * (locus_completeness**3)

        self._last_scores = final_locus_scores.copy()
        self._last_completeness = locus_completeness.copy()

        best_locus_idx = int(np.argmax(final_locus_scores))
        best_locus_name = self._db.loci.ids[best_locus_idx]

        # Reconstruction phase -----------------------------------------------------------------------------------------
        valid_alns = gene_alns

        # Cull alignments, prioritizing genes belonging to the best match locus
        valid_indices = valid_alns.q_names.astype(np.int32)
        priority_mask = self._db.gene_locus_indices[valid_indices] == best_locus_idx

        culled_alns = valid_alns.cull_overlaps(by_query=False, priority_mask=priority_mask, max_overlap_fraction=0.1)

        # Re-extract arrays for the culled batch
        culled_gene_indices = culled_alns.q_names.astype(np.int32)
        t_indices = np.array([genome.id_map[n] for n in culled_alns.t_names], dtype=np.uint32)
        # Cluster intervals by contig using max_locus_length as tolerance
        culled_intervals = culled_alns.to_intervals(by_query=False)
        piece_ids = culled_intervals.cluster_spatial(tolerance=self._db.max_locus_length, group_by=t_indices)

        # Identify expected genes and the clusters they fall into
        is_expected = (self._db.gene_locus_indices[culled_gene_indices] == best_locus_idx) & ~self._db.extra_genes[
            culled_gene_indices
        ]
        valid_cluster_ids = np.unique(piece_ids[is_expected])
        is_extra = self._db.extra_genes[culled_gene_indices]

        # Calculate gene alignment coverages based on total coverage across all non-overlapping fragments
        coverages = np.clip(total_q_covs[culled_gene_indices] * 100.0, 0.0, 100.0)

        # Identify the primary hit for each expected gene for bounding box calculation
        primary_expected = np.zeros(len(culled_alns), dtype=bool)
        is_expected_hits = np.where(is_expected)[0]
        if len(is_expected_hits) > 0:
            exp_gene_indices = culled_gene_indices[is_expected_hits]
            exp_scores = culled_alns.scores[is_expected_hits]
            order = np.lexsort((-exp_scores, exp_gene_indices))
            sorted_exp_gene_indices = exp_gene_indices[order]
            _, unique_indices = np.unique(sorted_exp_gene_indices, return_index=True)
            best_hits = is_expected_hits[order[unique_indices]]
            primary_expected[best_hits] = True

        # Construct bounding locus pieces from valid pieces
        l_ctg_indices, l_starts, l_ends, l_strands = [], [], [], []
        l_expected_means = []
        for c_id in valid_cluster_ids:
            piece_mask = piece_ids == c_id

            piece_primary = piece_mask & primary_expected
            if np.any(piece_primary):
                ctg_idx = t_indices[piece_mask][0]
                l_ctg_indices.append(ctg_idx)
                start = np.min(culled_intervals.starts[piece_primary])
                end = np.max(culled_intervals.ends[piece_primary])
                l_starts.append(start)
                l_ends.append(end)

                exp_genes = culled_gene_indices[piece_primary]
                l_expected_means.append(np.mean(self._db.gene_positions[exp_genes]))

                exp_strands = self._db.gene_intervals.strands[exp_genes]
                found_strands = culled_alns.strands[piece_primary]
                if np.sum(found_strands * exp_strands) < 0:
                    l_strands.append(-1)
                else:
                    l_strands.append(1)

        # Recompute is_inside using the strict expected boundaries
        is_inside = np.zeros(len(culled_alns), dtype=bool)
        for ctg_idx, start, end in zip(l_ctg_indices, l_starts, l_ends):
            on_ctg = t_indices == ctg_idx
            # An alignment is inside if it overlaps the bounding region
            inside_this = on_ctg & (culled_intervals.starts <= end) & (culled_intervals.ends >= start)
            is_inside |= inside_this

        # Sort pieces by expected mean position
        piece_order = np.argsort(l_expected_means)

        locus_pieces = LocusPieces(
            ctg_indices=np.array(l_ctg_indices, dtype=np.uint32)[piece_order],
            starts=np.array(l_starts, dtype=np.int32)[piece_order],
            ends=np.array(l_ends, dtype=np.int32)[piece_order],
            strands=np.array(l_strands, dtype=np.int8)[piece_order],
        )

        # Identify missing expected genes
        expected_genes_mask = (self._db.gene_locus_indices == best_locus_idx) & ~self._db.extra_genes
        expected_gene_indices = np.where(expected_genes_mask)[0]
        # Which ones did we find inside the locus?
        found_expected_gene_indices = culled_gene_indices[is_expected & is_inside]
        missing_indices = np.setdiff1d(expected_gene_indices, found_expected_gene_indices, assume_unique=True)
        missing_expected_genes = tuple(self._db.genes.ids[i] for i in missing_indices)

        # Calculate actual completeness based on reconstructed locus
        actual_locus_completeness = (
            1.0 - (len(missing_indices) / len(expected_gene_indices)) if len(expected_gene_indices) > 0 else 1.0
        )

        gene_hits = GeneHits(
            gene_indices=culled_gene_indices,
            q_starts=culled_alns.q_starts,
            q_ends=culled_alns.q_ends,
            t_indices=t_indices,
            t_starts=culled_alns.t_starts,
            t_ends=culled_alns.t_ends,
            strands=culled_alns.strands,
            is_expected=is_expected,
            is_inside=is_inside,
            is_extra=is_extra,
            expected_positions=self._db.gene_positions[culled_gene_indices].astype(np.int32),
            expected_strands=self._db.gene_intervals.strands[culled_gene_indices],
            gene_ids=np.array([self._db.genes.ids[i].encode("utf-8") for i in culled_gene_indices], dtype="S32"),
            cluster_names=np.array(
                [self._db.cluster_keys[self._db.gene_cluster_ids[i]].encode("utf-8") for i in culled_gene_indices],
                dtype="S10",
            ),
            product_descriptions=np.array(
                [
                    self._db.description_keys[self._db.gene_description_ids[i]].encode("utf-8")
                    for i in culled_gene_indices
                ],
                dtype="S64",
            ),
            coverages=coverages,
        )

        # Locus extraction phase ---------------------------------------------------------------------------------------
        if len(locus_pieces) > 0:  # Extract locus sequences using the batched SoA locus pieces
            locus_seqs = genome.contigs.extract(
                locus_pieces.ctg_indices,  # type: ignore
                locus_pieces.starts,
                locus_pieces.ends,
                locus_pieces.strands,
            )
        else:
            locus_seqs = Sequences.empty()

        # Calculate coverage and length discrepancy
        assem_len = np.sum(locus_pieces.ends - locus_pieces.starts)
        ref_len = self._db.loci.lengths[best_locus_idx]
        pcov = float(min(100.0, (assem_len / ref_len) * 100.0)) if ref_len > 0 else 0.0
        if len(locus_pieces) == 1:
            length_discrepancy = float(assem_len - ref_len)
        else:
            length_discrepancy = float("nan")

        # Gene state phase ---------------------------------------------------------------------------------------------
        gene_seqs = genome.contigs.extract_intervals(  # Extract gene nucleotides from their contigs
            gene_hits.t_indices,
            gene_hits.t_intervals,
            new_ids=tuple(self._db.genes.ids[i] for i in gene_hits.gene_indices),
        )
        # Translate nucleotides to amino acids, compensating for the reading frames of the alignments
        # Truncate at the first stop codon to match Old Kaptive's behavior, which prevents frameshifts
        # from pulling down the identity score of the valid upstream alignment.
        prot_seqs = gene_seqs.translate(frames=gene_hits.frames, to_stop=True)  # type: ignore

        # Initialize states
        gene_states = np.full(len(gene_hits), GeneState.NORMAL.value, dtype=np.int8)
        is_partial = culled_alns.is_partial(self.partial_edge_tolerance)
        db_gene_lengths = self._db.genes.lengths[gene_hits.gene_indices]

        # In Old Kaptive, truncation was calculated based on the *translated protein* length
        # (up to the first stop codon) divided by the reference protein length.
        prot_covs = (prot_seqs.lengths * 3.0) / db_gene_lengths

        # We must also update the coverage to reflect the protein coverage, so it prints correctly in the output
        gene_hits.coverages[:] = np.clip(prot_covs * 100.0, 0.0, 100.0)

        # A partial gene colliding with a contig edge is excluded from being truncated
        is_truncated = (~is_partial) & (prot_covs < 0.90)
        gene_states[is_partial] = GeneState.PARTIAL.value
        gene_states[is_truncated] = GeneState.TRUNCATED.value
        prot_alns = self._protein_aligner(prot_seqs, self._db.translations[gene_hits.gene_indices])  # type: ignore
        prot_idents = prot_alns.pidents.astype(np.float32)

        # Drop genes outside the locus that fall below the identity threshold,
        # mirroring Old Kaptive's behavior of ignoring spurious homologies
        is_spurious = (~gene_hits.is_inside) & (prot_idents < self._db.metadata.id_threshold)
        if np.any(is_spurious):
            keep_mask = ~is_spurious
            gene_hits = gene_hits[keep_mask]
            gene_seqs = gene_seqs[keep_mask]
            prot_seqs = prot_seqs[keep_mask]
            gene_states = gene_states[keep_mask]
            prot_idents = prot_idents[keep_mask]

        # Normal genes that fall below the identity threshold are considered NOVEL
        below_threshold = (gene_states == GeneState.NORMAL.value) & (prot_idents < self._db.metadata.id_threshold)
        gene_states[below_threshold] = GeneState.NOVEL.value
        valid_pidents = prot_idents[gene_states == GeneState.NORMAL.value]
        pident = float(np.mean(valid_pidents)) if valid_pidents.size > 0 else 0.0

        # Phenotype Evaluation phase -----------------------------------------------------------------------------------
        base_phenotype = self._db.serotypes[best_locus_idx]
        phenotypes = self._db.phenotypes

        if len(phenotypes) > 0:
            # A cluster is considered 'active' if it's found NORMAL or PARTIAL
            q_active = np.zeros(len(self._db.cluster_keys), dtype=bool)
            is_active = (gene_states == GeneState.NORMAL.value) | (gene_states == GeneState.PARTIAL.value)
            if np.any(is_active):
                active_clusters = self._db.gene_cluster_ids[gene_hits.gene_indices[is_active]]
                q_active[active_clusters] = True

            # Vectorized rule evaluation using BLAS matrix multiplication
            locus_match = phenotypes.locus_masks[:, best_locus_idx]
            q_active_int = q_active.astype(np.int8)
            extra_match = np.dot(phenotypes.extra_masks, q_active_int) == phenotypes.extra_counts

            has_inactive_rule = phenotypes.inactive_masks.sum(axis=1) > 0

            expected_mask = np.zeros(len(self._db.cluster_keys), dtype=np.int8)
            offset = self._db.locus_gene_offsets[best_locus_idx]
            length = self._db.locus_gene_lengths[best_locus_idx]
            expected_clusters = self._db.gene_cluster_ids[offset : offset + length]
            expected_mask[expected_clusters] = 1

            applicable_inactive_masks = phenotypes.inactive_masks & expected_mask
            has_applicable_inactive = applicable_inactive_masks.sum(axis=1) > 0

            q_inactive_int = (~q_active).astype(np.int8)
            inactive_hits = np.dot(applicable_inactive_masks, q_inactive_int)

            inactive_match = (~has_inactive_rule) | (has_applicable_inactive & (inactive_hits > 0))

            if np.any(valid_mask := locus_match & extra_match & inactive_match):
                valid_indices = np.where(valid_mask)[0]
                is_suffix = phenotypes.as_suffix[valid_indices]

                if len(replacements := valid_indices[~is_suffix]) > 0:
                    best_rep_idx = replacements[np.argmax(phenotypes.priorities[replacements])]
                    base_phenotype = phenotypes.ids[best_rep_idx].decode("utf-8")

                if len(suffixes := valid_indices[is_suffix]) > 0:
                    sorted_suffixes = suffixes[np.argsort(-phenotypes.priorities[suffixes])]
                    suffix_strs = [phenotypes.ids[i].decode("utf-8") for i in sorted_suffixes]
                    base_phenotype = f"{base_phenotype}{''.join(suffix_strs)}"

        # Confidence evaluation phase ----------------------------------------------------------------------------------
        typeable = True
        if actual_locus_completeness < self.min_completeness:
            typeable = False

        # Truncated unexpected genes do not count towards the limit, mirroring Old Kaptive
        is_unexpected = gene_hits.is_inside & ~gene_hits.is_expected & ~gene_hits.is_extra
        is_not_truncated = gene_states != GeneState.TRUNCATED.value
        unexpected_count = np.count_nonzero(is_unexpected & is_not_truncated)
        if unexpected_count > self.max_other_genes:
            typeable = False

        # 3. Check for any genes falling below the identity threshold
        if not self.allow_below_threshold:
            if np.any(gene_hits.is_inside & (gene_states == GeneState.NOVEL.value)):
                typeable = False

        # Return result object -----------------------------------------------------------------------------------------
        return SerotypingResult(
            kaptive_version=__version__,
            database_name=self._db.metadata.name,
            database_version=self._db.metadata.version,
            database_organism=self._db.metadata.organism,
            database_taxon=self._db.metadata.taxon,
            genome=genome.id,
            best_locus_idx=best_locus_idx,
            best_locus_name=best_locus_name,
            best_locus_score=locus_scores[best_locus_idx],
            best_locus_completeness=actual_locus_completeness,
            length_discrepancy=length_discrepancy,
            gene_hits=gene_hits,
            gene_states=gene_states,
            locus_pieces=locus_pieces,
            locus_seqs=locus_seqs,
            gene_seqs=gene_seqs,  # type: ignore
            translations=prot_seqs,  # type: ignore
            percent_identity=pident,
            percent_coverage=pcov,
            protein_identities=prot_idents,
            phenotype=base_phenotype,
            typeable=typeable,
            missing_expected_genes=missing_expected_genes,
        )
