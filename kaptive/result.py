"""
Copyright 2023 Tom Stanton (tomdstanton@gmail.com)
https://github.com/tomdstanton/kaptive-mapper

This file is part of kaptive-mapper. kaptive-mapper is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by the Free Software Foundation,
either version 3 of the License, or (at your option) any later version. kaptive-mapper is distributed
in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
details. You should have received a copy of the GNU General Public License along with kaptive-mapper.
If not, see <https://www.gnu.org/licenses/>.
"""
from typing import List, Tuple

from kaptive.database import Locus, Gene
from kaptive.minimap import PafLine, Minimap2Result, weighted_identity, get_overlapping_alignments, \
    target_ranges_covered_by_alignments
from kaptive.snps import VcfRecord, get_mutation
from kaptive.misc import range_overlap
from kaptive.log import warning


class Result:
    def __init__(self, locus: Locus, sample):
        self.sample = sample
        self.name = sample.name
        self.locus = locus
        self.db = locus.db


class MapResult(Result):
    """
    Result class to hold results for a single locus
    """

    def __init__(self, locus: Locus, sample, alignments: Minimap2Result):
        super().__init__(locus, sample)
        self.alignments = alignments
        self.aligned_ranges = self.alignments.target_ranges[self.locus.name]
        self.length = len(self.locus)
        self.aligned_bases = sum(i[1] - i[0] for i in self.aligned_ranges)
        self.coverage = self.aligned_bases / self.length * 100
        self.identity = weighted_identity(self.alignments.targets[self.locus.name], 0, self.length)
        self.ises = {}
        self.snps = {}
        self.expected_genes = []
        self.missing_genes = []

        for gene in self.locus.genes:
            if overlap := sum(
                    range_overlap(aligned_range, (int(gene.feature.location.start), int(gene.feature.location.end)))
                    for aligned_range in self.aligned_ranges):
                self.expected_genes.append(GeneResult(self, self.alignments.targets[self.locus.name], gene, overlap))
            else:
                self.missing_genes.append(gene.name)

    def __repr__(self):
        return f'{self.name} {self.locus.name}'

    def get_results(self, synonymous: bool = False):
        """
        Function to get results for a single locus
        """
        cols = [self.name,
                # str(self.locus.db),
                self.locus.name, self.locus.type, f"{self.identity:.1f}",
                f"{self.coverage:.1f}", ";".join(str(i) for i in self.aligned_ranges)]

        gene_column, mutations_column, ise_column, snp_column = [], [], [], []

        for gene in self.expected_genes:
            gene_column.append(f"{gene.name} {gene.identity:.1f}%")
            if gene.mutation:
                mutations_column.append(f"{gene.name} {gene.mutation}")

        for ise in self.ises.values():
            ise_column.append(f"{ise.name} {ise.identity:.1f}% {ise.coverage:.1f}% {ise.locus_region[0]} "
                              f"{ise.locus_region[1]} {' '.join(i.name for i in ise.genes)}")

        for snp in self.snps.values():
            if snp.consequence == "synonymous":
                if synonymous:
                    snp_column.append(f"{snp.biotype} {snp.dna_pos}")
            else:
                snp_column.append(f"{snp.biotype} {snp.dna_pos}")

        cols += [";".join(gene_column), ";".join(self.missing_genes), ";".join(ise_column),
                 ";".join(mutations_column), ";".join(snp_column)]

        return "\t".join(cols)

    def add_snps(self, vcf: str, all_csqs: bool):
        """
        Function to call SNPs in the locus using the snp_calling_pipeline function
        """
        # Get consensus fasta
        # self.consensus_fasta = seq_to_dict(vcf2consensus(vcf, self.locus.get_fasta()))
        # if len(self.consensus_fasta) > 1:
        #     raise ValueError("More than one sequence in consensus fasta")

        # Parse vcf output into VcfRecord objects
        snps = [VcfRecord(record) for record in vcf.splitlines() if record and not record.startswith("#")]
        # If SNP is part of same locus, add to snp dict
        self.snps |= {snp.pos: snp for snp in snps if snp.chrom == self.locus.name}
        # Add snps to genes
        for gene in self.expected_genes:
            gene.snps = [snp for snp in snps if snp.gene_id == gene.name]
            if gene.snps:
                gene.mutation = get_mutation(gene.snps, all_mutations=all_csqs)

    def add_ises(self, ise_paf: Minimap2Result, ise_range_merge_tolerance: int = 0):
        # Get reads that aligned to the best locus that are also in the ISE PAF
        locus_reads = {}
        for paf in self.alignments.targets[self.locus.name]:
            if paf.read in ise_paf.queries:
                if paf.read in locus_reads:
                    locus_reads[paf.read].append(paf)
                else:
                    locus_reads[paf.read] = [paf]

        if not locus_reads:
            return  # Failsafe 1

        # Get ISEs that have reads that aligned to the best locus
        ise_queries_in_locus = {}
        for read in locus_reads.keys():
            for paf in ise_paf.queries[read]:
                if paf.target_name in ise_queries_in_locus:
                    ise_queries_in_locus[paf.target_name].append(paf)
                else:
                    ise_queries_in_locus[paf.target_name] = [paf]

        if not ise_queries_in_locus:
            return  # Failsafe 2

        # Regions on the locus that are covered by alignments corresponding to ISEs alignments
        locus_paf = [i for paf in locus_reads.values() for i in paf]
        locus_ise_regions = target_ranges_covered_by_alignments(locus_paf, ise_range_merge_tolerance)[self.locus.name]

        # Iterate over each region and determine which ISE is in the region
        for n, region in enumerate(locus_ise_regions):
            best_ise = None
            best_overlap = 0
            # Iterate over each ISE query
            for ise, pafs in ise_queries_in_locus.items():
                # Get all alignments for the ISE that are in the locus
                locus_alignments_for_ise = [i for paf in pafs for i in locus_reads[paf.read]]
                
                ise_regions = target_ranges_covered_by_alignments(locus_alignments_for_ise, ise_range_merge_tolerance)[self.locus.name]

                for ise_region in ise_regions:
                    overlap = range_overlap(ise_region, region)
                    if overlap > best_overlap:
                        best_overlap = overlap

                        is_element = ISElement(ise, pafs, self, ise_region)
                        is_element.locus_region = region
                        if not best_ise or is_element.identity > best_ise.identity:
                            best_ise = is_element

            if best_ise:
                self.ises[f"ise_{n + 1}"] = best_ise
                for gene in self.expected_genes:
                    if range_overlap(region, (int(gene.feature.location.start), int(gene.feature.location.end))):
                        gene.ises.append(best_ise)
                        best_ise.genes.append(gene)
            else:
                warning(f"No ISE found for region {region} in {self.locus.name}")


class GeneResult:
    """
    Class to store results for a single gene including SNPs and IS elements based on the output of bedtools coverage
    """

    def __init__(self, result: MapResult, alignments: List[PafLine], gene: Gene, overlapping_bases: int):
        self.result = result
        self.alignments = get_overlapping_alignments(alignments, int(gene.feature.location.start),
                                                     int(gene.feature.location.end))
        self.gene = gene
        self.name = gene.name
        self.locus = gene.locus
        self.db = gene.locus.db
        self.overlapping_bases = overlapping_bases
        self.coverage = overlapping_bases / len(gene.feature) * 100
        self.feature = gene.feature
        self.identity = weighted_identity(self.alignments, int(gene.feature.location.start),
                                          int(gene.feature.location.end))
        self.snps = []
        self.ises = []
        self.mutation = ''

    def __repr__(self):
        return self.name


class ISElement:
    """
    Class that represents alignments of a single ISE corresponding to a region on the locus
    """
    def __init__(self, name: str, alignments: List[PafLine], result: MapResult, alignment_range: Tuple[int, int]):
        self.result = result
        self.locus = result.locus
        self.db = result.locus.db
        self.name = name
        self.alignments = alignments
        self.reads = set(paf.read for paf in alignments)
        self.length = alignments[0].target_length
        # self.alignment_ranges = target_ranges_covered_by_alignments(alignments)[name]
        self.alignment_range = alignment_range
        # self.overlapping_bases = sum(end - start for start, end in self.alignment_ranges)
        self.overlapping_bases = self.alignment_range[1] - self.alignment_range[0]
        self.coverage = self.overlapping_bases / self.length * 100
        self.identity = weighted_identity(alignments, 0, self.length)
        self.genes = []
        self.locus_region = None

    def __repr__(self):
        return f'{self.name} {self.coverage:.2f}% {self.identity:.2f}%'
#
#
# class AlignResult(Result):
#     def __init__(self, locus, sample):
#         super().__init__(locus, sample)
#         self.scores = None
#         self.blast_hits = []
#         self.hit_ranges = IntRange()
#         self.assembly_pieces = []
#         self.identity = 0.0
#         self.expected_hits_inside_locus = []
#         self.missing_expected_genes = []
#         self.expected_hits_outside_locus = []
#         self.other_hits_inside_locus = []
#         self.other_hits_outside_locus = []
#         # self.locus_fasta = Path(f'{self.args.out}_{self.assembly.name}.fasta')
#
#     def add_best_match_locus(self, locus: 'Locus', scores):
#         self.locus = locus
#         self.seq = locus.seq
#         self.name = locus.name
#         self.type = locus.type
#         self.scores = scores
#
#     def find_assembly_pieces(self):
#         """
#         This function uses the BLAST hits in the given locus type to find the corresponding pieces of
#         the given assembly. It saves its results in the Locus object.
#         """
#         if not self.blast_hits:
#             return
#         assembly_pieces = [x.get_assembly_piece(self.assembly) for x in self.blast_hits]
#         merged_pieces = merge_assembly_pieces(assembly_pieces)
#         length_filtered_pieces = [x for x in merged_pieces if len(x) >= self.args.min_assembly_piece]
#         if not length_filtered_pieces:
#             return
#         self.assembly_pieces = fill_assembly_piece_gaps(length_filtered_pieces, self.args.gap_fill_size, self.assembly)
#         # Now check to see if the biggest assembly piece seems to capture the whole locus. If so, this
#         # is an ideal match.
#         biggest_piece = sorted(self.assembly_pieces, key=lambda z: len(z), reverse=True)[0]
#         start = biggest_piece.earliest_hit_coordinate()
#         end = biggest_piece.latest_hit_coordinate()
#         if good_start_and_end(start, end, len(self.locus), self.args.start_end_margin):
#             self.assembly_pieces = [biggest_piece]
#
#         # If it isn't the ideal case, we still want to check if the start and end of the locus were
#         # found in the same contig. If so, fill all gaps in between so we include the entire
#         # intervening sequence.
#         else:
#             earliest, latest, same_contig_and_strand = self.get_earliest_and_latest_pieces()
#             start = earliest.earliest_hit_coordinate()
#             end = latest.latest_hit_coordinate()
#             if good_start_and_end(start, end, len(self.locus), self.args.start_end_margin) and \
#                     same_contig_and_strand:
#                 gap_filling_piece = AssemblyPiece(self.assembly, earliest.contig_name, earliest.start,
#                                                   latest.end, earliest.strand)
#                 self.assembly_pieces = merge_assembly_pieces(self.assembly_pieces + [gap_filling_piece])
#         pieces_across_contigs = len(set([piece.contig_name for piece in self.assembly_pieces]))
#         self.identity = get_mean_identity(self.assembly_pieces)
#
#     def add_blast_hit(self, hit):
#         """Adds a BLAST hit and updates the hit ranges."""
#         self.blast_hits.append(hit)
#         self.hit_ranges.add_range(hit.qstart, hit.qend)
#
#     def get_mean_blast_hit_identity(self):
#         """Returns the mean identity (weighted by hit length) for all BLAST hits in the locus."""
#         identity_sum = 0.0
#         length_sum = 0
#         for hit in self.blast_hits:
#             length_sum += hit.length
#             identity_sum += hit.length * hit.pident
#         if identity_sum == 0.0:
#             return 0.0
#         else:
#             return identity_sum / length_sum
#
#     def get_coverage(self, cumulative=False):
#         """
#         Returns the % of this locus which is covered by BLAST hits in the given assembly.
#         If cumulative==True, returns the % of the non-fragmented, expected genes this locus
#         Updated in Kaptive v3 to reflect the biology of the gene-based approach.
#         Missing and fragmented genes drag down the coverage score, and the two are correlated, so it has a snowball
#         effect on the confidence.
#         If we consider the cumulative coverage of the genes we KNOW are inside the locus, it gives us a more
#         realistic idea of the coverage, but we can still judge confidence based on missing genes.
#         """
#         try:
#             if not cumulative:
#                 return 100.0 * self.hit_ranges.get_total_length() / len(self.seq)
#             else:
#                 good_hits = [hit for hit in self.expected_hits_inside_locus if not hit.fragmented]
#                 return sum([i.query_cov for i in good_hits]) / len(good_hits)
#         except ZeroDivisionError:
#             return 0.0
#
#     def get_coverage_string(self):
#         return '%.2f' % self.get_coverage() + '%'
#
#     def get_identity_string(self):
#         return f"{'%.2f' % self.identity}%"
#
#     def clean_up_blast_hits(self):
#         """
#         This function removes unnecessary BLAST hits from self.blast_hits.
#         For each BLAST hit, we keep it if it offers new parts of the locus. If, on the other
#         hand, it lies entirely within an existing hit (in locus positions), we ignore it. Since
#         we first sort the BLAST hits longest to shortest, this strategy will prioritise long hits
#         over short ones.
#         """
#         self.blast_hits.sort(key=lambda x: x.length, reverse=True)
#         kept_hits = []
#         range_so_far = IntRange()
#         for hit in self.blast_hits:
#             hit_range = hit.get_query_range()
#             if not range_so_far.contains(hit_range):
#                 range_so_far.merge_in_range(hit_range)
#                 kept_hits.append(hit)
#         self.blast_hits = kept_hits
#
#     def get_match_uncertainty_chars(self):
#         """
#         Returns the character code which indicates uncertainty with how this locus was found in
#         the current assembly.
#         '?' means the locus was found in multiple discontinuous assembly pieces.
#         '-' means that one or more expected genes were missing.
#         '+' means that one or more additional genes were found in the locus assembly parts.
#         '*' means that at least one of the expected genes in the locus is low identity.
#         """
#         uncertainty_chars = ''
#         if len(self.assembly_pieces) > 1:
#             uncertainty_chars += '?'
#         if self.missing_expected_genes:
#             uncertainty_chars += '-'
#         if self.other_hits_inside_locus:
#             uncertainty_chars += '+'
#         if not all([x.over_identity_threshold for x in self.expected_hits_inside_locus]):
#             uncertainty_chars += '*'
#         return uncertainty_chars
#
#     def get_length_discrepancy(self):
#         """
#         Returns an integer of the base discrepancy between the locus in the assembly and the
#         reference locus sequence.
#         E.g. if the assembly match was 5 bases shorter than the reference, this returns -5.
#         This function only applies to cases where the locus was found in a single piece. In
#         other cases, it returns None.
#         """
#         if len(self.assembly_pieces) != 1:
#             return None
#         only_piece = self.assembly_pieces[0]
#         a_start = only_piece.start
#         a_end = only_piece.end
#         start = only_piece.earliest_hit_coordinate()
#         end = only_piece.latest_hit_coordinate()
#         expected_length = end - start
#         actual_length = a_end - a_start
#         return actual_length - expected_length
#
#     def get_length_discrepancy_string(self):
#         """
#         Returns the length discrepancy, not as an integer but as a string with a sign and units.
#         """
#         length_discrepancy = self.get_length_discrepancy()
#         if length_discrepancy is None:
#             return 'n/a'
#         length_discrepancy_string = str(length_discrepancy) + ' bp'
#         if length_discrepancy > 0:
#             length_discrepancy_string = '+' + length_discrepancy_string
#         return length_discrepancy_string
#
#     def get_earliest_and_latest_pieces(self):
#         """
#         Returns the AssemblyPiece with the earliest coordinate (closest to the locus start) and
#         the AssemblyPiece with the latest coordinate (closest to the locus end)
#         """
#         earliest_piece = sorted(self.assembly_pieces, key=lambda x: x.earliest_hit_coordinate())[0]
#         latest_piece = sorted(self.assembly_pieces, key=lambda x: x.latest_hit_coordinate())[-1]
#         same_contig_and_strand = earliest_piece.contig_name == latest_piece.contig_name and \
#                                  earliest_piece.strand == latest_piece.strand
#
#         # Even though the pieces are on the same contig and strand, we still need to check whether
#         # the earliest piece comes before the latest piece in that contig.
#         if same_contig_and_strand:
#             if earliest_piece.strand == '+' and earliest_piece.start > latest_piece.end:
#                 same_contig_and_strand = False
#             elif earliest_piece.strand == '-' and earliest_piece.start < latest_piece.end:
#                 same_contig_and_strand = False
#         return earliest_piece, latest_piece, same_contig_and_strand
#
#     def apply_special_logic(self):
#         """
#         This function has special logic for dealing with the locus -> type situations that depend on
#         other genes in the genome.
#         """
#         found_loci = sorted({x.locus.name.replace(  # Set comprehension to remove duplicates (fixes O2ac bug)
#             'Extra_genes_', '') for x in self.other_hits_outside_locus if x.qseqid.startswith('Extra_genes_')})
#         new_types = [
#             logic['type'] for logic in self.locus.special_logic if
#             self.name == logic['locus'] and found_loci == logic['extra_loci']
#         ]
#
#         if len(new_types) == 0:
#             self.type = 'unknown'
#         elif len(new_types) == 1:
#             self.type = new_types[0]
#         else:  # multiple matches - shouldn't happen!
#             quit_with_error('redundancy in special logic file')
#
#     def get_match_confidence(self):
#         """
#         These confidence thresholds match those specified in the paper supp. text, with the
#         addition of two new top-level categories: perfect and very high
#         """
#         single_piece = len(self.assembly_pieces) == 1
#         cov = self.get_coverage(cumulative=True)  # Use the cumulative coverage of non-fragmented expected genes
#         ident = self.identity
#         missing = len(self.missing_expected_genes)
#         extra = len(self.other_hits_inside_locus)
#
#         # Define thresholds for each confidence category
#         thresholds = {
#             "ident": [100, 95.0, 95.0, 95.0, 95.0],
#             "cov": [100, 95.0, 90.0, 85.0, 80.0],  # Relax the coverage thresholds
#             "missing": [0, 0, 2, 3, 4],
#             "extra": [0, 0, 0, 1, 2]
#         }
#         # Perfect logic
#         if single_piece and \
#                 cov == thresholds['cov'][0] and \
#                 ident == thresholds['ident'][0] and \
#                 missing == thresholds['missing'][0] and \
#                 extra == thresholds['extra'][0] and \
#                 self.get_length_discrepancy() == 0:
#             confidence = 'Perfect'
#         # Very high logic
#         elif single_piece and \
#                 cov >= thresholds['cov'][1] and \
#                 ident >= thresholds['ident'][1] and \
#                 missing == thresholds['missing'][1] and \
#                 extra == thresholds['extra'][1]:
#             confidence = 'Very high'
#         # High logic
#         elif single_piece and \
#                 cov >= thresholds['cov'][2] and \
#                 missing <= thresholds['missing'][2] and \
#                 extra == thresholds['extra'][2]:
#             confidence = 'High'
#         # Good logic
#         elif (single_piece or cov >= thresholds['cov'][3]) and \
#                 missing <= thresholds['missing'][3] and \
#                 extra <= thresholds['extra'][3]:
#             confidence = 'Good'
#         # Low logic
#         elif (single_piece or cov >= thresholds['cov'][4]) and \
#                 missing <= thresholds['missing'][4] and \
#                 extra <= thresholds['extra'][4]:
#             confidence = 'Low'
#         # None logic
#         else:
#             confidence = 'None'
#         return confidence
#
#     def add_to_table(self):
#         """
#         Writes a line to the output table describing all that we've learned about the given locus and
#         writes to stdout as well.
#         """
#         try:
#             expected_in_locus_per = 100.0 * len(self.expected_hits_inside_locus) / len(self.locus.gene_names)
#             expected_out_locus_per = 100.0 * len(self.expected_hits_outside_locus) / len(self.locus.gene_names)
#             expected_genes_in_locus_str = f'{len(self.expected_hits_inside_locus)} / ' \
#                                           f'{len(self.locus.gene_names)} ({expected_in_locus_per:.1f}%)'
#             expected_genes_out_locus_str = f'{len(self.expected_hits_outside_locus)} / ' \
#                                            f'{len(self.locus.gene_names)} ({expected_out_locus_per:.1f}%)'
#             missing_per = 100.0 * len(self.missing_expected_genes) / len(self.locus.gene_names)
#             missing_genes_str = f'{len(self.missing_expected_genes)} / {len(self.locus.gene_names)} ({missing_per}%)'
#         except ZeroDivisionError:
#             expected_genes_in_locus_str, expected_genes_out_locus_str, missing_genes_str = '', '', ''
#
#         line = [self.assembly.name,
#                 self.name,
#                 self.type,
#                 self.get_match_confidence(),
#                 self.get_match_uncertainty_chars(),
#                 self.get_coverage_string(),
#                 self.get_identity_string(),
#                 self.get_length_discrepancy_string(),
#                 expected_genes_in_locus_str,
#                 get_gene_info_string(self.expected_hits_inside_locus),
#                 ';'.join(self.missing_expected_genes),
#                 str(len(self.other_hits_inside_locus)),
#                 get_gene_info_string(self.other_hits_inside_locus),
#                 expected_genes_out_locus_str,
#                 get_gene_info_string(self.expected_hits_outside_locus),
#                 str(len(self.other_hits_outside_locus)),
#                 get_gene_info_string(self.other_hits_outside_locus)]
#
#         for gene_name in self.db.type_gene_names:
#             hit = self.type_gene_results[gene_name]
#             line.append('-' if not hit else hit.result)
#
#         with open(f'{self.args.out}_table.txt', 'at') as table:
#             table.write('\t'.join(line) + '\n')
#
#     def add_to_scores(self):
#         with open(f'{self.args.out}_scores.txt', 'at') as scores:
#             scores.write(
#                 '\n'.join([
#                     '\t'.join([self.assembly.name, i[0]] + [str(x) for x in i[1].values()]) for i in self.scores
#                 ]) + '\n')
#
#     def add_to_json(self):
#         expected_genes_in_locus = {x.qseqid: x for x in self.expected_hits_inside_locus}
#         expected_hits_outside_locus = {x.qseqid: x for x in self.expected_hits_outside_locus}
#         other_hits_inside_locus = {x.qseqid: x for x in self.other_hits_inside_locus}
#         other_hits_outside_locus = {x.qseqid: x for x in self.other_hits_outside_locus}
#         uncertainty_chars = self.get_match_uncertainty_chars()
#
#         json_record = OrderedDict({
#             'Assembly name': self.assembly.name,
#             'Best match': OrderedDict({
#                 'Locus name': self.name,
#                 'Type': self.locus.type,
#                 'Match confidence': self.get_match_confidence(),
#                 'Reference': OrderedDict({
#                     'Length': len(self.seq), 'Sequence': self.seq
#                 })
#             }),
#             'Problems': OrderedDict({
#                 'Locus assembled in multiple pieces': str('?' in uncertainty_chars),
#                 'Missing genes in locus': str('-' in uncertainty_chars),
#                 'Extra genes in locus': str('+' in uncertainty_chars),
#                 'At least one low identity gene': str('*' in uncertainty_chars)
#             }),
#             'blastn result': OrderedDict({
#                 'Coverage': self.get_coverage_string(),
#                 'Identity': self.get_identity_string(),
#                 'Length discrepancy': self.get_length_discrepancy_string(),
#                 'Locus assembly pieces': [
#                     OrderedDict({
#                         'Contig name': piece.contig_name,
#                         'Contig start position': piece.start + 1,
#                         'Contig end position': piece.end,
#                         'Contig strand': piece.strand,
#                         'Length': len(piece_seq := piece.get_sequence()),
#                         'Sequence': piece_seq
#                     }) for i, piece in enumerate(self.assembly_pieces)
#                 ]
#             }),
#             'Locus genes': [
#                 OrderedDict({
#                     'Name': gene.name,
#                     'Result': (
#                         'Found in locus' if gene.name in expected_genes_in_locus else 'Found outside locus' if
#                         gene.name in expected_hits_outside_locus else 'Not found'),
#                     'Reference': gene.get_reference_info_json_dict(),
#                     'blastn result': hit.get_blast_result_json_dict() if (
#                         hit := {**expected_genes_in_locus, **expected_hits_outside_locus}.get(
#                             gene.name)) else None,
#                     'Match confidence': hit.get_match_confidence() if hit else 'Not found'
#                 }) for gene in self.locus.genes],
#             'Other genes in locus': OrderedDict({
#                 gene_name: OrderedDict({
#                     'Reference': hit.gene.get_reference_info_json_dict(),
#                     'blastn result': hit.get_blast_result_json_dict()
#                 }) for gene_name, hit in other_hits_inside_locus.items()
#             }),
#             'Other genes outside locus': OrderedDict({
#                 gene_name: OrderedDict({
#                     'Reference': hit.gene.get_reference_info_json_dict(),
#                     'blastn result': hit.get_blast_result_json_dict()
#                 }) for gene_name, hit in other_hits_outside_locus.items()}),
#             'Scores': self.scores if self.args.scores else None
#         })
#
#         if self.db.type_gene_names:
#             allelic_typing = OrderedDict()
#             for gene_name in self.db.type_gene_names:
#                 allelic_type = OrderedDict()
#                 if not self.type_gene_results[gene_name]:
#                     allelic_type['Allele'] = 'Not found'
#                 else:
#                     blast_hit = self.type_gene_results[gene_name]
#                     allele = blast_hit.result
#                     if allele.endswith('*'):
#                         perfect_match = False
#                         allele = allele[:-1]
#                     else:
#                         perfect_match = True
#                     try:
#                         allele = int(allele)
#                     except ValueError:
#                         pass
#                     allelic_type['Allele'] = allele
#                     allelic_type['Perfect match'] = str(perfect_match)
#                     allelic_type['blastn result'] = blast_hit.get_blast_result_json_dict()
#                 allelic_typing[gene_name] = allelic_type
#             json_record['Allelic_typing'] = allelic_typing
#
#         write_json_file(self.args.out, json_record)
#
#     def write_fasta(self):
#         """
#         Creates a single FASTA file for all the assembly pieces.
#         Assumes all assembly pieces are from the same assembly.
#         """
#         if self.assembly_pieces:
#             with open(self.locus_fasta, 'wt') as f:
#                 f.write('\n'.join(
#                     [f'>{self.assembly.name}_{piece.get_header()}\n{fasta_wrap(piece.get_sequence())}' for piece in
#                      self.assembly_pieces]))
#                 f.write('\n')
#
#     def print_result(self):
#         out_string = self.assembly.name
#         if self.locus:
#             out_string += f': {self.name}{self.get_match_uncertainty_chars()}'
#         for gene_name in self.type_gene_results.items():
#             result = 'Not found' if not self.type_gene_results[gene_name] else self.type_gene_results[gene_name].result
#             out_string += f', {gene_name}({result})'
#         print(out_string)
