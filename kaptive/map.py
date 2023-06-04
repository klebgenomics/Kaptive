#!/usr/bin/env python
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
from tempfile import TemporaryDirectory
from pathlib import Path

from kaptive.log import warning, log
from kaptive.reads import load_read_files
from kaptive.database import Database
from kaptive.minimap import minimap2, paftools_sam2paf, Minimap2Result
from kaptive.result import MapResult
from kaptive.snps import snp_calling_pipeline


def run_map(args, db: Database):
    samples = load_read_files(args.reads)

    with TemporaryDirectory() as tmpdir:
    # tmpdir = Path.cwd() / "tmp"
    # tmpdir.mkdir(exist_ok=True)

        bed = Path(tmpdir) / "loci.bed"
        bed.write_text(db.get_bed())

        gff = Path(tmpdir) / "loci.gff"
        gff.write_text(db.get_gff())

        ref_fasta = Path(tmpdir) / "loci.fasta"
        ref_fasta.write_text(db.get_fasta())

        if args.header:
            print(
                'Assembly\tLocus\tType\tPercent_identity\tPercent_coverage\tAlignment_ranges\tExpected_genes\t'
                'Missing_genes\tISE_in_locus\tGene_mutations\tSNPs', end='\n'
            )

        for sample in samples:
            if args.verbose:
                log(f"Mapping {sample.name}...")

            # Choose best locus ----------------------------------------------------------------------------------------
            # This can be done in a number of different ways, but I've found the correct locus corresponds to the
            # highest number of alignments (i.e. depth) for either the whole locus or genes in each locus.
            # If we want to perform SNP calling, we need to output a sam file, so we can get the best locus like this:

            locus_sam = minimap2(
                query=sample.reads, target=db.get_fasta(), output=Path(tmpdir) / f"{sample.name}.sam",
                threads=args.threads, verbose=args.verbose, preset=args.preset,
            )
            if locus_sam is None:
                warning(f"No alignments found for {sample.name}")
                continue

            locus_paf = Minimap2Result(paftools_sam2paf(locus_sam))

            result = max(
                [MapResult(locus=db.loci[k], sample=sample, alignments=locus_paf) for k in db.loci if k in locus_paf.targets],
                key=lambda x: x.coverage)

            vcf = snp_calling_pipeline(ref_fasta, locus_sam.path, gff, args.threads, args.min_depth, args.min_qual,
                                       args.min_mq, args.csq_phase, [result.locus.name])
            if vcf:
                result.add_snps(vcf, args.all_csqs)
                if args.vcf:
                    vcf_file = args.vcf / f"{sample.name}_{result.locus.name}.vcf"
                    vcf_file.write_text(vcf)

            if args.keep_sam:
                locus_sam.path.rename(Path.cwd() / locus_sam.path.name)
            else:
                locus_sam.path.unlink()

            if args.ise:
                ise_paf = Minimap2Result(
                    minimap2(
                        query=sample.reads, target=args.ise, threads=args.threads, verbose=args.verbose,
                        preset=args.preset
                    )
                )
                if ise_paf:
                    result.add_ises(ise_paf, args.ise_range_merge_tolerance)

            print(result.get_results(args.synonymous), end="\n")
            # args.out.write(result.get_results(args.synonymous))
