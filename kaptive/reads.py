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
from re import compile, IGNORECASE
from gzip import open as gzopen
from tempfile import TemporaryDirectory
from pathlib import Path

from kaptive.log import log, warning, quit_with_error
from kaptive.database import Database
from kaptive.minimap import minimap2, paftools_sam2paf, Minimap2Result
from kaptive.result import MapResult
from kaptive.snps import snp_calling_pipeline


# File path regexes ---------------------------------------------------------------------------------------------------
FASTQ_REGEX = compile('|'.join([
    '_R[12]\.(f(?:ast)?q(?:\.gz)?)$',
    '_R[12]_[0-9]+?\.(f(?:ast)?q(?:\.gz)?)$',
    '_R[12].[0-9]+?\.(f(?:ast)?q(?:\.gz)?)$',
    '_[12]\.(f(?:ast)?q(?:\.gz)?)$',
    '_[12]_[0-9]+?\.(f(?:ast)?q(?:\.gz)?)$',
    '_[12].[0-9]+?\.(f(?:ast)?q(?:\.gz)?)$',
    '\.(f(?:ast)?q(?:\.gz)?)$']), IGNORECASE)


# Classes -------------------------------------------------------------------------------------------------------------
class ReadGroup:
    """
    Sample class to hold sample name and reads
    Acts like a ReadGroup but will also be used to hold results
    """

    def __init__(self, reads: list['ReadFile'], sample_name):
        self.reads = reads
        self.name = sample_name
        [setattr(read, 'sample', self) for read in self.reads]
        self.path = self.reads[0].path.parent
        self.prefix = self.path / self.name

    def __repr__(self):
        return self.name


class ReadFile:
    def __init__(self, filepath: Path, extension: str):
        self.path = filepath
        self.extension = extension
        self.sample = None
        self.is_gzipped = self.extension.endswith('.gz')
        self.open_mode = 'rt' if self.is_gzipped else 'r'
        self.open_function = open if not self.is_gzipped else gzopen

    def __repr__(self):
        return str(self.path)


def load_read_files(file_list: List[Path]):
    reads = {}
    for file in file_list:
        extension_match = FASTQ_REGEX.search(str(file))
        if extension_match:
            sample_name = file.name.replace(extension_match[0], '')
            if sample_name in reads.keys():
                reads[sample_name].append(ReadFile(file, extension_match[0]))
            else:
                reads[sample_name] = [ReadFile(file, extension_match[0])]
        else:
            warning(f'{file.name} does not match the fastq extension regex')
    samples = []
    if reads:
        for sample, files in reads.items():
            if len(files) != 2:
                warning(f'Files for {sample} files not paired: {" ".join(str(i) for i in files)}')
            else:
                samples.append(ReadGroup(files, sample))
    if not samples:
        quit_with_error('No files to analyse')
    return samples


def reads_pipeline(args, db: Database):
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
