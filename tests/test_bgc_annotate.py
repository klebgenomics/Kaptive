"""Unit tests for kaptive.bgc.annotate (Annotator, Genes, AnnotationResult)."""

from pathlib import Path
import numpy as np
import pytest

from kaptive.bgc.annotate import AnnotationResult, Annotator, Genes
from kaptive.core.genome import GenomeAssembly
from kaptive.core.interval import Interval, Intervals, Strand
from kaptive.core.kmers import FracMinHashIndex, RandstrobeIndex, Seeds
from kaptive.core.seq import SeqRecord, Sequences
from kaptive.db import Database, DatabaseMetadata, Phenotypes


def create_sample_genes() -> Genes:
    """Helper to create a sample Genes container with 2 genes."""
    ivs = Intervals.from_intervals(
        [
            Interval(10, 100, strand=Strand.FORWARD),
            Interval(200, 500, strand=Strand.REVERSE),
        ]
    )
    seqs = Sequences.from_records(
        [
            SeqRecord(id="gene1", seq=b"MKTLLILAV"),
            SeqRecord(id="gene2", seq=b"MSKGEELFTG"),
        ]
    )
    contig_indices = np.array([0, 0], dtype=np.uint32)
    return Genes(intervals=ivs, translations=seqs, contig_indices=contig_indices)


def create_sample_db() -> Database:
    """Helper to create a valid Database instance."""
    db_metadata = DatabaseMetadata(
        name="Test Database",
        keyword="test_db",
        genbank="test.gbk",
        organism="Klebsiella pneumoniae",
        taxon=573,
        antigen="K",
        pathway="Wzx/Wzy",
        version="1.0",
        id_threshold=80.0,
        doi=["10.1000/182"],
        owner="klebgenomics",
        repo="Kaptive",
        branch="main",
        contact={"Curator": "curator@example.com"},
        phenotype_logic={},
        antigenic_units={},
    )
    loci = Sequences.from_records([SeqRecord(id="K1", seq=b"ATGCGT")])
    genes = Sequences.from_records([SeqRecord(id="gene1", seq=b"ATGCGT")])
    translations = Sequences.from_records([SeqRecord(id="prot1", seq=b"MR")])
    sketches = FracMinHashIndex.build(loci)
    return Database(
        metadata=db_metadata,
        loci=loci,
        serotypes=("K1",),
        locus_gene_offsets=np.array([0], dtype=np.uint32),
        locus_gene_lengths=np.array([1], dtype=np.uint32),
        gene_intervals=Intervals.from_intervals([Interval(0, 6, strand=Strand.FORWARD)]),
        genes=genes,
        translations=translations,
        extra_genes=np.array([False], dtype=bool),
        gene_locus_indices=np.array([0], dtype=np.uint16),
        cluster_keys=("c1",),
        gene_cluster_ids=np.array([0], dtype=np.uint16),
        description_keys=("d1",),
        gene_description_ids=np.array([0], dtype=np.uint16),
        gene_positions=np.array([1], dtype=np.uint16),
        phenotypes=Phenotypes.empty(),
        loci_sketches=sketches,
    )


def test_genes_basic_operations():
    """Test len, indexing, empty, and concat on Genes container."""
    genes = create_sample_genes()
    assert len(genes) == 2

    # Integer indexing
    iv, seq_rec, contig_idx = genes[0]
    assert iv.start == 10
    assert iv.end == 100
    assert seq_rec.id == "gene1"
    assert contig_idx == 0

    # Slicing
    sliced = genes[1:]
    assert isinstance(sliced, Genes)
    assert len(sliced) == 1
    assert sliced[0][1].id == "gene2"

    # Empty
    empty_genes = Genes.empty()
    assert len(empty_genes) == 0

    # Concat
    concatenated = Genes.concat([genes, sliced])
    assert len(concatenated) == 3
    assert concatenated[0][1].id == "gene1"
    assert concatenated[2][1].id == "gene2"

    # Concat empty list
    empty_concat = Genes.concat([])
    assert len(empty_concat) == 0


def test_annotation_result_bed_export(tmp_path: Path):
    """Test BED file generation from AnnotationResult."""
    genes = create_sample_genes()
    seqs = Sequences.from_records([SeqRecord(id="g1", seq=b"MKTLLILAV")])
    rs_idx = RandstrobeIndex.build(seqs)
    seeds = Seeds.empty()

    hits_mask = np.array([True, False], dtype=bool)
    top_hit_names = np.array(["ref_prot_A", "none"], dtype=object)
    top_hit_scores = np.array([95.5, 0.0], dtype=np.float32)
    contig_names = ("contig_1",)

    result = AnnotationResult(
        genes=genes,
        translations_idx=rs_idx,
        seeds=seeds,
        hits_mask=hits_mask,
        top_hit_names=top_hit_names,
        top_hit_scores=top_hit_scores,
        contig_names=contig_names,
    )

    bed_file = tmp_path / "output_hits.bed"
    result.write_bed(bed_file, hits_only=True)
    lines = bed_file.read_text().strip().splitlines()
    assert len(lines) == 1
    assert "contig_1" in lines[0]
    assert "top_hit=ref_prot_A" in lines[0]

    bed_file_all = tmp_path / "output_all.bed"
    result.write_bed(bed_file_all, hits_only=False)
    lines_all = bed_file_all.read_text().strip().splitlines()
    assert len(lines_all) == 2


def test_annotator_mock_run():
    """Test Annotator initialization and annotation call on synthetic genome."""
    db = create_sample_db()
    annotator = Annotator(db=db, align=False, whole_genome=False)

    genome_seqs = Sequences.from_records([SeqRecord(id="chr1", seq=b"ATGAAAACTCTGCTGATCCTGGCTGTCTAA")])
    assembly = GenomeAssembly(id="assembly1", contigs=genome_seqs)

    result = annotator(assembly)
    assert isinstance(result, AnnotationResult)
    assert len(result.contig_names) == 1
    assert result.contig_names[0] == "chr1"
