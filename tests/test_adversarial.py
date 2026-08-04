import tempfile
from pathlib import Path

import numpy as np

from kaptive.bgc.annotate import AnnotationResult, Annotator, Genes
from kaptive.core.genome import GenomeAssembly
from kaptive.core.interval import Intervals
from kaptive.core.seq import SeqRecord, Sequences


def test_genes_soa_edge_cases() -> None:
    empty_genes = Genes(
        intervals=Intervals.empty(),
        translations=Sequences.empty(),
        contig_indices=np.empty(0, dtype=np.uint32),
    )
    assert len(empty_genes) == 0
    assert len(empty_genes[0:0]) == 0
    assert len(Genes.concat([empty_genes, empty_genes])) == 0

    g1 = Genes(
        intervals=Intervals(
            np.array([10], dtype=np.int32),
            np.array([100], dtype=np.int32),
            np.array([1], dtype=np.int8),
        ),
        translations=Sequences.from_bytes([b"ATGC"]),
        contig_indices=np.array([0], dtype=np.uint32),
    )
    assert len(g1) == 1
    c1 = Genes.concat([empty_genes, g1])
    assert len(c1) == 1

    c2 = Genes.concat([g1, empty_genes])
    assert len(c2) == 1

    slice_empty = g1[5:10]
    assert len(slice_empty) == 0


def test_annotator_short_and_empty_contigs() -> None:
    from kaptive.db.manager import DatabaseManager

    db = DatabaseManager.get("ab_k")
    annotator = Annotator(db=db)

    # Empty assembly
    empty_genome = GenomeAssembly("empty_asm", Sequences.empty())
    res_empty = annotator(empty_genome)
    assert len(res_empty.genes) == 0
    assert len(res_empty.hits_mask) == 0
    assert len(res_empty.top_hit_names) == 0
    assert len(res_empty.top_hit_scores) == 0

    with tempfile.NamedTemporaryFile(mode="w+", suffix=".bed", delete=False) as f:
        out_path = Path(f.name)

    res_empty.write_bed(out_path, hits_only=True)
    assert out_path.read_text() == ""

    res_empty.write_bed(out_path, hits_only=False)
    assert out_path.read_text() == ""

    # Contigs < 3 bp
    short_genome = GenomeAssembly(
        "short_asm",
        Sequences.from_records(
            [
                SeqRecord(id="short_contig_1", seq=b"A"),
                SeqRecord(id="short_contig_2", seq=b"AT"),
                SeqRecord(id="short_contig_3", seq=b"CG"),
                SeqRecord(id="short_contig_4", seq=b""),
            ]
        ),
    )
    res_short = annotator(short_genome)
    assert len(res_short.genes) == 0

    res_short.write_bed(out_path, hits_only=True)
    assert out_path.read_text() == ""
    out_path.unlink()


def test_bed_tag_formatting_and_zero_hits() -> None:
    syn_genes = Genes(
        intervals=Intervals(
            np.array([100, 500], dtype=np.int32),
            np.array([300, 800], dtype=np.int32),
            np.array([1, -1], dtype=np.int8),
        ),
        translations=Sequences.from_bytes([b"MKV", b"MKL"]),
        contig_indices=np.array([0, 0], dtype=np.uint32),
    )

    syn_res = AnnotationResult(
        genes=syn_genes,
        translations_idx=None,
        seeds=None,
        hits_mask=np.array([True, False], dtype=bool),
        top_hit_names=np.array(["wzi_1", ""], dtype=object),
        top_hit_scores=np.array([98.543, 0.0], dtype=np.float32),
        contig_names=("contig1",),
    )

    with tempfile.NamedTemporaryFile(mode="w+", suffix=".bed", delete=False) as f:
        out_path = Path(f.name)

    syn_res.write_bed(out_path, hits_only=True)
    lines = out_path.read_text().strip().split("\n")
    assert len(lines) == 1
    assert lines[0] == "contig1\t100\t300\t0\t0\t+\ttop_hit=wzi_1;score=98.54"

    syn_res.write_bed(out_path, hits_only=False)
    lines_all = out_path.read_text().strip().split("\n")
    assert len(lines_all) == 2
    assert lines_all[0] == "contig1\t100\t300\t0\t0\t+\ttop_hit=wzi_1;score=98.54"
    assert lines_all[1] == "contig1\t500\t800\t1\t0\t-\t."

    out_path.unlink()


def test_score_boundary_formatting() -> None:
    syn_genes = Genes(
        intervals=Intervals(
            np.array([100, 500], dtype=np.int32),
            np.array([300, 800], dtype=np.int32),
            np.array([1, -1], dtype=np.int8),
        ),
        translations=Sequences.from_bytes([b"MKV", b"MKL"]),
        contig_indices=np.array([0, 0], dtype=np.uint32),
    )

    syn_res_boundary = AnnotationResult(
        genes=syn_genes,
        translations_idx=None,
        seeds=None,
        hits_mask=np.array([True, True], dtype=bool),
        top_hit_names=np.array(["geneA", "geneB"], dtype=object),
        top_hit_scores=np.array([0.0, 100.0], dtype=np.float32),
        contig_names=("contig1",),
    )

    with tempfile.NamedTemporaryFile(mode="w+", suffix=".bed", delete=False) as f:
        out_path = Path(f.name)

    syn_res_boundary.write_bed(out_path, hits_only=True)
    b_lines = out_path.read_text().strip().split("\n")
    assert b_lines[0] == "contig1\t100\t300\t0\t0\t+\ttop_hit=geneA;score=0.00"
    assert b_lines[1] == "contig1\t500\t800\t1\t0\t-\ttop_hit=geneB;score=100.00"

    out_path.unlink()
