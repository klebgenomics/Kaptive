# /// script
# requires-python = ">=3.14"
# dependencies = [
#     "marimo>=0.23.8",
# ]
# ///

import marimo

__generated_with = "0.23.10"
app = marimo.App(width="medium")

with app.setup:
    from pathlib import Path
    from urllib.request import urlopen
    import pickle
    import json
    from gzip import open as gz_open

    import marimo as mo
    import numpy as np
    import pandas as pd

    from kaptive._version import __version__

    def _get_latest_version() -> str:
        with urlopen(f"https://pypi.org/pypi/kaptive/json") as response:
            if response.status == 200:
                return json.load(response)["info"]["version"]            
        return "Package not found"

    latest_version = _get_latest_version()
    assembly_dtype = pd.CategoricalDtype(
        categories=['10x', '20x', '30x', '40x', '50x', '60x', '70x', '80x', '90x', '100x', 'Unsubsampled'], ordered=True)


@app.cell(hide_code=True)
def _():
    mo.md(f"""
    # Kaptive Update Testing 🧪

    ## Introduction ✏️

    Here we are comparing `kaptive` v{__version__} to the current pypi version {latest_version} using the same tests
    we performed in the Kaptive 3 paper[^1].
    """)
    return


@app.cell(hide_code=True)
def _():
    mo.md(r"""
    ## Loading the Kaptive 3 Paper Results 📄
    We can load in the Kaptive 2 vs Kaptive 3 ground truth comparisons (Table S2) from the Kaptive 3 Zenodo repository.
    """)
    return


@app.function
def load_kaptive3_data(filepath: Path = Path('tests') / 'data' / 'kaptive3_s2.csv.gz'):    
    if not filepath.is_file():
        with urlopen("https://api.figshare.com/v2/file/download/52165283") as _response:
            with gz_open(filepath, 'wb') as _fh:
                _fh.write(_response.read())
    
    with gz_open(filepath, 'rb') as _fh:
        return (
            pd.read_csv(_fh)
            .loc[
                lambda i:
                (i.Species == 'KpSC') &
                (i.Polysaccharide == 'K') &
                (i.Assembly_type != 'Complete') &
                (i.Kaptive_version == 3) &
                (i.Ground_truth != 'novel_locus')
            ]
            .set_index('Assembly_name')
            .assign(Assembly_type = lambda i: i.Assembly_type.astype(assembly_dtype))
        )


@app.cell
def _():
    truths = load_kaptive3_data()
    mo.ui.table(truths)
    return (truths,)


@app.cell
def _(truths):
    (truths
    .groupby('Assembly_type')
    .agg(
        N_genomes=('Locus_outcome_to_ground_truth', 'count'),
        N_match=('Locus_outcome_to_ground_truth', lambda x: (x == 'correct').sum())
    )
    .assign(Prop_match=lambda i: i.N_match/i.N_genomes)
    )
    return


@app.cell(hide_code=True)
def _():
    mo.md("""
    ## Running the new Kaptive 🚀

    The newest version of Kaptive is an API-first implementation that is perfectly happy being run in a notebook or code chunk

    The running on 3802 KpSC (sub-sampled Illumina) assebmlies from the Kaptive 3 paper should take ~10 minutes - however we can dump the results to a compressed pickle so we don't have to keep re-running the serotyping.

    Dumping the `kaptive.seroyping.SerotypingResult` objects allows us to keep all information including the scoring metrics of each assembly.
    """)
    return


@app.cell
def get_kpsc_k_results(truths, zstd):
    def get_kpsc_k_results(filepath: Path = Path('tests') / 'data' / 'kpsc_k_results.pkl.gz'):
        if not filepath.is_file() or filepath.stat().st_size > 10_000:
            from kaptive.serotyping import Serotyper
            from kaptive.db import Database 
            with Serotyper(Database.load('kpsc_k'), indexing_threads=8) as st, gz_open(filepath, 'wb') as fh:
                pickle.dump([r for g in Path('../kpsc').glob('*.fasta.gz') if (r := st(g)) is not None], fh)

        from kaptive.serotyping import ConfidenceEvaluator, SerotypingProblem
        with zstd.open(filepath, 'rb') as _fh:    
            ce = ConfidenceEvaluator()
            _records = [{
                'Assembly_name': r.genome,
                'Best_match_locus': r.best_locus_name,
                'Score': r.best_locus_score,
                'Completeness': r.best_locus_completeness,
                'Identity': r.percent_identity,
                'Coverage': r.percent_coverage,
                'Problems': r.problems.value,
                'Typeable': ce(r)
            } for r in pickle.load(_fh)]
        
            return (
                pd.DataFrame.from_records(_records)
                .set_index('Assembly_name')
                .join(
                    truths[['Genome', 'Assembly_type', 'Ground_truth', 'Best_match_locus', 'Match_confidence', 'Notes']]
                        .rename(columns={
                            'Best_match_locus': 'Best_match_locus_prev', 
                            'Match_confidence': 'Typeable_prev'
                        }),
                    how='inner'
                )
                .assign(
                    Truth_match=lambda i: i.Best_match_locus == i.Ground_truth,
                    Truth_match_prev=lambda i: i.Best_match_locus_prev == i.Ground_truth,
                    Typeable_prev=lambda i: i.Typeable_prev == 'Typeable',
                    Missing_genes=lambda i: (i.Problems & SerotypingProblem.MISSING_GENES.value).astype(bool),
                    Truncated_genes=lambda i: (i.Problems & SerotypingProblem.TRUNCATED_GENES.value).astype(bool),
                    Unexpected_genes=lambda i: (i.Problems & SerotypingProblem.UNEXPECTED_GENES.value).astype(bool),
                    Novel_genes=lambda i: (i.Problems & SerotypingProblem.NOVEL_GENES.value).astype(bool),
                    Fragmented=lambda i: (i.Problems & SerotypingProblem.FRAGMENTED.value).astype(bool),
                )
            )

    return (get_kpsc_k_results,)


@app.cell(hide_code=True)
def _():
    mo.md(r"""
    ### Parsing the results
    Here we extract relevant information from the `kaptive.seroyping.SerotypingResult` objects into a pandas dataframe for vectorised analysis and plotting.
    """)
    return


@app.cell
def _(get_kpsc_k_results):
    kpsc_k_results = get_kpsc_k_results()
    mo.ui.table(kpsc_k_results)
    return (kpsc_k_results,)


@app.cell
def _(kpsc_k_results):
    kpsc_k_results.query('Truth_match == False and Truth_match_prev == True')
    return


@app.cell
def _(kpsc_k_results):
    kpsc_k_results.query('Truth_match == True and Truth_match_prev == False')
    return


@app.cell
def _(kpsc_k_results):
    _metrics = [
        'Truth_match', 'Typeable', 'Missing_genes', 
        'Truncated_genes', 'Extra_genes', 'Novel_genes', 'Fragmented'
    ]

    (
        kpsc_k_results
        .groupby('Assembly_type')
        .agg(
            N_genomes=('Truth_match', 'count'),
            **{f'N_{m}': (m, 'sum') for m in _metrics}
        )
        .assign(**{
            f'Prop_{m}': lambda i, col=f'N_{m}': i[col] / i['N_genomes'] 
            for m in _metrics
        })
    )
    return


@app.cell
def _():
    list(Path('../ab_mrsn/').glob('*.tsv'))
    return


@app.cell
def _(truths, zstd):
    def get_kpsc_k_results(filepath: Path = Path('tests') / 'data' / 'kpsc_k_results.pkl.gz'):
        if not filepath.is_file() or filepath.stat().st_size > 10_000:
            from kaptive.serotyping import Serotyper
            from kaptive.db import Database 
            with Serotyper(Database.load('kpsc_k'), indexing_threads=8) as st, gz_open(filepath, 'wb') as fh:
                pickle.dump([r for g in Path('../kpsc').glob('*.fasta.gz') if (r := st(g)) is not None], fh)

        from kaptive.serotyping import ConfidenceEvaluator, SerotypingProblem
        with zstd.open(filepath, 'rb') as _fh:    
            ce = ConfidenceEvaluator()
            _records = [{
                'Assembly_name': r.genome,
                'Best_match_locus': r.best_locus_name,
                'Score': r.best_locus_score,
                'Completeness': r.best_locus_completeness,
                'Identity': r.percent_identity,
                'Coverage': r.percent_coverage,
                'Problems': r.problems.value,
                'Typeable': ce(r)
            } for r in pickle.load(_fh)]
        
            return (
                pd.DataFrame.from_records(_records)
                .set_index('Assembly_name')
                .join(
                    truths[['Genome', 'Assembly_type', 'Ground_truth', 'Best_match_locus', 'Match_confidence', 'Notes']]
                        .rename(columns={
                            'Best_match_locus': 'Best_match_locus_prev', 
                            'Match_confidence': 'Typeable_prev'
                        }),
                    how='inner'
                )
                .assign(
                    Truth_match=lambda i: i.Best_match_locus == i.Ground_truth,
                    Truth_match_prev=lambda i: i.Best_match_locus_prev == i.Ground_truth,
                    Typeable_prev=lambda i: i.Typeable_prev == 'Typeable',
                    Missing_genes=lambda i: (i.Problems & SerotypingProblem.MISSING_GENES.value).astype(bool),
                    Truncated_genes=lambda i: (i.Problems & SerotypingProblem.TRUNCATED_GENES.value).astype(bool),
                    Unexpected_genes=lambda i: (i.Problems & SerotypingProblem.UNEXPECTED_GENES.value).astype(bool),
                    Novel_genes=lambda i: (i.Problems & SerotypingProblem.NOVEL_GENES.value).astype(bool),
                    Fragmented=lambda i: (i.Problems & SerotypingProblem.FRAGMENTED.value).astype(bool),
                )
            )

    return (get_kpsc_k_results,)


@app.cell
def _():
    pd.read_csv("../ab_mrsn/ab_mrsn_truths.tsv", delimiter='\t')
    return


@app.cell
def _():
    pd.read_csv("../ab_mrsn/data_summary.tsv", delimiter='\t')
    return


@app.cell
def _():
    ab_dataset = pd.read_csv("../ab_mrsn/ncbi_dataset.tsv", delimiter='\t')
    ab_dataset
    return


@app.cell(hide_code=True)
def _():
    mo.md(rf"""
    ## References 📚

    [^1]: Stanton TD, Hetland MAK, Löhr IH, Holt KE, Wyres KL. Fast and
        Accurate in silico Antigen Typing with Kaptive 3.
        2025 _Microbial Genomics_ 11(6):001428.
        <https://doi.org/10.1099/mgen.0.001428>
    """)
    return


if __name__ == "__main__":
    app.run()
