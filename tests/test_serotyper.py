import pytest
from pathlib import Path

from kaptive.serotyping import Serotyper
from kaptive.db import Database
from kaptive.core.genome import GenomeAssembly

genomes_dir = Path('/Users/tsta0015/Programming/kpsc')

# Define the test cases using pytest.mark.parametrize
@pytest.mark.parametrize("genome_file, expected_locus", [
    ('KP_NORM_BLD_111598_70.fasta.gz', 'KL48'),
    ('NK_H21_015_80.fasta.gz', 'KL103'),
    ('NK_H12_048_10.fasta.gz', 'KL30'),
    ('NK_H17_042_10.fasta.gz', 'KL10'),
])

def test_serotyper_problem_genomes(genome_file, expected_locus):
    """
    Tests the serotyper on problematic genomes to ensure the correct K-locus is identified.
    """
    genome_path = genomes_dir / genome_file
    assert genome_path.is_file(), f"Genome file not found: {genome_path}"

    genome = GenomeAssembly.from_file(genome_path)
    
    db_path = 'kpsc_k'  # Assuming 'kpsc_k' can be loaded by the database
    db = Database.load(db_path)

    with Serotyper(db) as serotyper:
        result = serotyper(genome)
        
        assert result is not None, f"Serotyper returned no result for {genome_file}"
        assert result.best_match_locus == expected_locus, \
            f"For {genome_file}, expected {expected_locus} but got {result.best_match_locus}"



def test():
    from pathlib import Path

    from kaptive.serotyping import Serotyper
    from kaptive.db import Database
    from kaptive.core.genome import GenomeAssembly

    genomes_dir = Path('/Users/tsta0015/Programming/kpsc')
    test_genomes = [
        'KP_NORM_BLD_111598_70.fasta.gz',
        'NK_H21_015_80.fasta.gz',
        'NK_H12_048_10.fasta.gz',
        'NK_H17_042_10.fasta.gz',
    ]

    with Serotyper(Database.load('kpsc_k')) as serotyper:
        results = [serotyper(GenomeAssembly.from_file(genomes_dir / i)) for i in test_genomes]