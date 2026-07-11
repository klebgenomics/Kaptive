import csv
import subprocess
import os
import time
import sys

def main():
    truth_file = "../ab_mrsn/ab_mrsn_k_locus_truths.tsv"
    assembly_dir = "../ab_mrsn"

    # read ground truth
    ground_truths = {}
    with open(truth_file, 'r') as f:
        reader = csv.DictReader(f, delimiter='\t')
        for row in reader:
            ground_truths[row['Assembly']] = row['Expected_locus']

    # Assemblies to test are the ones that had mismatches in the new run
    mismatched_assemblies = [
        "GCF_006658595.1_ASM665859v1_genomic",
        "GCF_006658565.1_ASM665856v1_genomic",
        "GCF_006494615.1_ASM649461v1_genomic",
        "GCF_006494265.1_ASM649426v1_genomic",
        "GCF_006494225.1_ASM649422v1_genomic",
        "GCF_006494155.1_ASM649415v1_genomic",
        "GCF_006494085.1_ASM649408v1_genomic",
        "GCF_006493685.1_ASM649368v1_genomic",
        "GCF_006493045.1_ASM649304v1_genomic",
        "GCF_006492945.1_ASM649294v1_genomic",
        "GCF_006492865.1_ASM649286v1_genomic",
        "GCF_006492805.1_ASM649280v1_genomic",
        "GCF_006492705.1_ASM649270v1_genomic",
        "GCF_006492675.1_ASM649267v1_genomic",
        "GCF_006492515.1_ASM649251v1_genomic",
        "GCF_006492495.1_ASM649249v1_genomic",
        "GCF_006492245.1_ASM649224v1_genomic",
        "GCF_006492215.1_ASM649221v1_genomic",
        "GCF_006492175.1_ASM649217v1_genomic",
        "GCF_006492125.1_ASM649212v1_genomic",
        "GCF_006492075.1_ASM649207v1_genomic",
        "GCF_006491875.1_ASM649187v1_genomic",
        "GCF_006491855.1_ASM649185v1_genomic"
    ]

    print("Assembly | Truth | Pred (New) | Pred (Old) | Old Conf")

    import textwrap

    script = textwrap.dedent("""
    import sys
    from kaptive.core.genome import GenomeAssembly
    from kaptive.db import Database
    from kaptive.serotyping_old import Serotyper
    db = Database.load('ab_k')
    typer = Serotyper(db)
    res = typer(sys.argv[1])
    print(res.best_locus_name + "\t" + ("Typeable" if res.typeable else "Untypeable"))
    """)

    with open("scratch/run_old_kaptive.py", "w") as f:
        f.write(script)

    for assembly in mismatched_assemblies:
        expected_locus = ground_truths[assembly]
        assembly_path = os.path.join(assembly_dir, f"{assembly}.fna.gz")
        
        result = subprocess.run(
            ["uv", "run", "python", "scratch/run_old_kaptive.py", assembly_path],
            capture_output=True, text=True
        )

        if result.returncode != 0:
            old_pred = "ERROR"
            old_conf = "ERROR"
        else:
            parts = result.stdout.strip().split("\t")
            old_pred = parts[0]
            old_conf = parts[1] if len(parts) > 1 else ""

        print(f"{assembly} | {expected_locus} | ... | {old_pred} | {old_conf}")

if __name__ == "__main__":
    main()
