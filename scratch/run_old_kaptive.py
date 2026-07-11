import sys
from kaptive.db import Database
from kaptive.core.genome import GenomeAssembly
from kaptive.serotyping_old import Serotyper

def main():
    db = Database.load("ab_k")
    assembly = GenomeAssembly.ensure("../ab_mrsn/GCF_006492945.1_ASM649294v1_genomic.fna.gz")
    with Serotyper(db) as serotyper:
        res = serotyper(assembly)
        print(f"Old Kaptive -> Locus: {res.best_locus_name}, Typeable: {res.typeable}")

if __name__ == "__main__":
    main()
