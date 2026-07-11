import pandas as pd
import re

truth_file = '../ab_mrsn/ab_mrsn_k_locus_truths.tsv'
kaptive_file = 'ab_mrsn_k_loci.tsv'

# Read truths
truth_df = pd.read_csv(truth_file, sep='\t')
truth_map = {}
for _, row in truth_df.iterrows():
    asm = str(row['Assembly']).strip()
    expected = str(row['Expected_locus']).strip()
    truth_map[asm] = expected

# Read kaptive results
kaptive_df = pd.read_csv(kaptive_file, sep='\t')

mismatches = []

for _, row in kaptive_df.iterrows():
    asm = str(row['Assembly']).strip()
    
    # Kaptive appends .fna.gz or similar sometimes? No, it usually strips extension or keeps it.
    # Let's match by prefix
    matched_key = None
    for key in truth_map:
        if asm.startswith(key) or key.startswith(asm):
            matched_key = key
            break
            
    if not matched_key:
        continue
        
    expected = truth_map[matched_key]
    predicted = str(row['Best match locus']).strip()
    
    # Simplify expected (e.g. "KL24 + GI-2" -> "KL24")
    expected_base = expected.split(' ')[0]
    
    if expected_base != predicted and not expected.startswith(predicted) and not predicted.startswith(expected):
        mismatches.append({
            'Assembly': asm,
            'Expected': expected,
            'Predicted': predicted,
            'Confidence': row['Match confidence'],
            'Coverage': row['Coverage'],
            'Expected_genes_missing': row['Missing expected genes'],
            'Other_genes': row['Other genes in locus'],
            'Other_genes_details': row['Other genes in locus, details']
        })

print(f"Found {len(mismatches)} mismatches.")
for m in mismatches:
    print(f"\nAssembly: {m['Assembly']}")
    print(f"Expected: {m['Expected']} | Predicted: {m['Predicted']} ({m['Confidence']})")
    print(f"Coverage: {m['Coverage']}")
    print(f"Missing: {m['Expected_genes_missing']}")
    print(f"Other genes: {m['Other_genes']} - {m['Other_genes_details']}")

