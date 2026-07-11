import csv
import subprocess
import os
import time

def main():
    truth_file = "../ab_mrsn/ab_mrsn_k_locus_truths.tsv"
    assembly_dir = "../ab_mrsn"

    # read ground truth
    ground_truths = {}
    with open(truth_file, 'r') as f:
        reader = csv.DictReader(f, delimiter='\t')
        for row in reader:
            ground_truths[row['Assembly']] = row['Expected_locus']

    print(f"Total assemblies to test: {len(ground_truths)}")

    mismatches = 0
    matches = 0
    missing = 0
    total_time = 0

    mismatch_details = []

    for idx, (assembly, expected_locus) in enumerate(ground_truths.items()):
        assembly_path = os.path.join(assembly_dir, f"{assembly}.fna.gz")
        if not os.path.exists(assembly_path):
            missing += 1
            continue

        start_time = time.time()
        # run kaptive
        result = subprocess.run(
            ["uv", "run", "kaptive", "type", "ab_k", assembly_path],
            capture_output=True, text=True
        )
        end_time = time.time()
        total_time += (end_time - start_time)

        if result.returncode != 0:
            print(f"Error running kaptive on {assembly}: {result.stderr}")
            continue

        # output is tab separated, second line has the data
        lines = result.stdout.strip().split('\n')
        if len(lines) < 2:
            print(f"Unexpected output for {assembly}: {result.stdout}")
            continue

        headers = lines[0].split('\t')
        data = lines[1].split('\t')
        
        # parse best match locus and match confidence
        try:
            locus_idx = headers.index('Best match locus')
            conf_idx = headers.index('Match confidence')
            pred_locus = data[locus_idx]
            pred_conf = data[conf_idx]
        except ValueError as e:
            print(f"Error parsing headers: {headers}")
            continue

        if pred_locus != expected_locus:
            mismatches += 1
            mismatch_details.append(f"{assembly} | Truth: {expected_locus} | Pred: {pred_locus} | Conf: {pred_conf}")
        else:
            matches += 1
        
        if (idx + 1) % 10 == 0:
            print(f"Processed {idx + 1}/{len(ground_truths)}...")

    print(f"\n--- Results ---")
    print(f"Total evaluated: {matches + mismatches}")
    print(f"Matches: {matches}")
    print(f"Mismatches: {mismatches}")
    print(f"Missing assemblies: {missing}")
    print(f"Average time per assembly: {total_time / (matches + mismatches):.2f}s")
    print(f"Total time: {total_time:.2f}s")
    
    if mismatches > 0:
        print("\n--- Mismatches ---")
        for m in mismatch_details:
            print(m)

if __name__ == "__main__":
    main()
