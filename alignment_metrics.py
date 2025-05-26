import argparse
import os
import glob
import math
from Bio import AlignIO
from Bio.Align.substitution_matrices import load as load_matrix
import pandas as pd
import matplotlib.pyplot as plt



def parse_alignment(filepath, fmt="fasta"):
    return AlignIO.read(filepath, fmt)


def compute_sp_score(alignment, seq_type='nt'):
    num_seqs = len(alignment)
    L = alignment.get_alignment_length()
    num_pairs = num_seqs * (num_seqs - 1) / 2

    if seq_type == 'aa':
        matrix = load_matrix("BLOSUM62")
        alphabet = matrix.alphabet
        max_score = max(matrix[a, b] for a in alphabet for b in alphabet)
    else:
        matrix = None
        max_score = 1

    total_score = 0
    for col_idx in range(L):
        col = [record.seq[col_idx] for record in alignment]
        for i in range(num_seqs):
            for j in range(i+1, num_seqs):
                a = col[i]
                b = col[j]
                if seq_type == 'nt':
                    score = 1 if a == b and a != '-' else 0
                else:
                    a_up = a.upper()
                    b_up = b.upper()
                    if a_up == '-' or b_up == '-':
                        score = 0
                    else:
                        try:
                            score = matrix[a_up, b_up]
                        except KeyError:
                            try:
                                score = matrix[b_up, a_up]
                            except KeyError:
                                score = 0
                total_score += score

    return total_score / (num_pairs * L * max_score)


def compute_tc_score(alignment):
    L = alignment.get_alignment_length()
    count = 0
    for col_idx in range(L):
        col = [res for res in alignment[:, col_idx] if res != '-']
        if col and len(set(col)) == 1:
            count += 1
    return count / L


def compute_cs_score(alignment, threshold=0.8):
    num_seqs = len(alignment)
    L = alignment.get_alignment_length()
    count = 0
    for col_idx in range(L):
        col = [res for res in alignment[:, col_idx] if res != '-']
        if not col:
            continue
        freqs = {}
        for r in col:
            freqs[r] = freqs.get(r, 0) + 1
        if max(freqs.values()) / num_seqs >= threshold:
            count += 1
    return count / L


def compute_shannon_entropy(alignment):
    L = alignment.get_alignment_length()
    entropies = []
    for col_idx in range(L):
        col = [res for res in alignment[:, col_idx] if res != '-']
        if not col:
            entropies.append(0)
            continue
        freqs = {}
        for r in col:
            freqs[r] = freqs.get(r, 0) + 1
        H = 0
        length = len(col)
        for count in freqs.values():
            p = count / length
            H -= p * math.log2(p)
        entropies.append(H)
    return sum(entropies) / L


def list_alignment_files(directory, extensions=(".fa", ".fna", ".fasta")):
    files = []
    for ext in extensions:
        files.extend(glob.glob(os.path.join(directory, f"*{ext}")))
    return sorted(files)


def extract_algorithm_name(filepath):
    name = os.path.basename(filepath).lower()
    for algo in ['mafft', 'muscle', 'clustalo', 'clustalw', 'kalign']:
        if algo in name:
            return algo.capitalize()
    return os.path.splitext(os.path.basename(filepath))[0]


def main():
    parser = argparse.ArgumentParser(description="Compute alignment metrics, save CSV and heatmaps.")
    parser.add_argument("-d", "--directory", required=True, help="Directory with alignment files.")
    parser.add_argument("-t", "--type", choices=['nt','aa'], required=True, help="Sequence type: 'nt' or 'aa'.")
    parser.add_argument("-x", "--threshold", type=float, default=0.8, help="Threshold for CS score.")
    parser.add_argument("-o", "--output", help="Output CSV file (optional). If set, saves pivot table as CSV.")
    parser.add_argument("--imgdir", help="Directory to save heatmaps (optional). Creates SP and TC heatmaps.")
    args = parser.parse_args()

    files = list_alignment_files(args.directory)
    results = []
    for filepath in files:
        alignment = parse_alignment(filepath)
        sp = compute_sp_score(alignment, seq_type=args.type)
        tc = compute_tc_score(alignment)
        cs = compute_cs_score(alignment, threshold=args.threshold)
        ent = compute_shannon_entropy(alignment)
        algo = extract_algorithm_name(filepath)
        results.append({'algorithm': algo, 'SP': sp, 'TC': tc, 'CS': cs, 'Entropy': ent})

    df = pd.DataFrame(results)
    pivot = df.set_index('algorithm').T

    # Console output
    print("Metrics vs Algorithms:\n")
    print(pivot.to_string(float_format="%.4f"))

    # Save CSV
    if args.output:
        pivot.to_csv(args.output, sep=';', float_format='%.4f')
        print(f"\nResults written to {args.output} (semicolon-separated for Excel compatibility)")

    # Generate heatmaps
    if args.imgdir:
        os.makedirs(args.imgdir, exist_ok=True)
        sample_name = os.path.basename(os.path.normpath(args.directory))
        # SP heatmap
        sp_vals = pivot.loc[['SP']]
        plt.figure()
        plt.imshow(sp_vals.values, aspect='auto')
        cbar = plt.colorbar() 
        plt.xticks(range(len(sp_vals.columns)), sp_vals.columns, rotation=45, ha='right')
        plt.yticks([0], ['SP'])
        plt.title(f"SP-score: {sample_name}")
        plt.tight_layout()
        sp_path = os.path.join(args.imgdir, f"{sample_name}_SP_heatmap.png")
        plt.savefig(sp_path)
        plt.close()
        print(f"SP heatmap saved to {sp_path}")

        # TC heatmap
        tc_vals = pivot.loc[['TC']]
        plt.figure()
        plt.imshow(tc_vals.values, aspect='auto')
        cbar = plt.colorbar() 
        plt.xticks(range(len(tc_vals.columns)), tc_vals.columns, rotation=45, ha='right')
        plt.yticks([0], ['TC'])
        plt.title(f"TC-score: {sample_name}")
        plt.tight_layout()
        tc_path = os.path.join(args.imgdir, f"{sample_name}_TC_heatmap.png")
        plt.savefig(tc_path)
        plt.close()
        print(f"TC heatmap saved to {tc_path}")

if __name__ == '__main__':
    main()
