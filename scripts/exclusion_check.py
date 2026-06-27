#!/usr/bin/env python3
"""
Exclusion List Checker
======================
Checks designed GFP sequences against the competition exclusion list
to ensure they differ by >6 amino acids from all excluded sequences.

Competition rule: submitted sequences must differ by more than 6 amino acids
from ALL sequences in the Exclusion_List.csv.

Usage:
    python exclusion_check.py --submission submission.csv --exclusion Exclusion_List.csv
    python exclusion_check.py --sequence "MSKGEEL..." --exclusion Exclusion_List.csv

Requires: pandas, numpy
"""

import argparse
import os
import sys
import pandas as pd
import numpy as np


def load_exclusion_list(exclusion_path):
    """Load exclusion list CSV and return array of sequences."""
    df = pd.read_csv(exclusion_path)
    # The first column contains sequences
    seq_col = df.columns[0]
    sequences = df[seq_col].astype(str).str.strip().values
    print(f"Exclusion list loaded: {len(sequences)} sequences")
    return sequences


def count_mismatches(seq1, seq2):
    """Count amino acid mismatches between two sequences of equal length."""
    if len(seq1) != len(seq2):
        return float('inf')
    return sum(1 for a, b in zip(seq1, seq2) if a != b)


def fast_exclusion_check(designed_seq, exclusion_sequences, min_mismatches=7):
    """
    Fast exclusion check using vectorized comparison.
    
    Args:
        designed_seq: designed protein sequence string
        exclusion_sequences: array of exclusion sequence strings
        min_mismatches: minimum required mismatches (default 7, i.e. >6)
    
    Returns:
        (passed: bool, min_mismatches_found: int)
    """
    designed = np.array(list(designed_seq))
    global_min = float('inf')
    
    for excl_seq in exclusion_sequences:
        if len(excl_seq) != len(designed_seq):
            continue
        excl_arr = np.array(list(excl_seq))
        mismatches = int(np.sum(designed != excl_arr))
        if mismatches < global_min:
            global_min = mismatches
        if mismatches < min_mismatches:
            return False, mismatches
    
    return True, int(global_min)


def check_single_sequence(seq, seq_name, exclusion_sequences, min_mismatches=7):
    """Check a single sequence against the exclusion list."""
    # Exact match check
    exact_match = seq in set(exclusion_sequences)
    if exact_match:
        print(f"  {seq_name}: EXCLUDED (exact match in exclusion list)")
        return False, 0
    
    # Minimum mismatch check (same-length sequences only)
    same_len = [s for s in exclusion_sequences if len(s) == len(seq)]
    if not same_len:
        print(f"  {seq_name}: PASS (no exclusion sequences of same length {len(seq)})")
        return True, float('inf')
    
    passed, min_mm = fast_exclusion_check(seq, same_len, min_mismatches)
    
    if passed:
        print(f"  {seq_name}: PASS (min mismatches = {min_mm}, threshold > {min_mismatches-1})")
    else:
        print(f"  {seq_name}: FAIL (min mismatches = {min_mm}, need > {min_mismatches-1})")
    
    return passed, min_mm


def check_submission_csv(submission_path, exclusion_path, min_mismatches=7):
    """
    Check all sequences in a submission CSV against the exclusion list.
    
    Submission CSV format: Seq_ID, Sequence, Description
    """
    excl_seqs = load_exclusion_list(exclusion_path)
    sub_df = pd.read_csv(submission_path)
    
    print(f"\n{'='*60}")
    print(f"EXCLUSION CHECK FOR {len(sub_df)} SEQUENCES")
    print(f"Threshold: > {min_mismatches-1} mismatches required")
    print(f"{'='*60}\n")
    
    results = []
    all_pass = True
    
    for _, row in sub_df.iterrows():
        seq_id = row.get('Seq_ID', row.iloc[0])
        seq = row.get('Sequence', row.iloc[1])
        
        passed, min_mm = check_single_sequence(seq, f"Seq {seq_id}", excl_seqs, min_mismatches)
        results.append({
            'seq_id': seq_id,
            'sequence': seq,
            'passed': passed,
            'min_mismatches': min_mm,
            'length': len(seq)
        })
        if not passed:
            all_pass = False
    
    print(f"\n{'='*60}")
    if all_pass:
        print("ALL SEQUENCES PASS EXCLUSION CHECK")
    else:
        print("WARNING: SOME SEQUENCES FAIL EXCLUSION CHECK")
    print(f"{'='*60}")
    
    return pd.DataFrame(results), all_pass


def find_closest_matches(designed_seq, exclusion_sequences, top_n=5):
    """Find the closest matching sequences in the exclusion list."""
    same_len = [(s, count_mismatches(designed_seq, s)) 
                for s in exclusion_sequences if len(s) == len(designed_seq)]
    same_len.sort(key=lambda x: x[1])
    
    print(f"\nClosest {top_n} exclusion matches:")
    for s, mm in same_len[:top_n]:
        # Show which positions differ
        diff_pos = [(i+1, designed_seq[i], s[i]) 
                    for i in range(len(designed_seq)) if designed_seq[i] != s[i]]
        print(f"  Mismatches: {mm}")
        print(f"  Positions: {diff_pos[:15]}{'...' if len(diff_pos) > 15 else ''}")
    
    return same_len[:top_n]


# ── Main ─────────────────────────────────────────────────────────────
def main():
    parser = argparse.ArgumentParser(description="Check sequences against competition exclusion list")
    parser.add_argument("--submission", help="Submission CSV file (Seq_ID, Sequence, ...)")
    parser.add_argument("--sequence", help="Single sequence string to check")
    parser.add_argument("--name", default="designed", help="Name for single sequence")
    parser.add_argument("--exclusion", required=True, help="Exclusion_List.csv path")
    parser.add_argument("--threshold", type=int, default=7,
                        help="Minimum mismatches required (default 7, meaning >6)")
    parser.add_argument("--show-closest", action="store_true",
                        help="Show closest matching exclusion sequences")
    parser.add_argument("--output", default="exclusion_results.csv", help="Output CSV path")
    args = parser.parse_args()
    
    if args.submission:
        results_df, all_pass = check_submission_csv(args.submission, args.exclusion, args.threshold)
        results_df.to_csv(args.output, index=False)
        print(f"\nResults saved to {args.output}")
    
    elif args.sequence:
        excl_seqs = load_exclusion_list(args.exclusion)
        passed, min_mm = check_single_sequence(args.sequence, args.name, excl_seqs, args.threshold)
        
        if args.show_closest:
            find_closest_matches(args.sequence, excl_seqs)
    
    else:
        print("Provide --submission or --sequence")
        sys.exit(1)


if __name__ == "__main__":
    main()
