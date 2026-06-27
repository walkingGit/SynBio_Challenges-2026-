#!/usr/bin/env python3
"""Validate submission.csv against competition format requirements.

Competition rules:
- 220-250 amino acids per sequence
- Starts with M
- Uppercase standard 20 amino acids only
- No stop codons
- Must differ from Exclusion_List.csv by >6 AA mismatches
"""

import pandas as pd
import os

VALID_AAS = set("ACDEFGHIKLMNPQRSTVWY")
MIN_LENGTH = 220
MAX_LENGTH = 250


def validate_sequence(seq, seq_id):
    """Validate a single sequence against format rules. Returns list of errors."""
    errors = []

    if len(seq) < MIN_LENGTH or len(seq) > MAX_LENGTH:
        errors.append(f"Length {len(seq)} outside range [{MIN_LENGTH}, {MAX_LENGTH}]")

    if not seq.startswith("M"):
        errors.append("Does not start with M")

    if seq != seq.upper():
        errors.append("Contains lowercase characters")

    invalid = set(seq) - VALID_AAS
    if invalid:
        errors.append(f"Invalid amino acids: {invalid}")

    if "*" in seq:
        errors.append("Contains stop codon (*)")

    return errors


def check_exclusion(submission_seqs, exclusion_path=None):
    """Check that each sequence differs from all exclusion list entries by >6 AAs."""
    if exclusion_path is None or not os.path.exists(exclusion_path):
        print("  Exclusion list not found — skipping exclusion check")
        return True

    exclusion = pd.read_csv(exclusion_path)
    all_pass = True

    for seq_id, seq in submission_seqs.items():
        min_mismatches = len(seq)  # worst case
        for _, row in exclusion.iterrows():
            exc_seq = row.iloc[0] if len(row) == 1 else row["Sequence"]
            if len(exc_seq) != len(seq):
                continue
            mismatches = sum(1 for a, b in zip(seq, exc_seq) if a != b)
            min_mismatches = min(min_mismatches, mismatches)

        status = "PASS" if min_mismatches > 6 else "FAIL"
        print(f"  Seq {seq_id}: min {min_mismatches} mismatches -> {status}")
        if min_mismatches <= 6:
            all_pass = False

    return all_pass


def main():
    submission = pd.read_csv("submission.csv")

    # Check required columns
    required_cols = {"Team_Name", "Seq_ID", "Sequence"}
    missing = required_cols - set(submission.columns)
    if missing:
        print(f"ERROR: Missing columns: {missing}")
        return False

    print(f"Submission: {len(submission)} sequences, team '{submission['Team_Name'].iloc[0]}'")
    print()

    all_pass = True
    for _, row in submission.iterrows():
        seq_id = row["Seq_ID"]
        seq = row["Sequence"]
        errors = validate_sequence(seq, seq_id)

        cys = seq.count("C")
        neg = seq.count("D") + seq.count("E")
        pos = seq.count("K") + seq.count("R") + seq.count("H")
        charge = pos - neg

        if errors:
            print(f"Seq {seq_id}: FAIL")
            for e in errors:
                print(f"  - {e}")
            all_pass = False
        else:
            print(f"Seq {seq_id}: PASS ({len(seq)}aa, {cys} Cys, charge {charge})")

    # Exclusion check
    print("\nExclusion list check:")
    seqs = {row["Seq_ID"]: row["Sequence"] for _, row in submission.iterrows()}
    exclusion_path = os.path.join("data", "Exclusion_List.csv")
    if not check_exclusion(seqs, exclusion_path):
        all_pass = False

    print(f"\n{'ALL CHECKS PASSED' if all_pass else 'SOME CHECKS FAILED'}")
    return all_pass


if __name__ == "__main__":
    main()
