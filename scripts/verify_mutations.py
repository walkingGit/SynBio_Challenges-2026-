#!/usr/bin/env python3
"""Verify that each designed sequence's mutations match the claimed mutation list.

For each sequence, aligns against its parent scaffold (sfGFP, StayGold, or mBaoJin)
and checks that every claimed substitution is present and no unclaimed substitutions exist.
"""

import pandas as pd
from Bio.Align import PairwiseAligner

# ── Parent scaffold sequences ──────────────────────────────────────────
SFGFP = (
    "MSKGEELFTGVVPILVELDGDVNGHKFSVRGEGEGDATNGKLTLKFICTTGKLPVPWPTLVTTLTYGVQCFSRYPDHMKR"
    "HDFFKSAMPEGYVQERTISFKDDGTYKTRAEVKFEGDTLVNRIELKGIDFKEDGNILGHKLEYNFNSHNVYITADKQKNGI"
    "KANFKIRHNVEDGSVQLADHYQQNTPIGDGPVLLPDNHYLSTQSVLSKDPNEKRDHMVLLEFVTAAGITHGMDELYK"
)

STAYGOLD = (
    "MASTPFKFQLKGTINGKSFTVEGEGEGNSHEGSHKGKYVCTSGKLPMSWAALGTSFGYGMKYYTKYPSGLKNWFHEVMPEG"
    "FTYDRHIQYKGDGSIHAKHQHFMKNGTYHNIVEFTGQDFKENSPVLTGDMNVSLPNEVQHIPRDDGVECPVTLLYPLLSDK"
    "SKCVEAHQNTICKPLHNQPAPDVPYHWIRKQYTQSKDDTEERDHICQSETLEAHL"
)

MBAOJIN = (
    "MVSKGEEENMASTPFKFQLKGTINGKSFTVEGEGEGNSHEGSHKGKYVCTSGKLPMSWAALGTTFGYGMKYYTKYPSGLKN"
    "WFREVMPGGFTYDRHIQYKGDGSIHAKHQHFMKNGTYHNIVEFTGQDFKENSPVLTGDMNVSLPNEVPQIPRDDGVECPVT"
    "LLYPLLSDKSKYVEAHQYTICKPLHNQPAPDVPYHWIRKQYTQSKDDAEERDHICQSETLEAHLKGMDELYK"
)

# ── Claimed mutations per sequence ─────────────────────────────────────
# Format: {seq_id: (parent_name, parent_seq, {position: (wt, mut)})}
CLAIMED = {
    1: ("StayGold", STAYGOLD, {
        7: ("K", "E"), 11: ("K", "E"), 17: ("K", "E"), 40: ("C", "S"),
        105: ("K", "E"), 121: ("K", "E"), 135: ("L", "I"), 142: ("I", "L"),
        152: ("V", "L"), 153: ("T", "R"), 162: ("K", "E"), 164: ("K", "E"),
        174: ("C", "S"), 204: ("R", "E"), 208: ("C", "S"),
    }),
    2: ("mBaoJin", MBAOJIN, {
        16: ("K", "E"), 20: ("K", "E"), 26: ("K", "E"), 49: ("C", "S"),
        100: ("K", "E"), 114: ("K", "E"), 130: ("K", "E"), 171: ("K", "E"),
        173: ("K", "E"), 183: ("C", "S"), 213: ("R", "E"), 217: ("C", "S"),
    }),
    3: ("sfGFP", SFGFP, {
        19: ("D", "E"), 46: ("F", "L"), 48: ("C", "S"), 69: ("Q", "L"),
        70: ("C", "A"), 80: ("R", "E"), 101: ("K", "E"), 107: ("K", "E"),
        148: ("H", "S"), 158: ("K", "G"), 164: ("N", "Y"), 166: ("K", "E"),
        171: ("V", "I"), 175: ("S", "T"), 193: ("V", "L"), 198: ("N", "P"),
        203: ("T", "I"), 205: ("S", "T"), 206: ("V", "F"), 214: ("K", "E"),
        215: ("R", "E"), 231: ("H", "P"),
    }),
    4: ("sfGFP", SFGFP, {
        19: ("D", "E"), 46: ("F", "L"), 48: ("C", "S"), 69: ("Q", "L"),
        70: ("C", "A"), 80: ("R", "E"), 101: ("K", "E"), 107: ("K", "E"),
        148: ("H", "S"), 158: ("K", "G"), 164: ("N", "Y"), 166: ("K", "E"),
        171: ("V", "I"), 175: ("S", "T"), 193: ("V", "L"), 198: ("N", "P"),
        203: ("T", "I"), 205: ("S", "T"), 214: ("K", "E"), 215: ("R", "E"),
        231: ("H", "P"),
    }),
    5: ("mBaoJin", MBAOJIN, {
        16: ("K", "E"), 20: ("K", "E"), 26: ("K", "E"), 49: ("C", "S"),
        100: ("K", "E"), 114: ("K", "E"), 115: ("N", "D"), 125: ("T", "P"),
        130: ("K", "E"), 132: ("N", "D"), 159: ("C", "S"), 171: ("K", "E"),
        173: ("K", "E"), 183: ("C", "S"), 188: ("N", "D"), 189: ("Q", "E"),
        193: ("D", "P"), 213: ("R", "E"), 217: ("C", "S"),
    }),
    6: ("sfGFP", SFGFP, {
        19: ("D", "E"), 46: ("F", "L"), 48: ("C", "S"), 69: ("Q", "L"),
        70: ("C", "A"), 80: ("R", "E"), 101: ("K", "E"), 107: ("K", "E"),
        158: ("K", "G"), 164: ("N", "Y"), 166: ("K", "E"), 171: ("V", "I"),
        175: ("S", "T"), 193: ("V", "L"), 198: ("N", "P"), 203: ("T", "I"),
        205: ("S", "T"), 206: ("V", "F"), 214: ("K", "E"), 215: ("R", "E"),
        231: ("H", "P"),
    }),
}


def find_mutations_direct(parent, child):
    """Find substitutions between same-length sequences by direct comparison."""
    assert len(parent) == len(child), f"Length mismatch: {len(parent)} vs {len(child)}"
    mutations = {}
    for i in range(len(parent)):
        if parent[i] != child[i]:
            mutations[i + 1] = (parent[i], child[i])  # 1-indexed
    return mutations


def find_mutations_aligned(parent, child):
    """Find substitutions between different-length sequences via alignment."""
    aligner = PairwiseAligner()
    aligner.mode = "global"
    aligner.match_score = 2
    aligner.mismatch_score = -1
    aligner.open_gap_score = -2
    aligner.extend_gap_score = -0.5

    alignments = aligner.align(parent, child)
    best = alignments[0]
    aln_parent = str(best[0])
    aln_child = str(best[1])

    mutations = {}
    pos = 0
    for i in range(len(aln_parent)):
        if aln_parent[i] != "-":
            pos += 1
        if aln_parent[i] != aln_child[i] and aln_parent[i] != "-" and aln_child[i] != "-":
            mutations[pos] = (aln_parent[i], aln_child[i])
    return mutations


def main():
    submission = pd.read_csv("submission.csv")
    seqs = {row["Seq_ID"]: row["Sequence"] for _, row in submission.iterrows()}

    all_pass = True
    for seq_id in sorted(seqs.keys()):
        parent_name, parent_seq, claimed_dict = CLAIMED[seq_id]
        child_seq = seqs[seq_id]

        if len(parent_seq) == len(child_seq):
            actual = find_mutations_direct(parent_seq, child_seq)
            method = "direct"
        else:
            actual = find_mutations_aligned(parent_seq, child_seq)
            method = "aligned"

        print(f"\nSeq {seq_id} vs {parent_name} ({method})")
        print(f"  Claimed: {len(claimed_dict)} mutations, Actual: {len(actual)} mutations")

        seq_ok = True
        # Check each claimed mutation
        for pos in sorted(claimed_dict.keys()):
            wt_c, mut_c = claimed_dict[pos]
            if pos in actual:
                wt_a, mut_a = actual[pos]
                if wt_c == wt_a and mut_c == mut_a:
                    print(f"  OK  {wt_c}{pos}{mut_c}")
                else:
                    print(f"  ERR {wt_c}{pos}{mut_c} -> actual {wt_a}{pos}{mut_a}")
                    seq_ok = False
            else:
                print(f"  MISS {wt_c}{pos}{mut_c} not found in actual")
                seq_ok = False

        # Check for unclaimed mutations
        unclaimed = set(actual.keys()) - set(claimed_dict.keys())
        if unclaimed:
            for pos in sorted(unclaimed):
                wt_a, mut_a = actual[pos]
                print(f"  UNCLAIMED {wt_a}{pos}{mut_a}")
            seq_ok = False

        # Check metadata
        actual_cys = child_seq.count("C")
        neg = child_seq.count("D") + child_seq.count("E")
        pos_charge = child_seq.count("K") + child_seq.count("R") + child_seq.count("H")
        net_charge = pos_charge - neg

        print(f"  Length: {len(child_seq)}, Cys: {actual_cys}, Net charge: {net_charge}")

        if seq_ok:
            print(f"  PASS")
        else:
            print(f"  FAIL")
            all_pass = False

    print(f"\n{'ALL SEQUENCES VERIFIED' if all_pass else 'DISCREPANCIES FOUND'}")
    return all_pass


if __name__ == "__main__":
    main()
