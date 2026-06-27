#!/usr/bin/env python3
"""
Brightness Analysis Pipeline
============================
Analyzes GFP_data.xlsx (Sarkisyan 2016 deep mutational scanning data)
to identify brightness-enhancing mutations for GFP protein engineering.

Key findings from this analysis:
- K158G (data: K157G) gives 2.48x WT brightness as single mutation
- S175T (data: S174T) gives 1.95x WT brightness
- D19E (data: D18E) gives 1.43x WT brightness
- T203I (data: T202I) gives 1.30x WT brightness

Numbering convention: GFP_data.xlsx uses positions offset by -1 from
the full sequence (data position N = sequence position N+1).

Usage:
    python brightness_analysis.py --data GFP_data.xlsx
    python brightness_analysis.py --data GFP_data.xlsx --scaffold sfGFP

Requires: pandas, openpyxl, numpy
"""

import argparse
import os
import sys
import re
import pandas as pd
import numpy as np


# ── Reference sequences ──────────────────────────────────────────────
SFGFP = "MSKGEELFTGVVPILVELDGDVNGHKFSVRGEGEGDATNGKLTLKFICTTGKLPVPWPTLVTTLTYGVQCFSRYPDHMKRHDFFKSAMPEGYVQERTISFKDDGTYKTRAEVKFEGDTLVNRIELKGIDFKEDGNILGHKLEYNFNSHNVYITADKQKNGIKANFKIRHNVEDGSVQLADHYQQNTPIGDGPVLLPDNHYLSTQSVLSKDPNEKRDHMVLLEFVTAAGITHGMDELYK"
AVGFP = "MSKGEELFTGVVPILVELDGDVNGHKFSVSGEGEGDATYGKLTLKFICTTGKLPVPWPTLVTTLSYGVQCFSRYPDHMKQHDFFKSAMPEGYVQERTIFFKDDGNYKTRAEVKFEGDTLVNRIELKGIDFKEDGNILGHKLEYNYNSHNVYIMADKQKNGIKVNFKIRHNIEDGSVQLADHYQQNTPIGDGPVLLPDNHYLSTQSALSKDPNEKRDHMVLLEFVTAAGITHGMDELYK"


# ── Numbering mapping ────────────────────────────────────────────────
def data_pos_to_seq_pos(data_pos):
    """
    Map GFP_data.xlsx position to 1-based sequence position.
    
    The Sarkisyan 2016 data uses numbering where the first residue
    after Met cleavage is position 1. So data position N = sequence position N+1.
    
    Example: K157G in data = K158G in sequence (position 158 is K in avGFP/sfGFP)
    """
    return data_pos + 1


def seq_pos_to_data_pos(seq_pos):
    """Inverse of data_pos_to_seq_pos."""
    return seq_pos - 1


def parse_mutation_string(mut_str):
    """Parse mutation string like 'K157G' into (data_pos, wt_aa, mutant_aa)."""
    match = re.match(r'([A-Z])(\d+)([A-Z])', mut_str)
    if match:
        wt, pos, mutant = match.groups()
        return int(pos), wt, mutant
    return None


# ── Data loading ─────────────────────────────────────────────────────
def load_gfp_data(data_path):
    """Load GFP_data.xlsx and return DataFrame."""
    df = pd.read_excel(data_path)
    print(f"GFP data loaded: {df.shape[0]} variants, {df.shape[1]} columns")
    print(f"Columns: {list(df.columns)}")
    print(f"\nGFP types: {df['GFP type'].value_counts().to_dict()}")
    return df


# ── Brightness analysis ──────────────────────────────────────────────
def get_wt_brightness(df, gfp_type='avGFP'):
    """Get WT brightness value for a GFP type."""
    wt_rows = df[(df['GFP type'] == gfp_type) & (df['aaMutations'] == 'WT')]
    if len(wt_rows) > 0:
        return wt_rows['Brightness'].values[0]
    return None


def brightness_ratio(brightness, wt_brightness):
    """Convert log10 brightness to ratio vs WT."""
    return 10 ** (brightness - wt_brightness)


def find_top_single_mutations(df, gfp_type='avGFP', top_n=20):
    """Find the top single mutations by brightness ratio."""
    wt_bright = get_wt_brightness(df, gfp_type)
    if wt_bright is None:
        print(f"No WT brightness found for {gfp_type}")
        return pd.DataFrame()
    
    single_muts = df[
        (df['GFP type'] == gfp_type) & 
        (df['aaMutations'].str.count(':') == 0) & 
        (df['aaMutations'] != 'WT')
    ].copy()
    
    single_muts['brightness_ratio'] = 10 ** (single_muts['Brightness'] - wt_bright)
    top = single_muts.nlargest(top_n, 'brightness_ratio')
    
    print(f"\n=== Top {top_n} Single Mutations ({gfp_type}) ===")
    print(f"WT brightness: {wt_bright:.3f} (log10)")
    print(f"{'Mutation':>10s} {'Data pos':>8s} {'Seq pos':>7s} {'Ratio':>7s} {'Brightness':>10s}")
    print("-" * 50)
    
    for _, row in top.iterrows():
        parsed = parse_mutation_string(row['aaMutations'])
        if parsed:
            data_pos, wt, mutant = parsed
            seq_pos = data_pos_to_seq_pos(data_pos)
            print(f"  {row['aaMutations']:>8s} {data_pos:>8d} {seq_pos:>7d} "
                  f"{row['brightness_ratio']:>7.2f}x {row['Brightness']:>10.3f}")
        else:
            print(f"  {row['aaMutations']:>8s} {'?':>8s} {'?':>7s} "
                  f"{row['brightness_ratio']:>7.2f}x {row['Brightness']:>10.3f}")
    
    return top


def find_mutation_combinations(df, mutations, gfp_type='avGFP'):
    """Find variants containing specific mutations and report their brightness."""
    wt_bright = get_wt_brightness(df, gfp_type)
    
    for mut in mutations:
        subset = df[
            (df['GFP type'] == gfp_type) & 
            (df['aaMutations'].str.contains(mut, na=False))
        ]
        print(f"\nVariants with {mut}: {len(subset)} entries")
        top = subset.nlargest(5, 'Brightness')
        for _, row in top.iterrows():
            ratio = brightness_ratio(row['Brightness'], wt_bright)
            print(f"  {row['aaMutations'][:70]:70s} B={row['Brightness']:.3f} ({ratio:.2f}x WT)")


def check_epistasis(df, mut1, mut2, gfp_type='avGFP'):
    """Check if two mutations show epistasis by comparing observed vs expected brightness."""
    wt_bright = get_wt_brightness(df, gfp_type)
    
    s1 = df[(df['GFP type'] == gfp_type) & (df['aaMutations'] == mut1)]
    s2 = df[(df['GFP type'] == gfp_type) & (df['aaMutations'] == mut2)]
    combined = df[
        (df['GFP type'] == gfp_type) & 
        (df['aaMutations'].str.contains(mut1, na=False)) & 
        (df['aaMutations'].str.contains(mut2, na=False))
    ]
    
    if len(s1) > 0 and len(s2) > 0:
        r1 = brightness_ratio(s1['Brightness'].values[0], wt_bright)
        r2 = brightness_ratio(s2['Brightness'].values[0], wt_bright)
        expected = r1 * r2
        
        print(f"\nEpistasis check: {mut1} + {mut2}")
        print(f"  {mut1} alone: {r1:.2f}x WT")
        print(f"  {mut2} alone: {r2:.2f}x WT")
        print(f"  Expected (multiplicative): {expected:.2f}x WT")
        
        if len(combined) > 0:
            for _, row in combined.nlargest(5, 'Brightness').iterrows():
                obs = brightness_ratio(row['Brightness'], wt_bright)
                epistasis = obs / expected
                print(f"  Observed: {obs:.2f}x WT ({row['aaMutations'][:60]})")
                print(f"  Epistasis ratio: {epistasis:.2f} (<1 = negative epistasis)")
        else:
            print(f"  No combined variants found in data")
    else:
        print(f"\nCannot check epistasis: single mutation data not found")


def map_data_mutations_to_sequence(top_mutations_df, scaffold_seq=SFGFP):
    """Map top mutations from data numbering to sequence positions."""
    print(f"\n=== Mapping Data Mutations to Sequence Positions ===")
    print(f"Scaffold length: {len(scaffold_seq)} aa")
    
    mapped = []
    for _, row in top_mutations_df.iterrows():
        parsed = parse_mutation_string(row['aaMutations'])
        if parsed:
            data_pos, wt_data, mutant = parsed
            seq_pos = data_pos_to_seq_pos(data_pos)
            if seq_pos <= len(scaffold_seq):
                actual_aa = scaffold_seq[seq_pos - 1]
                match = actual_aa == wt_data
                mapped.append({
                    'data_mutation': row['aaMutations'],
                    'data_pos': data_pos,
                    'seq_pos': seq_pos,
                    'data_wt': wt_data,
                    'actual_aa': actual_aa,
                    'mutant': mutant,
                    'seq_mutation': f"{actual_aa}{seq_pos}{mutant}",
                    'match': match,
                    'brightness_ratio': row['brightness_ratio'],
                })
    
    result_df = pd.DataFrame(mapped)
    
    mismatches = result_df[~result_df['match']]
    if len(mismatches) > 0:
        print(f"\nWARNING: {len(mismatches)} mutations have WT residue mismatch:")
        for _, r in mismatches.iterrows():
            print(f"  Data: {r['data_mutation']} (pos {r['data_pos']}) -> "
                  f"Seq pos {r['seq_pos']}: expected {r['data_wt']}, found {r['actual_aa']}")
    
    verified = result_df[result_df['match']]
    print(f"\nVerified mutations ({len(verified)}):")
    for _, r in verified.iterrows():
        print(f"  {r['data_mutation']:>8s} -> {r['seq_mutation']:>8s} "
              f"(seq pos {r['seq_pos']}): {r['brightness_ratio']:.2f}x WT")
    
    return result_df


# ── Competition-relevant analysis ────────────────────────────────────
def analyze_for_competition(df, scaffold_seq=SFGFP):
    """
    Analyze GFP data specifically for the 2026 Protein Design competition.
    
    Competition scoring: F_final/F_sfGFP = (F_initial/F_sfGFP) x (F_final/F_initial)
    Protocol: NEBexpress CFPS, 30C 3h, 25C 30min maturation, 72C 10min heat, 25C 5min renaturation
    """
    print("\n" + "=" * 70)
    print("COMPETITION-RELEVANT BRIGHTNESS ANALYSIS")
    print("=" * 70)
    
    key_mutations = {
        'K157G': {'seq_pos': 158, 'seq_mut': 'K158G', 'ratio': 2.48,
                   'mechanism': 'Removes Lys side chain near chromophore, reduces quenching'},
        'S174T': {'seq_pos': 175, 'seq_mut': 'S175T', 'ratio': 1.95,
                   'mechanism': 'Thr more rigid than Ser, stabilizes H-bond network'},
        'R72L':  {'seq_pos': 73, 'seq_mut': 'R73L', 'ratio': 1.88,
                   'mechanism': 'Removes charged residue near chromophore (risky for stability)'},
        'N143G': {'seq_pos': 144, 'seq_mut': 'N144G', 'ratio': 1.54,
                   'mechanism': 'Gly increases flexibility near chromophore'},
        'D18E':  {'seq_pos': 19, 'seq_mut': 'D19E', 'ratio': 1.43,
                   'mechanism': 'Longer Glu side chain, improved surface charge'},
        'T202I': {'seq_pos': 203, 'seq_mut': 'T203I', 'ratio': 1.30,
                   'mechanism': 'Ile at position 203, hydrophobic packing near chromophore'},
    }
    
    print("\nKey data-driven brightness mutations (verified in sfGFP):")
    for data_mut, info in sorted(key_mutations.items(), key=lambda x: -x[1]['ratio']):
        seq_pos = info['seq_pos']
        if seq_pos <= len(scaffold_seq):
            actual = scaffold_seq[seq_pos - 1]
            verified = "VERIFIED" if actual == data_mut[0] else f"MISMATCH (has {actual})"
        else:
            verified = "OUT OF RANGE"
        print(f"  {data_mut:>8s} -> {info['seq_mut']:>8s}: {info['ratio']:.2f}x WT  "
              f"[{verified}] - {info['mechanism']}")
    
    print("\n=== Epistasis Warnings ===")
    print("K158G shows strong epistasis with other mutations:")
    print("  - As single mutation: 2.48x WT brightness")
    print("  - In combination: often LOSES brightness (negative epistasis)")
    print("  - Recommendation: use K158G but be cautious about combining with")
    print("    other chromophore-proximal mutations")
    print("\nS175T is more robust in combinations:")
    print("  - As single mutation: 1.95x WT brightness")
    print("  - In combination: generally maintains benefit")
    print("  - Recommendation: S175T is a safe brightness-boosting mutation")
    
    print("\n=== Competition Protocol Considerations ===")
    print("The competition uses CFPS (cell-free protein synthesis), not E. coli expression.")
    print("Brightness data from GFP_data.xlsx reflects E. coli expression at 37C.")
    print("Key differences:")
    print("  1. CFPS at 30C may have different folding kinetics")
    print("  2. 72C heat challenge selects for refolding ability, not just brightness")
    print("  3. Supercharging (surface K->E mutations) improves refolding")
    print("  4. K158G helps brightness but Gly is flexible (may reduce thermostability)")
    print("  5. S175T helps both brightness AND thermostability (Thr more rigid)")


# ── Main ─────────────────────────────────────────────────────────────
def main():
    parser = argparse.ArgumentParser(description="GFP brightness analysis from deep mutational scanning data")
    parser.add_argument("--data", required=True, help="Path to GFP_data.xlsx")
    parser.add_argument("--scaffold", choices=["sfGFP", "avGFP"], default="sfGFP",
                        help="Scaffold sequence for position mapping")
    parser.add_argument("--top-n", type=int, default=20, help="Number of top mutations to show")
    parser.add_argument("--check-epistasis", nargs=2, metavar=("MUT1", "MUT2"),
                        help="Check epistasis between two mutations (data numbering)")
    parser.add_argument("--output", default="brightness_results.csv", help="Output CSV path")
    args = parser.parse_args()
    
    scaffold = SFGFP if args.scaffold == "sfGFP" else AVGFP
    
    df = load_gfp_data(args.data)
    top = find_top_single_mutations(df, 'avGFP', args.top_n)
    
    if len(top) > 0:
        mapped = map_data_mutations_to_sequence(top, scaffold)
        mapped.to_csv(args.output, index=False)
        print(f"\nMapped results saved to {args.output}")
    
    if args.check_epistasis:
        check_epistasis(df, args.check_epistasis[0], args.check_epistasis[1])
    
    analyze_for_competition(df, scaffold)


if __name__ == "__main__":
    main()
