#!/usr/bin/env python3
"""
Supercharging Design Pipeline
=============================
Designs surface charge mutations (K/R -> E) for GFP proteins to improve
refolding after heat challenge, based on SASA analysis and structural data.

Key principles:
- Surface-exposed K/R residues (RSA > 30%) are candidates for K->E or R->E
- Interface residues (in dimer designs) are excluded from mutation
- Chromophore-proximal residues are excluded
- Cysteines are removed (C->S or C->A) for CFPS compatibility
- Net charge target: -7 to -13 (more negative = better refolding)

Usage:
    python supercharging_design.py --pdb 2B3P.pdb --scaffold sfGFP
    python supercharging_design.py --pdb 8BXT.pdb --scaffold StayGold --dimer

Requires: biopython (Bio.PDB, ShrakeRupley), numpy, pandas
"""

import argparse
import os
import sys
import numpy as np


# ── Reference sequences ──────────────────────────────────────────────
SFGFP = "MSKGEELFTGVVPILVELDGDVNGHKFSVRGEGEGDATNGKLTLKFICTTGKLPVPWPTLVTTLTYGVQCFSRYPDHMKRHDFFKSAMPEGYVQERTISFKDDGTYKTRAEVKFEGDTLVNRIELKGIDFKEDGNILGHKLEYNFNSHNVYITADKQKNGIKANFKIRHNVEDGSVQLADHYQQNTPIGDGPVLLPDNHYLSTQSVLSKDPNEKRDHMVLLEFVTAAGITHGMDELYK"
STAYGOLD = "MASTPFKFQLKGTINGKSFTVEGEGEGNSHEGSHKGKYVCTSGKLPMSWAALGTSFGYGMKYYTKYPSGLKNWFHEVMPEGFTYDRHIQYKGDGSIHAKHQHFMKNGTYHNIVEFTGQDFKENSPVLTGDMNVSLPNEVQHIPRDDGVECPVTLLYPLLSDKSKCVEAHQNTICKPLHNQPAPDVPYHWIRKQYTQSKDDTEERDHICQSETLEAHL"
MBAOJIN = "MVSKGEEENMASTPFKFQLKGTINGKSFTVEGEGEGNSHEGSHKGKYVCTSGKLPMSWAALGTTFGYGMKYYTKYPSGLKNWFREVMPGGFTYDRHIQYKGDGSIHAKHQHFMKNGTYHNIVEFTGQDFKENSPVLTGDMNVSLPNEVPQIPRDDGVECPVTLLYPLLSDKSKYVEAHQYTICKPLHNQPAPDVPYHWIRKQYTQSKDDAEERDHICQSETLEAHLKGMDELYK"

SCAFFOLDS = {
    'sfGFP': SFGFP,
    'StayGold': STAYGOLD,
    'mBaoJin': MBAOJIN,
}


# ── SASA calculation ─────────────────────────────────────────────────
# Max SASA values for RSA calculation (Tien et al. 2013)
MAX_SASA = {
    'ALA': 129.0, 'ARG': 274.0, 'ASN': 195.0, 'ASP': 193.0,
    'CYS': 167.0, 'GLN': 225.0, 'GLU': 223.0, 'GLY': 104.0,
    'HIS': 224.0, 'ILE': 182.0, 'LEU': 180.0, 'LYS': 236.0,
    'MET': 204.0, 'PHE': 218.0, 'PRO': 159.0, 'SER': 155.0,
    'THR': 172.0, 'TRP': 259.0, 'TYR': 229.0, 'VAL': 164.0
}

THREE_TO_ONE = {
    'ALA':'A','CYS':'C','ASP':'D','GLU':'E','PHE':'F','GLY':'G',
    'HIS':'H','ILE':'I','LYS':'K','LEU':'L','MET':'M','ASN':'N',
    'PRO':'P','GLN':'Q','ARG':'R','SER':'S','THR':'T','VAL':'V',
    'TRP':'W','TYR':'Y'
}


def calc_sasa_from_pdb(pdb_path, chain_id='A'):
    """
    Calculate per-residue SASA using Shrake-Rupley algorithm.
    
    Returns list of dicts with: pos, resname, aa, sasa, rsa, is_surface
    """
    from Bio.PDB import PDBParser, ShrakeRupley
    
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure('protein', pdb_path)
    
    sr = ShrakeRupley()
    sr.compute(structure, level='R')
    
    model = structure[0]
    chain = model[chain_id]
    
    results = []
    for res in chain:
        if res.id[0] != ' ':
            continue
        sasa = res.sasa
        aa = res.resname
        pos = res.id[1]
        rsa = (sasa / MAX_SASA.get(aa, 200.0)) * 100 if aa in MAX_SASA else 0
        one_letter = THREE_TO_ONE.get(aa, 'X')
        
        results.append({
            'pos': pos,
            'resname': aa,
            'aa': one_letter,
            'sasa': sasa,
            'rsa': rsa,
            'is_surface': rsa > 30  # RSA > 30% = surface exposed
        })
    
    return results


# ── Mutation design ──────────────────────────────────────────────────
def apply_mut(seq, mutations):
    """Apply mutations to sequence. mutations = list of (position, new_aa) using 1-based indexing."""
    seq_list = list(seq)
    for pos, new_aa in mutations:
        if 1 <= pos <= len(seq_list):
            old_aa = seq_list[pos - 1]
            seq_list[pos - 1] = new_aa
        else:
            print(f"  WARNING: position {pos} out of range (seq length {len(seq)})")
    return ''.join(seq_list)


def net_charge(seq):
    """Calculate net charge at pH 7: K(+1), R(+1), H(+0.5), D(-1), E(-1)."""
    pos = seq.count('K') + seq.count('R') + 0.5 * seq.count('H')
    neg = seq.count('D') + seq.count('E')
    return pos - neg


def find_surface_charged_residues(sasa_results, scaffold_seq, 
                                   interface_positions=None,
                                   chromophore_proximal=None,
                                   exclude_positions=None):
    """
    Identify surface K/R residues that are candidates for supercharging.
    
    Args:
        sasa_results: list of SASA dicts from calc_sasa_from_pdb
        scaffold_seq: full sequence string
        interface_positions: set of positions to exclude (dimer interface)
        chromophore_proximal: set of positions near chromophore to exclude
        exclude_positions: additional positions to exclude
    
    Returns:
        list of candidate positions for K->E or R->E mutations
    """
    if interface_positions is None:
        interface_positions = set()
    if chromophore_proximal is None:
        chromophore_proximal = set()
    if exclude_positions is None:
        exclude_positions = set()
    
    all_excluded = interface_positions | chromophore_proximal | exclude_positions
    
    candidates = []
    for r in sasa_results:
        if r['is_surface'] and r['aa'] in ('K', 'R'):
            if r['pos'] not in all_excluded:
                candidates.append(r)
    
    return candidates


def find_cysteines(scaffold_seq):
    """Find all cysteine positions in a sequence."""
    return [i + 1 for i, aa in enumerate(scaffold_seq) if aa == 'C']


def design_supercharging_mutations(scaffold_seq, sasa_results=None,
                                    target_charge=-10,
                                    interface_positions=None,
                                    chromophore_proximal=None,
                                    exclude_positions=None,
                                    remove_cysteines=True,
                                    keep_cysteines_at=None):
    """
    Design supercharging mutations for a GFP scaffold.
    
    Strategy:
    1. Remove all surface cysteines (C->S) for CFPS compatibility
    2. Convert surface K/R to E, prioritizing most exposed residues
    3. Stop when target net charge is reached
    
    Args:
        scaffold_seq: full protein sequence
        sasa_results: SASA data (if None, uses sequence-based heuristic)
        target_charge: target net charge (negative)
        interface_positions: positions to exclude (dimer interface)
        chromophore_proximal: positions near chromophore to exclude
        exclude_positions: additional positions to exclude
        remove_cysteines: whether to remove cysteines
        keep_cysteines_at: cysteine positions to keep (e.g. near chromophore)
    
    Returns:
        list of (position, new_aa) mutations, designed sequence
    """
    if keep_cysteines_at is None:
        keep_cysteines_at = set()
    
    mutations = []
    
    # Step 1: Remove cysteines
    if remove_cysteines:
        for pos in find_cysteines(scaffold_seq):
            if pos not in keep_cysteines_at:
                mutations.append((pos, 'S'))
    
    # Step 2: Find surface K/R candidates
    if sasa_results:
        candidates = find_surface_charged_residues(
            sasa_results, scaffold_seq, interface_positions, 
            chromophore_proximal, exclude_positions
        )
        # Sort by RSA (most exposed first)
        candidates.sort(key=lambda x: -x['rsa'])
    else:
        # Heuristic: all K/R positions (less accurate without structure)
        candidates = []
        for i, aa in enumerate(scaffold_seq):
            pos = i + 1
            if aa in ('K', 'R') and pos not in (exclude_positions or set()):
                candidates.append({'pos': pos, 'aa': aa, 'rsa': 50})  # assume surface
    
    # Step 3: Apply K/R -> E mutations until target charge reached
    current_seq = apply_mut(scaffold_seq, mutations)
    current_charge = net_charge(current_seq)
    
    for cand in candidates:
        if current_charge <= target_charge:
            break
        pos = cand['pos']
        old_aa = cand['aa']
        if old_aa == 'K':
            mutations.append((pos, 'E'))
            current_charge -= 2  # K(+1) -> E(-1) = -2 change
        elif old_aa == 'R':
            mutations.append((pos, 'E'))
            current_charge -= 2  # R(+1) -> E(-1) = -2 change
    
    designed_seq = apply_mut(scaffold_seq, mutations)
    final_charge = net_charge(designed_seq)
    
    print(f"\nSupercharging design results:")
    print(f"  Original charge: {net_charge(scaffold_seq):.1f}")
    print(f"  Target charge: {target_charge}")
    print(f"  Achieved charge: {final_charge:.1f}")
    print(f"  Total mutations: {len(mutations)}")
    print(f"  Cysteines removed: {sum(1 for p, a in mutations if a == 'S' and scaffold_seq[p-1] == 'C')}")
    print(f"  K/R->E mutations: {sum(1 for p, a in mutations if a == 'E')}")
    
    return mutations, designed_seq


# ── Proline loop stabilization ───────────────────────────────────────
def find_loop_proline_candidates(thermompnn_csv, ddg_threshold=-0.3):
    """
    Find stabilizing proline substitutions from ThermoMPNN results.
    
    Proline restricts backbone phi angle, reducing unfolded-state entropy
    and improving refolding yield after heat challenge.
    
    Args:
        thermompnn_csv: path to ThermoMPNN inference CSV
        ddg_threshold: ddG threshold for "stabilizing" (default -0.3)
    
    Returns:
        DataFrame of stabilizing proline substitutions
    """
    import pandas as pd
    df = pd.read_csv(thermompnn_csv)
    proline_subs = df[df['mutation'] == 'P'].sort_values('ddG_pred')
    stabilizing = proline_subs[proline_subs['ddG_pred'] < ddg_threshold]
    
    print(f"\nStabilizing proline substitutions (ddG < {ddg_threshold}):")
    for _, row in stabilizing.head(20).iterrows():
        print(f"  {row['wildtype']}{row['position']}P: ddG = {row['ddG_pred']:.3f} kcal/mol")
    
    return stabilizing


# ── Main ─────────────────────────────────────────────────────────────
def main():
    parser = argparse.ArgumentParser(description="Supercharging design pipeline for GFP proteins")
    parser.add_argument("--scaffold", choices=list(SCAFFOLDS.keys()), required=True,
                        help="Scaffold protein")
    parser.add_argument("--pdb", help="PDB structure file for SASA calculation")
    parser.add_argument("--target-charge", type=float, default=-10,
                        help="Target net charge (default -10)")
    parser.add_argument("--dimer", action="store_true",
                        help="Exclude dimer interface positions")
    parser.add_argument("--thermompnn", help="ThermoMPNN CSV for proline analysis")
    parser.add_argument("--output", default="supercharging_design.csv", help="Output CSV path")
    args = parser.parse_args()
    
    scaffold_seq = SCAFFOLDS[args.scaffold]
    print(f"Scaffold: {args.scaffold} ({len(scaffold_seq)} aa)")
    print(f"Original net charge: {net_charge(scaffold_seq):.1f}")
    print(f"Cysteines: {find_cysteines(scaffold_seq)}")
    
    # SASA analysis if PDB provided
    sasa_results = None
    if args.pdb:
        sasa_results = calc_sasa_from_pdb(args.pdb)
        surface_kr = [r for r in sasa_results if r['is_surface'] and r['aa'] in ('K', 'R')]
        print(f"\nSurface K/R residues (RSA > 30%): {len(surface_kr)}")
        for r in surface_kr:
            print(f"  {r['aa']}{r['pos']}: RSA={r['rsa']:.0f}%")
    
    # Dimer interface positions (known from GFP structures)
    interface_positions = set()
    if args.dimer:
        if args.scaffold == 'sfGFP':
            # A206F dimer interface residues
            interface_positions = {148, 149, 202, 203, 204, 205, 206, 207, 208, 209, 210}
        elif args.scaffold in ('StayGold', 'mBaoJin'):
            interface_positions = {91, 135, 136, 137, 138, 140, 142, 144, 149, 151, 152, 153, 155, 165, 167, 187, 189, 191, 216, 217}
        print(f"Excluding {len(interface_positions)} interface positions")
    
    # Design supercharging
    mutations, designed_seq = design_supercharging_mutations(
        scaffold_seq, sasa_results,
        target_charge=args.target_charge,
        interface_positions=interface_positions,
    )
    
    # Proline analysis if ThermoMPNN data provided
    if args.thermompnn:
        find_loop_proline_candidates(args.thermompnn)
    
    # Save results
    import pandas as pd
    mut_records = []
    for pos, new_aa in mutations:
        old_aa = scaffold_seq[pos - 1]
        mut_records.append({
            'position': pos,
            'wildtype': old_aa,
            'mutant': new_aa,
            'mutation': f"{old_aa}{pos}{new_aa}",
            'type': 'cysteine_removal' if old_aa == 'C' else 'supercharging'
        })
    
    result_df = pd.DataFrame(mut_records)
    result_df.to_csv(args.output, index=False)
    print(f"\nMutation list saved to {args.output}")
    print(f"\nDesigned sequence ({len(designed_seq)} aa):")
    print(designed_seq)


if __name__ == "__main__":
    main()
