#!/usr/bin/env python3
"""
ThermoMPNN Thermostability Prediction Pipeline
===============================================
Submits PDB structures to ThermoMPNN via HPC and extracts ddG predictions
for designed GFP mutations.

Usage:
    python run_thermompnn.py --pdb <pdb_file_or_id> --mutations <mutations_csv>

The mutations CSV should have columns: mutation (e.g. "H148S")
If no mutations CSV is provided, all mutations from the ThermoMPNN output are returned.

Requires: biomni.tool (HPC access), pandas, numpy
"""

import argparse
import os
import sys
import json
import pandas as pd
import numpy as np


# ── Reference sequences ──────────────────────────────────────────────
SFGFP = "MSKGEELFTGVVPILVELDGDVNGHKFSVRGEGEGDATNGKLTLKFICTTGKLPVPWPTLVTTLTYGVQCFSRYPDHMKRHDFFKSAMPEGYVQERTISFKDDGTYKTRAEVKFEGDTLVNRIELKGIDFKEDGNILGHKLEYNFNSHNVYITADKQKNGIKANFKIRHNVEDGSVQLADHYQQNTPIGDGPVLLPDNHYLSTQSVLSKDPNEKRDHMVLLEFVTAAGITHGMDELYK"
STAYGOLD = "MASTPFKFQLKGTINGKSFTVEGEGEGNSHEGSHKGKYVCTSGKLPMSWAALGTSFGYGMKYYTKYPSGLKNWFHEVMPEGFTYDRHIQYKGDGSIHAKHQHFMKNGTYHNIVEFTGQDFKENSPVLTGDMNVSLPNEVQHIPRDDGVECPVTLLYPLLSDKSKCVEAHQNTICKPLHNQPAPDVPYHWIRKQYTQSKDDTEERDHICQSETLEAHL"
MBAOJIN = "MVSKGEEENMASTPFKFQLKGTINGKSFTVEGEGEGNSHEGSHKGKYVCTSGKLPMSWAALGTTFGYGMKYYTKYPSGLKNWFREVMPGGFTYDRHIQYKGDGSIHAKHQHFMKNGTYHNIVEFTGQDFKENSPVLTGDMNVSLPNEVPQIPRDDGVECPVTLLYPLLSDKSKYVEAHQYTICKPLHNQPAPDVPYHWIRKQYTQSKDDAEERDHICQSETLEAHLKGMDELYK"


# ── Position mapping ─────────────────────────────────────────────────
def thermompnn_pos_to_seq_pos(tp_pos, pdb_id):
    """
    Map ThermoMPNN output position to 1-based sequence position.
    
    For 2B3P (sfGFP): TP_pos = seq_pos - 2
      (PDB starts at S2, M1 absent; ThermoMPNN 0-indexes from first PDB residue)
    For 8BXT (StayGold): TP_pos = seq_pos
      (PDB starts at M1)
    For 1GFL (avGFP WT): TP_pos = seq_pos - 2
    For 1EMA (eGFP): TP_pos = seq_pos - 2
    """
    mapping = {
        '2B3P': lambda tp: tp + 2,
        '1GFL': lambda tp: tp + 2,
        '1EMA': lambda tp: tp + 2,
        '8BXT': lambda tp: tp,
    }
    pdb_upper = pdb_id.upper()
    if pdb_upper in mapping:
        return mapping[pdb_upper](tp_pos)
    # Default: assume PDB starts at residue 2 (common for GFP structures)
    return tp_pos + 2


def seq_pos_to_thermompnn_pos(seq_pos, pdb_id):
    """Inverse of thermompnn_pos_to_seq_pos."""
    mapping = {
        '2B3P': lambda sp: sp - 2,
        '1GFL': lambda sp: sp - 2,
        '1EMA': lambda sp: sp - 2,
        '8BXT': lambda sp: sp,
    }
    pdb_upper = pdb_id.upper()
    if pdb_upper in mapping:
        return mapping[pdb_upper](seq_pos)
    return seq_pos - 2


# ── PDB download ─────────────────────────────────────────────────────
def download_pdb(pdb_id, output_dir="."):
    """Download PDB structure from RCSB."""
    import urllib.request
    pdb_url = f"https://files.rcsb.org/download/{pdb_id}.pdb"
    pdb_path = os.path.join(output_dir, f"{pdb_id}.pdb")
    if not os.path.exists(pdb_path):
        urllib.request.urlretrieve(pdb_url, pdb_path)
        print(f"Downloaded {pdb_id} to {pdb_path}")
    else:
        print(f"Using cached {pdb_path}")
    return pdb_path


# ── HPC submission ───────────────────────────────────────────────────
def submit_thermompnn(pdb_path):
    """Submit ThermoMPNN job to HPC cluster."""
    from biomni.tool import hpc_run_tool
    
    pdb_filename = os.path.basename(pdb_path)
    job = hpc_run_tool(
        tool_id="thermompnn",
        command=f"thermompnn_predict --pdb_path /input/{pdb_filename} --out_dir /output",
        input_files={pdb_filename: pdb_path}
    )
    print(f"ThermoMPNN job submitted: {job['job_id']}")
    return job['job_id']


def get_thermompnn_results(job_id):
    """Retrieve ThermoMPNN results from completed HPC job."""
    from biomni.tool import hpc_get_job_results
    result = hpc_get_job_results(job_id)
    output_dir = result.get('output_dir', '')
    
    # Find the output CSV
    import glob
    csv_files = glob.glob(os.path.join(output_dir, "ThermoMPNN_inference_*.csv"))
    if not csv_files:
        raise FileNotFoundError(f"No ThermoMPNN output CSV found in {output_dir}")
    
    df = pd.read_csv(csv_files[0])
    print(f"ThermoMPNN results: {len(df)} predictions from {os.path.basename(csv_files[0])}")
    return df, csv_files[0]


# ── ddG extraction ───────────────────────────────────────────────────
def extract_ddg_for_mutations(thermompnn_df, mutations, pdb_id):
    """
    Extract ddG predictions for specific mutations.
    
    Args:
        thermompnn_df: DataFrame from ThermoMPNN output
        mutations: dict of {mutation_name: (seq_position, wt_aa, mutant_aa)}
                   e.g. {"H148S": (148, "H", "S")}
        pdb_id: PDB ID for position mapping
    
    Returns:
        dict of {mutation_name: ddG_value}
    """
    results = {}
    for mut_name, (seq_pos, wt, mutant) in sorted(mutations.items(), key=lambda x: x[1][0]):
        tp_pos = seq_pos_to_thermompnn_pos(seq_pos, pdb_id)
        match = thermompnn_df[
            (thermompnn_df['position'] == tp_pos) & 
            (thermompnn_df['wildtype'] == wt) & 
            (thermompnn_df['mutation'] == mutant)
        ]
        if len(match) > 0:
            ddg = match['ddG_pred'].values[0]
            results[mut_name] = ddg
            stability = classify_ddg(ddg)
            print(f"  {mut_name:10s}  ddG = {ddg:+.3f} kcal/mol  {stability}")
        else:
            print(f"  {mut_name:10s}  NOT FOUND (tp_pos={tp_pos})")
    
    return results


def classify_ddg(ddg):
    """Classify ddG value into stability category."""
    if ddg < -1.0:
        return "STRONGLY STABILIZING"
    elif ddg < -0.3:
        return "STABILIZING"
    elif ddg < 0.3:
        return "NEUTRAL"
    elif ddg < 1.0:
        return "DESTABILIZING"
    else:
        return "STRONGLY DESTABILIZING"


def compute_plan_ddg(mutation_ddgs, plans):
    """
    Compute total ddG for each design plan.
    
    Args:
        mutation_ddgs: dict of {mutation_name: ddG_value}
        plans: dict of {plan_name: [list_of_mutation_names]}
    
    Returns:
        dict of {plan_name: total_ddG}
    """
    plan_totals = {}
    print("\n" + "=" * 65)
    print("Total ddG per Plan (additive model - ignores epistasis)")
    print("=" * 65)
    
    for plan_name, mutations in plans.items():
        total_ddg = sum(mutation_ddgs.get(m, 0) for m in mutations)
        plan_totals[plan_name] = total_ddg
        n_stabilizing = sum(1 for m in mutations if mutation_ddgs.get(m, 0) < 0)
        n_destabilizing = sum(1 for m in mutations if mutation_ddgs.get(m, 0) > 0)
        print(f"\n  {plan_name}")
        print(f"    Total ddG = {total_ddg:+.3f} kcal/mol")
        print(f"    Stabilizing: {n_stabilizing}/{len(mutations)}, "
              f"Destabilizing: {n_destabilizing}/{len(mutations)}")
        for m in mutations:
            ddg = mutation_ddgs.get(m, 0)
            print(f"      {m:10s} {ddg:+.3f}")
    
    return plan_totals


# ── Proline analysis ─────────────────────────────────────────────────
def find_stabilizing_prolines(thermompnn_df, top_n=20):
    """Find proline substitutions that are stabilizing (negative ddG)."""
    proline_subs = thermompnn_df[thermompnn_df['mutation'] == 'P'].sort_values('ddG_pred')
    
    print("=== Stabilizing Proline Substitutions (top {}) ===".format(top_n))
    print("Negative ddG = stabilizing\n")
    for _, row in proline_subs.head(top_n).iterrows():
        print(f"  {row['wildtype']}{row['position']}P: ddG = {row['ddG_pred']:.3f} kcal/mol")
    
    return proline_subs.head(top_n)


# ── Predicted Tm estimation ──────────────────────────────────────────
def estimate_tm(parent_tm, total_ddg, ddg_per_degree=1.5):
    """
    Rough Tm estimation from cumulative ddG.
    
    Very approximate: ~1.5°C decrease per kcal/mol destabilization.
    This is a rough heuristic; actual Tm depends on heat capacity change.
    
    Args:
        parent_tm: Tm of parent protein (°C)
        total_ddg: cumulative ddG (kcal/mol, positive = destabilizing)
        ddg_per_degree: °C change per kcal/mol (default 1.5)
    
    Returns:
        Estimated Tm in °C
    """
    return parent_tm - total_ddg * ddg_per_degree


# ── Main ─────────────────────────────────────────────────────────────
def main():
    parser = argparse.ArgumentParser(description="ThermoMPNN thermostability prediction pipeline")
    parser.add_argument("--pdb", required=True, help="PDB ID (e.g. 2B3P) or path to PDB file")
    parser.add_argument("--mutations", help="CSV file with mutations (column: mutation, e.g. H148S)")
    parser.add_argument("--submit", action="store_true", help="Submit ThermoMPNN job to HPC")
    parser.add_argument("--job-id", help="Retrieve results from existing job ID")
    parser.add_argument("--output", default="thermompnn_results.csv", help="Output CSV path")
    args = parser.parse_args()
    
    pdb_id = args.pdb.upper() if not os.path.exists(args.pdb) else os.path.basename(args.pdb).replace('.pdb', '')
    
    if args.submit:
        # Download PDB if needed
        if not os.path.exists(args.pdb):
            pdb_path = download_pdb(pdb_id)
        else:
            pdb_path = args.pdb
        
        # Submit job
        job_id = submit_thermompnn(pdb_path)
        print(f"\nJob ID: {job_id}")
        print("Use --job-id <id> to retrieve results when complete.")
        return
    
    if args.job_id:
        # Retrieve results
        df, csv_path = get_thermompnn_results(args.job_id)
        
        # If mutations specified, extract ddG for those mutations
        if args.mutations:
            mut_df = pd.read_csv(args.mutations)
            mutations = {}
            for _, row in mut_df.iterrows():
                mut_str = row['mutation']  # e.g. "H148S"
                pos = int(mut_str[1:-1])
                wt = mut_str[0]
                mutant = mut_str[-1]
                mutations[mut_str] = (pos, wt, mutant)
            
            results = extract_ddg_for_mutations(df, mutations, pdb_id)
            
            # Save results
            result_df = pd.DataFrame([
                {'mutation': k, 'seq_pos': v[0], 'ddG': results.get(k, float('nan'))}
                for k, v in mutations.items()
            ])
            result_df.to_csv(args.output, index=False)
            print(f"\nResults saved to {args.output}")
        else:
            # Save full ThermoMPNN output
            print(f"Full results available at: {csv_path}")
            print(f"Columns: {df.columns.tolist()}")
            print(f"Total predictions: {len(df)}")
        
        # Proline analysis
        find_stabilizing_prolines(df)
        
        return
    
    # If neither submit nor job-id, just analyze local CSV
    if os.path.exists(args.output):
        df = pd.read_csv(args.output)
        print(f"Loaded {len(df)} predictions from {args.output}")
    else:
        print("Specify --submit to run ThermoMPNN or --job-id to retrieve results.")
        sys.exit(1)


if __name__ == "__main__":
    main()
