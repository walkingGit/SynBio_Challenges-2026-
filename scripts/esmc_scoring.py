#!/usr/bin/env python3
"""
ESMC Mutation Scoring Pipeline
===============================
Uses ESMC (6B parameter protein language model) to score mutations
and combines with ThermoMPNN ddG for composite stability assessment.

ESMC score interpretation:
  score = mutant_log_prob - wildtype_log_prob
  Negative = FAVORABLE (mutation more probable than WT at that position)
  Positive = UNFAVORABLE (mutation less probable than WT)
  NOTE: This is the OPPOSITE convention from ThermoMPNN ddG!

Usage:
    # Submit ESMC scoring job
    python esmc_scoring.py --fasta protein.fasta --mutations mutations.csv --submit

    # Analyze results
    python esmc_scoring.py --esmc-results mutation_scores.csv --thermompnn-results ThermoMPNN_inference_2B3P.csv

    # Composite scoring
    python esmc_scoring.py --esmc-results mutation_scores.csv --thermompnn-results ThermoMPNN_inference_2B3P.csv --plans plans.json

Requires: biomni.tool (HPC access), pandas, numpy
"""

import argparse
import os
import sys
import json
import pandas as pd
import numpy as np


# ── ESMC HPC submission ──────────────────────────────────────────────
def create_esmc_fasta(sequence, name="protein", output_path="protein.fasta"):
    """Create FASTA input for ESMC."""
    with open(output_path, 'w') as f:
        f.write(f">{name}\n{sequence}\n")
    print(f"ESMC FASTA written to {output_path}")
    return output_path


def create_esmc_mutations_csv(mutations, output_path="mutations.csv"):
    """
    Create mutations CSV for ESMC score-mutations.
    
    Args:
        mutations: list of mutation strings (e.g. ["H148S", "K158G"])
    """
    with open(output_path, 'w') as f:
        f.write("mutation\n")
        for mut in mutations:
            f.write(f"{mut}\n")
    print(f"ESMC mutations CSV written to {output_path} ({len(mutations)} mutations)")
    return output_path


def submit_esmc_scoring(fasta_path, mutations_path):
    """Submit ESMC mutation scoring job to HPC."""
    from biomni.tool import hpc_run_tool
    
    fasta_filename = os.path.basename(fasta_path)
    mutations_filename = os.path.basename(mutations_path)
    
    job = hpc_run_tool(
        tool_id="esmc",
        command=f"esmc score-mutations --input /input/{fasta_filename} "
                f"--mutations /input/{mutations_filename} --output /output",
        input_files={
            fasta_filename: fasta_path,
            mutations_filename: mutations_path,
        }
    )
    print(f"ESMC scoring job submitted: {job['job_id']}")
    return job['job_id']


def get_esmc_results(job_id):
    """Retrieve ESMC results from completed HPC job."""
    from biomni.tool import hpc_get_job_results
    import glob
    
    result = hpc_get_job_results(job_id)
    output_dir = result.get('output_dir', '')
    
    csv_files = glob.glob(os.path.join(output_dir, "mutation_scores*.csv"))
    if not csv_files:
        csv_files = glob.glob(os.path.join(output_dir, "*.csv"))
    
    if not csv_files:
        raise FileNotFoundError(f"No ESMC output CSV found in {output_dir}")
    
    df = pd.read_csv(csv_files[0])
    print(f"ESMC results: {len(df)} mutation scores")
    return df


# ── Score interpretation ─────────────────────────────────────────────
def classify_esmc_score(score):
    """Classify ESMC score into favorability category."""
    if score < -0.5:
        return "FAVORABLE"
    elif score < 0:
        return "MILDLY FAVORABLE"
    elif score < 0.5:
        return "MILDLY UNFAVORABLE"
    else:
        return "UNFAVORABLE"


def classify_ddg(ddg):
    """Classify ThermoMPNN ddG into stability category."""
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


def compare_esmc_thermompnn(esmc_df, thermompnn_df, pdb_id='2B3P'):
    """
    Compare ESMC and ThermoMPNN scores for overlapping mutations.
    
    Both negative = good, but different scales and methods.
    ESMC captures evolutionary likelihood; ThermoMPNN captures structural stability.
    """
    # Position mapping for ThermoMPNN
    def tp_pos_to_seq_pos(tp_pos):
        if pdb_id.upper() in ('2B3P', '1GFL', '1EMA'):
            return tp_pos + 2
        return tp_pos  # 8BXT
    
    # Build ThermoMPNN lookup: mutation_string -> ddG
    tm_lookup = {}
    for _, row in thermompnn_df.iterrows():
        seq_pos = tp_pos_to_seq_pos(row['position'])
        mut_str = f"{row['wildtype']}{seq_pos}{row['mutation']}"
        tm_lookup[mut_str] = row['ddG_pred']
    
    print("\n" + "=" * 75)
    print("ESMC vs ThermoMPNN Comparison")
    print("(Both negative = good, but different scales)")
    print("=" * 75)
    print(f"\n{'Mutation':>10s} {'ESMC':>10s} {'ddG':>10s} {'Agreement':>15s}")
    print("-" * 50)
    
    agreements = {'agree_good': 0, 'agree_bad': 0, 'disagree': 0, 'no_tm': 0}
    
    for _, row in esmc_df.sort_values('score').iterrows():
        mut = row['mutation']
        esmc_score = row['score']
        ddg = tm_lookup.get(mut, None)
        
        if ddg is not None:
            if esmc_score < 0 and ddg < 0:
                agreement = "AGREE (good)"
                agreements['agree_good'] += 1
            elif esmc_score > 0 and ddg > 0:
                agreement = "AGREE (bad)"
                agreements['agree_bad'] += 1
            else:
                agreement = "DISAGREE"
                agreements['disagree'] += 1
            print(f"  {mut:>8s} {esmc_score:>+10.4f} {ddg:>+10.3f} {agreement:>15s}")
        else:
            agreements['no_tm'] += 1
    
    print(f"\nSummary: {agreements}")
    return agreements


# ── Composite plan scoring ───────────────────────────────────────────
def compute_composite_scores(esmc_df, thermompnn_df, plans, pdb_id='2B3P'):
    """
    Compute composite ESMC + ThermoMPNN scores for design plans.
    
    Args:
        esmc_df: ESMC mutation scores DataFrame
        thermompnn_df: ThermoMPNN predictions DataFrame
        plans: dict of {plan_name: [list_of_mutation_strings]}
        pdb_id: PDB ID for position mapping
    
    Returns:
        DataFrame with plan rankings
    """
    def tp_pos_to_seq_pos(tp_pos):
        if pdb_id.upper() in ('2B3P', '1GFL', '1EMA'):
            return tp_pos + 2
        return tp_pos
    
    # Build lookups
    esmc_lookup = dict(zip(esmc_df['mutation'], esmc_df['score']))
    tm_lookup = {}
    for _, row in thermompnn_df.iterrows():
        seq_pos = tp_pos_to_seq_pos(row['position'])
        mut_str = f"{row['wildtype']}{seq_pos}{row['mutation']}"
        tm_lookup[mut_str] = row['ddG_pred']
    
    print("\n" + "=" * 75)
    print("COMPOSITE PLAN SCORING: ESMC (sequence) + ThermoMPNN (structure)")
    print("=" * 75)
    
    rank_data = []
    for plan_name, mutations in plans.items():
        esmc_total = sum(esmc_lookup.get(m, 0) for m in mutations)
        esmc_available = sum(1 for m in mutations if m in esmc_lookup)
        tm_total = sum(tm_lookup.get(m, 0) for m in mutations)
        tm_available = sum(1 for m in mutations if m in tm_lookup)
        missing = [m for m in mutations if m not in esmc_lookup and m not in tm_lookup]
        
        print(f"\n  {plan_name}")
        print(f"    ESMC total: {esmc_total:+.3f} ({esmc_available}/{len(mutations)} scored)")
        print(f"    ThermoMPNN ddG: {tm_total:+.3f} kcal/mol ({tm_available}/{len(mutations)} scored)")
        if missing:
            print(f"    Missing from both: {missing}")
        
        rank_data.append({
            'Plan': plan_name,
            'ESMC_total': esmc_total,
            'TM_ddG_total': tm_total,
            'ESMC_scored': esmc_available,
            'TM_scored': tm_available,
        })
    
    # Rank plans
    df_rank = pd.DataFrame(rank_data)
    df_rank['ESMC_rank'] = df_rank['ESMC_total'].rank(ascending=True)  # more negative = better
    df_rank['TM_rank'] = df_rank['TM_ddG_total'].rank(ascending=True)
    df_rank['Combined_rank'] = df_rank['ESMC_rank'] + df_rank['TM_rank']
    df_rank = df_rank.sort_values('Combined_rank')
    
    print("\n" + "=" * 75)
    print("FINAL COMPOSITE RANKING")
    print("=" * 75)
    for i, (_, row) in enumerate(df_rank.iterrows(), 1):
        print(f"  #{i} {row['Plan']:30s}  "
              f"ESMC={row['ESMC_total']:+.3f} (rank {row['ESMC_rank']:.0f})  "
              f"ddG={row['TM_ddG_total']:+.3f} (rank {row['TM_rank']:.0f})  "
              f"Combined={row['Combined_rank']:.0f}")
    
    return df_rank


# ── Main ─────────────────────────────────────────────────────────────
def main():
    parser = argparse.ArgumentParser(description="ESMC mutation scoring pipeline")
    parser.add_argument("--fasta", help="Protein FASTA file for ESMC input")
    parser.add_argument("--mutations", help="Mutations CSV (column: mutation)")
    parser.add_argument("--submit", action="store_true", help="Submit ESMC scoring job")
    parser.add_argument("--job-id", help="Retrieve results from existing job")
    parser.add_argument("--esmc-results", help="ESMC results CSV (local)")
    parser.add_argument("--thermompnn-results", help="ThermoMPNN results CSV for comparison")
    parser.add_argument("--pdb-id", default="2B3P", help="PDB ID for position mapping")
    parser.add_argument("--plans", help="JSON file with plan definitions")
    parser.add_argument("--output", default="esmc_composite_scores.csv", help="Output CSV path")
    args = parser.parse_args()
    
    if args.submit:
        if not args.fasta or not args.mutations:
            print("Need --fasta and --mutations for submission")
            sys.exit(1)
        job_id = submit_esmc_scoring(args.fasta, args.mutations)
        print(f"\nJob ID: {job_id}")
        print(f"Use --job-id {job_id} to retrieve results.")
        return
    
    # Load ESMC results
    esmc_df = None
    if args.job_id:
        esmc_df = get_esmc_results(args.job_id)
    elif args.esmc_results:
        esmc_df = pd.read_csv(args.esmc_results)
        print(f"Loaded ESMC results: {len(esmc_df)} scores")
    
    if esmc_df is None:
        print("Provide --submit, --job-id, or --esmc-results")
        sys.exit(1)
    
    # Print ESMC scores
    print("\nESMC Mutation Scores:")
    for _, row in esmc_df.sort_values('score').iterrows():
        verdict = classify_esmc_score(row['score'])
        print(f"  {row['mutation']:>10s}: {row['score']:+.4f}  ({verdict})")
    
    # Compare with ThermoMPNN if available
    if args.thermompnn_results:
        tm_df = pd.read_csv(args.thermompnn_results)
        compare_esmc_thermompnn(esmc_df, tm_df, args.pdb_id)
        
        # Composite scoring if plans provided
        if args.plans:
            with open(args.plans) as f:
                plans = json.load(f)
            rank_df = compute_composite_scores(esmc_df, tm_df, plans, args.pdb_id)
            rank_df.to_csv(args.output, index=False)
            print(f"\nComposite scores saved to {args.output}")


if __name__ == "__main__":
    main()
