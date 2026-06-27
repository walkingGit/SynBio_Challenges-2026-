#!/usr/bin/env python3
"""
Structure Prediction Pipeline (Boltz-2 / Chai-1)
=================================================
Submits GFP sequences to Boltz-2 or Chai-1 for structure prediction
and extracts confidence metrics (pTM, pLDDT, iPTM).

Usage:
    # Boltz-2 monomer
    python run_structure_prediction.py --tool boltz-2 --sequence "MSKGEEL..." --name Seq1

    # Boltz-2 dimer
    python run_structure_prediction.py --tool boltz-2 --sequence "MSKGEEL..." --dimer --name Seq1_dimer

    # Chai-1 monomer
    python run_structure_prediction.py --tool chai-1 --sequence "MSKGEEL..." --name Seq3

    # Retrieve results
    python run_structure_prediction.py --tool boltz-2 --job-id <job_id> --name Seq1

Requires: biomni.tool (HPC access), pyyaml, numpy
"""

import argparse
import os
import sys
import json
import yaml
import numpy as np


# ── Boltz-2 YAML input preparation ──────────────────────────────────
def create_boltz2_yaml_monomer(sequence, name="protein", output_path="input.yaml"):
    """Create Boltz-2 YAML input for monomer structure prediction."""
    yaml_content = {
        "sequences": [
            {
                "protein": {
                    "id": "chainA",
                    "sequence": sequence
                }
            }
        ]
    }
    with open(output_path, 'w') as f:
        yaml.dump(yaml_content, f, default_flow_style=False)
    print(f"Boltz-2 monomer YAML written to {output_path}")
    return output_path


def create_boltz2_yaml_dimer(sequence, name="protein", output_path="input.yaml"):
    """Create Boltz-2 YAML input for homodimer structure prediction."""
    yaml_content = {
        "sequences": [
            {
                "protein": {
                    "id": "chainA",
                    "sequence": sequence
                }
            },
            {
                "protein": {
                    "id": "chainB",
                    "sequence": sequence
                }
            }
        ]
    }
    with open(output_path, 'w') as f:
        yaml.dump(yaml_content, f, default_flow_style=False)
    print(f"Boltz-2 dimer YAML written to {output_path}")
    return output_path


# ── Chai-1 FASTA input preparation ──────────────────────────────────
def create_chai1_fasta(sequence, name="protein", output_path="input.fasta"):
    """Create Chai-1 FASTA input for structure prediction."""
    fasta_content = f">protein|name={name}\n{sequence}\n"
    with open(output_path, 'w') as f:
        f.write(fasta_content)
    print(f"Chai-1 FASTA written to {output_path}")
    return output_path


# ── HPC submission ───────────────────────────────────────────────────
def submit_boltz2(yaml_path):
    """Submit Boltz-2 structure prediction job."""
    from biomni.tool import hpc_run_tool
    
    yaml_filename = os.path.basename(yaml_path)
    job = hpc_run_tool(
        tool_id="boltz-2",
        command=f"HF_HUB_OFFLINE=1 boltz predict /input/{yaml_filename} "
                f"--out_dir /output --cache /mnt/fsx/dbs/boltz/cache "
                f"--num_workers 0 --use_msa_server",
        input_files={yaml_filename: yaml_path}
    )
    print(f"Boltz-2 job submitted: {job['job_id']}")
    return job['job_id']


def submit_chai1(fasta_path):
    """Submit Chai-1 structure prediction job."""
    from biomni.tool import hpc_run_tool
    
    fasta_filename = os.path.basename(fasta_path)
    job = hpc_run_tool(
        tool_id="chai-1",
        command=f"chai-lab fold --use-msa-server /input/{fasta_filename} /output",
        input_files={fasta_filename: fasta_path}
    )
    print(f"Chai-1 job submitted: {job['job_id']}")
    return job['job_id']


# ── Result parsing ───────────────────────────────────────────────────
def parse_boltz2_results(job_id):
    """
    Parse Boltz-2 results: confidence scores, pLDDT, PAE.
    
    Returns dict with: pTM, pLDDT_mean, pLDDT_per_residue, iptm (for dimers)
    """
    from biomni.tool import hpc_get_job_results
    import glob
    
    result = hpc_get_job_results(job_id)
    output_dir = result.get('output_dir', '')
    
    metrics = {}
    
    # Parse confidence JSON
    conf_files = glob.glob(os.path.join(output_dir, "**/confidence_*.json"), recursive=True)
    if conf_files:
        with open(conf_files[0]) as f:
            conf = json.load(f)
        metrics['pTM'] = conf.get('ptm', None)
        metrics['complex_plddt'] = conf.get('complex_plddt', None)
        if 'iptm' in conf:
            metrics['iPTM'] = conf['iptm']
        print(f"  pTM: {metrics['pTM']:.3f}" if metrics['pTM'] else "  pTM: N/A")
        print(f"  pLDDT: {metrics['complex_plddt']:.3f}" if metrics['complex_plddt'] else "  pLDDT: N/A")
        if 'iPTM' in metrics:
            print(f"  iPTM: {metrics['iPTM']:.3f}")
    
    # Parse pLDDT per-residue (NPZ format)
    plddt_files = glob.glob(os.path.join(output_dir, "**/plddt_*.npz"), recursive=True)
    if plddt_files:
        plddt = np.load(plddt_files[0])
        plddt_arr = plddt['plddt']
        # Boltz-2 may use 0-1 scale; convert to 0-100
        if plddt_arr.max() <= 1.0:
            plddt_100 = plddt_arr * 100
        else:
            plddt_100 = plddt_arr
        
        metrics['pLDDT_mean'] = float(plddt_100.mean())
        metrics['pLDDT_per_residue'] = plddt_100.tolist()
        metrics['pLDDT_below_70'] = int(np.sum(plddt_100 < 70))
        metrics['pLDDT_below_50'] = int(np.sum(plddt_100 < 50))
        
        print(f"  pLDDT mean: {metrics['pLDDT_mean']:.1f}")
        print(f"  Residues pLDDT < 70: {metrics['pLDDT_below_70']}")
        print(f"  Residues pLDDT < 50: {metrics['pLDDT_below_50']}")
    
    # Find structure file (CIF)
    cif_files = glob.glob(os.path.join(output_dir, "**/model_*.cif"), recursive=True)
    if cif_files:
        metrics['structure_file'] = cif_files[0]
        print(f"  Structure: {os.path.basename(cif_files[0])}")
    
    return metrics


def parse_chai1_results(job_id):
    """
    Parse Chai-1 results: pTM, iPTM, aggregate_score from NPZ.
    
    Returns dict with: pTM, iPTM, aggregate_score, per_chain_ptm
    """
    from biomni.tool import hpc_get_job_results
    import glob
    
    result = hpc_get_job_results(job_id)
    output_dir = result.get('output_dir', '')
    
    metrics = {}
    
    # Parse scores NPZ
    scores_files = glob.glob(os.path.join(output_dir, "**/scores.model_idx_*.npz"), recursive=True)
    if scores_files:
        scores = np.load(scores_files[0])
        metrics['pTM'] = float(scores['ptm'].mean())
        metrics['iPTM'] = float(scores['iptm'].mean())
        metrics['aggregate_score'] = float(scores['aggregate_score'].mean())
        
        if 'per_chain_ptm' in scores:
            metrics['per_chain_ptm'] = scores['per_chain_ptm'].tolist()
        if 'has_inter_chain_clashes' in scores:
            metrics['has_clashes'] = bool(scores['has_inter_chain_clashes'])
        
        print(f"  pTM: {metrics['pTM']:.4f}")
        print(f"  iPTM: {metrics['iPTM']:.4f}")
        print(f"  Aggregate score: {metrics['aggregate_score']:.4f}")
    else:
        print("  No scores NPZ found")
    
    # Find structure file
    cif_files = glob.glob(os.path.join(output_dir, "**/*.cif"), recursive=True)
    if cif_files:
        metrics['structure_file'] = cif_files[0]
        print(f"  Structure: {os.path.basename(cif_files[0])}")
    
    return metrics


# ── Batch submission ─────────────────────────────────────────────────
def submit_batch_sequences(sequences, tool="boltz-2", dimer=False):
    """
    Submit multiple sequences for structure prediction.
    
    Args:
        sequences: dict of {name: sequence_string}
        tool: "boltz-2" or "chai-1"
        dimer: if True, submit as homodimer (Boltz-2 only)
    
    Returns:
        dict of {name: job_id}
    """
    jobs = {}
    for name, seq in sequences.items():
        if tool == "boltz-2":
            yaml_path = f"/workspace/{name}_input.yaml"
            if dimer:
                create_boltz2_yaml_dimer(seq, name, yaml_path)
            else:
                create_boltz2_yaml_monomer(seq, name, yaml_path)
            job_id = submit_boltz2(yaml_path)
        elif tool == "chai-1":
            fasta_path = f"/workspace/{name}_input.fasta"
            create_chai1_fasta(seq, name, fasta_path)
            job_id = submit_chai1(fasta_path)
        else:
            raise ValueError(f"Unknown tool: {tool}")
        
        jobs[name] = job_id
    
    print(f"\nSubmitted {len(jobs)} jobs via {tool}")
    for name, jid in jobs.items():
        print(f"  {name}: {jid}")
    
    return jobs


# ── Confidence interpretation ────────────────────────────────────────
def interpret_confidence(metrics):
    """Interpret structure prediction confidence scores."""
    ptm = metrics.get('pTM', 0)
    plddt = metrics.get('pLDDT_mean', metrics.get('complex_plddt', 0) * 100 
                        if metrics.get('complex_plddt', 0) <= 1.0 
                        else metrics.get('complex_plddt', 0))
    
    print("\n=== Confidence Interpretation ===")
    
    # pTM interpretation
    if ptm >= 0.9:
        print(f"  pTM = {ptm:.3f}: HIGH confidence - predicted fold is reliable")
    elif ptm >= 0.7:
        print(f"  pTM = {ptm:.3f}: MODERATE confidence - overall fold likely correct, some uncertainty")
    elif ptm >= 0.5:
        print(f"  pTM = {ptm:.3f}: LOW confidence - fold topology uncertain")
    else:
        print(f"  pTM = {ptm:.3f}: VERY LOW confidence - prediction unreliable")
    
    # pLDDT interpretation
    if plddt >= 90:
        print(f"  pLDDT = {plddt:.1f}: HIGH - atomic-level accuracy expected")
    elif plddt >= 70:
        print(f"  pLDDT = {plddt:.1f}: MODERATE - backbone correct, side chains uncertain")
    elif plddt >= 50:
        print(f"  pLDDT = {plddt:.1f}: LOW - topology likely correct, details uncertain")
    else:
        print(f"  pLDDT = {plddt:.1f}: VERY LOW - prediction unreliable")
    
    # iPTM for dimers
    if 'iPTM' in metrics:
        iptm = metrics['iPTM']
        if iptm >= 0.7:
            print(f"  iPTM = {iptm:.3f}: STRONG dimer interface prediction")
        elif iptm >= 0.5:
            print(f"  iPTM = {iptm:.3f}: MODERATE dimer interface - orientation uncertain")
        else:
            print(f"  iPTM = {iptm:.3f}: WEAK dimer interface - dimer may not form as predicted")


# ── Main ─────────────────────────────────────────────────────────────
def main():
    parser = argparse.ArgumentParser(description="Structure prediction pipeline (Boltz-2 / Chai-1)")
    parser.add_argument("--tool", choices=["boltz-2", "chai-1"], default="boltz-2",
                        help="Structure prediction tool")
    parser.add_argument("--sequence", help="Protein sequence string")
    parser.add_argument("--fasta", help="FASTA file with sequence")
    parser.add_argument("--name", default="protein", help="Sequence name/label")
    parser.add_argument("--dimer", action="store_true", help="Submit as homodimer (Boltz-2 only)")
    parser.add_argument("--job-id", help="Retrieve results from existing job")
    parser.add_argument("--batch-csv", help="CSV with 'name,sequence' columns for batch submission")
    args = parser.parse_args()
    
    if args.job_id:
        # Retrieve results
        print(f"Retrieving {args.tool} results for job {args.job_id}...")
        if args.tool == "boltz-2":
            metrics = parse_boltz2_results(args.job_id)
        else:
            metrics = parse_chai1_results(args.job_id)
        interpret_confidence(metrics)
        return
    
    # Get sequence
    sequence = args.sequence
    if args.fasta and not sequence:
        from Bio import SeqIO
        record = next(SeqIO.parse(args.fasta, "fasta"))
        sequence = str(record.seq)
        args.name = record.id
    
    if not sequence:
        print("Provide --sequence or --fasta")
        sys.exit(1)
    
    # Submit
    if args.tool == "boltz-2":
        yaml_path = f"/workspace/{args.name}_input.yaml"
        if args.dimer:
            create_boltz2_yaml_dimer(sequence, args.name, yaml_path)
        else:
            create_boltz2_yaml_monomer(sequence, args.name, yaml_path)
        job_id = submit_boltz2(yaml_path)
    else:
        fasta_path = f"/workspace/{args.name}_input.fasta"
        create_chai1_fasta(sequence, args.name, fasta_path)
        job_id = submit_chai1(fasta_path)
    
    print(f"\nJob ID: {job_id}")
    print(f"Use --job-id {job_id} --tool {args.tool} to retrieve results.")


if __name__ == "__main__":
    main()
