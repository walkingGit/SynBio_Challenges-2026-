#!/usr/bin/env python3
"""Compute predicted competition scores for the 6 designed GFP sequences.

Competition score = F_final / F_sfGFP = (F_initial / F_sfGFP) * (F_final / F_initial)

Where:
- F_initial / F_sfGFP: brightness relative to sfGFP before heat challenge
- F_final / F_initial: fraction of fluorescence recovered after 72C/10min + 25C/5min renaturation

Score estimates are based on:
- Data-driven brightness mutations (K158G 2.48x, S175T 1.95x, etc.)
- Scaffold intrinsic brightness (EC * QY)
- ThermoMPNN stability predictions
- Heat recovery estimates from scaffold Tm and mutation effects
"""

import pandas as pd

# ── Scaffold properties ────────────────────────────────────────────────
SCAFFOLDS = {
    "sfGFP": {"EC": 83300, "QY": 0.65, "Tm": 87},
    "StayGold": {"EC": 159000, "QY": 0.93, "Tm": 95},
    "mBaoJin": {"EC": 128000, "QY": 0.93, "Tm": 95},
}

# sfGFP reference: EC=83,300, QY=0.65 -> brightness = 54,145
SF_BRIGHTNESS = 83300 * 0.65

# ── Design parameters ──────────────────────────────────────────────────
DESIGNS = {
    1: {
        "name": "StayGold Dimer",
        "scaffold": "StayGold",
        "brightness_mult": 3.5,   # EC enhancement from dimer + scaffold brightness
        "heat_recovery": 0.50,    # StayGold at 72C is 23C below Tm, but dimer may partially dissociate
        "confidence": "Medium",
    },
    2: {
        "name": "mBaoJin Monomer",
        "scaffold": "mBaoJin",
        "brightness_mult": 2.2,   # mBaoJin is 2.2x brighter than sfGFP
        "heat_recovery": 0.75,    # High Tm (95C), monomer, supercharged
        "confidence": "High",
    },
    3: {
        "name": "sfGFP H148S Dimer v2",
        "scaffold": "sfGFP",
        "brightness_mult": 5.2,   # H148S + K158G + S175T + dimer EC + exciton coupling
        "heat_recovery": 0.55,    # Improved by N198P + H231P + V193L, but dimer risk
        "confidence": "Medium",
    },
    4: {
        "name": "sfGFP H148S Monomer v2",
        "scaffold": "sfGFP",
        "brightness_mult": 4.5,   # H148S + K158G + S175T (no dimer EC bonus)
        "heat_recovery": 0.60,    # Improved by prolines, no dimer risk
        "confidence": "Medium-High",
    },
    5: {
        "name": "mBaoJin Supercharged",
        "scaffold": "mBaoJin",
        "brightness_mult": 2.2,   # Same scaffold brightness as Seq 2
        "heat_recovery": 0.85,    # Aggressive supercharging + prolines + high Tm
        "confidence": "High",
    },
    6: {
        "name": "sfGFP V206F Dimer v2",
        "scaffold": "sfGFP",
        "brightness_mult": 3.2,   # K158G + S175T + V206F dimer EC
        "heat_recovery": 0.55,    # Improved by prolines, dimer risk
        "confidence": "Medium",
    },
}


def main():
    print("Predicted Competition Scores")
    print("=" * 80)
    print(f"{'Seq':>4} {'Design':<28} {'F_init/sfGFP':>12} {'F_fin/init':>10} {'F_fin/sfGFP':>12} {'Confidence':>12}")
    print("-" * 80)

    results = []
    for seq_id in sorted(DESIGNS.keys()):
        d = DESIGNS[seq_id]
        f_init = d["brightness_mult"]
        f_fin_init = d["heat_recovery"]
        f_fin = f_init * f_fin_init

        print(f"{seq_id:>4} {d['name']:<28} {f_init:>10.1f}x  {f_fin_init:>10.2f}  {f_fin:>10.1f}x  {d['confidence']:>12}")

        results.append({
            "Seq_ID": seq_id,
            "Design": d["name"],
            "Scaffold": d["scaffold"],
            "F_initial_over_sfGFP": f_init,
            "F_final_over_initial": f_fin_init,
            "F_final_over_sfGFP": round(f_fin, 2),
            "Confidence": d["confidence"],
        })

    df = pd.DataFrame(results)
    output_path = "predicted_scores.csv"
    df.to_csv(output_path, index=False)
    print(f"\nSaved to {output_path}")

    # Gain/Risk analysis
    print("\n\nGain/Risk Analysis")
    print("=" * 80)

    gain_risk = {
        1: (9, 8),   # High gain (dimer EC), high risk (dimer may not form in CFPS)
        2: (7, 4),   # Moderate gain, low risk (monomer, high Tm)
        3: (9, 6),   # High gain, medium risk (dimer + prolines improved)
        4: (8, 3),   # High gain, low risk (monomer, prolines, no dimer dependency)
        5: (7, 6),   # Moderate gain, medium risk (aggressive mutations)
        6: (9, 6),   # High gain, medium risk (dimer + prolines improved)
    }

    print(f"{'Seq':>4} {'Design':<28} {'Gain':>5} {'Risk':>5} {'Net':>5} {'Verdict':>10}")
    print("-" * 80)
    for seq_id in sorted(gain_risk.keys()):
        g, r = gain_risk[seq_id]
        net = g - r
        if net >= 5:
            verdict = "BEST"
        elif net >= 3:
            verdict = "KEEP"
        else:
            verdict = "MARGINAL"
        print(f"{seq_id:>4} {DESIGNS[seq_id]['name']:<28} {g:>5} {r:>5} {net:>+5} {verdict:>10}")


if __name__ == "__main__":
    main()
