#!/usr/bin/env python3
"""
T5 — Cross-trajectory geometry comparison.

Merges torsion classification (T4), energies (T2), and RMSDs (T3) on
(system, frame) and computes group statistics + significance tests.

Usage:
    python scripts/T5_cross_compare.py
    python scripts/T5_cross_compare.py \
        --torsion analysis/torsion_classification.csv \
        --energy analysis/energies_all.csv \
        --rmsd analysis/rmsd_all.csv \
        --output-prefix analysis/
"""

import argparse
import os
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy import stats


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def cohens_d(x, y):
    """Cohen's d for two independent samples."""
    nx, ny = len(x), len(y)
    dof = nx + ny - 2
    pooled_std = np.sqrt(((nx - 1) * x.var(ddof=1) + (ny - 1) * y.var(ddof=1)) / dof)
    if pooled_std == 0:
        return 0.0
    return (x.mean() - y.mean()) / pooled_std


def rank_biserial(x, y):
    """Rank-biserial correlation effect size for Mann-Whitney U."""
    nx, ny = len(x), len(y)
    u, _ = stats.mannwhitneyu(x, y, alternative="two-sided")
    return 1 - (2 * u) / (nx * ny)


def choose_test(a, b, subsample_max=5000):
    """
    Choose between t-test and Mann-Whitney U.
    Shapiro-Wilk on a subsample if > subsample_max.
    Returns: test_name, statistic, pvalue, effect_size, mean_a, mean_b
    """
    a = np.asarray(a, dtype=float)
    b = np.asarray(b, dtype=float)
    a = a[np.isfinite(a)]
    b = b[np.isfinite(b)]

    if len(a) < 2 or len(b) < 2:
        return None, None, None, None, None, None

    # subsample for Shapiro
    def _shapiro(arr):
        if len(arr) > subsample_max:
            arr = np.random.choice(arr, size=subsample_max, replace=False)
        if len(arr) < 3:
            return None
        stat, p = stats.shapiro(arr)
        return p

    pa = _shapiro(a)
    pb = _shapiro(b)

    normal = (pa is None or pa > 0.05) and (pb is None or pb > 0.05)

    if normal:
        t_stat, p_val = stats.ttest_ind(a, b, equal_var=False)
        d = cohens_d(a, b)
        return "welch_t", t_stat, p_val, d, float(a.mean()), float(b.mean())
    else:
        u_stat, p_val = stats.mannwhitneyu(a, b, alternative="two-sided")
        rbc = rank_biserial(a, b)
        return "mann_whitney_u", u_stat, p_val, rbc, float(np.median(a)), float(np.median(b))


def compute_fractions(df, group_cols):
    """Compute SAP/TSAP/UNK fractions for group combinations."""
    rows = []
    for keys, sub in df.groupby(group_cols):
        total = len(sub)
        if total == 0:
            continue
        row = dict(zip(group_cols, keys if isinstance(keys, tuple) else [keys]))
        row["n_frames"] = total
        row["frac_sap"] = (sub["actual_class"] == "SAP").sum() / total
        row["frac_tsap"] = (sub["actual_class"] == "TSAP").sum() / total
        row["frac_unk"] = (sub["actual_class"] == "UNK").sum() / total
        rows.append(row)
    return pd.DataFrame(rows)


def compute_means(df, group_cols):
    """Compute mean + sem for RMSD and energy by group."""
    rows = []
    for keys, sub in df.groupby(group_cols):
        total = len(sub)
        if total == 0:
            continue
        row = dict(zip(group_cols, keys if isinstance(keys, tuple) else [keys]))
        row["n_frames"] = total
        for col in ["rmsd_heavy_A", "rmsd_core_A", "relative_kjmol"]:
            if col in sub.columns:
                vals = sub[col].dropna()
                row[f"mean_{col}"] = vals.mean()
                row[f"sem_{col}"] = stats.sem(vals) if len(vals) > 1 else 0.0
        rows.append(row)
    return pd.DataFrame(rows)


# ---------------------------------------------------------------------------
# Main analysis
# ---------------------------------------------------------------------------

def main(args):
    # Ensure output dir exists
    output_dir = Path(args.output_prefix)
    output_dir.mkdir(parents=True, exist_ok=True)

    # Load data
    print(f"[T5] Loading {args.torsion} ...")
    tc = pd.read_csv(args.torsion)
    print(f"[T5] Loading {args.energy} ...")
    energy = pd.read_csv(args.energy)
    print(f"[T5] Loading {args.rmsd} ...")
    rmsd = pd.read_csv(args.rmsd)

    # --- merging ---
    # Inner join on (system, frame) to retain shared frames only
    merged = energy.merge(tc, on=["system", "frame"], how="inner", suffixes=("", "_tc"))
    merged = merged.merge(rmsd, on=["system", "frame"], how="inner", suffixes=("", "_rmsd"))

    print(f"[T5] Merged shape: {merged.shape}")
    if merged.empty:
        raise ValueError("Merged DataFrame is empty — check shared (system, frame) keys.")

    # --- derive columns ---
    merged["actual_class"] = merged["classification"]
    merged["trajectory_type"] = merged["topology"].str.upper()  # sap -> SAP, tsap -> TSAP

    # --- group statistics ---
    print("[T5] Computing group fractions ...")
    frac_species_traj = compute_fractions(merged, ["species", "trajectory_type"])
    frac_isomer_traj = compute_fractions(merged, ["isomer", "trajectory_type"])
    frac_hand_traj = compute_fractions(merged, ["handedness", "trajectory_type"])

    print("[T5] Computing group means (RMSD + energy) ...")
    means_species_class = compute_means(merged, ["species", "actual_class"])
    means_species_traj = compute_means(merged, ["species", "trajectory_type"])

    # --- statistical tests ---
    print("[T5] Running significance tests ...")
    test_results = []

    # Helper to run pairwise tests
    def run_pairwise(label, col_name, cat_col, val_a, val_b):
        sub_a = merged[merged[cat_col] == val_a][col_name].dropna()
        sub_b = merged[merged[cat_col] == val_b][col_name].dropna()
        test_name, stat, pval, eff, ma, mb = choose_test(sub_a, sub_b)
        test_results.append({
            "comparison": label,
            "metric": col_name,
            "group_a": val_a,
            "group_b": val_b,
            "n_a": len(sub_a),
            "n_b": len(sub_b),
            "test": test_name,
            "statistic": stat,
            "p_value": pval,
            "effect_size": eff,
            "location_a": ma,
            "location_b": mb,
        })

    # 1. me vs phe
    run_pairwise("me_vs_phe", "relative_kjmol", "species", "me", "phe")
    run_pairwise("me_vs_phe", "rmsd_heavy_A", "species", "me", "phe")
    run_pairwise("me_vs_phe", "rmsd_core_A", "species", "me", "phe")

    # 2. rrr vs sss
    run_pairwise("rrr_vs_sss", "relative_kjmol", "isomer", "rrr", "sss")
    run_pairwise("rrr_vs_sss", "rmsd_heavy_A", "isomer", "rrr", "sss")
    run_pairwise("rrr_vs_sss", "rmsd_core_A", "isomer", "rrr", "sss")

    # 3. D vs L
    run_pairwise("D_vs_L", "relative_kjmol", "handedness", "D", "L")
    run_pairwise("D_vs_L", "rmsd_heavy_A", "handedness", "D", "L")
    run_pairwise("D_vs_L", "rmsd_core_A", "handedness", "D", "L")

    # 4. SAP vs TSAP (actual_class)
    run_pairwise("SAP_vs_TSAP", "relative_kjmol", "actual_class", "SAP", "TSAP")
    run_pairwise("SAP_vs_TSAP", "rmsd_heavy_A", "actual_class", "SAP", "TSAP")
    run_pairwise("SAP_vs_TSAP", "rmsd_core_A", "actual_class", "SAP", "TSAP")

    test_df = pd.DataFrame(test_results)

    # Print test summary
    print("\n" + "="*60)
    print("Statistical Test Results")
    print("="*60)
    for _, row in test_df.iterrows():
        sig = "***" if row["p_value"] < 0.001 else "**" if row["p_value"] < 0.01 else "*" if row["p_value"] < 0.05 else "ns"
        print(f"{row['comparison']:14s} | {row['metric']:18s} | {row['test']:12s} | "
              f"p={row['p_value']:.2e} {sig:3s} | eff={row['effect_size']:.3f}")
    print("="*60 + "\n")

    # --- Plots ---
    print("[T5] Generating plots ...")

    # 1. Grouped bar: SAP/TSAP/UNK fractions by species, separated by trajectory_type
    fig, axes = plt.subplots(1, 2, figsize=(10, 5), sharey=True)
    plot_data = frac_species_traj.sort_values(["trajectory_type", "species"])
    for idx, (traj_type, ax) in enumerate(zip(["SAP", "TSAP"], axes)):
        sub = plot_data[plot_data["trajectory_type"] == traj_type]
        x = np.arange(len(sub))
        width = 0.25
        ax.bar(x - width, sub["frac_sap"], width, label="SAP", color="#1f77b4")
        ax.bar(x, sub["frac_tsap"], width, label="TSAP", color="#ff7f0e")
        ax.bar(x + width, sub["frac_unk"], width, label="UNK", color="#7f7f7f")
        ax.set_xticks(x)
        ax.set_xticklabels(sub["species"].tolist())
        ax.set_ylabel("Fraction")
        ax.set_title(f"{traj_type}-started trajectories")
        ax.set_ylim(0, 1.05)
        if idx == 1:
            ax.legend(loc="upper right")
    fig.suptitle("SAP / TSAP / UNK fractions by species")
    plt.tight_layout()
    fig.savefig(output_dir / "frac_by_species.png", dpi=300)
    plt.close(fig)
    print("  -> frac_by_species.png")

    # 2. Violin/Box plots: RMSD distributions by group
    fig, axes = plt.subplots(2, 3, figsize=(12, 8), sharex="col")
    groups = [
        ("species", ["me", "phe"]),
        ("isomer", ["rrr", "sss"]),
        ("handedness", ["D", "L"]),
    ]
    for col_idx, (group_col, group_vals) in enumerate(groups):
        for row_idx, (rmsd_col, ylabel) in enumerate(
            [("rmsd_heavy_A", "RMSD heavy (Å)"), ("rmsd_core_A", "RMSD core (Å)")]
        ):
            ax = axes[row_idx, col_idx]
            parts = []
            for val in group_vals:
                data = merged[merged[group_col] == val][rmsd_col].dropna()
                parts.append(data.values)
            vp = ax.violinplot(parts, positions=range(len(group_vals)), showmeans=False, showmedians=True, widths=0.7)
            for pc, color in zip(vp["bodies"], ["#1f77b4", "#ff7f0e"]):
                pc.set_facecolor(color)
                pc.set_alpha(0.6)
            # overlay boxes
            bp = ax.boxplot(parts, positions=range(len(group_vals)), widths=0.15, showfliers=False, patch_artist=True)
            for patch, color in zip(bp["boxes"], ["#1f77b4", "#ff7f0e"]):
                patch.set_facecolor(color)
                patch.set_alpha(0.3)
            ax.set_xticks(range(len(group_vals)))
            ax.set_xticklabels(group_vals)
            ax.set_ylabel(ylabel)
            ax.set_title(group_col.capitalize())
    fig.suptitle("RMSD distributions by group")
    plt.tight_layout()
    fig.savefig(output_dir / "rmsd_by_group.png", dpi=300)
    plt.close(fig)
    print("  -> rmsd_by_group.png")

    # 3. Energy distributions by group
    fig, axes = plt.subplots(2, 2, figsize=(12, 10))
    plot_groups = [
        ("species", "me", "phe", "Species"),
        ("isomer", "rrr", "sss", "Isomer"),
        ("handedness", "D", "L", "Handedness"),
        ("actual_class", "SAP", "TSAP", "Actual class"),
    ]
    for ax, (col, val_a, val_b, title) in zip(axes.flat, plot_groups):
        data_a = merged[merged[col] == val_a]["relative_kjmol"].dropna()
        data_b = merged[merged[col] == val_b]["relative_kjmol"].dropna()
        bins = np.linspace(
            min(data_a.quantile(0.01), data_b.quantile(0.01)),
            max(data_a.quantile(0.99), data_b.quantile(0.99)),
            60,
        )
        ax.hist(data_a, bins=bins, alpha=0.6, label=f"{val_a} (n={len(data_a)})", density=True, color="#1f77b4")
        ax.hist(data_b, bins=bins, alpha=0.6, label=f"{val_b} (n={len(data_b)})", density=True, color="#ff7f0e")
        ax.set_xlabel("Relative energy (kJ/mol)")
        ax.set_ylabel("Density")
        ax.set_title(title)
        ax.legend(loc="upper right")
    fig.suptitle("Energy distributions by group")
    plt.tight_layout()
    fig.savefig(output_dir / "energy_by_group.png", dpi=300)
    plt.close(fig)
    print("  -> energy_by_group.png")

    # --- build summary CSV ---
    summary_rows = []

    # fraction summaries (long format per plan)
    for df_frac, group_col in [
        (frac_species_traj, "species"),
        (frac_isomer_traj, "isomer"),
        (frac_hand_traj, "handedness"),
    ]:
        for _, row in df_frac.iterrows():
            summary_rows.append({
                "group": group_col,
                "subgroup_value": row[group_col],
                "trajectory_type": row["trajectory_type"],
                "actual_class": "all",
                "n_frames": row["n_frames"],
                "frac_sap": row["frac_sap"],
                "frac_tsap": row["frac_tsap"],
                "frac_unk": row["frac_unk"],
            })

    # mean summaries
    for df_mean, group_cols in [
        (means_species_class, ["species", "actual_class"]),
        (means_species_traj, ["species", "trajectory_type"]),
    ]:
        for _, row in df_mean.iterrows():
            d = {
                "group": ",".join(group_cols),
                "subgroup_value": ",".join(str(row[c]) for c in group_cols),
                "trajectory_type": row.get("trajectory_type", ""),
                "actual_class": row.get("actual_class", ""),
                "n_frames": row["n_frames"],
                "mean_rmsd_heavy_A": row.get("mean_rmsd_heavy_A", np.nan),
                "sem_rmsd_heavy_A": row.get("sem_rmsd_heavy_A", np.nan),
                "mean_rmsd_core_A": row.get("mean_rmsd_core_A", np.nan),
                "sem_rmsd_core_A": row.get("sem_rmsd_core_A", np.nan),
                "mean_relative_kjmol": row.get("mean_relative_kjmol", np.nan),
                "sem_relative_kjmol": row.get("sem_relative_kjmol", np.nan),
                "frac_sap": np.nan,
                "frac_tsap": np.nan,
                "frac_unk": np.nan,
            }
            summary_rows.append(d)

    summary_df = pd.DataFrame(summary_rows)
    # Ensure canonical column order
    col_order = [
        "group", "subgroup_value", "trajectory_type", "actual_class",
        "n_frames", "frac_sap", "frac_tsap", "frac_unk",
        "mean_rmsd_heavy_A", "sem_rmsd_heavy_A",
        "mean_rmsd_core_A", "sem_rmsd_core_A",
        "mean_relative_kjmol", "sem_relative_kjmol",
    ]
    summary_df = summary_df[[c for c in col_order if c in summary_df.columns]]

    out_csv = output_dir / "cross_trajectory_summary.csv"
    summary_df.to_csv(out_csv, index=False)
    print(f"  -> {out_csv}")

    # Save test results too
    out_tests = output_dir / "cross_trajectory_tests.csv"
    test_df.to_csv(out_tests, index=False)
    print(f"  -> {out_tests}")

    print("\n[T5] Done.")

    # Return data for programmatic use
    return {
        "merged": merged,
        "fractions": {
            "species": frac_species_traj,
            "isomer": frac_isomer_traj,
            "handedness": frac_hand_traj,
        },
        "means": {
            "species_class": means_species_class,
            "species_traj": means_species_traj,
        },
        "tests": test_df,
        "summary": summary_df,
    }


# ---------------------------------------------------------------------------
# CLI entry point
# ---------------------------------------------------------------------------

def parse_args(argv=None):
    parser = argparse.ArgumentParser(description="Cross-trajectory geometry comparison (T5)")
    parser.add_argument("--torsion", default="analysis/torsion_classification.csv",
                        help="Path to torsion classification CSV")
    parser.add_argument("--energy", default="analysis/energies_all.csv",
                        help="Path to energies CSV")
    parser.add_argument("--rmsd", default="analysis/rmsd_all.csv",
                        help="Path to RMSD CSV")
    parser.add_argument("--output-prefix", default="analysis/",
                        help="Directory prefix for output files")
    parser.add_argument("--random-seed", type=int, default=42,
                        help="Random seed for Shapiro subsampling")
    return parser.parse_args(argv)


if __name__ == "__main__":
    args = parse_args()
    np.random.seed(args.random_seed)
    main(args)
