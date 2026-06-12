#!/usr/bin/env python3
"""
check_solv_replicate_consistency.py — Assess solv_md replicate convergence.

The solv_md data consists of 4 × 100 ns concatenated replicates. This script
assigns each frame to a replicate ID and computes per-replicate centroids and
standard deviations in PC1/PC2 space. The consistency ratio (CR) is defined as:

    CR = max pairwise centroid distance among 4 replicates / overall PC std

A low CR (< 0.5) indicates that replicates cluster tightly relative to the
overall spread, i.e. good convergence. A high CR (> 1) suggests that
replicates occupy distinct sub-regions of PC space, i.e. poor convergence.

Produces:
  - solv_replicate_consistency.csv  (per-replicate stats + CR)
  - plot_replicate_consistency.png  (4x4 grid of replicate centroids)

Usage:
  python check_solv_replicate_consistency.py \
    --input cross_dataset/analysis/joint_projection_eu8_nochrom_xtb_solv.csv \
    --outdir cross_dataset/analysis
"""

import argparse
from itertools import combinations
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

# Canonical 16-system ordering
SYSTEM_ORDER = [
    "me_rrrD_sap",
    "me_rrrD_tsap",
    "me_rrrL_sap",
    "me_rrrL_tsap",
    "me_sssD_sap",
    "me_sssD_tsap",
    "me_sssL_sap",
    "me_sssL_tsap",
    "phe_rrrD_sap",
    "phe_rrrD_tsap",
    "phe_rrrL_sap",
    "phe_rrrL_tsap",
    "phe_sssD_sap",
    "phe_sssD_tsap",
    "phe_sssL_sap",
    "phe_sssL_tsap",
]

REPLICATE_COLORS = ["#1f77b4", "#2ca02c", "#ff7f0e", "#d62728"]  # blue, green, orange, red


def assign_replicates(df_solv):
    """Assign replicate ID (0-3) to each solv_md frame per system.

    Uses chunk_size = N // 4 so that each replicate corresponds to one
    contiguous 100 ns block of the concatenated trajectory.

    Parameters
    ----------
    df_solv : DataFrame
        Filtered to method='solv' only.

    Returns
    -------
    df_solv : DataFrame
        With added 'replicate' column.
    """
    df_solv = df_solv.copy()
    df_solv["replicate"] = -1  # placeholder

    for sys_name in df_solv["system"].unique():
        mask = df_solv["system"] == sys_name
        n_frames = mask.sum()
        chunk_size = n_frames // 4
        frames = df_solv.loc[mask, "frame"].values
        # replicate = min(frame // chunk_size, 3)
        reps = np.minimum(frames // chunk_size, 3)
        df_solv.loc[mask, "replicate"] = reps

    return df_solv


def compute_replicate_stats(df_solv):
    """Compute per-replicate centroids, stds, and consistency ratios.

    Parameters
    ----------
    df_solv : DataFrame
        Must have columns: system, species, isomer, handedness, conformer,
        PC1, PC2, replicate.

    Returns
    -------
    stats_df : DataFrame
        One row per system with replicate stats and CR values.
    """
    records = []

    for sys_name in SYSTEM_ORDER:
        sub = df_solv[df_solv["system"] == sys_name]
        if sub.empty:
            continue

        # System metadata (constant across rows)
        meta = sub.iloc[0][["species", "isomer", "handedness", "conformer"]]

        row = {
            "system": sys_name,
            "species": meta["species"],
            "isomer": meta["isomer"],
            "handedness": meta["handedness"],
            "conformer": meta["conformer"],
        }

        # Overall stats
        overall_pc1_std = sub["PC1"].std()
        overall_pc2_std = sub["PC2"].std()
        row["overall_PC1_std"] = overall_pc1_std
        row["overall_PC2_std"] = overall_pc2_std

        # Per-replicate stats
        centroids = []  # (pc1_mean, pc2_mean) per replicate
        for rep in range(4):
            rep_sub = sub[sub["replicate"] == rep]
            if rep_sub.empty:
                row[f"rep{rep}_PC1_mean"] = np.nan
                row[f"rep{rep}_PC2_mean"] = np.nan
                row[f"rep{rep}_PC1_std"] = np.nan
                row[f"rep{rep}_PC2_std"] = np.nan
                centroids.append((np.nan, np.nan))
            else:
                pc1_mean = rep_sub["PC1"].mean()
                pc2_mean = rep_sub["PC2"].mean()
                pc1_std = rep_sub["PC1"].std()
                pc2_std = rep_sub["PC2"].std()
                row[f"rep{rep}_PC1_mean"] = pc1_mean
                row[f"rep{rep}_PC2_mean"] = pc2_mean
                row[f"rep{rep}_PC1_std"] = pc1_std
                row[f"rep{rep}_PC2_std"] = pc2_std
                centroids.append((pc1_mean, pc2_mean))

        # Max pairwise centroid distance
        max_dist_pc1 = 0.0
        max_dist_pc2 = 0.0
        for i, j in combinations(range(4), 2):
            c1_i, c2_i = centroids[i]
            c1_j, c2_j = centroids[j]
            if not (np.isnan(c1_i) or np.isnan(c1_j)):
                dist_pc1 = abs(c1_i - c1_j)
                dist_pc2 = abs(c2_i - c2_j)
                max_dist_pc1 = max(max_dist_pc1, dist_pc1)
                max_dist_pc2 = max(max_dist_pc2, dist_pc2)

        # Consistency ratio
        row["consistency_ratio_PC1"] = (
            max_dist_pc1 / overall_pc1_std if overall_pc1_std > 0 else np.nan
        )
        row["consistency_ratio_PC2"] = (
            max_dist_pc2 / overall_pc2_std if overall_pc2_std > 0 else np.nan
        )

        records.append(row)

    return pd.DataFrame(records)


def plot_replicate_consistency(stats_df, output_path):
    """Plot 4x4 grid of replicate centroids per system.

    Parameters
    ----------
    stats_df : DataFrame
        From compute_replicate_stats().
    output_path : Path
        Output PNG path.
    """
    fig, axes = plt.subplots(4, 4, figsize=(20, 20), constrained_layout=True)

    # Compute shared axis limits from all replicate centroids
    all_pc1_means = []
    all_pc2_means = []
    for rep in range(4):
        pc1_col = f"rep{rep}_PC1_mean"
        pc2_col = f"rep{rep}_PC2_mean"
        vals1 = stats_df[pc1_col].dropna().values
        vals2 = stats_df[pc2_col].dropna().values
        all_pc1_means.extend(vals1)
        all_pc2_means.extend(vals2)

    pc1_min, pc1_max = min(all_pc1_means), max(all_pc1_means)
    pc2_min, pc2_max = min(all_pc2_means), max(all_pc2_means)
    # Add 10% padding
    pc1_pad = (pc1_max - pc1_min) * 0.1 or 1.0
    pc2_pad = (pc2_max - pc2_min) * 0.1 or 1.0
    xlim = (pc1_min - pc1_pad, pc1_max + pc1_pad)
    ylim = (pc2_min - pc2_pad, pc2_max + pc2_pad)

    for idx, sys_name in enumerate(SYSTEM_ORDER):
        row, col = divmod(idx, 4)
        ax = axes[row, col]

        sys_row = stats_df[stats_df["system"] == sys_name]
        if sys_row.empty:
            ax.set_title(sys_name, fontsize=10)
            ax.text(0.5, 0.5, "No data", transform=ax.transAxes,
                    ha="center", va="center", fontsize=12, color="gray")
            continue

        cr = sys_row["consistency_ratio_PC1"].values[0]

        for rep in range(4):
            pc1_mean = sys_row[f"rep{rep}_PC1_mean"].values[0]
            pc2_mean = sys_row[f"rep{rep}_PC2_mean"].values[0]
            pc1_std = sys_row[f"rep{rep}_PC1_std"].values[0]
            pc2_std = sys_row[f"rep{rep}_PC2_std"].values[0]

            if np.isnan(pc1_mean):
                continue

            ax.errorbar(
                pc1_mean, pc2_mean,
                xerr=pc1_std if not np.isnan(pc1_std) else 0,
                yerr=pc2_std if not np.isnan(pc2_std) else 0,
                fmt="o",
                color=REPLICATE_COLORS[rep],
                markersize=8,
                capsize=3,
                label=f"Rep {rep}",
                elinewidth=1.5,
            )

        ax.set_xlim(xlim)
        ax.set_ylim(ylim)
        ax.set_title(sys_name, fontsize=10)
        # CR text annotation in lower-right corner
        ax.text(
            0.98, 0.02,
            f"CR: {cr:.2f}",
            transform=ax.transAxes,
            fontsize=8,
            ha="right", va="bottom",
            bbox=dict(boxstyle="round", facecolor="white", alpha=0.7),
        )

        if col == 0:
            ax.set_ylabel("PC2", fontsize=9)
        if row == 3:
            ax.set_xlabel("PC1", fontsize=9)

    # Shared legend — use first axis with all replicates
    handles, labels = axes[0, 0].get_legend_handles_labels()
    if handles:
        fig.legend(handles, labels, loc="upper right", fontsize=10,
                   title="Replicate", title_fontsize=11)

    fig.suptitle(
        "solv_md Replicate Consistency (eu8_nochrom PC space)",
        fontsize=14, fontweight="bold",
    )

    fig.savefig(output_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved: {output_path}")


def main():
    parser = argparse.ArgumentParser(
        description="Assess solv_md replicate convergence via PC1/PC2 centroids"
    )
    parser.add_argument(
        "--input", type=Path, required=True,
        help="Combined projection CSV (joint_projection_eu8_nochrom_xtb_solv.csv)",
    )
    parser.add_argument(
        "--outdir", type=Path, required=True,
        help="Output directory for CSV and PNG",
    )
    args = parser.parse_args()

    args.outdir.mkdir(parents=True, exist_ok=True)

    # Load data
    print(f"Loading: {args.input}")
    df = pd.read_csv(args.input)
    print(f"  Total rows: {len(df)}")

    # Filter solv only
    df_solv = df[df["method"] == "solv"].copy()
    print(f"  solv_md rows: {len(df_solv)}")
    print(f"  Systems: {sorted(df_solv['system'].unique())}")

    # Assign replicates
    print("Assigning replicate IDs...")
    df_solv = assign_replicates(df_solv)
    for sys_name in sorted(df_solv["system"].unique()):
        sub = df_solv[df_solv["system"] == sys_name]
        n = len(sub)
        chunk = n // 4
        print(f"  {sys_name}: {n} frames, chunk_size={chunk}, "
              f"reps: {sorted(sub['replicate'].unique())}")

    # Compute stats
    print("Computing replicate statistics...")
    stats_df = compute_replicate_stats(df_solv)

    # Save CSV
    csv_path = args.outdir / "solv_replicate_consistency.csv"
    stats_df.to_csv(csv_path, index=False)
    print(f"  Saved: {csv_path}")

    # Summary
    avg_cr1 = stats_df["consistency_ratio_PC1"].mean()
    avg_cr2 = stats_df["consistency_ratio_PC2"].mean()
    most_consistent = stats_df.loc[
        stats_df["consistency_ratio_PC1"].idxmin(), "system"
    ]
    least_consistent = stats_df.loc[
        stats_df["consistency_ratio_PC1"].idxmax(), "system"
    ]
    print(f"\n  Average CR(PC1): {avg_cr1:.3f}")
    print(f"  Average CR(PC2): {avg_cr2:.3f}")
    print(f"  Most consistent:  {most_consistent} "
          f"(CR={stats_df['consistency_ratio_PC1'].min():.3f})")
    print(f"  Least consistent: {least_consistent} "
          f"(CR={stats_df['consistency_ratio_PC1'].max():.3f})")

    # Plot
    print("Generating replicate consistency plot...")
    png_path = args.outdir / "plot_replicate_consistency.png"
    plot_replicate_consistency(stats_df, png_path)

    print("\nDone.")


if __name__ == "__main__":
    main()
