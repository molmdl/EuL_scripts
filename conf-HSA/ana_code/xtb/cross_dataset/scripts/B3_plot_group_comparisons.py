#!/usr/bin/env python3
"""
B3_plot_group_comparisons.py — Enhanced group comparison bar charts
with bootstrap 95% CI error bars for me_vs_phe and D_vs_L metrics.

Reads xtb_vs_solv_group_metrics.csv and produces:
  - plot_me_vs_phe_comparison.png
  - plot_D_vs_L_comparison.png

Each plot shows centroid distance (PC units) per group, with xtb (blue)
and solv (red) bars side-by-side, and CI error bars.

Usage:
  python B3_plot_group_comparisons.py \
    --input cross_dataset/analysis/xtb_vs_solv_group_metrics.csv \
    --outdir cross_dataset/analysis
"""

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


def plot_group_comparison(df, metric_type, output_path):
    """Plot side-by-side bar chart with CI error bars.

    Parameters
    ----------
    df : DataFrame
        Group metrics data.
    metric_type : str
        One of 'me_vs_phe' or 'D_vs_L'.
    output_path : Path
        Output PNG path.
    """
    sub = df[df["metric_type"] == metric_type].copy()
    if sub.empty:
        print(f"  No data for metric_type={metric_type}, skipping.")
        return

    groups = sub["group_key"].values
    n_groups = len(groups)
    x = np.arange(n_groups)
    bar_width = 0.35

    # Centroid distances
    d_solv = sub["d_solv"].values
    d_xtb = sub["d_xtb"].values

    # CI error bar half-widths
    ci_lo_solv = sub["d_ci_low_solv"].values
    ci_hi_solv = sub["d_ci_high_solv"].values
    ci_lo_xtb = sub["d_ci_low_xtb"].values
    ci_hi_xtb = sub["d_ci_high_xtb"].values

    err_solv = np.vstack([d_solv - ci_lo_solv, ci_hi_solv - d_solv])
    err_xtb = np.vstack([d_xtb - ci_lo_xtb, ci_hi_xtb - d_xtb])

    fig, ax = plt.subplots(figsize=(max(10, n_groups * 1.5), 6))

    # Use logarithmic scale for D_vs_L (extreme ratio) or linear for me_vs_phe
    if metric_type == "D_vs_L":
        ax.set_yscale("log")
        ylabel = "Centroid distance (PC units, log scale)"
    else:
        ylabel = "Centroid distance (PC units)"

    bars_solv = ax.bar(
        x - bar_width / 2, d_solv, bar_width,
        label="solv_md", color="#d62728", alpha=0.85,
        yerr=err_solv, capsize=3, error_kw={"linewidth": 1.2},
    )
    bars_xtb = ax.bar(
        x + bar_width / 2, d_xtb, bar_width,
        label="xTB", color="#1f77b4", alpha=0.85,
        yerr=err_xtb, capsize=3, error_kw={"linewidth": 1.2},
    )

    ax.set_xticks(x)
    ax.set_xticklabels(groups, rotation=30, ha="right", fontsize=10)
    ax.set_ylabel(ylabel, fontsize=11)
    ax.set_xlabel("Group", fontsize=11)
    ax.legend(fontsize=10)
    ax.grid(axis="y", alpha=0.3, linestyle="--")

    # Title
    title_map = {
        "me_vs_phe": "me vs phe: Centroid Distance Comparison (95% CI)",
        "D_vs_L": "D vs L: Centroid Distance Comparison (95% CI)",
    }
    ax.set_title(title_map.get(metric_type, metric_type), fontsize=13, fontweight="bold")

    fig.tight_layout()
    fig.savefig(output_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved: {output_path}")


def main():
    parser = argparse.ArgumentParser(
        description="Generate enhanced group comparison bar charts with CIs"
    )
    parser.add_argument(
        "--input", type=Path, required=True,
        help="Group metrics CSV (xtb_vs_solv_group_metrics.csv)",
    )
    parser.add_argument(
        "--outdir", type=Path, required=True,
        help="Output directory for PNGs",
    )
    parser.add_argument(
        "--metric-types", nargs="+",
        default=["me_vs_phe", "D_vs_L"],
        choices=["me_vs_phe", "D_vs_L", "sap_vs_tsap"],
        help="Which metric types to plot (default: me_vs_phe, D_vs_L)",
    )
    args = parser.parse_args()

    args.outdir.mkdir(parents=True, exist_ok=True)

    print(f"Loading: {args.input}")
    df = pd.read_csv(args.input)
    print(f"  Rows: {len(df)}")
    print(f"  Metric types: {df['metric_type'].unique()}")

    for mt in args.metric_types:
        print(f"\nPlotting {mt}...")
        output_path = args.outdir / f"plot_{mt}_comparison.png"
        plot_group_comparison(df, mt, output_path)

    print("\nDone.")


if __name__ == "__main__":
    main()
