#!/usr/bin/env python3
"""
B2_plot_fes_grids_zoomed.py — Generate 4x4 FES grid for xtb data with
zoomed axis limits (1st-99th percentile + 10% padding).

The existing B2_plot_fes_grids.py uses shared limits from both xtb + solv
data, which compresses xtb substructure into tiny dots. This script zooms
into the xtb-relevant PC range to reveal sub-basin structure.

Produces:
  - plot_fes_grid_xtb_zoomed.png

Usage:
  python B2_plot_fes_grids_zoomed.py \
    --input cross_dataset/analysis/joint_projection_eu8_nochrom_xtb_solv.csv \
    --outdir cross_dataset/analysis
"""

import argparse
import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.ndimage import gaussian_filter


# ---------- Constants ----------
RT = 2.479  # kJ/mol at 298 K
FES_CAP = 12.0  # kJ/mol
CMAP = "cubehelix"
BINS = 50
SIGMA = 1.0
FIGSIZE = (20, 20)
DPI = 150

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


def compute_neff(series):
    """Effective sample size via lag-1 autocorrelation.

    N_eff = N * (1 - rho) / (1 + rho).
    """
    arr = np.asarray(series, dtype=float)
    n = len(arr)
    if n < 2:
        return n
    rho = np.corrcoef(arr[:-1], arr[1:])[0, 1]
    if np.isnan(rho) or rho <= 0:
        return n
    return max(1, int(n * (1 - rho) / (1 + rho)))


def compute_fes(pc1, pc2, bins=BINS, xlim=None, ylim=None):
    """2D free energy surface: ΔG = -RT ln(P/P_max)."""
    range_2d = None
    if xlim is not None and ylim is not None:
        range_2d = [xlim, ylim]

    H, xedges, yedges = np.histogram2d(pc1, pc2, bins=bins, range=range_2d)
    H_smooth = gaussian_filter(H.astype(float), sigma=SIGMA)
    H_max = H_smooth.max()
    if H_max == 0:
        fes = np.full_like(H_smooth, np.nan)
    else:
        with np.errstate(divide="ignore"):
            fes = -RT * np.log(H_smooth / H_max)
        fes[fes > FES_CAP] = np.nan
        fes[H_smooth == 0] = np.nan

    return fes, xedges, yedges


def compute_zoomed_limits(df, method="xtb", pctile_low=1, pctile_high=99, padding=0.10):
    """Compute zoomed axis limits from 1st-99th percentile + padding.

    Parameters
    ----------
    df : DataFrame
        Combined projection data.
    method : str
        Method to compute limits from.
    pctile_low, pctile_high : int
        Percentile bounds (default 1-99).
    padding : float
        Fraction of range to add as padding (default 0.10).

    Returns
    -------
    dict with pc1_min, pc1_max, pc2_min, pc2_max.
    """
    sub = df[df["method"] == method]
    pc1_lo, pc1_hi = np.percentile(sub["PC1"], [pctile_low, pctile_high])
    pc2_lo, pc2_hi = np.percentile(sub["PC2"], [pctile_low, pctile_high])
    pc1_pad = (pc1_hi - pc1_lo) * padding
    pc2_pad = (pc2_hi - pc2_lo) * padding
    return {
        "pc1_min": float(pc1_lo - pc1_pad),
        "pc1_max": float(pc1_hi + pc1_pad),
        "pc2_min": float(pc2_lo - pc2_pad),
        "pc2_max": float(pc2_hi + pc2_pad),
    }


def plot_fes_grid_zoomed(df, systems, xlim, ylim, output_path):
    """Plot 4x4 FES grid with zoomed xtb limits."""
    df_xtb = df[df["method"] == "xtb"].copy()

    fig, axes = plt.subplots(4, 4, figsize=FIGSIZE, constrained_layout=True)
    im = None

    for idx, sys_name in enumerate(systems):
        row, col = divmod(idx, 4)
        ax = axes[row, col]

        sub = df_xtb[df_xtb["system"] == sys_name]
        if sub.empty:
            ax.set_title(sys_name, fontsize=10)
            ax.set_xlim(xlim)
            ax.set_ylim(ylim)
            ax.text(0.5, 0.5, "No data", transform=ax.transAxes,
                    ha="center", va="center", fontsize=12, color="gray")
            continue

        pc1 = sub["PC1"].values
        pc2 = sub["PC2"].values

        fes, xedges, yedges = compute_fes(pc1, pc2, bins=BINS, xlim=xlim, ylim=ylim)

        n_eff = compute_neff(pc1)

        extent = [xedges[0], xedges[-1], yedges[0], yedges[-1]]
        im = ax.imshow(
            fes.T,
            origin="lower",
            extent=extent,
            aspect="auto",
            cmap=CMAP,
            vmin=0,
            vmax=FES_CAP,
            interpolation="bilinear",
        )

        # N_eff annotation
        ax.text(
            0.02, 0.98,
            f"N_eff = {n_eff}",
            transform=ax.transAxes,
            fontsize=6,
            verticalalignment="top",
            bbox=dict(boxstyle="round", facecolor="white", alpha=0.7),
        )

        ax.set_title(sys_name, fontsize=10)
        ax.set_xlim(xlim)
        ax.set_ylim(ylim)

        if col == 0:
            ax.set_ylabel("PC2", fontsize=9)
        if row == 3:
            ax.set_xlabel("PC1", fontsize=9)

    # Shared colorbar
    if im is not None:
        cbar = fig.colorbar(im, ax=axes, shrink=0.6, label="ΔG (kJ/mol)")

    fig.suptitle(
        "FES Grid — xTB zoomed (eu8_nochrom PC space, 1st-99th pctile limits)",
        fontsize=14, fontweight="bold",
    )

    fig.savefig(output_path, dpi=DPI, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved: {output_path}")


def main():
    parser = argparse.ArgumentParser(
        description="Generate zoomed 4x4 FES grid for xtb data"
    )
    parser.add_argument(
        "--input", type=Path, required=True,
        help="Combined projection CSV",
    )
    parser.add_argument(
        "--outdir", type=Path, required=True,
        help="Output directory",
    )
    parser.add_argument(
        "--pctile-low", type=int, default=1,
        help="Lower percentile for axis limits (default: 1)",
    )
    parser.add_argument(
        "--pctile-high", type=int, default=99,
        help="Upper percentile for axis limits (default: 99)",
    )
    parser.add_argument(
        "--padding", type=float, default=0.10,
        help="Fractional padding on range (default: 0.10)",
    )
    args = parser.parse_args()

    args.outdir.mkdir(parents=True, exist_ok=True)

    print(f"Loading: {args.input}")
    df = pd.read_csv(args.input)
    print(f"  Total rows: {len(df)}")

    # Compute zoomed limits from xtb data
    limits = compute_zoomed_limits(
        df, method="xtb",
        pctile_low=args.pctile_low, pctile_high=args.pctile_high,
        padding=args.padding,
    )
    xlim = (limits["pc1_min"], limits["pc1_max"])
    ylim = (limits["pc2_min"], limits["pc2_max"])
    print(f"  Zoomed PC1: [{xlim[0]:.2f}, {xlim[1]:.2f}]")
    print(f"  Zoomed PC2: [{ylim[0]:.2f}, {ylim[1]:.2f}]")

    # Compare with full limits
    full_pc1 = (float(df["PC1"].min()), float(df["PC1"].max()))
    full_pc2 = (float(df["PC2"].min()), float(df["PC2"].max()))
    print(f"  Full   PC1: [{full_pc1[0]:.2f}, {full_pc1[1]:.2f}]")
    print(f"  Full   PC2: [{full_pc2[0]:.2f}, {full_pc2[1]:.2f}]")
    zoom_factor = (full_pc1[1] - full_pc1[0]) / (xlim[1] - xlim[0])
    print(f"  Zoom factor: {zoom_factor:.1f}x narrower")

    # Save limits JSON
    limits_path = args.outdir / "xtb_zoomed_limits.json"
    with open(limits_path, "w") as f:
        json.dump(limits, f, indent=2)
    print(f"  Saved: {limits_path}")

    # Generate zoomed FES grid
    print("\nGenerating zoomed xtb FES grid...")
    output_path = args.outdir / "plot_fes_grid_xtb_zoomed.png"
    plot_fes_grid_zoomed(df, SYSTEM_ORDER, xlim, ylim, output_path)

    print("\nDone.")


if __name__ == "__main__":
    main()
