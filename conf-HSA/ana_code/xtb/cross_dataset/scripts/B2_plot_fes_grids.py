#!/usr/bin/env python3
"""
B2_plot_fes_grids.py — Generate 4×4 FES grids for xtb and solv_md
projected into the same eu8_nochrom PC space.

Produces:
  - plot_fes_grid_xtb.png  (xTB 100ps, 2000 frames/system)
  - plot_fes_grid_solv.png  (solv_md 400ns, ~4000 frames/system)
  - pc_axis_limits.json     (shared PC1/PC2 min/max)

Each subplot includes an N_eff annotation showing the effective sample
size (computed via lag-1 autocorrelation on PC1).

Usage:
  python B2_plot_fes_grids.py \\
    --input cross_dataset/analysis/joint_projection_eu8_nochrom_xtb_solv.csv \\
    --outdir cross_dataset/analysis

Free energy: ΔG = -RT ln(P/P_max), RT = 2.479 kJ/mol at 298 K.
Color cap: 12 kJ/mol. Colormap: cubehelix.
"""

import argparse
import json
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.ndimage import gaussian_filter


# ---------- Constants ----------
RT = 2.479  # kJ/mol at 298 K
FES_CAP = 12.0  # kJ/mol
CMAP = "cubehelix"
BINS = 50  # histogram bins per axis
SIGMA = 1.0  # Gaussian smoothing sigma (in bin units)
FIGSIZE = (20, 20)
DPI = 150

# Canonical 16-system ordering (species × isomer × handedness × conformer)
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
    """Compute effective sample size using lag-1 autocorrelation.

    For a time series with N observations and lag-1 autocorrelation rho,
    the effective sample size is N_eff = N * (1 - rho) / (1 + rho).

    Parameters
    ----------
    series : array-like
        1D array of observations (e.g., PC1 trajectory).

    Returns
    -------
    int
        Effective sample size (minimum 1).
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
    """
    Compute 2D free energy surface from PC1/PC2 coordinates.

    Parameters
    ----------
    pc1, pc2 : array-like
        PC1 and PC2 coordinate arrays.
    bins : int
        Number of bins per axis for 2D histogram.
    xlim, ylim : tuple of (min, max) or None
        Axis limits for the histogram range.

    Returns
    -------
    fes : ndarray (bins×bins)
        Free energy surface in kJ/mol. Values > FES_CAP set to NaN.
        Zero = minimum free energy.
    xedges, yedges : ndarray
        Histogram bin edges.
    """
    range_2d = None
    if xlim is not None and ylim is not None:
        range_2d = [xlim, ylim]

    H, xedges, yedges = np.histogram2d(pc1, pc2, bins=bins, range=range_2d)

    # Gaussian smoothing
    H_smooth = gaussian_filter(H.astype(float), sigma=SIGMA)

    # Convert to probability, then free energy
    H_max = H_smooth.max()
    if H_max == 0:
        fes = np.full_like(H_smooth, np.nan)
    else:
        with np.errstate(divide="ignore"):
            fes = -RT * np.log(H_smooth / H_max)
        # Cap at FES_CAP
        fes[fes > FES_CAP] = np.nan
        # Zero-bin regions (no data) → NaN
        fes[H_smooth == 0] = np.nan

    return fes, xedges, yedges


def plot_fes_grid(df, method, systems, xlim, ylim, output_path):
    """
    Plot a 4×4 grid of FES (PC1 vs PC2) for all 16 systems.

    Parameters
    ----------
    df : DataFrame
        Combined projection data (must have 'method' column).
    method : str
        Filter: "xtb" or "solv".
    systems : list of str
        System names in grid order (row-major).
    xlim, ylim : tuple of (min, max)
        Shared axis limits for all subplots.
    output_path : Path
        Output PNG file path.
    """
    df_method = df[df["method"] == method].copy()
    n = len(systems)

    fig, axes = plt.subplots(4, 4, figsize=FIGSIZE, constrained_layout=True)

    im = None
    for idx, sys_name in enumerate(systems):
        row, col = divmod(idx, 4)
        ax = axes[row, col]

        sub = df_method[df_method["system"] == sys_name]
        if sub.empty:
            ax.set_title(sys_name, fontsize=10)
            ax.set_xlim(xlim)
            ax.set_ylim(ylim)
            ax.text(
                0.5, 0.5, "No data", transform=ax.transAxes,
                ha="center", va="center", fontsize=12, color="gray",
            )
            continue

        pc1 = sub["PC1"].values
        pc2 = sub["PC2"].values

        fes, xedges, yedges = compute_fes(
            pc1, pc2, bins=BINS, xlim=xlim, ylim=ylim
        )

        # Compute effective sample size for this subplot
        n_eff = compute_neff(pc1)

        # imshow: origin lower so PC2 increases upward
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

        # N_eff annotation in upper-left corner
        ax.text(
            0.02, 0.98,
            f"N_eff ≈ {n_eff}",
            transform=ax.transAxes,
            fontsize=6,
            verticalalignment="top",
            bbox=dict(boxstyle="round", facecolor="white", alpha=0.7),
        )

        ax.set_title(sys_name, fontsize=10)
        ax.set_xlim(xlim)
        ax.set_ylim(ylim)

        # Label outer axes only
        if col == 0:
            ax.set_ylabel("PC2", fontsize=9)
        if row == 3:
            ax.set_xlabel("PC1", fontsize=9)

    # Shared colorbar
    if im is not None:
        cbar = fig.colorbar(im, ax=axes, shrink=0.6, label="ΔG (kJ/mol)")

    method_label = "xTB" if method == "xtb" else "solv_md"
    fig.suptitle(
        f"FES Grid — {method_label} (eu8_nochrom PC space)",
        fontsize=14, fontweight="bold",
    )

    fig.savefig(output_path, dpi=DPI, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved: {output_path}")


def compute_axis_limits(df):
    """
    Compute shared PC1/PC2 limits from all data (xtb + solv).

    Returns
    -------
    dict with keys: pc1_min, pc1_max, pc2_min, pc2_max
    """
    return {
        "pc1_min": float(df["PC1"].min()),
        "pc1_max": float(df["PC1"].max()),
        "pc2_min": float(df["PC2"].min()),
        "pc2_max": float(df["PC2"].max()),
    }


def main():
    parser = argparse.ArgumentParser(
        description="Generate 4×4 FES grids for xtb and solv_md in shared PC space"
    )
    parser.add_argument(
        "--input", type=Path, required=True,
        help="Combined projection CSV (joint_projection_eu8_nochrom_xtb_solv.csv)",
    )
    parser.add_argument(
        "--outdir", type=Path, required=True,
        help="Output directory for PNGs and JSON",
    )
    parser.add_argument(
        "--methods", nargs="+", default=["xtb", "solv"],
        choices=["xtb", "solv"],
        help="Which methods to generate grids for (default: both)",
    )
    args = parser.parse_args()

    args.outdir.mkdir(parents=True, exist_ok=True)

    # Load data
    print(f"Loading: {args.input}")
    df = pd.read_csv(args.input)
    print(f"  Total rows: {len(df)}")
    print(f"  Methods: {df['method'].unique()}")
    print(f"  Systems: {sorted(df['system'].unique())}")

    # Compute shared axis limits from ALL data
    limits = compute_axis_limits(df)
    xlim = (limits["pc1_min"], limits["pc1_max"])
    ylim = (limits["pc2_min"], limits["pc2_max"])
    print(f"  PC1 range: [{xlim[0]:.2f}, {xlim[1]:.2f}]")
    print(f"  PC2 range: [{ylim[0]:.2f}, {ylim[1]:.2f}]")

    # Save limits JSON
    limits_path = args.outdir / "pc_axis_limits.json"
    with open(limits_path, "w") as f:
        json.dump(limits, f, indent=2)
    print(f"  Saved: {limits_path}")

    # Generate grids
    for method in args.methods:
        method_label = "xTB" if method == "xtb" else "solv_md"
        print(f"\nGenerating FES grid for {method_label}...")
        output_path = args.outdir / f"plot_fes_grid_{method}.png"
        plot_fes_grid(df, method, SYSTEM_ORDER, xlim, ylim, output_path)

    print("\nDone.")


if __name__ == "__main__":
    main()
