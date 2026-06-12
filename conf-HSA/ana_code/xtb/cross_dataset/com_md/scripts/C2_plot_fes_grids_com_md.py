#!/usr/bin/env python3
"""
C2_plot_fes_grids_com_md.py — Generate 4×4 FES grids for com_md only
and 3-way overlay (xtb + solv + com_md) in eu8_nochrom PC space.

Produces:
  - plot_fes_grid_com_md_only.png  (com_md HSA-bound, ~4000 frames/system)
  - plot_fes_grid_3way_overlay.png (3 subplots per system: xtb / solv / com_md)

Free energy: ΔG = -RT ln(P/P_max), RT = 2.479 kJ/mol at 298 K.
Color cap: 12 kJ/mol. Colormap: cubehelix.
N_eff annotations computed via lag-1 autocorrelation on PC1.

Usage:
  python cross_dataset/com_md/scripts/C2_plot_fes_grids_com_md.py \
      --com-input cross_dataset/com_md/analysis/com_md_projection_eu8_nochrom.csv \
      --joint-input cross_dataset/com_md/analysis/joint_projection_3way_eu8_nochrom.csv \
      --outdir cross_dataset/com_md/analysis
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
RT = 2.479   # kJ/mol at 298 K
FES_CAP = 12.0  # kJ/mol
CMAP = "cubehelix"
BINS = 50  # histogram bins per axis
SIGMA = 1.0  # Gaussian smoothing sigma (in bin units)
FIGSIZE = (20, 20)
DPI = 150

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


def compute_neff(series):
    """Compute effective sample size using lag-1 autocorrelation.

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
    """
    Compute 2D free energy surface from PC1/PC2 coordinates.

    Parameters
    ----------
    pc1, pc2 : array-like
    bins : int
    xlim, ylim : tuple of (min, max) or None

    Returns
    -------
    fes : ndarray (bins×bins)
    xedges, yedges : ndarray
    """
    range_2d = None
    if xlim is not None and ylim is not None:
        range_2d = [xlim, ylim]

    H, xedges, yedges = np.histogram2d(pc1, pc2, bins=bins, range=range_2d)

    # Gaussian smoothing
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


def compute_axis_limits_3way(df):
    """Compute shared PC1/PC2 limits from ALL data (xtb + solv + com_md)."""
    return {
        "pc1_min": float(df["PC1"].min()),
        "pc1_max": float(df["PC1"].max()),
        "pc2_min": float(df["PC2"].min()),
        "pc2_max": float(df["PC2"].max()),
    }


# ---------- Plot 1: com_md only 4×4 FES grid ----------

def plot_fes_grid_com_md(df, systems, xlim, ylim, output_path):
    """Plot 4×4 FES grid for com_md data only."""
    df_method = df[df["method"] == "com_md"].copy()

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

        fes, xedges, yedges = compute_fes(pc1, pc2, bins=BINS, xlim=xlim, ylim=ylim)
        n_eff = compute_neff(pc1)

        extent = [xedges[0], xedges[-1], yedges[0], yedges[-1]]
        im = ax.imshow(
            fes.T, origin="lower", extent=extent, aspect="auto",
            cmap=CMAP, vmin=0, vmax=FES_CAP, interpolation="bilinear",
        )

        # N_eff annotation
        ax.text(
            0.02, 0.98,
            f"N_eff ≈ {n_eff}",
            transform=ax.transAxes, fontsize=6, verticalalignment="top",
            bbox=dict(boxstyle="round", facecolor="white", alpha=0.7),
        )

        ax.set_title(sys_name, fontsize=10)
        ax.set_xlim(xlim)
        ax.set_ylim(ylim)

        if col == 0:
            ax.set_ylabel("PC2", fontsize=9)
        if row == 3:
            ax.set_xlabel("PC1", fontsize=9)

    if im is not None:
        cbar = fig.colorbar(im, ax=axes, shrink=0.6, label="ΔG (kJ/mol)")

    fig.suptitle(
        "FES Grid — com_md (HSA-bound, eu8_nochrom PC space)",
        fontsize=14, fontweight="bold",
    )

    fig.savefig(output_path, dpi=DPI, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved: {output_path}")


# ---------- Plot 2: 3-way overlay (3 subplots per system) ----------

def plot_fes_grid_3way_overlay(df, systems, xlim, ylim, output_path):
    """Plot 4×4 grid where each cell has 3 sub-rows (xtb/solv/com_md).

    Each system gets a column-pair of 3 small subplots showing the three
    methods separately, with shared axis limits and N_eff annotations.
    """
    methods = ["xtb", "solv", "com_md"]
    method_labels = {"xtb": "xTB", "solv": "solv_md", "com_md": "com_md"}
    method_colors = {"xtb": "C0", "solv": "C2", "com_md": "C3"}

    # Layout: 4 columns of systems × 4 rows of systems
    # But each cell has 3 mini-rows → effective 12 rows × 4 cols
    # That's too tall. Instead, use a 4×4 grid with 3 overlaid contours
    # OR: 4×4 grid where each cell shows 3 separate small FES panels

    # Alternative: 4×4 grid of subplots, each subplot shows all 3 methods
    # as contour overlays with transparency
    fig, axes = plt.subplots(4, 4, figsize=FIGSIZE, constrained_layout=True)

    for idx, sys_name in enumerate(systems):
        row, col = divmod(idx, 4)
        ax = axes[row, col]
        ax.set_xlim(xlim)
        ax.set_ylim(ylim)

        has_data = False

        for method in methods:
            sub = df[(df["system"] == sys_name) & (df["method"] == method)]
            if sub.empty:
                continue

            pc1 = sub["PC1"].values
            pc2 = sub["PC2"].values
            n_eff = compute_neff(pc1)

            fes, xedges, yedges = compute_fes(pc1, pc2, bins=BINS, xlim=xlim, ylim=ylim)

            extent = [xedges[0], xedges[-1], yedges[0], yedges[-1]]
            # Use contour lines for overlay clarity
            X = np.linspace(xedges[0], xedges[-1], BINS)
            Y = np.linspace(yedges[0], yedges[-1], BINS)

            # Mask NaN for contour
            fes_masked = np.ma.array(fes.T, mask=np.isnan(fes.T))

            if fes_masked.count() > 0:
                has_data = True
                # Choose contour levels
                levels = [2, 4, 6, 8, 10, 12]
                cs = ax.contour(
                    X, Y, fes_masked,
                    levels=levels,
                    colors=method_colors[method],
                    linewidths=0.8,
                    linestyles="solid",
                    alpha=0.8,
                )

        # N_eff annotations for all methods
        neff_parts = []
        for method in methods:
            sub = df[(df["system"] == sys_name) & (df["method"] == method)]
            if not sub.empty:
                ne = compute_neff(sub["PC1"].values)
                short = method_labels[method]
                neff_parts.append(f"{short}:{ne}")
        neff_text = "  ".join(neff_parts)

        ax.text(
            0.02, 0.98, neff_text,
            transform=ax.transAxes, fontsize=5, verticalalignment="top",
            bbox=dict(boxstyle="round", facecolor="white", alpha=0.7),
        )

        ax.set_title(sys_name, fontsize=10)

        if col == 0:
            ax.set_ylabel("PC2", fontsize=9)
        if row == 3:
            ax.set_xlabel("PC1", fontsize=9)

    # Add legend for methods
    from matplotlib.lines import Line2D
    legend_elements = [
        Line2D([0], [0], color=method_colors[m], lw=2, label=method_labels[m])
        for m in methods
    ]
    fig.legend(
        handles=legend_elements, loc="upper center", ncol=3,
        fontsize=11, frameon=True, bbox_to_anchor=(0.5, 1.02),
    )

    fig.suptitle(
        "FES Grid — 3-Way Overlay (eu8_nochrom PC space)\n"
        "Contours: xTB (blue) | solv_md (green) | com_md (red)",
        fontsize=14, fontweight="bold", y=1.06,
    )

    fig.savefig(output_path, dpi=DPI, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved: {output_path}")


# ---------- Main ----------

def main():
    parser = argparse.ArgumentParser(
        description="Generate FES grids for com_md and 3-way overlay"
    )
    parser.add_argument(
        "--com-input", type=Path, required=True,
        help="com_md projection CSV",
    )
    parser.add_argument(
        "--joint-input", type=Path, required=True,
        help="3-way joint projection CSV",
    )
    parser.add_argument(
        "--outdir", type=Path, required=True,
        help="Output directory for PNGs",
    )
    args = parser.parse_args()

    args.outdir.mkdir(parents=True, exist_ok=True)

    # Load com_md data
    print(f"Loading com_md: {args.com_input}")
    df_com = pd.read_csv(args.com_input)
    print(f"  {len(df_com)} rows, {df_com['system'].nunique()} systems")

    # Load 3-way joint data
    print(f"Loading 3-way joint: {args.joint_input}")
    df_joint = pd.read_csv(args.joint_input)
    print(f"  {len(df_joint)} rows")
    print(f"  Methods: {df_joint['method'].unique().tolist()}")
    print(f"  Systems: {df_joint['system'].nunique()}")

    # Compute shared axis limits from 3-way data
    limits = compute_axis_limits_3way(df_joint)
    xlim = (limits["pc1_min"], limits["pc1_max"])
    ylim = (limits["pc2_min"], limits["pc2_max"])
    print(f"  PC1 range: [{xlim[0]:.2f}, {xlim[1]:.2f}]")
    print(f"  PC2 range: [{ylim[0]:.2f}, {ylim[1]:.2f}]")

    # Save limits
    limits_path = args.outdir / "pc_axis_limits_3way.json"
    with open(limits_path, "w") as f:
        json.dump(limits, f, indent=2)
    print(f"  Saved: {limits_path}")

    # Grid 1: com_md only
    print("\nGenerating com_md-only FES grid...")
    plot_fes_grid_com_md(
        df_com, SYSTEM_ORDER, xlim, ylim,
        args.outdir / "plot_fes_grid_com_md_only.png",
    )

    # Grid 2: 3-way overlay
    print("\nGenerating 3-way overlay FES grid...")
    plot_fes_grid_3way_overlay(
        df_joint, SYSTEM_ORDER, xlim, ylim,
        args.outdir / "plot_fes_grid_3way_overlay.png",
    )

    print("\nDone.")


if __name__ == "__main__":
    main()
