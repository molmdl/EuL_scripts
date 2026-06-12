#!/usr/bin/env python3
"""
HA5_plot_fes_grids.py -- Generate 4x4 FES grid plots for xtb, solv, com_md
in the allheavy PCA space.

Produces:
  - plot_fes_grid_xtb_allheavy.png
  - plot_fes_grid_solv_allheavy.png
  - plot_fes_grid_com_md_allheavy.png
  - pc_axis_limits_allheavy.json

Usage:
  python scripts/HA5_plot_fes_grids.py \
      --input analysis/joint_projection_3way_allheavy_scaffold.csv \
      --outdir analysis
"""

import argparse
import json
from pathlib import Path

import matplotlib as mpl
mpl.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.ndimage import gaussian_filter

RT = 2.479
FES_CAP = 12.0
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
    arr = np.asarray(series, dtype=float)
    n = len(arr)
    if n < 2:
        return n
    rho = np.corrcoef(arr[:-1], arr[1:])[0, 1]
    if np.isnan(rho) or rho <= 0:
        return n
    return max(1, int(n * (1 - rho) / (1 + rho)))


def compute_fes(pc1, pc2, bins=BINS, xlim=None, ylim=None):
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


def plot_fes_grid(df, method, systems, xlim, ylim, output_path, tag):
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

    if im is not None:
        cbar = fig.colorbar(im, ax=axes, shrink=0.6, label="$\Delta G$ (kJ/mol)")

    method_labels = {"xtb": "xTB", "solv": "solv_md", "com_md": "com_md (HSA-bound)"}
    method_label = method_labels.get(method, method)
    fig.suptitle(
        f"FES Grid -- {method_label} ({tag} PC space)",
        fontsize=14, fontweight="bold",
    )

    fig.savefig(output_path, dpi=DPI, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved: {output_path}")


def compute_axis_limits(df):
    return {
        "pc1_min": float(df["PC1"].min()),
        "pc1_max": float(df["PC1"].max()),
        "pc2_min": float(df["PC2"].min()),
        "pc2_max": float(df["PC2"].max()),
    }


def main():
    parser = argparse.ArgumentParser(
        description="Generate 4x4 FES grids for xtb, solv, com_md in allheavy PC space"
    )
    parser.add_argument(
        "--input", type=Path, required=True,
        help="3-way joint projection CSV from HA4",
    )
    parser.add_argument(
        "--outdir", type=Path, default=Path("analysis"),
        help="Output directory for PNGs and JSON",
    )
    parser.add_argument(
        "--methods", nargs="+", default=["xtb", "solv", "com_md"],
        choices=["xtb", "solv", "com_md"],
        help="Which methods to generate grids for (default: all three)",
    )
    parser.add_argument(
        "--tag", default="allheavy",
        help="Tag for output filenames",
    )
    args = parser.parse_args()

    args.outdir.mkdir(parents=True, exist_ok=True)

    print(f"Loading: {args.input}")
    df = pd.read_csv(args.input)
    print(f"  Total rows: {len(df)}")
    print(f"  Methods: {df['method'].unique().tolist()}")
    print(f"  Systems: {sorted(df['system'].unique())}")

    limits = compute_axis_limits(df)
    xlim = (limits["pc1_min"], limits["pc1_max"])
    ylim = (limits["pc2_min"], limits["pc2_max"])
    print(f"  PC1 range: [{xlim[0]:.2f}, {xlim[1]:.2f}]")
    print(f"  PC2 range: [{ylim[0]:.2f}, {ylim[1]:.2f}]")

    limits_path = args.outdir / f"pc_axis_limits_{args.tag}.json"
    with open(limits_path, "w") as f:
        json.dump(limits, f, indent=2)
    print(f"  Saved: {limits_path}")

    for method in args.methods:
        method_labels = {"xtb": "xTB", "solv": "solv_md", "com_md": "com_md (HSA-bound)"}
        method_label = method_labels.get(method, method)
        print(f"\nGenerating FES grid for {method_label}...")
        output_path = args.outdir / f"plot_fes_grid_{method}_{args.tag}.png"
        plot_fes_grid(df, method, SYSTEM_ORDER, xlim, ylim, output_path, args.tag)

    print("\nDone.")


if __name__ == "__main__":
    main()
