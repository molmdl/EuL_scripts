#!/usr/bin/env python3
"""
rmsd_plot.py — Read RMSD CSVs and produce histogram plots as SVG files.

Four plot types:
1. Per-ligand Unbound vs HSA-bound overlay histograms (16 SVGs)
2. xtb triple-alignment histograms (16 SVGs)
3. 1×4 grid for 4 selected systems (1 SVG)
4. Cross-method overlay histograms (16 SVGs)

Usage:
  python rmsd_analysis/scripts/rmsd_plot.py \\
      --csv-dir rmsd_analysis/csv \\
      --out-dir rmsd_analysis/svg

  python rmsd_analysis/scripts/rmsd_plot.py \\
      --csv-dir rmsd_analysis/csv \\
      --out-dir rmsd_analysis/svg \\
      --plot-type per_ligand

  python rmsd_analysis/scripts/rmsd_plot.py \\
      --csv-dir rmsd_analysis/csv \\
      --out-dir rmsd_analysis/svg \\
      --plot-type xtb_triple

  python rmsd_analysis/scripts/rmsd_plot.py \\
      --csv-dir rmsd_analysis/csv \\
      --out-dir rmsd_analysis/svg \\
      --plot-type selected_4

  python rmsd_analysis/scripts/rmsd_plot.py \\
      --csv-dir rmsd_analysis/csv \\
      --out-dir rmsd_analysis/svg \\
      --plot-type cross_method
"""

import argparse
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")  # headless backend
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.interpolate import CubicSpline


# ───────────────────────────────────────────────────────────────
# Constants
# ───────────────────────────────────────────────────────────────

SYSTEMS = [
    "me_rrrD_sap", "me_rrrD_tsap", "me_rrrL_sap", "me_rrrL_tsap",
    "me_sssD_sap", "me_sssD_tsap", "me_sssL_sap", "me_sssL_tsap",
    "phe_rrrD_sap", "phe_rrrD_tsap", "phe_rrrL_sap", "phe_rrrL_tsap",
    "phe_sssD_sap", "phe_sssD_tsap", "phe_sssL_sap", "phe_sssL_tsap",
]

# 4 selected systems for the 1×4 grid
SELECTED_4 = ["phe_sssD_sap", "phe_sssD_tsap", "phe_rrrL_sap", "phe_rrrL_tsap"]

# Colors — Tableau 10 palette (colorblind-safe, high-contrast)
COLOR_SOLV_MD = "#4e79a7"     # Tableau blue
COLOR_COM_MD = "#b07aa1"     # Tableau purple (replaces red)
COLOR_XTB_HEAVY = "#4e79a7"   # Tableau blue (consistent with solv_md)
COLOR_XTB_EU9 = "#59a14f"    # Tableau green (high contrast vs blue)
COLOR_XTB_COMREF = "#e15759"  # Tableau red (for xtb comref in cross-method)


# ───────────────────────────────────────────────────────────────
# Helper: adaptive bin edges
# ───────────────────────────────────────────────────────────────

def compute_bin_edges(rmsd_values, min_bins=40, max_bins=100, target_bin_width=0.03):
    """
    Compute adaptive bin edges from RMSD data.

    Targets ~0.03 Å bin width with floor of 40 and ceiling of 100 bins.
    """
    data = np.asarray(rmsd_values)
    data = data[~np.isnan(data)]  # remove NaN
    if len(data) == 0:
        return np.linspace(0, 1, 21)

    data_min = data.min()
    data_max = data.max()
    data_range = data_max - data_min

    if data_range < 1e-10:
        return np.linspace(data_min - 0.5, data_max + 0.5, 21)

    n_bins = max(min_bins, min(max_bins, int(data_range / target_bin_width)))
    return np.linspace(data_min, data_max, n_bins + 1)


# ───────────────────────────────────────────────────────────────
# Helper: CubicSpline smooth density curves
# ───────────────────────────────────────────────────────────────

def spline_density(data, bin_edges=None, min_bins=40, max_bins=100,
                   target_bin_width=0.03, fine_points=500):
    """
    Compute smooth percentage curve via CubicSpline through histogram bin centers.

    Unlike KDE, this preserves peak heights at bin midpoints while interpolating
    smoothly (Bezier-like) between them. Values are normalized as percentage of
    total samples per bin. The curve returns to y=0 at both edges of the data
    range, anchoring the distribution to the baseline.

    Parameters
    ----------
    data : array-like
        1-D RMSD values (may contain NaN — will be dropped).
    bin_edges : np.ndarray or None
        Pre-computed bin edges. If None, compute adaptively from data.
        Shared bin edges should be passed for overlay plots so that
        both curves use the same binning and are directly comparable.
    min_bins : int
        Minimum number of histogram bins (default 40).
    max_bins : int
        Maximum number of histogram bins (default 100).
    target_bin_width : float
        Target bin width in same units as data (default 0.03 Å).
    fine_points : int
        Number of evaluation points for the smooth curve (default 500).

    Returns
    -------
    x_fine : np.ndarray
        Evenly spaced evaluation points from data min to max.
    y_fine : np.ndarray
        Smooth percentage values (0–100%), clamped to >= 0.
    """
    arr = np.asarray(data, dtype=float)
    arr = arr[~np.isnan(arr)]
    if len(arr) < 2:
        return np.array([0.0, 1.0]), np.array([0.0, 0.0])

    # Compute adaptive bin edges if not provided
    if bin_edges is None:
        bin_edges = compute_bin_edges(arr, min_bins, max_bins, target_bin_width)

    # Histogram with raw counts, then normalize to percentage
    counts, edges = np.histogram(arr, bins=bin_edges, density=False)
    # Convert raw counts to percentage of total samples
    total = counts.sum()
    if total > 0:
        counts = counts / total * 100.0

    # Bin centers — these are the spline knot points
    centers = (edges[:-1] + edges[1:]) / 2

    # Prepend and append zero endpoints so curve returns to baseline
    x_knots = np.concatenate([[edges[0]], centers, [edges[-1]]])
    y_knots = np.concatenate([[0.0], counts, [0.0]])

    # CubicSpline (piecewise cubic = Bezier-like smooth interpolation)
    spline = CubicSpline(x_knots, y_knots)

    # Evaluate on fine grid
    x_fine = np.linspace(edges[0], edges[-1], fine_points)
    y_fine = spline(x_fine)

    # Clamp negative values (density cannot be < 0)
    y_fine = np.maximum(y_fine, 0.0)

    return x_fine, y_fine


# ───────────────────────────────────────────────────────────────
# Helper: standalone legend SVG
# ───────────────────────────────────────────────────────────────

def create_legend_svg(entries, output_path, linewidth=3.0, line_length=1.5,
                      font_size=12, row_height=0.5, width_inches=2.5):
    """
    Create a standalone SVG legend file.

    Each entry is drawn as a short colored line segment (not a box/rectangle)
    followed by the label text.

    Parameters
    ----------
    entries : list of (color_hex, label_text)
        Legend entries to draw.
    output_path : Path
        Where to save the SVG.
    linewidth : float
        Width of the line segment in points (thicker than plot lines for visibility).
    line_length : float
        Length of the line segment in inches.
    font_size : float
        Label text font size in points.
    row_height : float
        Vertical spacing between entries in inches.
    width_inches : float
        Total width of the legend figure in inches.
    """
    n = len(entries)
    height = n * row_height + 0.3  # top/bottom padding

    fig, ax = plt.subplots(figsize=(width_inches, height))
    ax.set_xlim(0, width_inches)
    ax.set_ylim(0, height)
    ax.axis("off")

    y_start = height - 0.3 - row_height / 2  # top entry vertically centered

    for i, (color, label) in enumerate(entries):
        y = y_start - i * row_height
        # Draw a short colored line segment
        x_start = 0.2
        x_end = x_start + line_length
        ax.plot([x_start, x_end], [y, y], color=color,
                linewidth=linewidth, linestyle="-", solid_capstyle="round")

        # Label text to the right of the line
        ax.text(x_end + 0.15, y, label, fontsize=font_size,
                va="center", ha="left")

    fig.savefig(output_path, format="svg", bbox_inches="tight")
    plt.close(fig)


# ───────────────────────────────────────────────────────────────
# Plot type 1: Per-ligand solv_md vs com_md overlay
# ───────────────────────────────────────────────────────────────

def plot_per_ligand(solv_df: pd.DataFrame, com_df: pd.DataFrame, out_dir: Path,
                    rmsd_col: str = "rmsd_heavy_comref_A"):
    """
    Generate 16 per-ligand SVGs: Unbound vs HSA-bound histogram overlay.

    Parameters
    ----------
    rmsd_col : str
        RMSD column to plot (default: rmsd_heavy_comref_A).
    """
    per_ligand_dir = out_dir / "per_ligand"
    per_ligand_dir.mkdir(parents=True, exist_ok=True)

    n_created = 0
    for sys_name in SYSTEMS:
        solv_sub = solv_df[solv_df["system"] == sys_name][rmsd_col].values
        com_sub = com_df[com_df["system"] == sys_name][rmsd_col].values

        if len(solv_sub) == 0 or len(com_sub) == 0:
            print(f"  SKIP {sys_name}: no data for one or both methods",
                  file=sys.stderr)
            continue

        # Shared xlim from combined range
        combined = np.concatenate([solv_sub, com_sub])
        xlim = (combined.min(), combined.max())

        # Shared bin edges from combined data for direct shape comparability
        shared_bins = compute_bin_edges(combined)

        fig, ax = plt.subplots(figsize=(6, 4))

        x_solv, y_solv = spline_density(solv_sub, bin_edges=shared_bins)
        x_com, y_com = spline_density(com_sub, bin_edges=shared_bins)

        ax.plot(x_solv, y_solv, color=COLOR_SOLV_MD, linewidth=0.5, alpha=0.7, label="Unbound")
        ax.plot(x_com, y_com, color=COLOR_COM_MD, linewidth=0.5, alpha=0.7, label="HSA-bound")

        ax.set_xlabel("RMSD (Å)", fontsize=14)
        ax.set_ylabel("Probability (%)", fontsize=14)
        ax.set_title(sys_name, fontsize=13)
        ax.tick_params(labelsize=12)
        # Add small left gap between y-axis and data start
        x_range = xlim[1] - xlim[0]
        ax.set_xlim(left=xlim[0] - 0.05 * x_range, right=xlim[1])
        ax.set_ylim(bottom=0)
        ax.margins(y=0)
        ax.locator_params(axis='x', nbins=5)
        ax.locator_params(axis='y', nbins=5)
        ax.minorticks_on()
        ax.tick_params(axis='y', which='minor', length=3, width=0.5)

        svg_path = per_ligand_dir / f"{sys_name}.svg"
        fig.savefig(svg_path, format="svg", bbox_inches="tight")
        plt.close(fig)
        n_created += 1

    print(f"  Created {n_created}/16 per-ligand SVGs in {per_ligand_dir}")


# ───────────────────────────────────────────────────────────────
# Plot type 2: xtb triple-alignment histograms
# ───────────────────────────────────────────────────────────────

def plot_xtb_triple(xtb_df: pd.DataFrame, out_dir: Path,
                    rmsd_col: str = "rmsd_heavy_comref_A"):
    """
    Generate 16 xtb triple-alignment SVGs: all-heavy vs Eu+9 vs comref RMSD overlay.

    Parameters
    ----------
    rmsd_col : str
        RMSD column for general use (not used here — this plot always shows
        both rmsd_heavy_A and rmsd_heavy_comref_A for alignment comparison).
    """
    xtb_triple_dir = out_dir / "xtb_triple"
    xtb_triple_dir.mkdir(parents=True, exist_ok=True)

    n_created = 0
    for sys_name in SYSTEMS:
        sub = xtb_df[xtb_df["system"] == sys_name]
        heavy = sub["rmsd_heavy_A"].values
        eu9 = sub["rmsd_eu9_A"].values
        comref = sub["rmsd_heavy_comref_A"].values

        if len(heavy) == 0:
            print(f"  SKIP {sys_name}: no xtb data", file=sys.stderr)
            continue

        combined = np.concatenate([heavy, eu9, comref])
        xlim = (combined.min(), combined.max())

        shared_bins = compute_bin_edges(combined)

        fig, ax = plt.subplots(figsize=(6, 4))

        x_heavy, y_heavy = spline_density(heavy, bin_edges=shared_bins)
        x_eu9, y_eu9 = spline_density(eu9, bin_edges=shared_bins)
        x_comref, y_comref = spline_density(comref, bin_edges=shared_bins)

        ax.plot(x_heavy, y_heavy, color=COLOR_XTB_HEAVY, linewidth=0.5, alpha=0.7,
                label="Heavy atoms (all-heavy align)")
        ax.plot(x_eu9, y_eu9, color=COLOR_XTB_EU9, linewidth=0.5, alpha=0.7,
                label="Eu+9 (Eu+9 align)")
        ax.plot(x_comref, y_comref, color=COLOR_XTB_COMREF, linewidth=0.5, alpha=0.7,
                label="Comref (com_md ref)")

        ax.set_xlabel("RMSD (Å)", fontsize=14)
        ax.set_ylabel("Probability (%)", fontsize=14)
        ax.set_title(sys_name, fontsize=13)
        ax.tick_params(labelsize=12)
        x_range = xlim[1] - xlim[0]
        ax.set_xlim(left=xlim[0] - 0.05 * x_range, right=xlim[1])
        ax.set_ylim(bottom=0)
        ax.margins(y=0)
        ax.locator_params(axis='x', nbins=5)
        ax.locator_params(axis='y', nbins=5)
        ax.minorticks_on()
        ax.tick_params(axis='y', which='minor', length=3, width=0.5)

        svg_path = xtb_triple_dir / f"{sys_name}.svg"
        fig.savefig(svg_path, format="svg", bbox_inches="tight")
        plt.close(fig)
        n_created += 1

    print(f"  Created {n_created}/16 xtb_triple SVGs in {xtb_triple_dir}")


# ───────────────────────────────────────────────────────────────
# Plot type 3: 1×4 grid for 4 selected systems
# ───────────────────────────────────────────────────────────────

def plot_selected_4(solv_df: pd.DataFrame, com_df: pd.DataFrame, out_dir: Path,
                    rmsd_col: str = "rmsd_heavy_comref_A"):
    """
    Generate 1 SVG: 1×4 grid of Unbound vs HSA-bound overlay for
    4 selected systems.

    Parameters
    ----------
    rmsd_col : str
        RMSD column to plot (default: rmsd_heavy_comref_A).
    """
    selected_dir = out_dir / "selected_4"
    selected_dir.mkdir(parents=True, exist_ok=True)

    global_max = 0
    for sys_name in SELECTED_4:
        solv_sub = solv_df[solv_df["system"] == sys_name][rmsd_col].values
        com_sub = com_df[com_df["system"] == sys_name][rmsd_col].values
        if len(solv_sub) > 0:
            global_max = max(global_max, solv_sub.max())
        if len(com_sub) > 0:
            global_max = max(global_max, com_sub.max())

    global_xlim = (0, global_max)
    # Add small left gap between y-axis and x=0
    x_range = global_xlim[1] - global_xlim[0]
    global_xlim_padded = (global_xlim[0] - 0.05 * x_range, global_xlim[1])

    fig, axes = plt.subplots(1, 4, figsize=(6, 2), sharex=True,
                              gridspec_kw={'wspace': 0.05})

    global_ymax = 0
    for i, sys_name in enumerate(SELECTED_4):
        ax = axes[i]
        solv_sub = solv_df[solv_df["system"] == sys_name][rmsd_col].values
        com_sub = com_df[com_df["system"] == sys_name][rmsd_col].values

        if len(solv_sub) == 0 or len(com_sub) == 0:
            ax.text(0.5, 0.5, f"{sys_name}\n(no data)", transform=ax.transAxes,
                    ha="center", va="center")
            continue

        # Shared bin edges from combined data for direct shape comparability
        combined_sub = np.concatenate([solv_sub, com_sub])
        shared_bins = compute_bin_edges(combined_sub)

        x_solv, y_solv = spline_density(solv_sub, bin_edges=shared_bins)
        x_com, y_com = spline_density(com_sub, bin_edges=shared_bins)

        # Track global y-max for shared y-axis scale
        global_ymax = max(global_ymax, y_solv.max(), y_com.max())

        ax.plot(x_solv, y_solv, color=COLOR_SOLV_MD, linewidth=0.5, alpha=0.7)
        ax.plot(x_com, y_com, color=COLOR_COM_MD, linewidth=0.5, alpha=0.7)

        ax.set_xlabel("RMSD (Å)", fontsize=9)
        ax.set_title(sys_name, fontsize=10)
        ax.tick_params(axis="both", labelsize=8)
        ax.locator_params(axis='x', nbins=5)

    # Set shared x limits across all subplots
    axes[0].set_xlim(global_xlim_padded)

    # Set shared y limits on ALL 4 panels (same scale)
    # Leftmost panel keeps y-label + ticks; inner panels remove them
    for i, ax in enumerate(axes):
        if global_ymax > 0:
            ax.set_ylim(0, global_ymax * 1.05)  # 5% top margin
        else:
            ax.set_ylim(bottom=0)

        if i == 0:
            # Leftmost panel: full y-axis labels and ticks
            ax.locator_params(axis='y', nbins=5)
            ax.minorticks_on()
            ax.tick_params(axis='y', which='minor', length=3, width=0.5)
            ax.set_ylabel("Probability (%)", fontsize=10)
        else:
            # Inner panels: no y-axis label, no y-tick marks, no y-tick labels
            ax.set_ylabel("")
            ax.tick_params(axis='y', left=False, labelleft=False)
            ax.tick_params(axis='y', which='minor', left=False)

    fig.tight_layout()

    svg_path = selected_dir / "selected_4_grid.svg"
    fig.savefig(svg_path, format="svg", bbox_inches="tight")
    plt.close(fig)

    print(f"  Created {svg_path}")


# ───────────────────────────────────────────────────────────────
# Plot type 4: Cross-method overlay (xtb vs solv_md vs com_md)
# ───────────────────────────────────────────────────────────────

def plot_cross_method(xtb_df: pd.DataFrame, solv_df: pd.DataFrame,
                      com_df: pd.DataFrame, out_dir: Path):
    """
    Generate 16 cross-method SVGs: xtb vs solv_md vs com_md overlay
    using rmsd_heavy_comref_A (common reference) for direct comparability.
    """
    cross_method_dir = out_dir / "cross_method"
    cross_method_dir.mkdir(parents=True, exist_ok=True)

    n_created = 0
    for sys_name in SYSTEMS:
        xtb_sub = xtb_df[xtb_df["system"] == sys_name]["rmsd_heavy_comref_A"].values
        solv_sub = solv_df[solv_df["system"] == sys_name]["rmsd_heavy_comref_A"].values
        com_sub = com_df[com_df["system"] == sys_name]["rmsd_heavy_comref_A"].values

        if len(xtb_sub) == 0 or len(solv_sub) == 0 or len(com_sub) == 0:
            print(f"  SKIP {sys_name}: no data for one or more methods",
                  file=sys.stderr)
            continue

        combined = np.concatenate([xtb_sub, solv_sub, com_sub])
        xlim = (combined.min(), combined.max())

        shared_bins = compute_bin_edges(combined)

        fig, ax = plt.subplots(figsize=(6, 4))

        x_xtb, y_xtb = spline_density(xtb_sub, bin_edges=shared_bins)
        x_solv, y_solv = spline_density(solv_sub, bin_edges=shared_bins)
        x_com, y_com = spline_density(com_sub, bin_edges=shared_bins)

        ax.plot(x_xtb, y_xtb, color=COLOR_XTB_COMREF, linewidth=0.5, alpha=0.7,
                label="xTB")
        ax.plot(x_solv, y_solv, color=COLOR_SOLV_MD, linewidth=0.5, alpha=0.7,
                label="solv_md")
        ax.plot(x_com, y_com, color=COLOR_COM_MD, linewidth=0.5, alpha=0.7,
                label="com_md")

        ax.set_xlabel("RMSD (Å)", fontsize=14)
        ax.set_ylabel("Probability (%)", fontsize=14)
        ax.set_title(sys_name, fontsize=13)
        ax.tick_params(labelsize=12)
        x_range = xlim[1] - xlim[0]
        ax.set_xlim(left=xlim[0] - 0.05 * x_range, right=xlim[1])
        ax.set_ylim(bottom=0)
        ax.margins(y=0)
        ax.locator_params(axis='x', nbins=5)
        ax.locator_params(axis='y', nbins=5)
        ax.minorticks_on()
        ax.tick_params(axis='y', which='minor', length=3, width=0.5)

        svg_path = cross_method_dir / f"{sys_name}.svg"
        fig.savefig(svg_path, format="svg", bbox_inches="tight")
        plt.close(fig)
        n_created += 1

    print(f"  Created {n_created}/16 cross_method SVGs in {cross_method_dir}")


# ───────────────────────────────────────────────────────────────
# Main
# ───────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(
        description="Produce RMSD histogram SVG plots from computed CSV data"
    )
    parser.add_argument("--csv-dir", type=Path, default=Path("rmsd_analysis/csv"),
                        help="Directory with RMSD CSV files (default: rmsd_analysis/csv)")
    parser.add_argument("--out-dir", type=Path, default=Path("rmsd_analysis/svg"),
                        help="Output directory for SVG files (default: rmsd_analysis/svg)")
    parser.add_argument("--plot-type",
                        choices=["per_ligand", "xtb_triple", "selected_4", "cross_method", "all"],
                        default="all", help="Which plot type to generate (default: all)")
    parser.add_argument("--rmsd-col",
                        choices=["rmsd_heavy_A", "rmsd_heavy_comref_A"],
                        default="rmsd_heavy_comref_A",
                        help="RMSD column for per_ligand and selected_4 plots (default: rmsd_heavy_comref_A)")

    args = parser.parse_args()

    xtb_csv = args.csv_dir / "rmsd_xtb.csv"
    solv_csv = args.csv_dir / "rmsd_solv_md.csv"
    com_csv = args.csv_dir / "rmsd_com_md.csv"

    xtb_df = None
    solv_df = None
    com_df = None

    if args.plot_type in ("xtb_triple", "cross_method", "all"):
        if xtb_csv.exists():
            xtb_df = pd.read_csv(xtb_csv)
            print(f"Loaded {xtb_csv}: {len(xtb_df)} rows")
        else:
            print(f"WARNING: {xtb_csv} not found", file=sys.stderr)

    if args.plot_type in ("per_ligand", "selected_4", "cross_method", "all"):
        if solv_csv.exists():
            solv_df = pd.read_csv(solv_csv)
            print(f"Loaded {solv_csv}: {len(solv_df)} rows")
        else:
            print(f"WARNING: {solv_csv} not found", file=sys.stderr)

        if com_csv.exists():
            com_df = pd.read_csv(com_csv)
            print(f"Loaded {com_csv}: {len(com_df)} rows")
        else:
            print(f"WARNING: {com_csv} not found", file=sys.stderr)

    if args.plot_type in ("per_ligand", "all") and solv_df is not None and com_df is not None:
        print("\n--- Per-ligand Unbound vs HSA-bound overlays ---")
        plot_per_ligand(solv_df, com_df, args.out_dir, rmsd_col=args.rmsd_col)

    if args.plot_type in ("xtb_triple", "all") and xtb_df is not None:
        print("\n--- xtb triple-alignment histograms ---")
        plot_xtb_triple(xtb_df, args.out_dir, rmsd_col=args.rmsd_col)

    if args.plot_type in ("selected_4", "all") and solv_df is not None and com_df is not None:
        print("\n--- 1×4 grid for selected 4 systems ---")
        plot_selected_4(solv_df, com_df, args.out_dir, rmsd_col=args.rmsd_col)

    if args.plot_type in ("cross_method", "all") and xtb_df is not None and solv_df is not None and com_df is not None:
        print("\n--- Cross-method overlay histograms ---")
        plot_cross_method(xtb_df, solv_df, com_df, args.out_dir)

    print("\n--- Creating standalone legend SVGs ---")

    create_legend_svg(
        entries=[
            (COLOR_SOLV_MD, "Unbound"),
            (COLOR_COM_MD, "HSA-bound"),
        ],
        output_path=args.out_dir / "legend_per_ligand.svg",
    )

    create_legend_svg(
        entries=[
            (COLOR_XTB_HEAVY, "Heavy atoms (all-heavy align)"),
            (COLOR_XTB_EU9, "Eu+9 (Eu+9 align)"),
            (COLOR_XTB_COMREF, "Comref (com_md ref)"),
        ],
        output_path=args.out_dir / "legend_xtb_triple.svg",
    )

    create_legend_svg(
        entries=[
            (COLOR_SOLV_MD, "Unbound"),
            (COLOR_COM_MD, "HSA-bound"),
        ],
        output_path=args.out_dir / "legend_selected_4.svg",
    )

    create_legend_svg(
        entries=[
            (COLOR_XTB_COMREF, "xTB"),
            (COLOR_SOLV_MD, "solv_md"),
            (COLOR_COM_MD, "com_md"),
        ],
        output_path=args.out_dir / "legend_cross_method.svg",
    )

    print("  Created legend SVGs in", args.out_dir)

    print("\n=== Validation ===")
    per_ligand_dir = args.out_dir / "per_ligand"
    xtb_triple_dir = args.out_dir / "xtb_triple"
    selected_dir = args.out_dir / "selected_4"
    cross_method_dir = args.out_dir / "cross_method"

    n_per_ligand = len(list(per_ligand_dir.glob("*.svg"))) if per_ligand_dir.exists() else 0
    n_xtb_triple = len(list(xtb_triple_dir.glob("*.svg"))) if xtb_triple_dir.exists() else 0
    n_selected = len(list(selected_dir.glob("*.svg"))) if selected_dir.exists() else 0
    n_cross_method = len(list(cross_method_dir.glob("*.svg"))) if cross_method_dir.exists() else 0

    print(f"  per_ligand:   {n_per_ligand}/16 SVGs")
    print(f"  xtb_triple:   {n_xtb_triple}/16 SVGs")
    print(f"  selected_4:   {n_selected}/1 SVGs")
    print(f"  cross_method: {n_cross_method}/16 SVGs")

    total_svg = n_per_ligand + n_xtb_triple + n_selected + n_cross_method
    all_ok = True
    for svg_dir in [per_ligand_dir, xtb_triple_dir, selected_dir, cross_method_dir]:
        if not svg_dir.exists():
            continue
        for svg_path in svg_dir.glob("*.svg"):
            size = svg_path.stat().st_size
            if size < 1024:
                print(f"  WARN: {svg_path.name} is {size} bytes (< 1 KB)")
                all_ok = False
            content = svg_path.read_text()[:200]
            if "<svg" not in content:
                print(f"  FAIL: {svg_path.name} missing <svg root tag")
                all_ok = False

    if all_ok and total_svg == 49:
        print("  All SVGs validated: size > 1 KB, <svg root present")
    elif not all_ok:
        print("  Some SVGs failed validation")

    for legend_name in ["legend_per_ligand.svg", "legend_xtb_triple.svg",
                        "legend_selected_4.svg", "legend_cross_method.svg"]:
        legend_path = args.out_dir / legend_name
        if legend_path.exists():
            size = legend_path.stat().st_size
            content = legend_path.read_text()[:200]
            if "<svg" in content and size > 500:
                print(f"  {legend_name}: OK ({size} bytes)")
            else:
                print(f"  WARN: {legend_name} validation failed")
                all_ok = False
        else:
            print(f"  FAIL: {legend_name} not found")
            all_ok = False

    print("\nDone.")


if __name__ == "__main__":
    main()
