#!/usr/bin/env python3
"""
A6_variance_decomposition.py — Decompose PCA variance by atom group.

Reads the eu8_nochrom loadings CSV and computes per-group fractional
contributions to each PC's variance. Produces:
  - analysis/variance_decomposition_eu8_nochrom.csv
  - analysis/plot_variance_decomposition_eu8_nochrom.png (PC1-3)
  - analysis/plot_variance_decomposition_eu8_nochrom_full.png (PC1-10)

Usage:
    python scripts/A6_variance_decomposition.py \
        --loadings analysis/joint_pca_loadings_eu8_nochrom.csv \
        --out-dir analysis \
        --tag eu8_nochrom \
        --n-pcs 3
"""

import argparse
import os
import shutil

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


# ───────────────────────────────────────────────────────────────
# Constants
# ───────────────────────────────────────────────────────────────
ATOM_GROUPS = {
    "Core (Eu+8)":  {0, 1, 2, 3, 4, 5, 6, 7, 54},
    "Macrocycle":    set(range(8, 28)),
    "Linker":        {55, 56, 57, 58, 59, 63, 64, 65},
    "Pendant":       set(range(66, 73)),
    "Ring A/B":      {73, 74, 75, 76, 80, 84, 91, 92, 93, 94, 96, 98},
}

GROUP_ORDER = ["Core (Eu+8)", "Macrocycle", "Linker", "Pendant", "Ring A/B"]

COLORS = {
    "Core (Eu+8)":   "#1f77b4",  # blue
    "Macrocycle":    "#ff7f0e",  # orange
    "Linker":        "#2ca02c",  # green
    "Pendant":       "#d62728",  # red
    "Ring A/B":      "#9467bd",  # purple
}


def backup_if_exists(path):
    """Back up an existing file with .bak suffix before overwriting."""
    if os.path.exists(path):
        backup_path = path + ".bak"
        shutil.copy2(path, backup_path)
        print(f"  Backed up existing file to {backup_path}")


def parse_args():
    parser = argparse.ArgumentParser(
        description="Variance decomposition by atom group for eu8_nochrom PCA")
    parser.add_argument("--loadings",
                        default="analysis/joint_pca_loadings_eu8_nochrom.csv",
                        help="Path to loadings CSV")
    parser.add_argument("--out-dir", default="analysis",
                        help="Output directory")
    parser.add_argument("--tag", default="eu8_nochrom",
                        help="Output file tag")
    parser.add_argument("--n-pcs", type=int, default=3,
                        help="Number of PCs to include in plot (default: 3)")
    parser.add_argument("--projection",
                        default=None,
                        help="Path to projection CSV for computing cumulative variance. "
                             "Default: analysis/joint_pca_projection_{tag}.csv")
    return parser.parse_args()


def main():
    args = parse_args()

    print("=" * 65)
    print(f"A6: Variance Decomposition by Atom Group ({args.tag})")
    print("=" * 65)

    # ── 1. Load loadings CSV ──────────────────────────────────
    print(f"\n[1/4] Loading loadings CSV: {args.loadings} ...")
    df = pd.read_csv(args.loadings)
    atom_indices = df["atom_index"].values
    n_atoms = len(atom_indices)
    print(f"  Loaded {n_atoms} feature atoms")

    # Identify PC columns
    pc_cols = [c for c in df.columns if c.startswith("pc")]
    n_pcs_total = len(pc_cols)
    print(f"  PC columns: {pc_cols}")

    # ── 2. Verify atom group assignment ───────────────────────
    print(f"\n[2/4] Verifying atom group assignments ...")

    # Build mapping: atom_index -> group name
    atom_to_group = {}
    all_group_atoms = set()
    for gname, gatoms in ATOM_GROUPS.items():
        all_group_atoms.update(gatoms)
        for aidx in gatoms:
            atom_to_group[aidx] = gname

    # Check which atoms in the loadings CSV are in groups
    in_group = [aidx for aidx in atom_indices if aidx in atom_to_group]
    unmapped = [aidx for aidx in atom_indices if aidx not in atom_to_group]

    print(f"  Atoms in groups: {len(in_group)}")
    print(f"  Unmapped atoms:  {len(unmapped)}")
    if unmapped:
        print(f"  [WARN] Unmapped atom indices: {unmapped}")

    # Per-group atom count among feature atoms
    print(f"\n  Per-group atom counts (in feature set):")
    for gname in GROUP_ORDER:
        gatoms = ATOM_GROUPS[gname]
        in_feat = gatoms & set(atom_indices)
        print(f"    {gname:15s}: {len(gatoms):3d} total, "
              f"{len(in_feat):3d} in features")

    # ── 3. Compute per-group variance contributions ────────────
    print(f"\n[3/4] Computing per-group variance contributions ...")

    # Loadings are magnitudes: sqrt(dx^2 + dy^2 + dz^2)
    # Squaring gives dx^2 + dy^2 + dz^2 = per-atom contribution to unit eigenvector norm
    # Since components are unit vectors, sum of squared loadings = 1.0 per PC

    results = {}
    for gname in GROUP_ORDER:
        gatoms = ATOM_GROUPS[gname]
        gfeat = [aidx for aidx in atom_indices if aidx in gatoms]
        gmask = df["atom_index"].isin(gatoms).values
        results[gname] = {}
        for pc_col in pc_cols:
            # Sum of squared loadings for this group
            group_sq = (df.loc[gmask, pc_col] ** 2).sum()
            # Total sum of squared loadings (should be ~1.0)
            total_sq = (df[pc_col] ** 2).sum()
            results[gname][pc_col] = group_sq / total_sq

    # Also compute total for verification
    total_row = {}
    for pc_col in pc_cols:
        total_row[pc_col] = sum(results[gname][pc_col] for gname in GROUP_ORDER)

    # Build output DataFrame
    rows = []
    for gname in GROUP_ORDER:
        row = {"group": gname}
        for pc_col in pc_cols:
            row[pc_col] = round(results[gname][pc_col], 4)
        rows.append(row)
    # Total row
    row = {"group": "total"}
    for pc_col in pc_cols:
        row[pc_col] = round(total_row[pc_col], 4)
    rows.append(row)

    out_df = pd.DataFrame(rows)

    # Print per-PC group fractions table
    print(f"\n  Per-group fractional contributions:")
    print(f"  {'Group':15s}", end="")
    for pc_col in pc_cols[:5]:
        print(f"  {pc_col:>8s}", end="")
    print()
    for _, row in out_df.iterrows():
        print(f"  {row['group']:15s}", end="")
        for pc_col in pc_cols[:5]:
            print(f"  {row[pc_col]:>8.4f}", end="")
        print()

    # Verify sums
    print(f"\n  Total per PC (should be ~1.0000):")
    for pc_col in pc_cols[:5]:
        print(f"    {pc_col}: {total_row[pc_col]:.4f}")

    # Top-contributing group for PC1-3
    for pc_col in pc_cols[:3]:
        best_group = max(GROUP_ORDER, key=lambda g: results[g][pc_col])
        best_frac = results[best_group][pc_col]
        print(f"  {pc_col} top group: {best_group} ({best_frac*100:.1f}%)")

    # Save CSV
    out_dir = args.out_dir
    os.makedirs(out_dir, exist_ok=True)
    csv_path = os.path.join(out_dir, f"variance_decomposition_{args.tag}.csv")
    backup_if_exists(csv_path)
    out_df.to_csv(csv_path, index=False)
    print(f"\n  Saved {csv_path}")

    # ── 4. Generate stacked bar plot ───────────────────────────
    print(f"\n[4/4] Generating stacked bar plots ...")

    n_pcs = args.n_pcs
    plot_pcs = pc_cols[:n_pcs]
    pc_labels = [f"PC{i+1}" for i in range(n_pcs)]

    # Extract fractions for plotting (exclude total row)
    plot_data = {}
    for gname in GROUP_ORDER:
        plot_data[gname] = [results[gname][pc] for pc in plot_pcs]

    # --- Main plot (PC1 to n_pcs) ---
    fig, ax = plt.subplots(figsize=(8, 5))

    y_positions = np.arange(n_pcs)
    left = np.zeros(n_pcs)

    for gname in GROUP_ORDER:
        values = np.array(plot_data[gname])
        ax.barh(y_positions, values, left=left, height=0.6,
                label=gname, color=COLORS[gname], edgecolor="white", linewidth=0.5)
        # Annotate segments > 5%
        for i, v in enumerate(values):
            if v > 0.05:
                ax.text(left[i] + v / 2, y_positions[i],
                        f"{v*100:.1f}%", ha="center", va="center",
                        fontsize=9, fontweight="bold", color="white")
        left += values

    ax.set_yticks(y_positions)
    ax.set_yticklabels(pc_labels)
    ax.set_xlabel("Fractional contribution")
    ax.set_xlim(0, 1.0)
    ax.invert_yaxis()  # PC1 at top

    # Compute cumulative variance for subtitle from projection CSV
    proj_path = args.projection
    if proj_path is None:
        # Derive from tag
        proj_path = os.path.join(args.out_dir, f"joint_pca_projection_{args.tag}.csv")
    cum_var_str = ""
    if os.path.exists(proj_path):
        try:
            proj_df = pd.read_csv(proj_path)
            pc_cols_proj = [c for c in proj_df.columns if c.startswith("pc")]
            if len(pc_cols_proj) >= 2:
                # Column variances = eigenvalues; cumulative = sum(PC1+PC2) / sum(all)
                eigs = [proj_df[c].var(ddof=1) for c in pc_cols_proj]
                cum_var = (eigs[0] + eigs[1]) / sum(eigs)
                cum_var_str = f"PC1+PC2 = {cum_var*100:.1f}% cumulative variance"
            else:
                cum_var_str = f"PC1+PC2 (variance unknown)"
        except Exception:
            cum_var_str = f"PC1+PC2 (variance unknown)"
    else:
        cum_var_str = f"PC1+PC2 (projection CSV not found)"
    ax.set_title(f"Variance Decomposition by Atom Group ({args.tag})\n"
                 f"{cum_var_str}",
                 fontsize=12)

    ax.legend(loc="lower right", fontsize=8, framealpha=0.9)
    fig.tight_layout()
    plot_path = os.path.join(out_dir, f"plot_variance_decomposition_{args.tag}.png")
    backup_if_exists(plot_path)
    fig.savefig(plot_path, dpi=300)
    plt.close(fig)
    print(f"  Saved {plot_path}")

    # --- Full plot (PC1-10) ---
    n_full = min(n_pcs_total, 10)
    full_pcs = pc_cols[:n_full]
    full_labels = [f"PC{i+1}" for i in range(n_full)]

    full_data = {}
    for gname in GROUP_ORDER:
        full_data[gname] = [results[gname][pc] for pc in full_pcs]

    fig2, ax2 = plt.subplots(figsize=(10, 7))
    y_pos = np.arange(n_full)
    left = np.zeros(n_full)

    for gname in GROUP_ORDER:
        values = np.array(full_data[gname])
        ax2.barh(y_pos, values, left=left, height=0.6,
                 label=gname, color=COLORS[gname], edgecolor="white", linewidth=0.5)
        for i, v in enumerate(values):
            if v > 0.05:
                ax2.text(left[i] + v / 2, y_pos[i],
                         f"{v*100:.1f}%", ha="center", va="center",
                         fontsize=8, fontweight="bold", color="white")
        left += values

    ax2.set_yticks(y_pos)
    ax2.set_yticklabels(full_labels)
    ax2.set_xlabel("Fractional contribution")
    ax2.set_xlim(0, 1.0)
    ax2.invert_yaxis()
    ax2.set_title(f"Variance Decomposition by Atom Group ({args.tag})\n"
                  f"PC1-PC10", fontsize=12)
    ax2.legend(loc="lower right", fontsize=8, framealpha=0.9)
    fig2.tight_layout()
    full_plot_path = os.path.join(out_dir,
                                   f"plot_variance_decomposition_{args.tag}_full.png")
    backup_if_exists(full_plot_path)
    fig2.savefig(full_plot_path, dpi=300)
    plt.close(fig2)
    print(f"  Saved {full_plot_path}")

    # ── Summary ────────────────────────────────────────────────
    print("\n" + "=" * 65)
    print("VARIANCE DECOMPOSITION SUMMARY")
    print("=" * 65)
    print(f"  Feature atoms:     {n_atoms}")
    print(f"  Atoms in groups:   {len(in_group)}")
    print(f"  Unmapped atoms:    {len(unmapped)}")
    print(f"\n  Per-group atom counts:")
    for gname in GROUP_ORDER:
        gatoms = ATOM_GROUPS[gname]
        in_feat = gatoms & set(atom_indices)
        print(f"    {gname:15s}: {len(in_feat):2d} / {len(gatoms):2d}")
    print(f"\n  PC1-3 group contributions:")
    for pc_col in pc_cols[:3]:
        print(f"    {pc_col}:")
        for gname in GROUP_ORDER:
            pct = results[gname][pc_col] * 100
            print(f"      {gname:15s}: {pct:5.1f}%")
    print(f"\n  Outputs:")
    print(f"    {csv_path}")
    print(f"    {plot_path}")
    print(f"    {full_plot_path}")
    print("=" * 65)


if __name__ == "__main__":
    main()
