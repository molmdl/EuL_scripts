#!/usr/bin/env python3
"""
T6_pca_fes.py
=============
PCA + Free Energy Surfaces for Eu(III) complex MD trajectories.

Performs PCA on heavy-atom coordinates, aligned on Eu+donor atoms,
separately for `me` (8 trajectories, 135 atoms) and `phe` (8 trajectories, 149 atoms).

Outputs:
  analysis/me_pca_projection.csv
  analysis/phe_pca_projection.csv
  analysis/plot_me_pca_scatter.png
  analysis/plot_phe_pca_scatter.png
  analysis/plot_me_fes.png
  analysis/plot_phe_fes.png

Usage:
  python scripts/T6_pca_fes.py
"""

import argparse
import json
import os
import sys
import warnings
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd
import MDAnalysis as mda
from sklearn.decomposition import PCA

warnings.filterwarnings("ignore", category=UserWarning)

# ───────────────────────────────────────────────────────────────
# Constants
# ───────────────────────────────────────────────────────────────
KT_KJ_PER_MOL = 2.479  # kJ/mol at 298 K
N_PCS = 10
COLORS = {"SAP": "#1f77b4", "TSAP": "#ff7f0e", "UNK": "#2ca02c"}
MARKERS = {"sap": "o", "tsap": "s"}


def kabsch(P, Q):
    """
    Kabsch algorithm.
    P, Q: (N, 3) numpy arrays of corresponding points.
    Returns optimal rotation matrix R (3, 3).
    """
    # Center
    P_centroid = P.mean(axis=0)
    Q_centroid = Q.mean(axis=0)
    Pc = P - P_centroid
    Qc = Q - Q_centroid
    # Covariance matrix
    H = Pc.T @ Qc
    U, S, Vt = np.linalg.svd(H)
    d = np.sign(np.linalg.det(U @ Vt))
    D = np.diag([1, 1, d])
    R = U @ D @ Vt
    return R, P_centroid, Q_centroid


def load_and_align(system_dir, align_indices, heavy_sel):
    """
    Load an XYZ trajectory via MDAnalysis, align frames to the first frame
    on `align_indices`, return ndarray of shape (n_frames, n_heavy*3).

    Parameters
    ----------
    system_dir : str
        Path to system directory containing traj.xyz.
    align_indices : list of int
        Atom indices (0-based) used for alignment.
    heavy_sel : MDAnalysis AtomGroup
        The heavy-atom selection for this universe.

    Returns
    -------
    coords : ndarray, shape (n_frames, n_heavy * 3)
    """
    traj_path = os.path.join(system_dir, "traj.xyz")
    u = mda.Universe(traj_path, format="XYZ")

    # Ensure align_indices are valid
    n_atoms = u.atoms.n_atoms
    align_idx = [int(i) for i in align_indices if int(i) < n_atoms]
    if len(align_idx) != len(align_indices):
        missing = set(align_indices) - set(align_idx)
        print(f"  [WARN] Ignoring out-of-range align indices: {missing}",
              file=sys.stderr)

    ref_positions = u.trajectory[0].positions.copy()
    ref_align = ref_positions[align_idx]

    n_frames = len(u.trajectory)
    n_heavy = heavy_sel.n_atoms
    coords = np.empty((n_frames, n_heavy * 3), dtype=np.float64)

    for i, ts in enumerate(u.trajectory):
        pos = ts.positions.copy()
        mobile_align = pos[align_idx]
        R, P_cent, Q_cent = kabsch(mobile_align, ref_align)
        # Apply rotation + translation to ALL heavy atoms
        heavy_pos = pos[heavy_sel.indices]
        aligned = (heavy_pos - P_cent) @ R + Q_cent
        coords[i] = aligned.ravel()

    return coords


def run_pca(system_names, data_dir, output_dir, ligand_label,
            indices_map, class_df):
    """
    Load, align, fit PCA, save projections, and generate plots
    for a set of system names belonging to one ligand type (me or phe).

    Parameters
    ----------
    system_names : list of str
    data_dir : str
    output_dir : str
    ligand_label : str ('me' or 'phe')
    indices_map : dict from JSON (keys 'sap', 'tsap')
    class_df : DataFrame from torsion_classification.csv
    """
    print(f"\n{'='*60}")
    print(f"PCA + FES for ligand: {ligand_label}")
    print(f"{'='*60}")

    all_coords = []
    metadata = []  # list of (system_name, frame)

    for sys_name in system_names:
        sys_dir = os.path.join(data_dir, sys_name)
        if not os.path.isdir(sys_dir):
            print(f"  [SKIP] Directory not found: {sys_dir}", file=sys.stderr)
            continue

        # Determine conformer type from system name
        if sys_name.endswith("_sap"):
            conf_type = "sap"
        elif sys_name.endswith("_tsap"):
            conf_type = "tsap"
        else:
            print(f"  [WARN] Cannot infer conformer type for {sys_name}",
                  file=sys.stderr)
            conf_type = "sap"  # fallback

        align_indices = indices_map[conf_type]["coord"]
        # Also include metal if present
        metal_idx = indices_map[conf_type].get("metal")
        if metal_idx is not None:
            align_indices = list(np.unique([metal_idx] + align_indices))

        # Load a scratch universe just to get heavy-atom selection
        traj_path = os.path.join(sys_dir, "traj.xyz")
        tmp_u = mda.Universe(traj_path, format="XYZ")
        heavy_sel = tmp_u.select_atoms("not element H")
        print(f"  {sys_name}: {len(tmp_u.trajectory)} frames, "
              f"{heavy_sel.n_atoms} heavy atoms, align_atoms={len(align_indices)}")

        coords = load_and_align(sys_dir, align_indices, heavy_sel)
        all_coords.append(coords)
        for f in range(coords.shape[0]):
            metadata.append((sys_name, f))

    # Concatenate
    X = np.vstack(all_coords)
    n_frames, n_features = X.shape
    print(f"  Total dataset: {n_frames} frames x {n_features} features")

    # ── PCA ────────────────────────────────────────────────────
    pca = PCA(n_components=N_PCS, svd_solver="randomized", random_state=42)
    proj = pca.fit_transform(X)
    evr = pca.explained_variance_ratio_
    print(f"  Explained variance ratio (PC1-PC3): {evr[0]:.4f}, "
          f"{evr[1]:.4f}, {evr[2]:.4f}")
    cumsum = evr[:2].sum()
    print(f"  PC1+PC2 cumulative: {cumsum:.4f} ({cumsum*100:.1f}%)")
    if cumsum < 0.20:
        print(f"  [WARN] PC1+PC2 < 20%  ({cumsum*100:.1f}%). "
              "FES may be less meaningful.", file=sys.stderr)

    # ── Save projection CSV ───────────────────────────────────
    df_proj = pd.DataFrame({
        "system": [m[0] for m in metadata],
        "frame": [m[1] for m in metadata],
    })
    pc_cols = [f"pc{i+1}" for i in range(N_PCS)]
    df_proj = pd.concat([df_proj, pd.DataFrame(proj, columns=pc_cols)], axis=1)
    csv_path = os.path.join(output_dir, f"{ligand_label}_pca_projection.csv")
    df_proj.to_csv(csv_path, index=False)
    print(f"  Saved: {csv_path}")

    # ── Merge classification ───────────────────────────────────
    df_merged = df_proj.merge(
        class_df[["system", "frame", "classification"]],
        on=["system", "frame"],
        how="left",
    )
    df_merged["classification"] = df_merged["classification"].fillna("UNK")
    df_merged["start_conf"] = df_merged["system"].apply(
        lambda s: "sap" if s.endswith("_sap") else "tsap"
    )

    pc1 = df_merged["pc1"].values
    pc2 = df_merged["pc2"].values

    # ── Scatter plot PC1 vs PC2 ────────────────────────────────
    fig, ax = plt.subplots(figsize=(6, 5), dpi=150)
    for cls in ["SAP", "TSAP", "UNK"]:
        for conf in ["sap", "tsap"]:
            mask = (
                (df_merged["classification"] == cls)
                & (df_merged["start_conf"] == conf)
            )
            if mask.sum() == 0:
                continue
            ax.scatter(
                pc1[mask],
                pc2[mask],
                c=COLORS[cls],
                marker=MARKERS[conf],
                alpha=0.3,
                s=10,
                label=f"{cls} ({conf})",
                edgecolors="none",
            )
    ax.set_xlabel("PC1")
    ax.set_ylabel("PC2")
    ax.set_title(f"PCA: {ligand_label} systems (PC1 vs PC2)")
    # Deduplicate legend (matplotlib will group by label)
    handles, labels = ax.get_legend_handles_labels()
    by_label = dict(zip(labels, handles))
    ax.legend(by_label.values(), by_label.keys(),
              loc="best", fontsize=7, ncol=1, markerscale=1.5)
    ax.set_xlim(pc1.min() - 1, pc1.max() + 1)
    ax.set_ylim(pc2.min() - 1, pc2.max() + 1)
    fig.tight_layout()
    scatter_path = os.path.join(output_dir, f"plot_{ligand_label}_pca_scatter.png")
    fig.savefig(scatter_path, dpi=300)
    plt.close(fig)
    print(f"  Saved: {scatter_path}")

    # ── FES on PC1/PC2 ─────────────────────────────────────────
    H, xedges, yedges = np.histogram2d(pc1, pc2, bins=50)
    H = H + 1e-10
    G = -KT_KJ_PER_MOL * np.log(H)
    G -= G.min()  # shift so minimum is 0 kJ/mol
    G = np.clip(G, None, 15)  # cap at 15 kJ/mol per AGENTS.md convention
    # Center of bins
    X = 0.5 * (xedges[:-1] + xedges[1:])
    Y = 0.5 * (yedges[:-1] + yedges[1:])
    X, Y = np.meshgrid(X, Y)

    fig, ax = plt.subplots(figsize=(6, 5), dpi=150)
    cf = ax.contourf(X, Y, G.T, levels=20, cmap="cubehelix")
    cbar = fig.colorbar(cf, ax=ax)
    cbar.set_label("Free Energy (kJ/mol)")
    ax.set_xlabel("PC1")
    ax.set_ylabel("PC2")
    ax.set_title(f"Free Energy Surface: {ligand_label} (PC1 vs PC2)")
    fig.tight_layout()
    fes_path = os.path.join(output_dir, f"plot_{ligand_label}_fes.png")
    fig.savefig(fes_path, dpi=300)
    plt.close(fig)
    print(f"  Saved: {fes_path}")

    return evr


def main():
    parser = argparse.ArgumentParser(
        description="PCA + FES for Eu(III) complex MD trajectories",
    )
    parser.add_argument(
        "--data-dir", default="data",
        help="Directory containing system subfolders (default: data)",
    )
    parser.add_argument(
        "--out-dir", default="analysis",
        help="Directory for outputs (default: analysis)",
    )
    parser.add_argument(
        "--indices", default="analysis/indices.json",
        help="Path to indices.json (default: analysis/indices.json)",
    )
    parser.add_argument(
        "--classification", default="analysis/torsion_classification.csv",
        help="Path to torsion classification CSV (default: analysis/torsion_classification.csv)",
    )
    parser.add_argument(
        "--ligand", choices=["me", "phe", "all"], default="all",
        help="Which ligand subset to process (default: all)",
    )
    args = parser.parse_args()

    os.makedirs(args.out_dir, exist_ok=True)

    # ── Load indices ───────────────────────────────────────────
    with open(args.indices, "r") as f:
        indices_map = json.load(f)
    print(f"Loaded indices from {args.indices}")

    # ── Load classification CSV ────────────────────────────────
    class_df = pd.read_csv(args.classification)
    print(f"Loaded classification: {len(class_df)} rows")

    # ── Enumerate systems ──────────────────────────────────────
    systems = sorted([
        d for d in os.listdir(args.data_dir)
        if os.path.isdir(os.path.join(args.data_dir, d))
        and not d.startswith(".")
    ])
    me_systems = [s for s in systems if s.startswith("me_")]
    phe_systems = [s for s in systems if s.startswith("phe_")]
    print(f"Found {len(me_systems)} me systems, {len(phe_systems)} phe systems")

    results = {}
    if args.ligand in ("me", "all"):
        evr_me = run_pca(
            me_systems, args.data_dir, args.out_dir,
            "me", indices_map, class_df,
        )
        results["me"] = evr_me

    if args.ligand in ("phe", "all"):
        evr_phe = run_pca(
            phe_systems, args.data_dir, args.out_dir,
            "phe", indices_map, class_df,
        )
        results["phe"] = evr_phe

    # ── Summary ────────────────────────────────────────────────
    print(f"\n{'='*60}")
    print("SUMMARY")
    print(f"{'='*60}")
    for ligand, evr in results.items():
        print(f"  {ligand}: PC1={evr[0]:.4f}, PC2={evr[1]:.4f}, "
              f"PC3={evr[2]:.4f}, PC1+PC2={evr[:2].sum()*100:.1f}%")


if __name__ == "__main__":
    main()
