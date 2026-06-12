#!/usr/bin/env python3
"""
T3_joint_pca.py — Joint PCA on all 16 xTB trajectories in one PC space.
[DEPRECATED] This script uses Eu+9 alignment (10 atoms incl. cap N63).
It has been superseded by T3_joint_pca_eu8.py which uses Eu+8 alignment
(9 atoms: indices [0,1,2,3,4,5,6,7,54], excluding N63). The Eu+8
revision was motivated by N63 asymmetry between SAP and TSAP conformers
(see Phase 02 analysis). This script is retained for reproducibility only;
do NOT use its outputs as the primary PCA results.

Uses common scaffold heavy atoms (68 atoms × 3 = 204 features) after
TSAP→SAP atom reordering and Eu+9 donor alignment to fit a single PCA
that enables direct cross-species and cross-conformer comparison.

Outputs:
  analysis/joint_pca_projection.csv         — PC projections for all 32,000 frames
  analysis/joint_pca_loadings.csv           — per-atom loading magnitudes
  analysis/joint_pca_projection_scaffold_all.csv — sub-analysis with H atoms (381 feats)
  analysis/joint_pca_projection_peripheral.csv  — sub-analysis without core (174 feats)
  analysis/plot_joint_pca_scatter.png       — PC1 vs PC2 by SAP/TSAP/UNK
  analysis/plot_joint_pca_species.png       — PC1 vs PC2 by me/phe
  analysis/plot_joint_pca_conformer.png     — PC1 vs PC2 by sap/tsap start
  analysis/plot_joint_fes.png               — Free energy surface

Usage:
  python scripts/T3_joint_pca.py \
      --data-dir data \
      --out-dir analysis \
      --mappings analysis/atom_mappings.json \
      --scaffold analysis/scaffold_mapping.json \
      --indices analysis/indices.json \
      --classification analysis/torsion_classification.csv
"""

import argparse
import json
import os
import shutil
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
KT_KJ_PER_MOL = 2.479   # kJ/mol at 298 K
N_PCS = 10
COLORS_CLASS = {"SAP": "#1f77b4", "TSAP": "#ff7f0e", "UNK": "#2ca02c"}
COLORS_SPECIES = {"me": "#008080", "phe": "#9467bd"}
COLORS_CONF = {"sap": "#1f77b4", "tsap": "#ff7f0e"}
MARKERS_SPECIES = {"me": "o", "phe": "s"}
MARKERS_CONF = {"sap": "o", "tsap": "s"}


# ───────────────────────────────────────────────────────────────
# Helper functions
# ───────────────────────────────────────────────────────────────

def kabsch(P, Q):
    """
    Kabsch algorithm for optimal rotation.

    Parameters
    ----------
    P, Q : ndarray, shape (N, 3)
        Mobile and reference coordinates.

    Returns
    -------
    R : ndarray (3, 3)
        Optimal rotation matrix.
    P_centroid, Q_centroid : ndarray (3,)
        Centroids of P and Q.
    """
    P_centroid = P.mean(axis=0)
    Q_centroid = Q.mean(axis=0)
    Pc = P - P_centroid
    Qc = Q - Q_centroid
    H = Pc.T @ Qc
    U, S, Vt = np.linalg.svd(H)
    d = np.sign(np.linalg.det(U @ Vt))
    D = np.diag([1, 1, d])
    R = U @ D @ Vt
    return R, P_centroid, Q_centroid


def backup_if_exists(path):
    """Back up an existing file with .bak suffix before overwriting."""
    if os.path.exists(path):
        backup_path = path + ".bak"
        shutil.copy2(path, backup_path)
        print(f"  Backed up existing file to {backup_path}")


def load_xyz_trajectory(traj_path):
    """
    Load XYZ trajectory via MDAnalysis.

    Returns
    -------
    coords : ndarray, shape (n_frames, n_atoms, 3)
    elements : list of str
    """
    u = mda.Universe(traj_path, format="XYZ")
    n_frames = len(u.trajectory)
    n_atoms = u.atoms.n_atoms
    coords = np.empty((n_frames, n_atoms, 3), dtype=np.float64)
    for i, ts in enumerate(u.trajectory):
        coords[i] = ts.positions
    elements = [a.name for a in u.atoms]
    return coords, elements


def apply_tsap_mapping(coords, elements, mapping):
    """
    Reorder TSAP coords/elements to SAP ordering.

    mapping[i] = SAP index that TSAP atom i corresponds to.
    Apply: mapped[:, mapping[j], :] = raw[:, j, :] for each j
    This places TSAP atom j at SAP position mapping[j].

    Returns
    -------
    mapped_coords : ndarray, same shape as coords
    mapped_elements : list of str
    """
    mapped_coords = np.empty_like(coords)
    mapped_coords[:, mapping, :] = coords  # mapped[mapping[j]] = raw[j]
    # Same permutation for elements: mapped_elements[mapping[j]] = elements[j]
    mapped_elements_arr = [None] * len(elements)
    for j, sap_idx in enumerate(mapping):
        mapped_elements_arr[sap_idx] = elements[j]
    return mapped_coords, mapped_elements_arr


def validate_scaffold_elements(mapped_elements, ref_elements, label):
    """
    Assert scaffold elements (0–126) match reference after TSAP→SAP mapping.
    """
    scaffold_elems = mapped_elements[:127]
    if scaffold_elems == ref_elements:
        print(f"  [PASS] {label}: scaffold elements match reference SAP")
    else:
        mismatches = [(i, scaffold_elems[i], ref_elements[i])
                      for i in range(127) if scaffold_elems[i] != ref_elements[i]]
        print(f"  [FAIL] {label}: {len(mismatches)} scaffold element mismatches:")
        for idx, got, expected in mismatches[:10]:
            print(f"    Index {idx}: got {got}, expected {expected}")
        raise ValueError(f"{label}: scaffold element mismatch after mapping")


# ───────────────────────────────────────────────────────────────
# Main pipeline
# ───────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(
        description="Joint PCA on all 16 xTB trajectories in one PC space")
    parser.add_argument("--data-dir", default="data",
                        help="Directory containing system subdirectories")
    parser.add_argument("--out-dir", default="analysis",
                        help="Output directory for CSVs and plots")
    parser.add_argument("--mappings", default="analysis/atom_mappings.json",
                        help="Path to atom_mappings.json")
    parser.add_argument("--scaffold", default="analysis/scaffold_mapping.json",
                        help="Path to scaffold_mapping.json")
    parser.add_argument("--indices", default="analysis/indices.json",
                        help="Path to indices.json")
    parser.add_argument("--classification", default="analysis/torsion_classification.csv",
                        help="Path to torsion_classification.csv")
    args = parser.parse_args()

    data_dir = args.data_dir
    out_dir = args.out_dir

    os.makedirs(out_dir, exist_ok=True)

    # ── 1. Load mapping files ───────────────────────────────────
    print("=" * 65)
    print("T3 Joint PCA — All 16 xTB Trajectories in One PC Space")
    print("=" * 65)

    print("\n[1/7] Loading mapping and classification files ...")

    if not os.path.exists(args.mappings):
        sys.exit(f"ERROR: mappings file not found: {args.mappings}")
    atom_mappings = json.load(open(args.mappings))

    if not os.path.exists(args.scaffold):
        sys.exit(f"ERROR: scaffold mapping file not found: {args.scaffold}")
    scaffold_mapping = json.load(open(args.scaffold))

    if not os.path.exists(args.indices):
        sys.exit(f"ERROR: indices file not found: {args.indices}")
    indices_data = json.load(open(args.indices))

    if not os.path.exists(args.classification):
        sys.exit(f"ERROR: classification file not found: {args.classification}")
    class_df = pd.read_csv(args.classification)

    # Alignment indices: Eu (54) + 9 donors from indices.json sap.coord + sap.metal
    # sap.coord = [0,1,2,3,4,5,6,7,63], sap.metal = 54
    sap_coord = indices_data["sap"]["coord"]
    sap_metal = indices_data["sap"]["metal"]
    ALIGN_INDICES = sorted(set(sap_coord + [sap_metal]))
    # = [0, 1, 2, 3, 4, 5, 6, 7, 54, 63]
    print(f"  Alignment indices (SAP): {ALIGN_INDICES}")

    # Scaffold atoms in SAP ordering (same for me and phe: 0–126)
    n_scaffold = scaffold_mapping["n_scaffold"]
    print(f"  Scaffold atoms (SAP): 0–{n_scaffold - 1} ({n_scaffold} atoms)")

    # ── 2. Identify systems ─────────────────────────────────────
    print("\n[2/7] Identifying systems ...")

    systems = sorted([d for d in os.listdir(data_dir)
                      if os.path.isdir(os.path.join(data_dir, d))
                      and not d.startswith(".")
                      and os.path.exists(os.path.join(data_dir, d, "traj.xyz"))])
    print(f"  Found {len(systems)} systems: {systems}")

    if len(systems) != 16:
        print(f"  [WARN] Expected 16 systems, found {len(systems)}")

    # ── 3. Load, map, validate, and collect all trajectories ────
    print("\n[3/7] Loading trajectories and applying TSAP→SAP mapping ...")

    ref_elements = None          # set from first SAP system
    all_systems = []             # list of dicts with coords, elements, metadata

    for sys_name in systems:
        traj_path = os.path.join(data_dir, sys_name, "traj.xyz")
        if not os.path.exists(traj_path):
            print(f"  [SKIP] {sys_name}: traj.xyz not found")
            continue

        species = "me" if sys_name.startswith("me_") else "phe"
        start_conf = "sap" if sys_name.endswith("_sap") else "tsap"

        print(f"  Loading {sys_name} ({species}, start={start_conf}) ...", end=" ")
        coords, elements = load_xyz_trajectory(traj_path)
        n_frames = coords.shape[0]
        n_atoms = coords.shape[1]
        print(f"{n_frames} frames, {n_atoms} atoms")

        # TSAP→SAP reordering
        if start_conf == "tsap":
            mapping_key = f"{species}_tsap_to_sap"
            if mapping_key not in atom_mappings:
                print(f"  [ERROR] Mapping key '{mapping_key}' not found in atom_mappings.json")
                sys.exit(1)
            mapping = atom_mappings[mapping_key]
            if len(mapping) != n_atoms:
                print(f"  [ERROR] Mapping length {len(mapping)} != n_atoms {n_atoms} for {sys_name}")
                sys.exit(1)
            coords, elements = apply_tsap_mapping(coords, elements, mapping)

            # Validate scaffold elements
            if ref_elements is not None:
                validate_scaffold_elements(elements, ref_elements, sys_name)
        else:
            # SAP system — set reference if first
            if ref_elements is None:
                ref_elements = elements[:127]
                print(f"  Set reference scaffold elements ({len(ref_elements)} atoms)")
            else:
                # Still validate
                validate_scaffold_elements(elements, ref_elements, sys_name)

        all_systems.append({
            "name": sys_name,
            "species": species,
            "start_conf": start_conf,
            "coords": coords,       # (n_frames, n_atoms, 3) in SAP order
            "elements": elements,    # element list in SAP order
            "n_frames": n_frames,
            "n_atoms": n_atoms,
        })

    total_frames = sum(s["n_frames"] for s in all_systems)
    print(f"\n  Total frames loaded: {total_frames}")

    # ── 4. Identify scaffold heavy atom indices ──────────────────
    print("\n[4/7] Identifying scaffold heavy atoms ...")

    SCAFFOLD_HEAVY_IDX = [i for i in range(n_scaffold) if ref_elements[i] != "H"]
    print(f"  Scaffold heavy atoms: {len(SCAFFOLD_HEAVY_IDX)} "
          f"(expecting 68; 68 × 3 = {len(SCAFFOLD_HEAVY_IDX) * 3} features)")

    if len(SCAFFOLD_HEAVY_IDX) != 68:
        print(f"  [WARN] Expected 68 heavy scaffold atoms, got {len(SCAFFOLD_HEAVY_IDX)}")

    # Count element types for verification
    heavy_elems = [ref_elements[i] for i in SCAFFOLD_HEAVY_IDX]
    elem_counts = {}
    for e in heavy_elems:
        elem_counts[e] = elem_counts.get(e, 0) + 1
    print(f"  Heavy element breakdown: {elem_counts}")

    # Sub-analysis indices
    ALL_SCAFFOLD_IDX = list(range(n_scaffold))  # 127 atoms including H
    CORE_SET = set(ALIGN_INDICES)
    PERIPHERAL_IDX = [i for i in SCAFFOLD_HEAVY_IDX if i not in CORE_SET]
    print(f"  All-scaffold features: {n_scaffold * 3}")
    print(f"  Peripheral heavy atoms: {len(PERIPHERAL_IDX)} "
          f"({len(PERIPHERAL_IDX) * 3} features)")

    # ── 5. Align all frames to reference ────────────────────────
    print("\n[5/7] Kabsch-aligning all frames to reference (me_rrrD_sap frame 0) ...")

    # Reference = first system (sorted: me_rrrD_sap), frame 0
    ref_frame = all_systems[0]["coords"][0]       # (n_atoms, 3)
    ref_align_pos = ref_frame[ALIGN_INDICES]       # (10, 3)

    for sys_data in all_systems:
        n_frames = sys_data["n_frames"]
        aligned = np.empty_like(sys_data["coords"])
        for f in range(n_frames):
            mobile = sys_data["coords"][f]
            mobile_align = mobile[ALIGN_INDICES]
            R, P_cent, Q_cent = kabsch(mobile_align, ref_align_pos)
            aligned[f] = (mobile - P_cent) @ R + Q_cent
        sys_data["aligned"] = aligned
        # Free raw coords to save memory
        del sys_data["coords"]
        print(f"  Aligned {sys_data['name']} ({n_frames} frames)")

    # ── 6. Extract features, concatenate, fit PCA ───────────────
    print("\n[6/7] Extracting features and fitting joint PCA ...")

    # Primary: scaffold heavy atoms
    all_features = []
    all_features_scaffold_all = []
    all_features_peripheral = []
    all_metadata = []

    for sys_data in all_systems:
        # Heavy scaffold features (68 × 3 = 204)
        heavy_coords = sys_data["aligned"][:, SCAFFOLD_HEAVY_IDX, :]  # (n_frames, 68, 3)
        features = heavy_coords.reshape(sys_data["n_frames"], -1)     # (n_frames, 204)
        all_features.append(features)

        # All-scaffold features (127 × 3 = 381)
        all_scaff_coords = sys_data["aligned"][:, ALL_SCAFFOLD_IDX, :]
        features_scaffold_all = all_scaff_coords.reshape(sys_data["n_frames"], -1)
        all_features_scaffold_all.append(features_scaffold_all)

        # Peripheral heavy features (~58 × 3 = ~174)
        periph_coords = sys_data["aligned"][:, PERIPHERAL_IDX, :]
        features_periph = periph_coords.reshape(sys_data["n_frames"], -1)
        all_features_peripheral.append(features_periph)

        # Metadata
        for f in range(sys_data["n_frames"]):
            all_metadata.append((sys_data["name"], f,
                                 sys_data["species"], sys_data["start_conf"]))

    X = np.vstack(all_features)
    X_all = np.vstack(all_features_scaffold_all)
    X_periph = np.vstack(all_features_peripheral)

    print(f"  Primary feature matrix: {X.shape}")
    print(f"  All-scaffold feature matrix: {X_all.shape}")
    print(f"  Peripheral feature matrix: {X_periph.shape}")

    # ── Fit PCA ──────────────────────────────────────────────────
    pca = PCA(n_components=N_PCS, svd_solver="randomized", random_state=42)
    proj = pca.fit_transform(X)

    evr = pca.explained_variance_ratio_
    cumulative = np.cumsum(evr)

    print(f"\n  Explained variance ratios:")
    for i in range(N_PCS):
        print(f"    PC{i+1:2d}: {evr[i]:.4f}  (cumulative: {cumulative[i]:.4f})")
    print(f"  PC1+PC2 cumulative: {cumulative[1]:.4f} ({cumulative[1]*100:.1f}%)")
    if cumulative[1] < 0.20:
        print(f"  [WARN] PC1+PC2 cumulative < 20% — consider more components")

    # ── 7. Save outputs ─────────────────────────────────────────
    print("\n[7/7] Saving outputs ...")

    # ── Projection CSV ───────────────────────────────────────────
    df = pd.DataFrame(all_metadata, columns=["system", "frame", "species", "start_conf"])
    for i in range(N_PCS):
        df[f"pc{i+1}"] = proj[:, i]

    # Merge classification from torsion_classification.csv
    df = df.merge(class_df[["system", "frame", "classification"]],
                  on=["system", "frame"], how="left")
    df["classification"] = df["classification"].fillna("UNK")
    df["conformer"] = df["classification"]  # alias per plan

    # Reorder columns
    col_order = ["system", "frame"] + [f"pc{i+1}" for i in range(N_PCS)] + \
                ["species", "conformer", "start_conf", "classification"]
    df = df[col_order]

    proj_path = os.path.join(out_dir, "joint_pca_projection.csv")
    backup_if_exists(proj_path)
    df.to_csv(proj_path, index=False)
    print(f"  Saved {proj_path} ({len(df)} rows)")

    # Check for NaN
    nan_count = df.isnull().sum().sum()
    if nan_count > 0:
        print(f"  [WARN] {nan_count} NaN values in projection CSV")
    else:
        print(f"  [OK] No NaN values in projection CSV")

    # ── Loadings CSV ─────────────────────────────────────────────
    # pca.components_ shape: (10, 204) → reshape to (10, 68, 3)
    components_3d = pca.components_.reshape(N_PCS, len(SCAFFOLD_HEAVY_IDX), 3)
    # Per-atom loading magnitude per PC
    atom_loading_mag = np.sqrt((components_3d ** 2).sum(axis=2))  # (10, 68)

    loadings_df = pd.DataFrame({
        "atom_index": SCAFFOLD_HEAVY_IDX,
        "element": [ref_elements[i] for i in SCAFFOLD_HEAVY_IDX],
    })
    for pc in range(N_PCS):
        loadings_df[f"pc{pc+1}"] = atom_loading_mag[pc]

    loadings_path = os.path.join(out_dir, "joint_pca_loadings.csv")
    backup_if_exists(loadings_path)
    loadings_df.to_csv(loadings_path, index=False)
    print(f"  Saved {loadings_path} ({len(loadings_df)} rows)")

    # ── Sub-analysis: all-scaffold PCA ───────────────────────────
    pca_all = PCA(n_components=N_PCS, svd_solver="randomized", random_state=42)
    proj_all = pca_all.fit_transform(X_all)
    evr_all = pca_all.explained_variance_ratio_

    df_all = pd.DataFrame(all_metadata, columns=["system", "frame", "species", "start_conf"])
    for i in range(N_PCS):
        df_all[f"pc{i+1}"] = proj_all[:, i]
    df_all = df_all.merge(class_df[["system", "frame", "classification"]],
                          on=["system", "frame"], how="left")
    df_all["classification"] = df_all["classification"].fillna("UNK")
    df_all["conformer"] = df_all["classification"]

    scaffold_all_path = os.path.join(out_dir, "joint_pca_projection_scaffold_all.csv")
    backup_if_exists(scaffold_all_path)
    df_all.to_csv(scaffold_all_path, index=False)
    print(f"  Saved {scaffold_all_path} ({len(df_all)} rows)")
    print(f"    All-scaffold PC1+PC2 cumulative: "
          f"{evr_all[0]+evr_all[1]:.4f} ({(evr_all[0]+evr_all[1])*100:.1f}%)")

    # ── Sub-analysis: peripheral PCA ─────────────────────────────
    pca_periph = PCA(n_components=N_PCS, svd_solver="randomized", random_state=42)
    proj_periph = pca_periph.fit_transform(X_periph)
    evr_periph = pca_periph.explained_variance_ratio_

    df_periph = pd.DataFrame(all_metadata, columns=["system", "frame", "species", "start_conf"])
    for i in range(N_PCS):
        df_periph[f"pc{i+1}"] = proj_periph[:, i]
    df_periph = df_periph.merge(class_df[["system", "frame", "classification"]],
                                on=["system", "frame"], how="left")
    df_periph["classification"] = df_periph["classification"].fillna("UNK")
    df_periph["conformer"] = df_periph["classification"]

    peripheral_path = os.path.join(out_dir, "joint_pca_projection_peripheral.csv")
    backup_if_exists(peripheral_path)
    df_periph.to_csv(peripheral_path, index=False)
    print(f"  Saved {peripheral_path} ({len(df_periph)} rows)")
    print(f"    Peripheral PC1+PC2 cumulative: "
          f"{evr_periph[0]+evr_periph[1]:.4f} ({(evr_periph[0]+evr_periph[1])*100:.1f}%)")

    # ── Plots ────────────────────────────────────────────────────

    # Helper for deduplicated legend
    def dedup_legend(ax):
        handles, labels = ax.get_legend_handles_labels()
        by_label = dict(zip(labels, handles))
        ax.legend(by_label.values(), by_label.keys(), fontsize=9, framealpha=0.9)

    # ── Plot 1: PC1 vs PC2 colored by classification ────────────
    fig, ax = plt.subplots(figsize=(7, 5.5))
    for cls_name, color in COLORS_CLASS.items():
        mask = df["classification"] == cls_name
        if mask.sum() == 0:
            continue
        ax.scatter(df.loc[mask, "pc1"], df.loc[mask, "pc2"],
                   c=color, label=cls_name, alpha=0.3, s=10, edgecolors="none")
    ax.set_xlabel(f"PC1 ({evr[0]*100:.1f}% variance)")
    ax.set_ylabel(f"PC2 ({evr[1]*100:.1f}% variance)")
    ax.set_title("Joint PCA — All 16 Systems (colored by SAP/TSAP/UNK)")
    dedup_legend(ax)
    fig.tight_layout()
    plot1_path = os.path.join(out_dir, "plot_joint_pca_scatter.png")
    backup_if_exists(plot1_path)
    fig.savefig(plot1_path, dpi=300)
    plt.close(fig)
    print(f"  Saved {plot1_path}")

    # ── Plot 2: PC1 vs PC2 colored by species ────────────────────
    fig, ax = plt.subplots(figsize=(7, 5.5))
    for sp_name, color in COLORS_SPECIES.items():
        mask = df["species"] == sp_name
        if mask.sum() == 0:
            continue
        marker = MARKERS_SPECIES[sp_name]
        ax.scatter(df.loc[mask, "pc1"], df.loc[mask, "pc2"],
                   c=color, label=sp_name, marker=marker,
                   alpha=0.3, s=10, edgecolors="none")
    ax.set_xlabel(f"PC1 ({evr[0]*100:.1f}% variance)")
    ax.set_ylabel(f"PC2 ({evr[1]*100:.1f}% variance)")
    ax.set_title("Joint PCA — All 16 Systems (colored by me vs phe)")
    dedup_legend(ax)
    fig.tight_layout()
    plot2_path = os.path.join(out_dir, "plot_joint_pca_species.png")
    backup_if_exists(plot2_path)
    fig.savefig(plot2_path, dpi=300)
    plt.close(fig)
    print(f"  Saved {plot2_path}")

    # ── Plot 3: PC1 vs PC2 colored by starting conformer ────────
    fig, ax = plt.subplots(figsize=(7, 5.5))
    for conf_name, color in COLORS_CONF.items():
        mask = df["start_conf"] == conf_name
        if mask.sum() == 0:
            continue
        ax.scatter(df.loc[mask, "pc1"], df.loc[mask, "pc2"],
                   c=color, label=f"{conf_name}-start", alpha=0.3, s=10, edgecolors="none")
    ax.set_xlabel(f"PC1 ({evr[0]*100:.1f}% variance)")
    ax.set_ylabel(f"PC2 ({evr[1]*100:.1f}% variance)")
    ax.set_title("Joint PCA — All 16 Systems (colored by starting conformer)")
    dedup_legend(ax)
    fig.tight_layout()
    plot3_path = os.path.join(out_dir, "plot_joint_pca_conformer.png")
    backup_if_exists(plot3_path)
    fig.savefig(plot3_path, dpi=300)
    plt.close(fig)
    print(f"  Saved {plot3_path}")

    # ── Plot 4: Free Energy Surface ─────────────────────────────
    fig, ax = plt.subplots(figsize=(7, 5.5))
    H, xedges, yedges = np.histogram2d(df["pc1"], df["pc2"], bins=50)
    H += 1e-10  # avoid log(0)
    G = -KT_KJ_PER_MOL * np.log(H)
    G -= G.min()  # shift so minimum = 0

    # Create meshgrid for contourf
    xcenters = 0.5 * (xedges[:-1] + xedges[1:])
    ycenters = 0.5 * (yedges[:-1] + yedges[1:])
    Xg, Yg = np.meshgrid(xcenters, ycenters)

    cf = ax.contourf(Xg, Yg, G.T, levels=20, cmap="cubehelix")
    cbar = fig.colorbar(cf, ax=ax, label="G (kJ/mol)")
    ax.set_xlabel(f"PC1 ({evr[0]*100:.1f}% variance)")
    ax.set_ylabel(f"PC2 ({evr[1]*100:.1f}% variance)")
    ax.set_title("Joint FES — All 16 Systems (PC1 vs PC2)")
    fig.tight_layout()
    plot4_path = os.path.join(out_dir, "plot_joint_fes.png")
    backup_if_exists(plot4_path)
    fig.savefig(plot4_path, dpi=300)
    plt.close(fig)
    print(f"  Saved {plot4_path}")

    # ── Plot 5: 4x4 grid — FES per system (shared color scale) ──
    systems_sorted = sorted(df["system"].unique())
    n_systems = len(systems_sorted)
    ncols = 4
    nrows = (n_systems + ncols - 1) // ncols

    # Pre-compute global G range for shared color scale
    g_min, g_max = np.inf, -np.inf
    for sys_name in systems_sorted:
        sys_df = df[df["system"] == sys_name]
        Hs, xe, ye = np.histogram2d(sys_df["pc1"], sys_df["pc2"], bins=30)
        Hs = Hs + 1e-10
        Gs = -KT_KJ_PER_MOL * np.log(Hs)
        Gs -= Gs.min()
        g_min = min(g_min, Gs.min())
        g_max = max(g_max, np.percentile(Gs, 99))  # clip extreme single-bin spikes

    fig, axes = plt.subplots(nrows, ncols, figsize=(16, 16),
                             sharex=True, sharey=True)
    axes_flat = axes.flatten() if nrows > 1 else axes.flatten()
    pc1_min, pc1_max = df["pc1"].min(), df["pc1"].max()
    pc2_min, pc2_max = df["pc2"].min(), df["pc2"].max()
    pc1_pad = 0.05 * (pc1_max - pc1_min)
    pc2_pad = 0.05 * (pc2_max - pc2_min)

    for idx, sys_name in enumerate(systems_sorted):
        ax = axes_flat[idx]
        sys_df = df[df["system"] == sys_name]
        Hs, xe, ye = np.histogram2d(sys_df["pc1"], sys_df["pc2"], bins=30)
        Hs = Hs + 1e-10
        Gs = -KT_KJ_PER_MOL * np.log(Hs)
        Gs -= Gs.min()
        xc = 0.5 * (xe[:-1] + xe[1:])
        yc = 0.5 * (ye[:-1] + ye[1:])
        Xs, Ys = np.meshgrid(xc, yc)
        cf = ax.contourf(Xs, Ys, Gs.T, levels=15, cmap="cubehelix",
                         vmin=g_min, vmax=g_max)
        ax.set_title(sys_name, fontsize=9)
        ax.set_xlim(pc1_min - pc1_pad, pc1_max + pc1_pad)
        ax.set_ylim(pc2_min - pc2_pad, pc2_max + pc2_pad)
        ax.tick_params(labelsize=7)

    for idx in range(n_systems, len(axes_flat)):
        axes_flat[idx].axis("off")

    cbar_ax = fig.add_axes([0.92, 0.15, 0.015, 0.7])
    fig.colorbar(cf, cax=cbar_ax, label="G (kJ/mol)")
    fig.text(0.5, 0.02, f"PC1 ({evr[0]*100:.1f}% variance)",
             ha="center", fontsize=12)
    fig.text(0.02, 0.5, f"PC2 ({evr[1]*100:.1f}% variance)",
             va="center", rotation="vertical", fontsize=12)
    fig.suptitle("Joint FES — Individual System Energy Landscapes (shared axes & color scale)",
                 fontsize=14, y=0.98)
    fig.tight_layout(rect=[0.03, 0.03, 0.90, 0.96])
    plot5_path = os.path.join(out_dir, "plot_joint_pca_4x4_fes_grid.png")
    backup_if_exists(plot5_path)
    fig.savefig(plot5_path, dpi=300)
    plt.close(fig)
    print(f"  Saved {plot5_path}")

    # ── Summary ──────────────────────────────────────────────────
    print("\n" + "=" * 65)
    print("SUMMARY")
    print("=" * 65)
    print(f"  Total frames:         {total_frames}")
    print(f"  Feature matrix:       {X.shape}")
    print(f"  PCA components:       {N_PCS}")
    print(f"  PC1 explained var:    {evr[0]:.4f} ({evr[0]*100:.1f}%)")
    print(f"  PC2 explained var:    {evr[1]:.4f} ({evr[1]*100:.1f}%)")
    print(f"  PC3 explained var:    {evr[2]:.4f} ({evr[2]*100:.1f}%)")
    print(f"  PC1+PC2 cumulative:   {cumulative[1]:.4f} ({cumulative[1]*100:.1f}%)")
    print(f"  Top 5 loading atoms (PC1):")
    top5_pc1 = np.argsort(-atom_loading_mag[0])[:5]
    for rank, idx in enumerate(top5_pc1):
        print(f"    {rank+1}. Atom {SCAFFOLD_HEAVY_IDX[idx]} "
              f"({ref_elements[SCAFFOLD_HEAVY_IDX[idx]]}): "
              f"loading = {atom_loading_mag[0, idx]:.4f}")
    print(f"  Top 5 loading atoms (PC2):")
    top5_pc2 = np.argsort(-atom_loading_mag[1])[:5]
    for rank, idx in enumerate(top5_pc2):
        print(f"    {rank+1}. Atom {SCAFFOLD_HEAVY_IDX[idx]} "
              f"({ref_elements[SCAFFOLD_HEAVY_IDX[idx]]}): "
              f"loading = {atom_loading_mag[1, idx]:.4f}")
    print(f"\n  Outputs:")
    print(f"    {proj_path}")
    print(f"    {loadings_path}")
    print(f"    {scaffold_all_path}")
    print(f"    {peripheral_path}")
    print(f"    {plot1_path}")
    print(f"    {plot2_path}")
    print(f"    {plot3_path}")
    print(f"    {plot4_path}")
    print("=" * 65)


if __name__ == "__main__":
    main()
