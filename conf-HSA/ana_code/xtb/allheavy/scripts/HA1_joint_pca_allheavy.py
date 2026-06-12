#!/usr/bin/env python3
"""
HA1_joint_pca_allheavy.py — Joint PCA on all 16 xTB trajectories with
ALL-HEAVY-ATOM alignment (68 scaffold heavy atoms instead of 9 Eu+8 atoms).

Instead of Eu+8 alignment [0,1,2,3,4,5,6,7,54] (9 atoms), this script uses
all 68 scaffold heavy atoms for Kabsch alignment, discovered dynamically from
the reference frame element list.

Variants:
  --tag allheavy_scaffold (default):
    68-atom alignment, all 68 scaffold heavy atoms as features (204 feats)
  --tag allheavy_scaffold_nochrom --exclude-features 102,103,104,105,107,109,113,114,115,116,118,120:
    68-atom alignment, exclude Ring C+D distal chromophore (56 atoms, 168 feats)

Outputs (tag-substituted):
  analysis/joint_pca_projection_{tag}.csv
  analysis/joint_pca_loadings_{tag}.csv
  analysis/allheavy_alignment_indices_{tag}.json
  analysis/plot_joint_pca_{tag}_scatter.png
  analysis/plot_joint_pca_{tag}_species.png
  analysis/plot_joint_pca_{tag}_conformer.png
  analysis/plot_joint_pca_{tag}_fes.png
  analysis/plot_joint_pca_{tag}_4x4_fes_grid.png

Usage:
  python scripts/HA1_joint_pca_allheavy.py \
      --data-dir data \
      --out-dir analysis \
      --mappings ../analysis/atom_mappings.json \
      --scaffold ../analysis/scaffold_mapping.json \
      --indices ../analysis/indices.json \
      --classification ../analysis/torsion_classification.csv \
      --tag allheavy_scaffold

  # Nochrom variant:
  python scripts/HA1_joint_pca_allheavy.py \
      --tag allheavy_scaffold_nochrom \
      --exclude-features 102,103,104,105,107,109,113,114,115,116,118,120
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

KT_KJ_PER_MOL = 2.479
N_PCS = 10

CORE_EU8_IDX = sorted([0, 1, 2, 3, 4, 5, 6, 7, 54])

COLORS_CLASS = {"SAP": "#1f77b4", "TSAP": "#ff7f0e", "UNK": "#2ca02c"}
COLORS_SPECIES = {"me": "#008080", "phe": "#9467bd"}
COLORS_CONF = {"sap": "#1f77b4", "tsap": "#ff7f0e"}
MARKERS_SPECIES = {"me": "o", "phe": "s"}


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
    if os.path.exists(path):
        backup_path = path + ".bak"
        shutil.copy2(path, backup_path)
        print(f"  Backed up existing file to {backup_path}")


def load_xyz_trajectory(traj_path):
    u = mda.Universe(traj_path, format="XYZ")
    n_frames = len(u.trajectory)
    n_atoms = u.atoms.n_atoms
    coords = np.empty((n_frames, n_atoms, 3), dtype=np.float64)
    for i, ts in enumerate(u.trajectory):
        coords[i] = ts.positions
    elements = [a.name for a in u.atoms]
    return coords, elements


def apply_tsap_mapping(coords, elements, mapping):
    mapped_coords = np.empty_like(coords)
    mapped_coords[:, mapping, :] = coords
    mapped_elements_arr = [None] * len(elements)
    for j, sap_idx in enumerate(mapping):
        mapped_elements_arr[sap_idx] = elements[j]
    return mapped_coords, mapped_elements_arr


def validate_scaffold_elements(mapped_elements, ref_elements, label):
    scaffold_elems = mapped_elements[:len(ref_elements)]
    if scaffold_elems == ref_elements:
        print(f"  [PASS] {label}: scaffold elements match reference SAP")
    else:
        mismatches = [(i, scaffold_elems[i], ref_elements[i])
                      for i in range(len(ref_elements)) if scaffold_elems[i] != ref_elements[i]]
        print(f"  [FAIL] {label}: {len(mismatches)} scaffold element mismatches:")
        for idx, got, expected in mismatches[:10]:
            print(f"    Index {idx}: got {got}, expected {expected}")
        raise ValueError(f"{label}: scaffold element mismatch after mapping")


def parse_comma_ints(s):
    if s is None or s.strip() == "":
        return None
    return sorted([int(x.strip()) for x in s.split(",") if x.strip()])


def compute_alignment_rmsd(aligned_coords, ref_align_pos, align_indices):
    n_frames = aligned_coords.shape[0]
    per_frame = np.empty(n_frames)
    for f in range(n_frames):
        diff = aligned_coords[f, align_indices] - ref_align_pos
        per_frame[f] = np.sqrt(np.mean(np.sum(diff ** 2, axis=1)))
    return per_frame.mean(), per_frame


def compute_alignment_rmsd_subset(aligned_coords, ref_frame, atom_indices):
    n_frames = aligned_coords.shape[0]
    ref_pos = ref_frame[np.array(atom_indices)]
    rmsds = np.empty(n_frames)
    for f in range(n_frames):
        diff = aligned_coords[f, atom_indices] - ref_pos
        rmsds[f] = np.sqrt(np.mean(np.sum(diff ** 2, axis=1)))
    return rmsds.mean(), rmsds


def run_pca(args):
    data_dir = args.data_dir
    out_dir = args.out_dir
    tag = args.tag

    os.makedirs(out_dir, exist_ok=True)

    # ── 1. Load mapping files ───────────────────────────────────
    print("=" * 65)
    print(f"HA1 Joint PCA (all-heavy-atom alignment) — tag={tag}")
    print("=" * 65)

    print("\n[1/8] Loading mapping and classification files ...")

    if not os.path.exists(args.mappings):
        sys.exit(f"ERROR: mappings file not found: {args.mappings}")
    with open(args.mappings) as f:
        atom_mappings = json.load(f)

    if not os.path.exists(args.scaffold):
        sys.exit(f"ERROR: scaffold mapping file not found: {args.scaffold}")
    with open(args.scaffold) as f:
        scaffold_mapping = json.load(f)

    if not os.path.exists(args.indices):
        sys.exit(f"ERROR: indices file not found: {args.indices}")
    with open(args.indices) as f:
        indices_data = json.load(f)

    if not os.path.exists(args.classification):
        sys.exit(f"ERROR: classification file not found: {args.classification}")
    class_df = pd.read_csv(args.classification)

    n_scaffold = scaffold_mapping["n_scaffold"]
    print(f"  Scaffold atoms (SAP): 0-{n_scaffold - 1} ({n_scaffold} atoms)")

    exclude_set = set()
    exclude_features = parse_comma_ints(args.exclude_features)
    if exclude_features is not None:
        exclude_set = set(exclude_features)
        print(f"  Excluded feature atoms: {sorted(exclude_set)} ({len(exclude_set)} atoms)")

    # ── 2. Identify systems ─────────────────────────────────────
    print("\n[2/8] Identifying systems ...")

    systems = sorted([d for d in os.listdir(data_dir)
                      if os.path.isdir(os.path.join(data_dir, d))
                      and not d.startswith(".")
                      and os.path.exists(os.path.join(data_dir, d, "traj.xyz"))])
    print(f"  Found {len(systems)} systems: {systems}")

    if len(systems) != 16:
        print(f"  [WARN] Expected 16 systems, found {len(systems)}")

    # ── 3. Load, map, validate, and collect all trajectories ────
    print("\n[3/8] Loading trajectories and applying TSAP->SAP mapping ...")

    ref_elements = None
    all_systems = []

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

            if ref_elements is not None:
                validate_scaffold_elements(elements, ref_elements, sys_name)
        else:
            if ref_elements is None:
                ref_elements = elements[:n_scaffold]
                print(f"  Set reference scaffold elements ({len(ref_elements)} atoms)")
            else:
                validate_scaffold_elements(elements, ref_elements, sys_name)

        all_systems.append({
            "name": sys_name,
            "species": species,
            "start_conf": start_conf,
            "coords": coords,
            "elements": elements,
            "n_frames": n_frames,
            "n_atoms": n_atoms,
        })

    total_frames = sum(s["n_frames"] for s in all_systems)
    print(f"\n  Total frames loaded: {total_frames}")

    # ── 4. Discover scaffold heavy atoms for alignment ───────────
    print("\n[4/8] Discovering scaffold heavy atoms for alignment ...")

    SCAFFOLD_HEAVY_IDX = [i for i in range(n_scaffold) if ref_elements[i] != "H"]
    print(f"  Scaffold heavy atoms: {len(SCAFFOLD_HEAVY_IDX)} "
          f"(expecting 68; 68 x 3 = {len(SCAFFOLD_HEAVY_IDX) * 3} features)")

    if len(SCAFFOLD_HEAVY_IDX) != 68:
        print(f"  [WARN] Expected 68 heavy scaffold atoms, got {len(SCAFFOLD_HEAVY_IDX)}")

    heavy_elems = [ref_elements[i] for i in SCAFFOLD_HEAVY_IDX]
    elem_counts = {}
    for e in heavy_elems:
        elem_counts[e] = elem_counts.get(e, 0) + 1
    print(f"  Heavy element breakdown: {elem_counts}")

    # Resolve alignment indices
    align_indices_cli = parse_comma_ints(args.align_indices)
    if align_indices_cli is None:
        align_indices = SCAFFOLD_HEAVY_IDX
        print(f"  Alignment indices (SAP): auto-discovered {len(align_indices)} scaffold heavy atoms")
    else:
        align_indices = align_indices_cli
        print(f"  Alignment indices (SAP): user-specified {len(align_indices)} atoms")
    print(f"  Alignment indices (SAP): {align_indices} ({len(align_indices)} atoms)")

    # Apply feature exclusions
    if exclude_set:
        FEATURE_IDX = [i for i in SCAFFOLD_HEAVY_IDX if i not in exclude_set]
        n_excluded = len(SCAFFOLD_HEAVY_IDX) - len(FEATURE_IDX)
        print(f"  Excluded {n_excluded} heavy atoms from features: {sorted(exclude_set)}")
        print(f"  Retained feature atoms: {len(FEATURE_IDX)} "
              f"({len(FEATURE_IDX)} x 3 = {len(FEATURE_IDX) * 3} features)")
    else:
        FEATURE_IDX = SCAFFOLD_HEAVY_IDX
        print(f"  No exclusions -- all {len(FEATURE_IDX)} heavy atoms in feature set "
              f"({len(FEATURE_IDX) * 3} features)")

    # ── 5. Align all frames to reference ────────────────────────
    print(f"\n[5/8] Kabsch-aligning all frames to reference "
          f"(me_rrrD_sap frame 0) on {len(align_indices)} atoms ...")

    ref_frame = all_systems[0]["coords"][0]
    ref_align_pos = ref_frame[np.array(align_indices)]

    align_rmsd_summary = {}
    core9_rmsd_summary = {}

    for sys_data in all_systems:
        n_frames = sys_data["n_frames"]
        aligned = np.empty_like(sys_data["coords"])
        for f in range(n_frames):
            mobile = sys_data["coords"][f]
            mobile_align = mobile[np.array(align_indices)]
            R, P_cent, Q_cent = kabsch(mobile_align, ref_align_pos)
            aligned[f] = (mobile - P_cent) @ R + Q_cent
        sys_data["aligned"] = aligned

        mean_rmsd, _ = compute_alignment_rmsd(aligned, ref_align_pos, align_indices)
        align_rmsd_summary[sys_data["name"]] = mean_rmsd

        # Core-9 RMSD diagnostic (Eu+8 alignment quality after all-heavy alignment)
        mean_core9, _ = compute_alignment_rmsd_subset(aligned, ref_frame, CORE_EU8_IDX)
        core9_rmsd_summary[sys_data["name"]] = mean_core9

        del sys_data["coords"]
        print(f"  Aligned {sys_data['name']} ({n_frames} frames, "
              f"align RMSD={mean_rmsd:.4f} A, core-9 RMSD={mean_core9:.4f} A)")

    print(f"\n  Alignment RMSD summary (tag={tag}):")
    for sys_name in sorted(align_rmsd_summary.keys()):
        rmsd = align_rmsd_summary[sys_name]
        core9 = core9_rmsd_summary[sys_name]
        warn_flag = " *** WARN > 0.5 A" if core9 > 0.5 else ""
        print(f"    {sys_name:25s}  align={rmsd:.4f} A  core-9={core9:.4f} A{warn_flag}")

    n_warn = sum(1 for v in core9_rmsd_summary.values() if v > 0.5)
    if n_warn > 0:
        print(f"  [WARN] {n_warn}/16 systems have core-9 RMSD > 0.5 A after all-heavy alignment")
    else:
        print(f"  [OK] All 16 systems have core-9 RMSD < 0.5 A after all-heavy alignment")

    # ── 6. Save alignment indices JSON ──────────────────────────
    print(f"\n[6/8] Saving alignment indices JSON ...")

    indices_json = {
        "tag": tag,
        "align_indices": align_indices,
        "feature_indices": FEATURE_IDX,
        "n_align": len(align_indices),
        "n_features": len(FEATURE_IDX),
        "n_feature_dims": len(FEATURE_IDX) * 3,
        "excluded_features": sorted(exclude_set) if exclude_set else [],
        "n_excluded": len(exclude_set),
        "core_eu8_indices": CORE_EU8_IDX,
        "scaffold_heavy_indices": SCAFFOLD_HEAVY_IDX,
        "n_scaffold": n_scaffold,
        "n_scaffold_heavy": len(SCAFFOLD_HEAVY_IDX),
        "reference_system": all_systems[0]["name"],
        "reference_frame": 0,
    }
    indices_json_path = os.path.join(out_dir, f"allheavy_alignment_indices_{tag}.json")
    backup_if_exists(indices_json_path)
    with open(indices_json_path, "w") as f:
        json.dump(indices_json, f, indent=2)
    print(f"  Saved {indices_json_path}")

    # ── 7. Extract features, concatenate, fit PCA ──────────────
    print(f"\n[7/8] Extracting features and fitting joint PCA "
          f"({len(FEATURE_IDX)} atoms x 3 = {len(FEATURE_IDX) * 3} features) ...")

    feature_idx_arr = np.array(FEATURE_IDX)
    all_features = []
    all_metadata = []

    for sys_data in all_systems:
        feat_coords = sys_data["aligned"][:, feature_idx_arr, :]
        features = feat_coords.reshape(sys_data["n_frames"], -1)
        all_features.append(features)

        for f in range(sys_data["n_frames"]):
            all_metadata.append((sys_data["name"], f,
                                 sys_data["species"], sys_data["start_conf"]))

    X = np.vstack(all_features)
    print(f"  Feature matrix: {X.shape}")

    pca = PCA(n_components=N_PCS, svd_solver="randomized", random_state=42)
    proj = pca.fit_transform(X)

    evr = pca.explained_variance_ratio_
    cumulative = np.cumsum(evr)

    print(f"\n  Explained variance ratios:")
    for i in range(N_PCS):
        print(f"    PC{i+1:2d}: {evr[i]:.4f}  (cumulative: {cumulative[i]:.4f})")
    print(f"  PC1+PC2 cumulative: {cumulative[1]:.4f} ({cumulative[1]*100:.1f}%)")
    if cumulative[1] < 0.20:
        print(f"  [WARN] PC1+PC2 cumulative < 20% -- consider more components")

    # ── 8. Save outputs ─────────────────────────────────────────
    print(f"\n[8/8] Saving outputs (tag={tag}) ...")

    df = pd.DataFrame(all_metadata, columns=["system", "frame", "species", "start_conf"])
    for i in range(N_PCS):
        df[f"pc{i+1}"] = proj[:, i]

    df = df.merge(class_df[["system", "frame", "classification"]],
                  on=["system", "frame"], how="left")
    df["classification"] = df["classification"].fillna("UNK")
    df["conformer"] = df["classification"]

    col_order = ["system", "frame"] + [f"pc{i+1}" for i in range(N_PCS)] + \
                ["species", "conformer", "start_conf", "classification"]
    df = df[col_order]

    proj_path = os.path.join(out_dir, f"joint_pca_projection_{tag}.csv")
    backup_if_exists(proj_path)
    df.to_csv(proj_path, index=False)
    print(f"  Saved {proj_path} ({len(df)} rows)")

    nan_count = df.isnull().sum().sum()
    if nan_count > 0:
        print(f"  [WARN] {nan_count} NaN values in projection CSV")
    else:
        print(f"  [OK] No NaN values in projection CSV")

    # ── Loadings CSV ────────────────────────────────────────────
    n_feat_atoms = len(FEATURE_IDX)
    components_3d = pca.components_.reshape(N_PCS, n_feat_atoms, 3)
    atom_loading_mag = np.sqrt((components_3d ** 2).sum(axis=2))

    loadings_df = pd.DataFrame({
        "atom_index": FEATURE_IDX,
        "element": [ref_elements[i] for i in FEATURE_IDX],
    })
    for pc in range(N_PCS):
        loadings_df[f"pc{pc+1}"] = atom_loading_mag[pc]

    loadings_path = os.path.join(out_dir, f"joint_pca_loadings_{tag}.csv")
    backup_if_exists(loadings_path)
    loadings_df.to_csv(loadings_path, index=False)
    print(f"  Saved {loadings_path} ({len(loadings_df)} rows)")

    # ── Plots ───────────────────────────────────────────────────

    def dedup_legend(ax):
        handles, labels = ax.get_legend_handles_labels()
        by_label = dict(zip(labels, handles))
        ax.legend(by_label.values(), by_label.keys(), fontsize=9, framealpha=0.9)

    # ── Plot 1: PC1 vs PC2 colored by classification ───────────
    fig, ax = plt.subplots(figsize=(7, 5.5))
    for cls_name, color in COLORS_CLASS.items():
        mask = df["classification"] == cls_name
        if mask.sum() == 0:
            continue
        ax.scatter(df.loc[mask, "pc1"], df.loc[mask, "pc2"],
                   c=color, label=cls_name, alpha=0.3, s=10, edgecolors="none")
    ax.set_xlabel(f"PC1 ({evr[0]*100:.1f}% variance)")
    ax.set_ylabel(f"PC2 ({evr[1]*100:.1f}% variance)")
    ax.set_title(f"Joint PCA ({tag}) — colored by SAP/TSAP/UNK")
    dedup_legend(ax)
    fig.tight_layout()
    plot1_path = os.path.join(out_dir, f"plot_joint_pca_{tag}_scatter.png")
    backup_if_exists(plot1_path)
    fig.savefig(plot1_path, dpi=300)
    plt.close(fig)
    print(f"  Saved {plot1_path}")

    # ── Plot 2: PC1 vs PC2 colored by species ───────────────────
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
    ax.set_title(f"Joint PCA ({tag}) — colored by me vs phe")
    dedup_legend(ax)
    fig.tight_layout()
    plot2_path = os.path.join(out_dir, f"plot_joint_pca_{tag}_species.png")
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
                    c=color, label=f"{conf_name}-start", alpha=0.3, s=10,
                    edgecolors="none")
    ax.set_xlabel(f"PC1 ({evr[0]*100:.1f}% variance)")
    ax.set_ylabel(f"PC2 ({evr[1]*100:.1f}% variance)")
    ax.set_title(f"Joint PCA ({tag}) — colored by starting conformer")
    dedup_legend(ax)
    fig.tight_layout()
    plot3_path = os.path.join(out_dir, f"plot_joint_pca_{tag}_conformer.png")
    backup_if_exists(plot3_path)
    fig.savefig(plot3_path, dpi=300)
    plt.close(fig)
    print(f"  Saved {plot3_path}")

    # ── Plot 4: Free Energy Surface (cubehelix colormap) ────────
    fig, ax = plt.subplots(figsize=(7, 5.5))
    H, xedges, yedges = np.histogram2d(df["pc1"], df["pc2"], bins=50)
    H += 1e-10
    G = -KT_KJ_PER_MOL * np.log(H)
    G -= G.min()
    G = np.clip(G, 0, 15)

    xcenters = 0.5 * (xedges[:-1] + xedges[1:])
    ycenters = 0.5 * (yedges[:-1] + yedges[1:])
    Xg, Yg = np.meshgrid(xcenters, ycenters)

    levels = np.linspace(0, 15, 16)
    cf = ax.contourf(Xg, Yg, G.T, levels=levels, cmap="cubehelix", vmin=0, vmax=15)
    cbar = fig.colorbar(cf, ax=ax, label="G (kJ/mol)")
    ax.set_xlabel(f"PC1 ({evr[0]*100:.1f}% variance)")
    ax.set_ylabel(f"PC2 ({evr[1]*100:.1f}% variance)")
    ax.set_title(f"Joint FES ({tag}) — PC1 vs PC2")
    fig.tight_layout()
    plot4_path = os.path.join(out_dir, f"plot_joint_pca_{tag}_fes.png")
    backup_if_exists(plot4_path)
    fig.savefig(plot4_path, dpi=300)
    plt.close(fig)
    print(f"  Saved {plot4_path}")

    # ── Plot 5: 4x4 grid — FES per system (shared color scale) ──
    systems_sorted = sorted(df["system"].unique())
    n_systems = len(systems_sorted)
    ncols = 4
    nrows = (n_systems + ncols - 1) // ncols

    G_MAX_CLIP = 12
    levels_shared = np.linspace(0, G_MAX_CLIP, 13)

    fig, axes = plt.subplots(nrows, ncols, figsize=(16, 16),
                              sharex=True, sharey=True)
    axes_flat = axes.flatten() if nrows > 1 else axes.flatten()
    pc1_min, pc1_max = df["pc1"].min(), df["pc1"].max()
    pc2_min, pc2_max = df["pc2"].min(), df["pc2"].max()
    pc1_pad = 0.05 * (pc1_max - pc1_min)
    pc2_pad = 0.05 * (pc2_max - pc2_min)

    cf_last = None
    for idx, sys_name in enumerate(systems_sorted):
        ax = axes_flat[idx]
        sys_df = df[df["system"] == sys_name]
        Hs, xe, ye = np.histogram2d(sys_df["pc1"], sys_df["pc2"], bins=30)
        Hs = Hs + 1e-10
        Gs = -KT_KJ_PER_MOL * np.log(Hs)
        Gs -= Gs.min()
        Gs = np.clip(Gs, 0, G_MAX_CLIP)
        xc = 0.5 * (xe[:-1] + xe[1:])
        yc = 0.5 * (ye[:-1] + ye[1:])
        Xs, Ys = np.meshgrid(xc, yc)
        cf_last = ax.contourf(Xs, Ys, Gs.T, levels=levels_shared, cmap="cubehelix",
                         vmin=0, vmax=G_MAX_CLIP)
        ax.set_title(sys_name, fontsize=9)
        ax.set_xlim(pc1_min - pc1_pad, pc1_max + pc1_pad)
        ax.set_ylim(pc2_min - pc2_pad, pc2_max + pc2_pad)
        ax.tick_params(labelsize=7)

    for idx in range(n_systems, len(axes_flat)):
        axes_flat[idx].axis("off")

    cbar_ax = fig.add_axes([0.92, 0.15, 0.015, 0.7])
    if cf_last is not None:
        fig.colorbar(cf_last, cax=cbar_ax, label="G (kJ/mol)")
    fig.text(0.5, 0.02, f"PC1 ({evr[0]*100:.1f}% variance)",
             ha="center", fontsize=12)
    fig.text(0.02, 0.5, f"PC2 ({evr[1]*100:.1f}% variance)",
             va="center", rotation="vertical", fontsize=12)
    fig.suptitle(f"Joint FES ({tag}) — Individual System Energy Landscapes "
                 f"(shared axes & color scale)",
                 fontsize=14, y=0.98)
    fig.tight_layout(rect=[0.03, 0.03, 0.90, 0.96])
    plot5_path = os.path.join(out_dir, f"plot_joint_pca_{tag}_4x4_fes_grid.png")
    backup_if_exists(plot5_path)
    fig.savefig(plot5_path, dpi=300)
    plt.close(fig)
    print(f"  Saved {plot5_path}")

    # ── Summary ─────────────────────────────────────────────────
    print("\n" + "=" * 65)
    print(f"SUMMARY (tag={tag})")
    print("=" * 65)
    print(f"  Total frames:         {total_frames}")
    print(f"  Feature matrix:       {X.shape}")
    print(f"  Alignment indices:    {len(align_indices)} atoms (all scaffold heavy)")
    print(f"  Feature atoms:        {len(FEATURE_IDX)} "
          f"({len(FEATURE_IDX) * 3} features)")
    if exclude_set:
        print(f"  Excluded atoms:       {sorted(exclude_set)} ({len(exclude_set)} atoms)")
    else:
        print(f"  Excluded atoms:       none")
    print(f"  PCA components:       {N_PCS}")
    print(f"  PC1 explained var:    {evr[0]:.4f} ({evr[0]*100:.1f}%)")
    print(f"  PC2 explained var:    {evr[1]:.4f} ({evr[1]*100:.1f}%)")
    print(f"  PC3 explained var:    {evr[2]:.4f} ({evr[2]*100:.1f}%)")
    print(f"  PC1+PC2 cumulative:  {cumulative[1]:.4f} ({cumulative[1]*100:.1f}%)")
    print(f"  Core-9 RMSD range:   {min(core9_rmsd_summary.values()):.4f} – "
          f"{max(core9_rmsd_summary.values()):.4f} A")
    print(f"  Top 5 loading atoms (PC1):")
    top5_pc1 = np.argsort(-atom_loading_mag[0])[:5]
    for rank, idx in enumerate(top5_pc1):
        print(f"    {rank+1}. Atom {FEATURE_IDX[idx]} "
              f"({ref_elements[FEATURE_IDX[idx]]}): "
              f"loading = {atom_loading_mag[0, idx]:.4f}")
    print(f"  Top 5 loading atoms (PC2):")
    top5_pc2 = np.argsort(-atom_loading_mag[1])[:5]
    for rank, idx in enumerate(top5_pc2):
        print(f"    {rank+1}. Atom {FEATURE_IDX[idx]} "
              f"({ref_elements[FEATURE_IDX[idx]]}): "
              f"loading = {atom_loading_mag[1, idx]:.4f}")
    print(f"\n  Outputs:")
    print(f"    {proj_path}")
    print(f"    {loadings_path}")
    print(f"    {indices_json_path}")
    print(f"    {plot1_path}")
    print(f"    {plot2_path}")
    print(f"    {plot3_path}")
    print(f"    {plot4_path}")
    print(f"    {plot5_path}")
    print("=" * 65)

    return {
        "pca": pca,
        "projection_df": df,
        "loadings_df": loadings_df,
        "feature_idx": FEATURE_IDX,
        "ref_elements": ref_elements,
        "align_rmsd_summary": align_rmsd_summary,
        "core9_rmsd_summary": core9_rmsd_summary,
    }


def build_parser():
    parser = argparse.ArgumentParser(
        description="Joint PCA on all 16 xTB trajectories with all-heavy-atom "
                    "alignment (68 scaffold heavy atoms instead of Eu+8).")
    parser.add_argument("--data-dir", default="data",
                        help="Directory containing system subdirectories")
    parser.add_argument("--out-dir", default="analysis",
                        help="Output directory for CSVs and plots")
    parser.add_argument("--mappings", default="../analysis/atom_mappings.json",
                        help="Path to atom_mappings.json")
    parser.add_argument("--scaffold", default="../analysis/scaffold_mapping.json",
                        help="Path to scaffold_mapping.json")
    parser.add_argument("--indices", default="../analysis/indices.json",
                        help="Path to indices.json")
    parser.add_argument("--classification", default="../analysis/torsion_classification.csv",
                        help="Path to torsion_classification.csv")
    parser.add_argument("--align-indices", type=str, default=None,
                        help="Comma-separated SAP atom indices for Kabsch alignment "
                             "(default: auto-discover 68 scaffold heavy atoms)")
    parser.add_argument("--exclude-features", type=str, default=None,
                        help="Comma-separated SAP atom indices to EXCLUDE from "
                             "feature set (e.g., Ring C+D: 102,103,104,105,107,109,"
                             "113,114,115,116,118,120)")
    parser.add_argument("--tag", type=str, default="allheavy_scaffold",
                        help="Output file tag")
    parser.add_argument("--indices-json", type=str, default=None,
                        help="Path to save alignment indices JSON "
                             "(default: analysis/allheavy_alignment_indices_{tag}.json)")
    return parser


if __name__ == "__main__":
    parser = build_parser()
    args = parser.parse_args()
    if args.indices_json is None:
        args.indices_json = os.path.join(args.out_dir,
                                         f"allheavy_alignment_indices_{args.tag}.json")
    run_pca(args)
