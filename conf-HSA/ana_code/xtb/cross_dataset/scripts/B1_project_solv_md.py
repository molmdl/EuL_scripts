#!/usr/bin/env python3
"""
B1_project_solv_md.py — Project solv_md trajectories into existing xTB eu8_nochrom PCA space.

Pipeline:
1. Re-fit PCA on xTB data (all 16 systems, deterministic random_state=42)
2. Save PCA mean and components as .npy files
3. For each solv_md system: extract MOL, apply mappings, align, project
4. Save combined CSV with xtb + solv projections

Usage:
  python cross_dataset/scripts/B1_project_solv_md.py \
      --data-dir data \
      --solv-dir solv_md \
      --out-dir cross_dataset/analysis \
      --mappings analysis/atom_mappings.json \
      --scaffold analysis/scaffold_mapping.json

  # Validate against existing projections:
  python cross_dataset/scripts/B1_project_solv_md.py --validate analysis/joint_pca_projection_eu8_nochrom.csv
"""

import argparse
import json
import os
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
import MDAnalysis as mda
from sklearn.decomposition import PCA

# ───────────────────────────────────────────────────────────────
# Constants
# ───────────────────────────────────────────────────────────────
N_PCS = 10
ALIGN_EU8 = sorted([0, 1, 2, 3, 4, 5, 6, 7, 54])  # Eu + 8 donor atoms

# Ring C+D distal chromophore atoms to exclude (eu8_nochrom variant)
EXCLUDE_FEATURES = {102, 103, 104, 105, 107, 109, 113, 114, 115, 116, 118, 120}

# Canonical 16 systems
SYSTEMS = [
    "me_rrrD_sap", "me_rrrD_tsap", "me_rrrL_sap", "me_rrrL_tsap",
    "me_sssD_sap", "me_sssD_tsap", "me_sssL_sap", "me_sssL_tsap",
    "phe_rrrD_sap", "phe_rrrD_tsap", "phe_rrrL_sap", "phe_rrrL_tsap",
    "phe_sssD_sap", "phe_sssD_tsap", "phe_sssL_sap", "phe_sssL_tsap",
]


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


def apply_tsap_mapping(coords, mapping):
    """
    Reorder TSAP coords to SAP ordering.

    mapping[i] = SAP index that TSAP atom i corresponds to.
    Apply: mapped[:, mapping[j], :] = raw[:, j, :] for each j
    This places TSAP atom j at SAP position mapping[j].

    Parameters
    ----------
    coords : ndarray, shape (n_frames, n_atoms, 3) or (n_atoms, 3)
    mapping : list of int

    Returns
    -------
    mapped_coords : ndarray, same shape as coords
    """
    mapped_coords = np.empty_like(coords)
    if coords.ndim == 3:
        # (n_frames, n_atoms, 3)
        mapped_coords[:, mapping, :] = coords
    else:
        # (n_atoms, 3)
        mapped_coords[mapping, :] = coords
    return mapped_coords


def parse_system_name(sys_name):
    """
    Parse system name into components.

    Examples: 'me_rrrD_sap' -> (me, rrr, D, sap)
              'phe_sssL_tsap' -> (phe, sss, L, tsap)
    """
    parts = sys_name.split("_")
    species = parts[0]  # me or phe
    isomer = parts[1][:3]  # rrr or sss
    handedness = parts[1][3]  # D or L
    conformer = parts[2]  # sap or tsap
    return species, isomer, handedness, conformer


def load_xyz_trajectory(traj_path):
    """
    Load XYZ trajectory via MDAnalysis.

    Returns
    -------
    coords : ndarray, shape (n_frames, n_atoms, 3)
    elements : list of str
    """
    u = mda.Universe(str(traj_path), format="XYZ")
    n_frames = len(u.trajectory)
    n_atoms = u.atoms.n_atoms
    coords = np.empty((n_frames, n_atoms, 3), dtype=np.float64)
    for i, ts in enumerate(u.trajectory):
        coords[i] = ts.positions
    elements = [a.name for a in u.atoms]
    return coords, elements


def get_reference_frame(data_dir):
    """
    Load reference frame from me_rrrD_sap xtbopt.xyz (frame 0).

    Returns
    -------
    ref_frame : ndarray (n_atoms, 3)
    ref_elements : list of str
    """
    ref_path = Path(data_dir) / "me_rrrD_sap" / "xtbopt.xyz"
    if not ref_path.exists():
        # Fallback: use traj.xyz frame 0
        ref_path = Path(data_dir) / "me_rrrD_sap" / "traj.xyz"
    u = mda.Universe(str(ref_path), format="XYZ")
    ref_frame = u.trajectory[0].positions.copy()
    ref_elements = [a.name for a in u.atoms]
    # NOTE: xtbopt.xyz is the GFN2-xTB optimized geometry, which is identical
    # to frame 0 of traj.xyz (the first saved frame of the MD). B1 uses xtbopt.xyz
    # as the reference; T3_joint_pca_eu8.py uses traj.xyz frame 0. Both sources
    # should produce the same reference coordinates. An optional assertion can be
    # added: np.allclose(ref_xtbopt, ref_traj0, atol=1e-6).
    return ref_frame, ref_elements


# ───────────────────────────────────────────────────────────────
# Step 1: Re-fit PCA on xTB data
# ───────────────────────────────────────────────────────────────

def refit_pca_xtb(data_dir, mappings_path, scaffold_path):
    """
    Re-fit PCA on all 16 xTB trajectories using eu8_nochrom variant.

    Returns
    -------
    pca : fitted PCA object
    feature_idx : list of int — feature atom indices (SAP ordering)
    ref_elements : list of str — reference element list (SAP)
    ref_frame : ndarray (n_atoms, 3) — reference frame
    xtb_proj : ndarray (32000, N_PCS) — xTB projections
    xtb_metadata : list of tuples — (system, species, isomer, handedness, conformer, frame)
    """
    # Load mappings
    with open(mappings_path) as f:
        atom_mappings = json.load(f)
    with open(scaffold_path) as f:
        scaffold_mapping = json.load(f)

    n_scaffold = scaffold_mapping["n_scaffold"]

    # Get reference
    ref_frame, ref_elements = get_reference_frame(data_dir)
    ref_align_pos = ref_frame[np.array(ALIGN_EU8)]

    # Identify scaffold heavy atoms
    SCAFFOLD_HEAVY_IDX = [i for i in range(n_scaffold) if ref_elements[i] != "H"]
    assert len(SCAFFOLD_HEAVY_IDX) == 68, f"Expected 68 scaffold heavy atoms, got {len(SCAFFOLD_HEAVY_IDX)}"

    # Apply feature exclusions (eu8_nochrom)
    FEATURE_IDX = [i for i in SCAFFOLD_HEAVY_IDX if i not in EXCLUDE_FEATURES]
    assert len(FEATURE_IDX) == 56, f"Expected 56 feature atoms, got {len(FEATURE_IDX)}"

    # Load and process all 16 xTB systems
    all_features = []
    all_metadata = []

    for sys_name in SYSTEMS:
        traj_path = Path(data_dir) / sys_name / "traj.xyz"
        if not traj_path.exists():
            print(f"  [SKIP] {sys_name}: traj.xyz not found")
            continue

        species, isomer, handedness, conformer = parse_system_name(sys_name)
        print(f"  Loading {sys_name} ({species}, {conformer}) ...", end=" ", flush=True)

        coords, elements = load_xyz_trajectory(traj_path)
        n_frames = coords.shape[0]
        print(f"{n_frames} frames")

        # TSAP->SAP reordering
        if conformer == "tsap":
            mapping_key = f"{species}_tsap_to_sap"
            mapping = atom_mappings[mapping_key]
            coords = apply_tsap_mapping(coords, mapping)

        # Align each frame
        aligned = np.empty_like(coords)
        for f in range(n_frames):
            mobile = coords[f]
            mobile_align = mobile[np.array(ALIGN_EU8)]
            R, P_cent, Q_cent = kabsch(mobile_align, ref_align_pos)
            aligned[f] = (mobile - P_cent) @ R + Q_cent

        # Extract feature atoms
        feature_idx_arr = np.array(FEATURE_IDX)
        feat_coords = aligned[:, feature_idx_arr, :]  # (n_frames, 56, 3)
        features = feat_coords.reshape(n_frames, -1)   # (n_frames, 168)
        all_features.append(features)

        # Metadata
        for f in range(n_frames):
            all_metadata.append((sys_name, species, isomer, handedness, conformer, f))

        # Free memory
        del coords, aligned

    # Stack features
    X = np.vstack(all_features)
    print(f"  xTB feature matrix: {X.shape}")

    # Fit PCA
    pca = PCA(n_components=N_PCS, svd_solver="randomized", random_state=42)
    xtb_proj = pca.fit_transform(X)

    evr = pca.explained_variance_ratio_
    cumulative = np.cumsum(evr)
    print(f"  PC1+PC2 cumulative: {cumulative[1]:.4f} ({cumulative[1]*100:.1f}%)")

    return pca, FEATURE_IDX, ref_elements, ref_frame, xtb_proj, all_metadata


# ───────────────────────────────────────────────────────────────
# Step 2: Save PCA state
# ───────────────────────────────────────────────────────────────

def save_pca_state(pca, out_dir):
    """Save PCA mean_ and components_ as .npy files."""
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    mean_path = out_dir / "pca_mean_eu8_nochrom.npy"
    comp_path = out_dir / "pca_components_eu8_nochrom.npy"

    np.save(str(mean_path), pca.mean_)
    np.save(str(comp_path), pca.components_)

    print(f"  Saved PCA mean: {mean_path} (shape={pca.mean_.shape})")
    print(f"  Saved PCA components: {comp_path} (shape={pca.components_.shape})")

    return mean_path, comp_path


# ───────────────────────────────────────────────────────────────
# Step 3: Project solv_md trajectories
# ───────────────────────────────────────────────────────────────

def project_solv_md(solv_dir, pca, feature_idx, ref_frame, mappings_path, scaffold_path):
    """
    Project solv_md trajectories into PCA space.

    For each system: load TPR+XTC, extract MOL, apply mappings, align, project.

    Returns
    -------
    solv_proj : ndarray (total_frames, N_PCS)
    solv_metadata : list of tuples — (system, species, isomer, handedness, conformer, frame)
    frame_counts : dict — {system: n_frames}
    """
    # Load mappings
    with open(mappings_path) as f:
        atom_mappings = json.load(f)
    with open(scaffold_path) as f:
        scaffold_mapping = json.load(f)

    n_scaffold = scaffold_mapping["n_scaffold"]
    feature_idx_arr = np.array(feature_idx)
    ref_align_pos = ref_frame[np.array(ALIGN_EU8)]

    # Build phe_to_me reordering from scaffold_mapping
    # phe_to_me_sap: maps phe SAP index -> me SAP index
    phe_to_me_sap = scaffold_mapping["phe_to_me_sap"]


    all_proj = []
    all_metadata = []
    frame_counts = {}

    for sys_name in SYSTEMS:
        tpr_path = Path(solv_dir) / sys_name / "prod_0.tpr"
        xtc_path = Path(solv_dir) / sys_name / "solv_all.xtc"

        if not tpr_path.exists() or not xtc_path.exists():
            print(f"  [SKIP] {sys_name}: missing TPR or XTC")
            continue

        # Check XTC is non-empty
        if xtc_path.stat().st_size == 0:
            print(f"  [SKIP] {sys_name}: empty XTC")
            continue

        species, isomer, handedness, conformer = parse_system_name(sys_name)
        print(f"  Loading {sys_name} ({species}, {conformer}) ...", end=" ", flush=True)

        try:
            u = mda.Universe(str(tpr_path), str(xtc_path))
        except Exception as e:
            print(f"[ERROR] Cannot load: {e}")
            continue

        mol = u.select_atoms("resname MOL")
        if mol.n_atoms == 0:
            print(f"[SKIP] No MOL atoms found")
            continue

        n_frames = len(u.trajectory)
        print(f"{n_frames} frames, {mol.n_atoms} MOL atoms")

        # Process frame by frame to save memory
        frame_features = []
        for frame_idx, ts in enumerate(u.trajectory):
            mol_coords = mol.positions.copy()  # (n_mol_atoms, 3)

            # TSAP->SAP reordering
            if conformer == "tsap":
                mapping_key = f"{species}_tsap_to_sap"
                mapping = atom_mappings[mapping_key]
                mol_coords = apply_tsap_mapping(mol_coords, mapping)

            # phe->me scaffold reordering (after TSAP->SAP)
            # For phe systems: need to reorder scaffold atoms to me ordering
            # so that feature indices (defined in me SAP space) select the correct atoms
            if species == "phe":
                # phe_to_me_sap maps phe_SAP_index -> me_SAP_index
                # We need: reorder phe coords so that atom at me_SAP_index j
                # gets coords from phe atom that maps to j
                # i.e., for each me index j, find phe index i where phe_to_me_sap[i] == j
                # Build inverse: me_to_phe_sap
                # Actually, we just need the scaffold portion reordered
                # The inverse mapping: for each me_sap position j (0..126),
                # the phe_sap atom that maps to j
                me_to_phe_sap = {}
                for phe_idx_str, me_idx in phe_to_me_sap.items():
                    me_to_phe_sap[int(me_idx)] = int(phe_idx_str)

                # Build a reorder array for scaffold portion
                # reordered[j] = phe_coords[me_to_phe_sap[j]] for j in 0..n_scaffold-1
                reorder = np.array([me_to_phe_sap[j] for j in range(n_scaffold)])
                mol_coords[:n_scaffold] = mol_coords[reorder]

            # Kabsch align on Eu+8 alignment atoms
            mobile_align = mol_coords[np.array(ALIGN_EU8)]
            R, P_cent, Q_cent = kabsch(mobile_align, ref_align_pos)
            mol_coords = (mol_coords - P_cent) @ R + Q_cent

            # Extract feature atoms and flatten
            feat = mol_coords[feature_idx_arr].flatten()  # (168,)
            frame_features.append(feat)

        if len(frame_features) == 0:
            print(f"  [WARN] No frames processed for {sys_name}")
            continue

        X_solv = np.array(frame_features)  # (n_frames, 168)
        proj = pca.transform(X_solv)

        all_proj.append(proj)
        for f in range(n_frames):
            all_metadata.append((sys_name, species, isomer, handedness, conformer, f))
        frame_counts[sys_name] = n_frames

        print(f"    Projected {n_frames} frames")

        # Free memory
        del u, frame_features, X_solv, proj

    if len(all_proj) == 0:
        print("  [WARN] No solv_md projections collected")
        return np.empty((0, N_PCS)), [], {}, {}

    solv_proj = np.vstack(all_proj)
    print(f"  Total solv_md frames: {solv_proj.shape[0]}")

    return solv_proj, all_metadata, frame_counts


# ───────────────────────────────────────────────────────────────
# Step 4: Save combined CSV
# ───────────────────────────────────────────────────────────────

def save_combined_csv(xtb_proj, xtb_metadata, solv_proj, solv_metadata, out_path):
    """
    Save combined xtb + solv projection CSV.

    Columns: system, species, isomer, handedness, conformer, frame, method, PC1..PC10
    """
    # Build xtb dataframe
    xtb_rows = []
    for i, (sys_name, species, isomer, handedness, conformer, frame) in enumerate(xtb_metadata):
        row = {
            "system": sys_name,
            "species": species,
            "isomer": isomer,
            "handedness": handedness,
            "conformer": conformer,
            "frame": frame,
            "method": "xtb",
        }
        for pc in range(N_PCS):
            row[f"PC{pc+1}"] = xtb_proj[i, pc]
        xtb_rows.append(row)

    # Build solv dataframe
    solv_rows = []
    for i, (sys_name, species, isomer, handedness, conformer, frame) in enumerate(solv_metadata):
        row = {
            "system": sys_name,
            "species": species,
            "isomer": isomer,
            "handedness": handedness,
            "conformer": conformer,
            "frame": frame,
            "method": "solv",
        }
        for pc in range(N_PCS):
            row[f"PC{pc+1}"] = solv_proj[i, pc]
        solv_rows.append(row)

    df_xtb = pd.DataFrame(xtb_rows)
    df_solv = pd.DataFrame(solv_rows)
    df = pd.concat([df_xtb, df_solv], ignore_index=True)

    # Column order
    col_order = ["system", "species", "isomer", "handedness", "conformer", "frame", "method"]
    col_order += [f"PC{i+1}" for i in range(N_PCS)]
    df = df[col_order]

    out_path = Path(out_path)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(str(out_path), index=False)
    print(f"  Saved combined CSV: {out_path} ({len(df)} rows, {len(df.columns)} cols)")
    print(f"    xtb rows: {len(df_xtb)}, solv rows: {len(df_solv)}")

    # Check for NaN
    nan_count = df.isnull().sum().sum()
    if nan_count > 0:
        print(f"  [WARN] {nan_count} NaN values in combined CSV")
    else:
        print(f"  [OK] No NaN values in combined CSV")

    return df


# ───────────────────────────────────────────────────────────────
# Validation
# ───────────────────────────────────────────────────────────────

def validate_xtb_projection(xtb_proj, xtb_metadata, existing_csv_path):
    """
    Validate that re-fitted xTB projections match existing CSV.

    Checks Pearson correlation for PC1 > 0.999.
    """
    existing_df = pd.read_csv(existing_csv_path)

    # Build our xTB projection dataframe for comparison
    our_rows = []
    for i, (sys_name, species, isomer, handedness, conformer, frame) in enumerate(xtb_metadata):
        our_rows.append({
            "system": sys_name,
            "frame": frame,
            "PC1_new": xtb_proj[i, 0],
            "PC2_new": xtb_proj[i, 1],
        })
    our_df = pd.DataFrame(our_rows)

    # Merge on (system, frame)
    merged = existing_df.merge(our_df, on=["system", "frame"], how="inner")

    # Check correlation for PC1
    from scipy.stats import pearsonr
    r_pc1, p_pc1 = pearsonr(merged["pc1"], merged["PC1_new"])
    r_pc2, p_pc2 = pearsonr(merged["pc2"], merged["PC2_new"])

    print(f"  PC1 Pearson r = {r_pc1:.6f} (p = {p_pc1:.2e})")
    print(f"  PC2 Pearson r = {r_pc2:.6f} (p = {p_pc2:.2e})")

    if r_pc1 > 0.999:
        print(f"  [PASS] PC1 correlation > 0.999 — re-fit matches existing projections")
    else:
        print(f"  [WARN] PC1 correlation = {r_pc1:.6f} < 0.999 — possible numerical discrepancy")

    return r_pc1, r_pc2


def validate_solv_projection(solv_proj, solv_metadata, frame_counts):
    """
    Validate solv_md projections: check for NaN, extreme outliers, report per-system stats.
    """
    if solv_proj.shape[0] == 0:
        print("  [WARN] No solv_md projections to validate")
        return

    # Check NaN
    nan_count = np.isnan(solv_proj).sum()
    if nan_count > 0:
        print(f"  [WARN] {nan_count} NaN values in solv_md projections")
    else:
        print(f"  [OK] No NaN values in solv_md projections")

    # Per-PC range
    for pc in range(N_PCS):
        pc_vals = solv_proj[:, pc]
        print(f"  PC{pc+1}: min={pc_vals.min():.2f}, max={pc_vals.max():.2f}, "
              f"mean={pc_vals.mean():.2f}, std={pc_vals.std():.2f}")

    # Per-system frame counts
    n_systems_with_data = len(frame_counts)
    print(f"  Systems with solv_md data: {n_systems_with_data}/16")
    for sys_name, n_frames in sorted(frame_counts.items()):
        print(f"    {sys_name}: {n_frames} frames")


# ───────────────────────────────────────────────────────────────
# Main pipeline
# ───────────────────────────────────────────────────────────────

def run_pipeline(args):
    """Execute the full B1 pipeline."""
    t_start = time.time()

    data_dir = args.data_dir
    solv_dir = args.solv_dir
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    print("=" * 70)
    print("B1: Project solv_md into xTB eu8_nochrom PCA space")
    print("=" * 70)

    # ── Step 1: Re-fit PCA on xTB ────────────────────────────
    print("\n[Step 1] Re-fitting PCA on xTB data (eu8_nochrom, 56 atoms, 168 feats) ...")
    t1 = time.time()

    pca, feature_idx, ref_elements, ref_frame, xtb_proj, xtb_metadata = \
        refit_pca_xtb(data_dir, args.mappings, args.scaffold)

    t1_end = time.time()
    print(f"  Step 1 runtime: {t1_end - t1:.1f} s")

    # ── Step 2: Save PCA state ───────────────────────────────
    print("\n[Step 2] Saving PCA state ...")
    t2 = time.time()

    save_pca_state(pca, out_dir)

    t2_end = time.time()
    print(f"  Step 2 runtime: {t2_end - t2:.1f} s")

    # ── Step 3: Validate xTB re-fit ──────────────────────────
    print("\n[Step 3] Validating xTB re-fit against existing projections ...")
    t3 = time.time()

    existing_csv = Path(args.validate) if args.validate else \
        Path("analysis/joint_pca_projection_eu8_nochrom.csv")
    if existing_csv.exists():
        r_pc1, r_pc2 = validate_xtb_projection(xtb_proj, xtb_metadata, str(existing_csv))
    else:
        print(f"  [SKIP] Existing projection CSV not found: {existing_csv}")
        r_pc1, r_pc2 = None, None

    t3_end = time.time()
    print(f"  Step 3 runtime: {t3_end - t3:.1f} s")

    # ── Step 4: Project solv_md ──────────────────────────────
    print("\n[Step 4] Projecting solv_md trajectories into PCA space ...")
    t4 = time.time()

    solv_proj, solv_metadata, frame_counts = \
        project_solv_md(solv_dir, pca, feature_idx, ref_frame,
                        args.mappings, args.scaffold)

    t4_end = time.time()
    print(f"  Step 4 runtime: {t4_end - t4:.1f} s")

    # ── Step 5: Save combined CSV ────────────────────────────
    print("\n[Step 5] Saving combined CSV ...")
    t5 = time.time()

    csv_path = out_dir / "joint_projection_eu8_nochrom_xtb_solv.csv"
    df = save_combined_csv(xtb_proj, xtb_metadata, solv_proj, solv_metadata, csv_path)

    t5_end = time.time()
    print(f"  Step 5 runtime: {t5_end - t5:.1f} s")

    # ── Step 6: Validation ───────────────────────────────────
    print("\n[Step 6] Final validation ...")
    t6 = time.time()

    if existing_csv.exists() and r_pc1 is not None:
        print("\n  xTB re-fit validation:")
        # Already done in Step 3

    print("\n  solv_md projection validation:")
    validate_solv_projection(solv_proj, solv_metadata, frame_counts)

    t6_end = time.time()
    print(f"  Step 6 runtime: {t6_end - t6:.1f} s")

    # ── Summary ──────────────────────────────────────────────
    total_time = time.time() - t_start
    print("\n" + "=" * 70)
    print("B1 SUMMARY")
    print("=" * 70)
    print(f"  Total runtime: {total_time:.1f} s")
    print(f"  xTB frames: {xtb_proj.shape[0]}")
    print(f"  solv_md frames: {solv_proj.shape[0]}")
    print(f"  solv_md systems with data: {len(frame_counts)}/16")
    print(f"  Combined CSV: {csv_path} ({len(df)} rows)")
    if r_pc1 is not None:
        print(f"  xTB re-fit PC1 Pearson r: {r_pc1:.6f}")
    evr = pca.explained_variance_ratio_
    print(f"  PC1+PC2 variance: {evr[0]*100:.1f}% + {evr[1]*100:.1f}% = {(evr[0]+evr[1])*100:.1f}%")
    print(f"  Outputs:")
    print(f"    {out_dir / 'pca_mean_eu8_nochrom.npy'}")
    print(f"    {out_dir / 'pca_components_eu8_nochrom.npy'}")
    print(f"    {csv_path}")
    print("=" * 70)

    return {
        "pca": pca,
        "xtb_proj": xtb_proj,
        "solv_proj": solv_proj,
        "combined_df": df,
        "frame_counts": frame_counts,
        "r_pc1": r_pc1,
    }


def build_parser():
    """Build argument parser."""
    parser = argparse.ArgumentParser(
        description="B1: Project solv_md trajectories into existing xTB eu8_nochrom PCA space")
    parser.add_argument("--data-dir", default="data",
                        help="Directory with xTB system subdirectories")
    parser.add_argument("--solv-dir", default="solv_md",
                        help="Directory with solv_md system subdirectories")
    parser.add_argument("--out-dir", default="cross_dataset/analysis",
                        help="Output directory for .npy and CSV files")
    parser.add_argument("--mappings", default="analysis/atom_mappings.json",
                        help="Path to atom_mappings.json")
    parser.add_argument("--scaffold", default="analysis/scaffold_mapping.json",
                        help="Path to scaffold_mapping.json")
    parser.add_argument("--validate", default=None,
                        help="Path to existing eu8_nochrom projection CSV for validation "
                             "(default: analysis/joint_pca_projection_eu8_nochrom.csv)")
    return parser


if __name__ == "__main__":
    parser = build_parser()
    args = parser.parse_args()
    # Default validation path
    if args.validate is None:
        args.validate = "analysis/joint_pca_projection_eu8_nochrom.csv"
    run_pipeline(args)
