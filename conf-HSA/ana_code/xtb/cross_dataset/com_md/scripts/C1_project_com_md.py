#!/usr/bin/env python3
"""
C1_project_com_md.py — Project com_md (HSA-bound) trajectories into
existing xTB eu8_nochrom PCA space.

Pipeline:
1. Load existing PCA model (mean + components from cross_dataset/analysis/)
2. Load xTB reference frame for Kabsch alignment
3. For each com_md system: extract MOL, apply mappings, align, project
4. Save com_md-only projection CSV
5. Build joint 3-way CSV (xtb + solv + com_md)

Usage:
  python cross_dataset/com_md/scripts/C1_project_com_md.py \
      --data-dir data \
      --com-dir com_md \
      --out-dir cross_dataset/com_md/analysis \
      --mappings analysis/atom_mappings.json \
      --scaffold analysis/scaffold_mapping.json \
      --pca-dir cross_dataset/analysis \
      --existing-joint cross_dataset/analysis/joint_projection_eu8_nochrom_xtb_solv.csv
"""

import argparse
import json
import sys
import time
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
import MDAnalysis as mda
from sklearn.decomposition import PCA

# Suppress MDAnalysis warnings about EU3 element
warnings.filterwarnings("ignore", message="Unknown element")
warnings.filterwarnings("ignore", message="Failed to guess the mass")

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
# Helper functions (from B1, reused)
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
        Optimal rotation matrix (rotate P onto Q).
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
        mapped_coords[:, mapping, :] = coords
    else:
        mapped_coords[mapping, :] = coords
    return mapped_coords


def parse_system_name(sys_name):
    """
    Parse system name into components.

    Examples: 'me_rrrD_sap' -> (me, rrr, D, sap)
              'phe_sssL_tsap' -> (phe, sss, L, tsap)
    """
    parts = sys_name.split("_")
    species = parts[0]       # me or phe
    isomer = parts[1][:3]    # rrr or sss
    handedness = parts[1][3] # D or L
    conformer = parts[2]     # sap or tsap
    return species, isomer, handedness, conformer


def load_reference_frame(data_dir):
    """
    Load reference frame from me_rrrD_sap xtbopt.xyz (frame 0).

    Returns
    -------
    ref_frame : ndarray (n_atoms, 3)
    ref_elements : list of str
    """
    ref_path = Path(data_dir) / "me_rrrD_sap" / "xtbopt.xyz"
    if not ref_path.exists():
        ref_path = Path(data_dir) / "me_rrrD_sap" / "traj.xyz"
    u = mda.Universe(str(ref_path), format="XYZ")
    ref_frame = u.trajectory[0].positions.copy()
    ref_elements = [a.name for a in u.atoms]
    return ref_frame, ref_elements


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


# ───────────────────────────────────────────────────────────────
# Step 1: Load PCA model
# ───────────────────────────────────────────────────────────────

def load_pca_model(pca_dir):
    """
    Load existing PCA model (mean + components) from .npy files.

    Returns
    -------
    mean : ndarray (168,)
    components : ndarray (10, 168)
    """
    pca_dir = Path(pca_dir)
    mean_path = pca_dir / "pca_mean_eu8_nochrom.npy"
    comp_path = pca_dir / "pca_components_eu8_nochrom.npy"

    mean = np.load(str(mean_path))
    components = np.load(str(comp_path))

    print(f"  Loaded PCA mean: {mean_path} (shape={mean.shape})")
    print(f"  Loaded PCA components: {comp_path} (shape={components.shape})")

    return mean, components


# ───────────────────────────────────────────────────────────────
# Step 2: Project com_md trajectories
# ───────────────────────────────────────────────────────────────

def project_com_md(com_dir, pca_mean, pca_components, ref_frame, ref_elements,
                   mappings_path, scaffold_path):
    """
    Project com_md trajectories into PCA space.

    For each system: load PDB+XTC, extract MOL, apply mappings, align, project.

    Parameters
    ----------
    com_dir : str or Path
        Base directory with com_md system subdirectories.
    pca_mean : ndarray (168,)
    pca_components : ndarray (10, 168)
    ref_frame : ndarray (n_atoms, 3)
    ref_elements : list of str
    mappings_path : str or Path
    scaffold_path : str or Path

    Returns
    -------
    com_proj : ndarray (total_frames, N_PCS)
    com_metadata : list of tuples — (system, species, isomer, handedness, conformer, frame)
    frame_counts : dict — {system: n_frames}
    """
    # Load mappings
    with open(mappings_path) as f:
        atom_mappings = json.load(f)
    with open(scaffold_path) as f:
        scaffold_mapping = json.load(f)

    n_scaffold = scaffold_mapping["n_scaffold"]

    # Identify feature atoms from reference (SAP ordering)
    SCAFFOLD_HEAVY_IDX = [i for i in range(n_scaffold) if ref_elements[i] != "H"]
    assert len(SCAFFOLD_HEAVY_IDX) == 68, \
        f"Expected 68 scaffold heavy atoms, got {len(SCAFFOLD_HEAVY_IDX)}"

    FEATURE_IDX = [i for i in SCAFFOLD_HEAVY_IDX if i not in EXCLUDE_FEATURES]
    assert len(FEATURE_IDX) == 56, \
        f"Expected 56 feature atoms, got {len(FEATURE_IDX)}"

    feature_idx_arr = np.array(FEATURE_IDX)
    ref_align_pos = ref_frame[np.array(ALIGN_EU8)]

    # Build phe_to_me reordering from scaffold_mapping
    phe_to_me_sap = scaffold_mapping["phe_to_me_sap"]

    all_proj = []
    all_metadata = []
    frame_counts = {}

    for sys_name in SYSTEMS:
        pdb_path = Path(com_dir) / f"{sys_name}_fp" / "v1.pdb"
        xtc_path = Path(com_dir) / f"{sys_name}_fp" / "v1.xtc"

        if not pdb_path.exists() or not xtc_path.exists():
            print(f"  [SKIP] {sys_name}: missing PDB or XTC")
            continue

        # Check XTC is non-empty
        if xtc_path.stat().st_size == 0:
            print(f"  [SKIP] {sys_name}: empty XTC")
            continue

        species, isomer, handedness, conformer = parse_system_name(sys_name)
        print(f"  Loading {sys_name} ({species}, {conformer}) ...", end=" ", flush=True)

        try:
            u = mda.Universe(str(pdb_path), str(xtc_path))
        except Exception as e:
            print(f"[ERROR] Cannot load: {e}")
            continue

        mol = u.select_atoms("resname MOL")
        if mol.n_atoms == 0:
            print(f"[SKIP] No MOL atoms found")
            continue

        # Verify expected atom count for ligand
        expected_atoms = 135 if species == "me" else 149
        assert mol.n_atoms == expected_atoms, \
            f"Expected {expected_atoms} atoms for {sys_name}, got {mol.n_atoms}"

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
            if species == "phe":
                # Build inverse: for each me_sap position j (0..126),
                # find the phe_sap atom that maps to j
                me_to_phe_sap = {}
                for phe_idx_str, me_idx in phe_to_me_sap.items():
                    me_to_phe_sap[int(me_idx)] = int(phe_idx_str)

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

        X_com = np.array(frame_features)  # (n_frames, 168)

        # Project using PCA model (manual: center + dot with components)
        X_centered = X_com - pca_mean
        proj = X_centered @ pca_components.T  # (n_frames, 10)

        all_proj.append(proj)
        for f in range(n_frames):
            all_metadata.append(
                (sys_name, species, isomer, handedness, conformer, f)
            )
        frame_counts[sys_name] = n_frames

        print(f"    Projected {n_frames} frames")

        # Free memory
        del u, frame_features, X_com, proj

    if len(all_proj) == 0:
        print("  [WARN] No com_md projections collected")
        return np.empty((0, N_PCS)), [], {}

    com_proj = np.vstack(all_proj)
    print(f"  Total com_md frames: {com_proj.shape[0]}")

    return com_proj, all_metadata, frame_counts


# ───────────────────────────────────────────────────────────────
# Step 3: Save CSVs
# ───────────────────────────────────────────────────────────────

def save_com_md_csv(com_proj, com_metadata, out_path):
    """
    Save com_md-only projection CSV.

    Columns: system, species, isomer, handedness, conformer, method, frame,
             PC1..PC10, Neff_pt1, Neff_pt2
    """
    rows = []
    for i, (sys_name, species, isomer, handedness, conformer, frame) in enumerate(com_metadata):
        row = {
            "system": sys_name,
            "species": species,
            "isomer": isomer,
            "handedness": handedness,
            "conformer": conformer,
            "method": "com_md",
            "frame": frame,
        }
        for pc in range(N_PCS):
            row[f"PC{pc+1}"] = com_proj[i, pc]
        rows.append(row)

    df = pd.DataFrame(rows)

    # Compute N_eff per system and merge
    neff_data = []
    for sys_name in sorted(df["system"].unique()):
        sub = df[df["system"] == sys_name]
        neff1 = compute_neff(sub["PC1"].values)
        neff2 = compute_neff(sub["PC2"].values)
        neff_data.append({
            "system": sys_name,
            "Neff_pt1": neff1,
            "Neff_pt2": neff2,
        })
    neff_df = pd.DataFrame(neff_data)
    df = df.merge(neff_df, on="system", how="left")

    # Column order
    col_order = ["system", "species", "isomer", "handedness", "conformer",
                 "method", "frame", "PC1", "PC2", "Neff_pt1", "Neff_pt2"]
    # Add remaining PCs
    for pc in range(3, N_PCS + 1):
        col_order.append(f"PC{pc}")

    df = df[col_order]

    out_path = Path(out_path)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(str(out_path), index=False)
    print(f"  Saved com_md CSV: {out_path} ({len(df)} rows)")

    # Report per-system N_eff
    for _, row in neff_df.iterrows():
        print(f"    {row['system']}: Neff_pt1={row['Neff_pt1']}, Neff_pt2={row['Neff_pt2']}")

    return df


def build_joint_3way_csv(com_df, existing_joint_path, out_path):
    """
    Build 3-way joint projection CSV by combining existing xtb+solv
    with new com_md data.

    Parameters
    ----------
    com_df : DataFrame
        com_md-only projection DataFrame.
    existing_joint_path : str or Path
        Path to existing xtb+solv joint projection CSV (read-only).
    out_path : str or Path
        Output path for 3-way joint CSV.
    """
    existing_df = pd.read_csv(existing_joint_path)
    print(f"  Existing joint CSV: {len(existing_df)} rows")
    print(f"    Methods: {existing_df['method'].unique().tolist()}")
    print(f"    Systems: {existing_df['system'].nunique()}")

    # Ensure compatible columns
    # Existing: system, species, isomer, handedness, conformer, frame, method, PC1..PC10
    # New: system, species, isomer, handedness, conformer, method, frame, PC1, PC2, Neff_pt1, Neff_pt2, PC3..PC10

    # Add Neff_pt1, Neff_pt2 columns to existing data (will be NaN for xtb/solv)
    existing_df["Neff_pt1"] = np.nan
    existing_df["Neff_pt2"] = np.nan

    # Reorder com_df columns to match
    common_cols = ["system", "species", "isomer", "handedness", "conformer",
                   "method", "frame", "PC1", "PC2", "Neff_pt1", "Neff_pt2"]
    for pc in range(3, N_PCS + 1):
        common_cols.append(f"PC{pc}")

    com_aligned = com_df[common_cols]
    existing_aligned = existing_df[common_cols]

    # Concatenate
    joint_df = pd.concat([existing_aligned, com_aligned], ignore_index=True)

    out_path = Path(out_path)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    joint_df.to_csv(str(out_path), index=False)

    print(f"  Saved 3-way joint CSV: {out_path} ({len(joint_df)} rows)")
    print(f"    xtb: {(joint_df['method'] == 'xtb').sum()} rows")
    print(f"    solv: {(joint_df['method'] == 'solv').sum()} rows")
    print(f"    com_md: {(joint_df['method'] == 'com_md').sum()} rows")

    # Check for NaN
    nan_count = joint_df[["PC1", "PC2"]].isnull().sum().sum()
    if nan_count > 0:
        print(f"  [WARN] {nan_count} NaN values in PC1/PC2 columns")
    else:
        print(f"  [OK] No NaN values in PC1/PC2 columns")

    return joint_df


# ───────────────────────────────────────────────────────────────
# Main pipeline
# ───────────────────────────────────────────────────────────────

def run_pipeline(args):
    """Execute the full C1 pipeline."""
    t_start = time.time()

    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    print("=" * 70)
    print("C1: Project com_md into xTB eu8_nochrom PCA space")
    print("=" * 70)

    # ── Step 1: Load PCA model ────────────────────────────────
    print("\n[Step 1] Loading existing PCA model ...")
    t1 = time.time()

    pca_mean, pca_components = load_pca_model(args.pca_dir)

    t1_end = time.time()
    print(f"  Step 1 runtime: {t1_end - t1:.1f} s")

    # ── Step 2: Load xTB reference frame ─────────────────────
    print("\n[Step 2] Loading xTB reference frame for Kabsch alignment ...")
    t2 = time.time()

    ref_frame, ref_elements = load_reference_frame(args.data_dir)
    print(f"  Reference: {ref_frame.shape[0]} atoms")

    t2_end = time.time()
    print(f"  Step 2 runtime: {t2_end - t2:.1f} s")

    # ── Step 3: Project com_md ────────────────────────────────
    print("\n[Step 3] Projecting com_md trajectories into PCA space ...")
    t3 = time.time()

    com_proj, com_metadata, frame_counts = project_com_md(
        args.com_dir, pca_mean, pca_components, ref_frame, ref_elements,
        args.mappings, args.scaffold
    )

    t3_end = time.time()
    print(f"  Step 3 runtime: {t3_end - t3:.1f} s")

    # ── Step 4: Validate ──────────────────────────────────────
    print("\n[Step 4] Validation ...")
    t4 = time.time()

    if com_proj.shape[0] == 0:
        print("  [ERROR] No com_md projections — cannot continue")
        return None

    # NaN check
    nan_count = np.isnan(com_proj).sum()
    if nan_count > 0:
        print(f"  [WARN] {nan_count} NaN values in com_md projections")
    else:
        print(f"  [OK] No NaN values in com_md projections")

    # Per-PC range
    for pc in range(N_PCS):
        pc_vals = com_proj[:, pc]
        print(f"  PC{pc+1}: min={pc_vals.min():.2f}, max={pc_vals.max():.2f}, "
              f"mean={pc_vals.mean():.2f}, std={pc_vals.std(ddof=1):.2f}")

    # Per-system frame counts
    n_systems_with_data = len(frame_counts)
    print(f"  Systems with com_md data: {n_systems_with_data}/16")
    for sys_name, n_frames in sorted(frame_counts.items()):
        print(f"    {sys_name}: {n_frames} frames")

    t4_end = time.time()
    print(f"  Step 4 runtime: {t4_end - t4:.1f} s")

    # ── Step 5: Save com_md CSV ──────────────────────────────
    print("\n[Step 5] Saving com_md projection CSV ...")
    t5 = time.time()

    com_csv_path = out_dir / "com_md_projection_eu8_nochrom.csv"
    com_df = save_com_md_csv(com_proj, com_metadata, com_csv_path)

    t5_end = time.time()
    print(f"  Step 5 runtime: {t5_end - t5:.1f} s")

    # ── Step 6: Build 3-way joint CSV ─────────────────────────
    print("\n[Step 6] Building 3-way joint projection CSV ...")
    t6 = time.time()

    joint_csv_path = out_dir / "joint_projection_3way_eu8_nochrom.csv"
    joint_df = build_joint_3way_csv(com_df, args.existing_joint, joint_csv_path)

    t6_end = time.time()
    print(f"  Step 6 runtime: {t6_end - t6:.1f} s")

    # ── Summary ──────────────────────────────────────────────
    total_time = time.time() - t_start
    print("\n" + "=" * 70)
    print("C1 SUMMARY")
    print("=" * 70)
    print(f"  Total runtime: {total_time:.1f} s")
    print(f"  com_md frames: {com_proj.shape[0]}")
    print(f"  com_md systems with data: {len(frame_counts)}/16")
    print(f"  com_md CSV: {com_csv_path} ({len(com_df)} rows)")
    print(f"  3-way joint CSV: {joint_csv_path} ({len(joint_df)} rows)")
    print(f"  Outputs:")
    print(f"    {com_csv_path}")
    print(f"    {joint_csv_path}")
    print("=" * 70)

    return {
        "com_proj": com_proj,
        "com_df": com_df,
        "joint_df": joint_df,
        "frame_counts": frame_counts,
    }


def build_parser():
    """Build argument parser."""
    parser = argparse.ArgumentParser(
        description="C1: Project com_md trajectories into xTB eu8_nochrom PCA space"
    )
    parser.add_argument("--data-dir", default="data",
                        help="Directory with xTB system subdirectories (for reference frame)")
    parser.add_argument("--com-dir", default="com_md",
                        help="Directory with com_md system subdirectories")
    parser.add_argument("--out-dir", default="cross_dataset/com_md/analysis",
                        help="Output directory for CSV files")
    parser.add_argument("--mappings", default="analysis/atom_mappings.json",
                        help="Path to atom_mappings.json")
    parser.add_argument("--scaffold", default="analysis/scaffold_mapping.json",
                        help="Path to scaffold_mapping.json")
    parser.add_argument("--pca-dir", default="cross_dataset/analysis",
                        help="Directory with PCA model .npy files")
    parser.add_argument("--existing-joint",
                        default="cross_dataset/analysis/joint_projection_eu8_nochrom_xtb_solv.csv",
                        help="Path to existing xtb+solv joint projection CSV")
    return parser


if __name__ == "__main__":
    parser = build_parser()
    args = parser.parse_args()
    run_pipeline(args)
