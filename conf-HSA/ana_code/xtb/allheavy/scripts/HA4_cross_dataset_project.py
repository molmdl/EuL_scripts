#!/usr/bin/env python3
"""
HA4_cross_dataset_project.py -- Re-fit allheavy PCA on xTB, project solv_md + com_md.

Pipeline:
1. Re-fit PCA on xTB data (all 16 systems, all-heavy alignment on 68 scaffold atoms)
2. Validate re-fit against existing HA1 projections (Pearson r > 0.999)
3. Save PCA state as .npy files
4. Project solv_md trajectories (TPR+XTC, resname MOL)
5. Project com_md trajectories (PDB+XTC, resname MOL, manual projection)
6. Save combined CSVs (xtb+solv, xtb+solv+com_md)

Usage:
  python scripts/HA4_cross_dataset_project.py \
      --data-dir data \
      --solv-dir solv_md \
      --com-dir com_md \
      --out-dir analysis \
      --tag allheavy_scaffold \
      --mappings ../analysis/atom_mappings.json \
      --scaffold ../analysis/scaffold_mapping.json \
      --validate analysis/joint_pca_projection_allheavy_scaffold.csv
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
from scipy.stats import pearsonr

warnings.filterwarnings("ignore", message="Unknown element")
warnings.filterwarnings("ignore", message="Failed to guess the mass")

N_PCS = 10

SYSTEMS = [
    "me_rrrD_sap", "me_rrrD_tsap", "me_rrrL_sap", "me_rrrL_tsap",
    "me_sssD_sap", "me_sssD_tsap", "me_sssL_sap", "me_sssL_tsap",
    "phe_rrrD_sap", "phe_rrrD_tsap", "phe_rrrL_sap", "phe_rrrL_tsap",
    "phe_sssD_sap", "phe_sssD_tsap", "phe_sssL_sap", "phe_sssL_tsap",
]


def kabsch(P, Q):
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
    mapped_coords = np.empty_like(coords)
    if coords.ndim == 3:
        mapped_coords[:, mapping, :] = coords
    else:
        mapped_coords[mapping, :] = coords
    return mapped_coords


def parse_system_name(sys_name):
    parts = sys_name.split("_")
    species = parts[0]
    isomer = parts[1][:3]
    handedness = parts[1][3]
    conformer = parts[2]
    return species, isomer, handedness, conformer


def load_xyz_trajectory(traj_path):
    u = mda.Universe(str(traj_path), format="XYZ")
    n_frames = len(u.trajectory)
    n_atoms = u.atoms.n_atoms
    coords = np.empty((n_frames, n_atoms, 3), dtype=np.float64)
    for i, ts in enumerate(u.trajectory):
        coords[i] = ts.positions
    elements = [a.name for a in u.atoms]
    return coords, elements


def get_reference_frame(data_dir):
    ref_path = Path(data_dir) / "me_rrrD_sap" / "xtbopt.xyz"
    if not ref_path.exists():
        ref_path = Path(data_dir) / "me_rrrD_sap" / "traj.xyz"
    u = mda.Universe(str(ref_path), format="XYZ")
    ref_frame = u.trajectory[0].positions.copy()
    ref_elements = [a.name for a in u.atoms]
    return ref_frame, ref_elements


def compute_neff(series):
    arr = np.asarray(series, dtype=float)
    n = len(arr)
    if n < 2:
        return n
    rho = np.corrcoef(arr[:-1], arr[1:])[0, 1]
    if np.isnan(rho) or rho <= 0:
        return n
    return max(1, int(n * (1 - rho) / (1 + rho)))


def refit_pca_xtb(data_dir, mappings_path, scaffold_path):
    with open(mappings_path) as f:
        atom_mappings = json.load(f)
    with open(scaffold_path) as f:
        scaffold_mapping = json.load(f)

    n_scaffold = scaffold_mapping["n_scaffold"]

    ref_frame, ref_elements = get_reference_frame(data_dir)
    SCAFFOLD_HEAVY_IDX = [i for i in range(n_scaffold) if ref_elements[i] != "H"]
    assert len(SCAFFOLD_HEAVY_IDX) == 68, f"Expected 68 scaffold heavy atoms, got {len(SCAFFOLD_HEAVY_IDX)}"

    FEATURE_IDX = SCAFFOLD_HEAVY_IDX
    ALIGN_IDX = SCAFFOLD_HEAVY_IDX
    ref_align_pos = ref_frame[np.array(ALIGN_IDX)]

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

        if conformer == "tsap":
            mapping_key = f"{species}_tsap_to_sap"
            mapping = atom_mappings[mapping_key]
            coords = apply_tsap_mapping(coords, mapping)

        aligned = np.empty_like(coords)
        for f in range(n_frames):
            mobile = coords[f]
            mobile_align = mobile[np.array(ALIGN_IDX)]
            R, P_cent, Q_cent = kabsch(mobile_align, ref_align_pos)
            aligned[f] = (mobile - P_cent) @ R + Q_cent

        feature_idx_arr = np.array(FEATURE_IDX)
        feat_coords = aligned[:, feature_idx_arr, :]
        features = feat_coords.reshape(n_frames, -1)
        all_features.append(features)

        for f in range(n_frames):
            all_metadata.append((sys_name, species, isomer, handedness, conformer, f))

        del coords, aligned

    X = np.vstack(all_features)
    print(f"  xTB feature matrix: {X.shape}")

    pca = PCA(n_components=N_PCS, svd_solver="randomized", random_state=42)
    xtb_proj = pca.fit_transform(X)

    evr = pca.explained_variance_ratio_
    cumulative = np.cumsum(evr)
    print(f"  PC1+PC2 cumulative: {cumulative[1]:.4f} ({cumulative[1]*100:.1f}%)")

    return pca, FEATURE_IDX, ref_elements, ref_frame, xtb_proj, all_metadata


def validate_xtb_projection(xtb_proj, xtb_metadata, existing_csv_path, tag):
    existing_df = pd.read_csv(existing_csv_path)

    our_rows = []
    for i, (sys_name, species, isomer, handedness, conformer, frame) in enumerate(xtb_metadata):
        our_rows.append({
            "system": sys_name,
            "frame": frame,
            "PC1_new": xtb_proj[i, 0],
            "PC2_new": xtb_proj[i, 1],
        })
    our_df = pd.DataFrame(our_rows)

    merged = existing_df.merge(our_df, on=["system", "frame"], how="inner")

    pc1_col = "pc1" if "pc1" in merged.columns else "PC1"
    pc2_col = "pc2" if "pc2" in merged.columns else "PC2"

    r_pc1, p_pc1 = pearsonr(merged[pc1_col], merged["PC1_new"])
    r_pc2, p_pc2 = pearsonr(merged[pc2_col], merged["PC2_new"])

    print(f"  PC1 Pearson r = {r_pc1:.6f} (p = {p_pc1:.2e})")
    print(f"  PC2 Pearson r = {r_pc2:.6f} (p = {p_pc2:.2e})")

    if r_pc1 > 0.999:
        print(f"  [PASS] PC1 correlation > 0.999 -- re-fit matches existing projections")
    else:
        print(f"  [WARN] PC1 correlation = {r_pc1:.6f} < 0.999 -- possible discrepancy")

    return r_pc1, r_pc2


def save_pca_state(pca, out_dir, tag):
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    mean_path = out_dir / f"pca_mean_{tag}.npy"
    comp_path = out_dir / f"pca_components_{tag}.npy"

    np.save(str(mean_path), pca.mean_)
    np.save(str(comp_path), pca.components_)

    print(f"  Saved PCA mean: {mean_path} (shape={pca.mean_.shape})")
    print(f"  Saved PCA components: {comp_path} (shape={pca.components_.shape})")

    return mean_path, comp_path


def build_phe_reorder(scaffold_mapping):
    phe_to_me_sap = scaffold_mapping["phe_to_me_sap"]
    n_scaffold = scaffold_mapping["n_scaffold"]
    me_to_phe_sap = {}
    for phe_idx_str, me_idx in phe_to_me_sap.items():
        me_to_phe_sap[int(me_idx)] = int(phe_idx_str)
    reorder = np.array([me_to_phe_sap[j] for j in range(n_scaffold)])
    return reorder, n_scaffold


def project_solv_md(solv_dir, pca, feature_idx, ref_frame, mappings_path, scaffold_path, tag):
    with open(mappings_path) as f:
        atom_mappings = json.load(f)
    with open(scaffold_path) as f:
        scaffold_mapping = json.load(f)

    n_scaffold = scaffold_mapping["n_scaffold"]
    feature_idx_arr = np.array(feature_idx)
    SCAFFOLD_HEAVY_IDX = feature_idx
    ref_align_pos = ref_frame[np.array(SCAFFOLD_HEAVY_IDX)]

    phe_reorder, _ = build_phe_reorder(scaffold_mapping)

    all_proj = []
    all_metadata = []
    frame_counts = {}

    for sys_name in SYSTEMS:
        tpr_path = Path(solv_dir) / sys_name / "prod_0.tpr"
        xtc_path = Path(solv_dir) / sys_name / "solv_all.xtc"

        if not tpr_path.exists() or not xtc_path.exists():
            print(f"  [SKIP] {sys_name}: missing TPR or XTC")
            continue

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

        frame_features = []
        for frame_idx, ts in enumerate(u.trajectory):
            mol_coords = mol.positions.copy()

            if conformer == "tsap":
                mapping_key = f"{species}_tsap_to_sap"
                mapping = atom_mappings[mapping_key]
                mol_coords = apply_tsap_mapping(mol_coords, mapping)

            if species == "phe":
                mol_coords[:n_scaffold] = mol_coords[phe_reorder]

            mobile_align = mol_coords[np.array(SCAFFOLD_HEAVY_IDX)]
            R, P_cent, Q_cent = kabsch(mobile_align, ref_align_pos)
            mol_coords = (mol_coords - P_cent) @ R + Q_cent

            feat = mol_coords[feature_idx_arr].flatten()
            frame_features.append(feat)

        if len(frame_features) == 0:
            print(f"  [WARN] No frames processed for {sys_name}")
            continue

        X_solv = np.array(frame_features)
        proj = pca.transform(X_solv)

        all_proj.append(proj)
        for f in range(n_frames):
            all_metadata.append((sys_name, species, isomer, handedness, conformer, f))
        frame_counts[sys_name] = n_frames

        print(f"    Projected {n_frames} frames")

        del u, frame_features, X_solv, proj

    if len(all_proj) == 0:
        print("  [WARN] No solv_md projections collected")
        return np.empty((0, N_PCS)), [], {}

    solv_proj = np.vstack(all_proj)
    print(f"  Total solv_md frames: {solv_proj.shape[0]}")

    return solv_proj, all_metadata, frame_counts


def project_com_md(com_dir, pca_mean, pca_components, ref_frame, ref_elements,
                   mappings_path, scaffold_path, tag):
    with open(mappings_path) as f:
        atom_mappings = json.load(f)
    with open(scaffold_path) as f:
        scaffold_mapping = json.load(f)

    n_scaffold = scaffold_mapping["n_scaffold"]

    SCAFFOLD_HEAVY_IDX = [i for i in range(n_scaffold) if ref_elements[i] != "H"]
    assert len(SCAFFOLD_HEAVY_IDX) == 68, f"Expected 68 scaffold heavy atoms, got {len(SCAFFOLD_HEAVY_IDX)}"

    FEATURE_IDX = SCAFFOLD_HEAVY_IDX
    feature_idx_arr = np.array(FEATURE_IDX)
    ref_align_pos = ref_frame[np.array(SCAFFOLD_HEAVY_IDX)]

    phe_reorder, _ = build_phe_reorder(scaffold_mapping)

    all_proj = []
    all_metadata = []
    frame_counts = {}

    for sys_name in SYSTEMS:
        pdb_path = Path(com_dir) / sys_name / "fp" / "v1.pdb"
        xtc_path = Path(com_dir) / sys_name / "fp" / "v1.xtc"

        if not pdb_path.exists() or not xtc_path.exists():
            print(f"  [SKIP] {sys_name}: missing PDB or XTC")
            continue

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

        expected_atoms = 135 if species == "me" else 149
        assert mol.n_atoms == expected_atoms, \
            f"Expected {expected_atoms} atoms for {sys_name}, got {mol.n_atoms}"

        n_frames = len(u.trajectory)
        print(f"{n_frames} frames, {mol.n_atoms} MOL atoms")

        frame_features = []
        for frame_idx, ts in enumerate(u.trajectory):
            mol_coords = mol.positions.copy()

            if conformer == "tsap":
                mapping_key = f"{species}_tsap_to_sap"
                mapping = atom_mappings[mapping_key]
                mol_coords = apply_tsap_mapping(mol_coords, mapping)

            if species == "phe":
                mol_coords[:n_scaffold] = mol_coords[phe_reorder]

            mobile_align = mol_coords[np.array(SCAFFOLD_HEAVY_IDX)]
            R, P_cent, Q_cent = kabsch(mobile_align, ref_align_pos)
            mol_coords = (mol_coords - P_cent) @ R + Q_cent

            feat = mol_coords[feature_idx_arr].flatten()
            frame_features.append(feat)

        if len(frame_features) == 0:
            print(f"  [WARN] No frames processed for {sys_name}")
            continue

        X_com = np.array(frame_features)

        X_centered = X_com - pca_mean
        proj = X_centered @ pca_components.T

        all_proj.append(proj)
        for f in range(n_frames):
            all_metadata.append((sys_name, species, isomer, handedness, conformer, f))
        frame_counts[sys_name] = n_frames

        print(f"    Projected {n_frames} frames")

        del u, frame_features, X_com, proj

    if len(all_proj) == 0:
        print("  [WARN] No com_md projections collected")
        return np.empty((0, N_PCS)), [], {}

    com_proj = np.vstack(all_proj)
    print(f"  Total com_md frames: {com_proj.shape[0]}")

    return com_proj, all_metadata, frame_counts


def build_df(proj, metadata, method_name):
    rows = []
    for i, (sys_name, species, isomer, handedness, conformer, frame) in enumerate(metadata):
        row = {
            "system": sys_name,
            "species": species,
            "isomer": isomer,
            "handedness": handedness,
            "conformer": conformer,
            "frame": frame,
            "method": method_name,
        }
        for pc in range(N_PCS):
            row[f"PC{pc+1}"] = proj[i, pc]
        rows.append(row)
    return pd.DataFrame(rows)


def save_combined_csvs(xtb_proj, xtb_metadata, solv_proj, solv_metadata,
                       com_proj, com_metadata, out_dir, tag):
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    col_order = ["system", "species", "isomer", "handedness", "conformer", "frame", "method"]
    col_order += [f"PC{i+1}" for i in range(N_PCS)]

    df_xtb = build_df(xtb_proj, xtb_metadata, "xtb")
    df_solv = build_df(solv_proj, solv_metadata, "solv")
    df_com = build_df(com_proj, com_metadata, "com_md")

    df_xtb_solv = pd.concat([df_xtb, df_solv], ignore_index=True)[col_order]
    csv_2way = out_dir / f"joint_projection_{tag}_xtb_solv.csv"
    df_xtb_solv.to_csv(str(csv_2way), index=False)
    print(f"  Saved 2-way CSV: {csv_2way} ({len(df_xtb_solv)} rows)")
    print(f"    xtb rows: {len(df_xtb)}, solv rows: {len(df_solv)}")

    nan_2way = df_xtb_solv.isnull().sum().sum()
    if nan_2way > 0:
        print(f"  [WARN] {nan_2way} NaN values in 2-way CSV")
    else:
        print(f"  [OK] No NaN values in 2-way CSV")

    df_3way = pd.concat([df_xtb, df_solv, df_com], ignore_index=True)[col_order]
    csv_3way = out_dir / f"joint_projection_3way_{tag}.csv"
    df_3way.to_csv(str(csv_3way), index=False)
    print(f"  Saved 3-way CSV: {csv_3way} ({len(df_3way)} rows)")
    print(f"    xtb rows: {len(df_xtb)}, solv rows: {len(df_solv)}, com_md rows: {len(df_com)}")

    nan_3way = df_3way.isnull().sum().sum()
    if nan_3way > 0:
        print(f"  [WARN] {nan_3way} NaN values in 3-way CSV")
    else:
        print(f"  [OK] No NaN values in 3-way CSV")

    return df_xtb_solv, df_3way


def validate_projection(proj, metadata, frame_counts, label):
    if proj.shape[0] == 0:
        print(f"  [WARN] No {label} projections to validate")
        return

    nan_count = np.isnan(proj).sum()
    if nan_count > 0:
        print(f"  [WARN] {nan_count} NaN values in {label} projections")
    else:
        print(f"  [OK] No NaN values in {label} projections")

    for pc in range(min(3, N_PCS)):
        pc_vals = proj[:, pc]
        print(f"  PC{pc+1}: min={pc_vals.min():.2f}, max={pc_vals.max():.2f}, "
              f"mean={pc_vals.mean():.2f}, std={pc_vals.std(ddof=1):.2f}")

    n_systems_with_data = len(frame_counts)
    print(f"  Systems with {label} data: {n_systems_with_data}/16")
    for sys_name, n_frames in sorted(frame_counts.items()):
        print(f"    {sys_name}: {n_frames} frames")


def run_pipeline(args):
    t_start = time.time()

    data_dir = args.data_dir
    solv_dir = args.solv_dir
    com_dir = args.com_dir
    out_dir = Path(args.out_dir)
    tag = args.tag
    out_dir.mkdir(parents=True, exist_ok=True)

    print("=" * 70)
    print(f"HA4: Cross-dataset projection (allheavy alignment, tag={tag})")
    print("=" * 70)

    # Step 1: Re-fit PCA on xTB
    print("\n[Step 1] Re-fitting PCA on xTB data (68 scaffold heavy atoms, 204 feats) ...")
    t1 = time.time()

    pca, feature_idx, ref_elements, ref_frame, xtb_proj, xtb_metadata = \
        refit_pca_xtb(data_dir, args.mappings, args.scaffold)

    t1_end = time.time()
    print(f"  Step 1 runtime: {t1_end - t1:.1f} s")

    # Step 2: Validate xTB re-fit
    print("\n[Step 2] Validating xTB re-fit against existing projections ...")
    t2 = time.time()

    existing_csv = Path(args.validate) if args.validate else \
        out_dir / f"joint_pca_projection_{tag}.csv"
    r_pc1, r_pc2 = None, None
    if existing_csv.exists():
        r_pc1, r_pc2 = validate_xtb_projection(xtb_proj, xtb_metadata, str(existing_csv), tag)
    else:
        print(f"  [SKIP] Existing projection CSV not found: {existing_csv}")

    t2_end = time.time()
    print(f"  Step 2 runtime: {t2_end - t2:.1f} s")

    # Step 3: Save PCA state
    print("\n[Step 3] Saving PCA state ...")
    t3 = time.time()

    save_pca_state(pca, out_dir, tag)

    t3_end = time.time()
    print(f"  Step 3 runtime: {t3_end - t3:.1f} s")

    # Step 4: Project solv_md
    print("\n[Step 4] Projecting solv_md trajectories into PCA space ...")
    t4 = time.time()

    solv_proj, solv_metadata, solv_frame_counts = \
        project_solv_md(solv_dir, pca, feature_idx, ref_frame, args.mappings, args.scaffold, tag)

    t4_end = time.time()
    print(f"  Step 4 runtime: {t4_end - t4:.1f} s")

    # Step 5: Project com_md
    print("\n[Step 5] Projecting com_md trajectories into PCA space ...")
    t5 = time.time()

    com_proj, com_metadata, com_frame_counts = \
        project_com_md(com_dir, pca.mean_, pca.components_, ref_frame, ref_elements,
                       args.mappings, args.scaffold, tag)

    t5_end = time.time()
    print(f"  Step 5 runtime: {t5_end - t5:.1f} s")

    # Step 6: Save combined CSVs
    print("\n[Step 6] Saving combined CSVs ...")
    t6 = time.time()

    df_2way, df_3way = save_combined_csvs(
        xtb_proj, xtb_metadata, solv_proj, solv_metadata,
        com_proj, com_metadata, out_dir, tag
    )

    t6_end = time.time()
    print(f"  Step 6 runtime: {t6_end - t6:.1f} s")

    # Step 7: Validation
    print("\n[Step 7] Final validation ...")
    t7 = time.time()

    print("\n  xTB re-fit validation:")
    if r_pc1 is not None:
        print(f"    PC1 Pearson r: {r_pc1:.6f}")
    else:
        print(f"    [SKIP] No validation performed")

    print("\n  solv_md projection validation:")
    validate_projection(solv_proj, solv_metadata, solv_frame_counts, "solv_md")

    print("\n  com_md projection validation:")
    validate_projection(com_proj, com_metadata, com_frame_counts, "com_md")

    t7_end = time.time()
    print(f"  Step 7 runtime: {t7_end - t7:.1f} s")

    # Summary
    total_time = time.time() - t_start
    print("\n" + "=" * 70)
    print("HA4 SUMMARY")
    print("=" * 70)
    print(f"  Total runtime: {total_time:.1f} s")
    print(f"  xTB frames: {xtb_proj.shape[0]}")
    print(f"  solv_md frames: {solv_proj.shape[0]}")
    print(f"  com_md frames: {com_proj.shape[0]}")
    print(f"  solv_md systems: {len(solv_frame_counts)}/16")
    print(f"  com_md systems: {len(com_frame_counts)}/16")
    print(f"  2-way CSV: {out_dir / f'joint_projection_{tag}_xtb_solv.csv'} ({len(df_2way)} rows)")
    print(f"  3-way CSV: {out_dir / f'joint_projection_3way_{tag}.csv'} ({len(df_3way)} rows)")
    if r_pc1 is not None:
        print(f"  xTB re-fit PC1 Pearson r: {r_pc1:.6f}")
    evr = pca.explained_variance_ratio_
    print(f"  PC1+PC2 variance: {evr[0]*100:.1f}% + {evr[1]*100:.1f}% = {(evr[0]+evr[1])*100:.1f}%")
    print(f"  Outputs:")
    print(f"    {out_dir / f'pca_mean_{tag}.npy'}")
    print(f"    {out_dir / f'pca_components_{tag}.npy'}")
    print(f"    {out_dir / f'joint_projection_{tag}_xtb_solv.csv'}")
    print(f"    {out_dir / f'joint_projection_3way_{tag}.csv'}")
    print("=" * 70)

    return {
        "pca": pca,
        "xtb_proj": xtb_proj,
        "solv_proj": solv_proj,
        "com_proj": com_proj,
        "r_pc1": r_pc1,
    }


def build_parser():
    parser = argparse.ArgumentParser(
        description="HA4: Cross-dataset projection with all-heavy alignment")
    parser.add_argument("--data-dir", default="data",
                        help="Directory with xTB system subdirectories")
    parser.add_argument("--solv-dir", default="solv_md",
                        help="Directory with solv_md system subdirectories")
    parser.add_argument("--com-dir", default="com_md",
                        help="Directory with com_md system subdirectories")
    parser.add_argument("--out-dir", default="analysis",
                        help="Output directory for .npy and CSV files")
    parser.add_argument("--tag", default="allheavy_scaffold",
                        help="Output file tag")
    parser.add_argument("--mappings", default="../analysis/atom_mappings.json",
                        help="Path to atom_mappings.json")
    parser.add_argument("--scaffold", default="../analysis/scaffold_mapping.json",
                        help="Path to scaffold_mapping.json")
    parser.add_argument("--validate", default="analysis/joint_pca_projection_allheavy_scaffold.csv",
                        help="Path to existing projection CSV for validation")
    return parser


if __name__ == "__main__":
    parser = build_parser()
    args = parser.parse_args()
    run_pipeline(args)
