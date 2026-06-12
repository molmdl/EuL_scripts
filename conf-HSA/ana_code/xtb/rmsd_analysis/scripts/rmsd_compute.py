#!/usr/bin/env python3
"""
rmsd_compute.py — Compute per-frame RMSD for all 16 chelate complex systems
across three MD methods (xtb, solv_md, com_md).

For each method, align each frame to its own frame-0 reference and compute
heavy-atom RMSD. For xtb, also compute Eu+9-coordinated-atom RMSD.
Additionally, all methods compute rmsd_heavy_comref_A using the com_md
frame-0 MOL heavy-atom positions as a common reference structure, enabling
direct cross-method RMSD comparison.

Key design: alignment indices for Eu+8/Eu+9 are determined from the
conformer type (SAP vs TSAP). SAP uses canonical indices [0..7,54] and
[0..7,54,63]; TSAP uses the inverse atom mapping from atom_mappings.json
to find the corresponding indices [0..6,84,106] and [0..6,84,106,107].

The common reference (rmsd_heavy_comref_A) uses com_md frame-0 as the
universal reference structure. For com_md, this equals rmsd_heavy_A
within float precision (self-reference). For xtb/solv_md, values > 0
reflect structural offset from the com_md reference.

Output: 4 CSV files in out-dir:
  - rmsd_xtb.csv     (per-frame xtb RMSD; columns: rmsd_heavy_A,
                      rmsd_eu9_A, rmsd_heavy_comref_A)
  - rmsd_solv_md.csv (per-frame solv_md RMSD; columns: rmsd_heavy_A,
                      rmsd_eu8_A, rmsd_heavy_comref_A)
  - rmsd_com_md.csv  (per-frame com_md RMSD; columns: rmsd_heavy_A,
                      rmsd_eu8_A, rmsd_heavy_comref_A)
  - rmsd_summary.csv (summary stats per system per method)

Usage:
  python rmsd_analysis/scripts/rmsd_compute.py \\
      --data-dir data --solv-dir solv_md --com-dir com_md \\
      --out-dir rmsd_analysis/csv --method all

  python rmsd_analysis/scripts/rmsd_compute.py \\
      --data-dir data --solv-dir solv_md --com-dir com_md \\
      --out-dir rmsd_analysis/csv --method xtb

  python rmsd_analysis/scripts/rmsd_compute.py \\
      --data-dir data --solv-dir solv_md --com-dir com_md \\
      --out-dir rmsd_analysis/csv --method solv_md --system phe_sssD_sap
"""

import argparse
import gc
import json
import sys
import time
from pathlib import Path

import MDAnalysis as mda
import numpy as np
import pandas as pd

# ───────────────────────────────────────────────────────────────
# Constants
# ───────────────────────────────────────────────────────────────

SYSTEMS = [
    "me_rrrD_sap", "me_rrrD_tsap", "me_rrrL_sap", "me_rrrL_tsap",
    "me_sssD_sap", "me_sssD_tsap", "me_sssL_sap", "me_sssL_tsap",
    "phe_rrrD_sap", "phe_rrrD_tsap", "phe_rrrL_sap", "phe_rrrL_tsap",
    "phe_sssD_sap", "phe_sssD_tsap", "phe_sssL_sap", "phe_sssL_tsap",
]

# Eu+8 alignment indices in SAP ordering: 8 donors (3O+5N) + Eu
ALIGN_EU8_SAP = sorted([0, 1, 2, 3, 4, 5, 6, 7, 54])

# Eu+9 alignment indices in SAP ordering: 9 donors (3O+5N+cap N63) + Eu
ALIGN_EU9_SAP = sorted([0, 1, 2, 3, 4, 5, 6, 7, 54, 63])


# ───────────────────────────────────────────────────────────────
# Kabsch alignment
# ───────────────────────────────────────────────────────────────

def kabsch(P, Q):
    """
    Kabsch algorithm for optimal rotation. R rotates P onto Q.

    Parameters
    ----------
    P, Q : ndarray, shape (N, 3)
        Mobile (P) and reference (Q) coordinates.

    Returns
    -------
    R : ndarray (3, 3)
        Optimal rotation matrix. Apply as: (coords - P_cent) @ R + Q_cent
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


# ───────────────────────────────────────────────────────────────
# Helper: parse system name into components
# ───────────────────────────────────────────────────────────────

def parse_system_name(system: str):
    """Parse system name like 'me_rrrD_sap' into components."""
    parts = system.split("_")
    return {
        "species": parts[0],       # me or phe
        "isomer": parts[1][:3],    # rrr or sss
        "handedness": parts[1][3], # D or L
        "conformer": parts[2],     # sap or tsap
    }


# ───────────────────────────────────────────────────────────────
# Helper: TSAP→SAP index mapping
# ───────────────────────────────────────────────────────────────

def load_tsap_mappings(mappings_path: Path):
    """
    Load TSAP→SAP atom index mappings and build inverse (SAP→TSAP).

    The atom_mappings.json format:
      me_tsap_to_sap[tsap_index] = sap_index
      phe_tsap_to_sap[tsap_index] = sap_index

    Returns
    -------
    sap_to_tsap : dict
        {species: list} where list[sap_index] = tsap_index
    """
    with open(mappings_path) as f:
        data = json.load(f)

    sap_to_tsap = {}
    for species in ("me", "phe"):
        key = f"{species}_tsap_to_sap"
        if key not in data:
            continue
        tsap_to_sap = data[key]
        n = len(tsap_to_sap)
        inv = [0] * n
        for tsap_idx, sap_idx in enumerate(tsap_to_sap):
            inv[sap_idx] = tsap_idx
        sap_to_tsap[species] = inv

    return sap_to_tsap


def get_eu8_align_idx(conformer: str, species: str, sap_to_tsap: dict):
    """
    Get Eu+8 alignment indices for the given conformer and species.

    For SAP: uses canonical indices ALIGN_EU8_SAP.
    For TSAP: uses inverse mapping to find corresponding TSAP indices.

    Returns
    -------
    align_idx : list of int
        Sorted alignment atom indices in the native ordering.
    """
    if conformer == "sap":
        return ALIGN_EU8_SAP
    else:
        inv = sap_to_tsap[species]
        return sorted([inv[i] for i in ALIGN_EU8_SAP])


def get_eu9_align_idx(conformer: str, species: str, sap_to_tsap: dict):
    """
    Get Eu+9 alignment indices for the given conformer and species.

    For SAP: uses canonical indices ALIGN_EU9_SAP.
    For TSAP: uses inverse mapping to find corresponding TSAP indices.

    Returns
    -------
    align_idx : list of int
        Sorted alignment atom indices in the native ordering.
    """
    if conformer == "sap":
        return ALIGN_EU9_SAP
    else:
        inv = sap_to_tsap[species]
        return sorted([inv[i] for i in ALIGN_EU9_SAP])


# ───────────────────────────────────────────────────────────────
# Helper: heavy-atom indices
# ───────────────────────────────────────────────────────────────

def get_heavy_indices(atom_group):
    """
    Return indices of non-hydrogen atoms within an MDAnalysis AtomGroup.

    Uses atom.element when available; falls back to name-based detection
    for atoms with empty element records (common in GROMACS PDB files).
    """
    heavy_idx = []
    for i, atom in enumerate(atom_group):
        # Skip if element is definitely hydrogen
        if atom.element == "H":
            continue
        # Skip if element is empty but name suggests hydrogen
        if atom.element in ("", None) and atom.name.startswith("H") and len(atom.name) <= 2:
            continue
        # Everything else is heavy
        heavy_idx.append(i)
    return heavy_idx


def get_heavy_indices_from_names(names):
    """Return indices of non-hydrogen atoms from a list of atom names (xtb XYZ)."""
    heavy_idx = []
    for i, name in enumerate(names):
        if not name.startswith("H"):
            heavy_idx.append(i)
    return heavy_idx


# ───────────────────────────────────────────────────────────────
# Helper: load com_md frame-0 as common reference
# ───────────────────────────────────────────────────────────────

def load_comref_heavy(com_dir: Path, sys_name: str):
    """
    Load com_md frame-0 MOL heavy-atom positions as common reference.

    Loads the PDB, selects MOL residue, extracts heavy-atom positions
    at frame 0, then deletes the Universe and calls gc.collect() to
    free memory before the caller loads a trajectory Universe.

    Parameters
    ----------
    com_dir : Path
        Path to com_md root directory (contains {sys_name}_fp/v1.pdb).
    sys_name : str
        System name, e.g. 'me_rrrD_sap'.

    Returns
    -------
    (ref_heavy_pos, heavy_idx) : tuple of (ndarray, list[int])
        ref_heavy_pos: shape (N_heavy, 3), heavy-atom positions at frame 0.
        heavy_idx: indices of heavy atoms within MOL selection.
    None : if PDB file does not exist or MOL selection fails.
    """
    pdb_path = com_dir / f"{sys_name}_fp" / "v1.pdb"
    if not pdb_path.exists():
        print(f"  WARN: comref PDB not found: {pdb_path}", file=sys.stderr)
        return None
    u_ref = mda.Universe(str(pdb_path))
    mol_ref = u_ref.select_atoms("resname MOL")
    if mol_ref.n_atoms == 0:
        print(f"  WARN: no MOL atoms in comref {pdb_path}", file=sys.stderr)
        del u_ref
        gc.collect()
        return None
    u_ref.trajectory[0]
    ref_heavy_idx = get_heavy_indices(mol_ref)
    ref_heavy_pos = mol_ref.positions[ref_heavy_idx].copy()
    del u_ref
    gc.collect()
    return ref_heavy_pos, ref_heavy_idx


# ───────────────────────────────────────────────────────────────
# Task 1a: xtb RMSD computation
# ───────────────────────────────────────────────────────────────

def compute_xtb_rmsd(data_dir: Path, out_dir: Path, sap_to_tsap: dict,
                     com_dir: Path, systems=None):
    """
    Compute per-frame RMSD for xtb trajectories.

    Three RMSD columns per frame:
    1. rmsd_heavy_A: all-heavy alignment to own frame-0 reference
    2. rmsd_eu9_A: Eu+9 alignment to own frame-0 reference
    3. rmsd_heavy_comref_A: all-heavy alignment to com_md frame-0 common reference

    Reference for columns 1-2 = frame 0 of own trajectory.
    Common reference for column 3 = com_md frame-0 MOL heavy atoms.
    Eu+9 alignment indices determined from conformer type (SAP/TSAP).
    """
    if systems is None:
        systems = SYSTEMS

    rows = []
    for sys_name in systems:
        t0 = time.time()
        traj_path = data_dir / sys_name / "traj.xyz"

        if not traj_path.exists():
            print(f"  SKIP {sys_name}: traj.xyz not found", file=sys.stderr)
            continue

        comp = parse_system_name(sys_name)
        conformer = comp["conformer"]
        species = comp["species"]

        comref_result = load_comref_heavy(com_dir, sys_name)
        if comref_result is None:
            print(f"  SKIP {sys_name}: comref not available", file=sys.stderr)
            continue
        comref_heavy_pos, comref_heavy_idx = comref_result

        print(f"  xtb: {sys_name} ...", end="", flush=True)

        u = mda.Universe(str(traj_path), in_memory=True)
        n_atoms = u.atoms.n_atoms
        n_frames = len(u.trajectory)

        atom_names = [atom.name for atom in u.atoms]
        heavy_idx = get_heavy_indices_from_names(atom_names)
        n_heavy = len(heavy_idx)

        if len(heavy_idx) != len(comref_heavy_idx):
            print(f"  SKIP {sys_name}: heavy count mismatch "
                  f"(xtb={len(heavy_idx)}, comref={len(comref_heavy_idx)})",
                  file=sys.stderr)
            del u
            gc.collect()
            continue

        eu9_align_idx = get_eu9_align_idx(conformer, species, sap_to_tsap)
        n_eu9 = len(eu9_align_idx)

        u.trajectory[0]
        ref_positions = u.atoms.positions.copy()
        ref_heavy = ref_positions[heavy_idx]
        ref_eu9 = ref_positions[eu9_align_idx]

        for frame_idx in range(n_frames):
            u.trajectory[frame_idx]
            frame_pos = u.atoms.positions.copy()

            R_h, P_cent_h, Q_cent_h = kabsch(frame_pos[heavy_idx], ref_heavy)
            aligned_heavy = (frame_pos[heavy_idx] - P_cent_h) @ R_h + Q_cent_h
            rmsd_heavy = np.sqrt(np.mean((aligned_heavy - ref_heavy) ** 2))

            R_e, P_cent_e, Q_cent_e = kabsch(frame_pos[eu9_align_idx], ref_eu9)
            aligned_eu9 = (frame_pos[eu9_align_idx] - P_cent_e) @ R_e + Q_cent_e
            rmsd_eu9 = np.sqrt(np.mean((aligned_eu9 - ref_eu9) ** 2))

            R_c, P_cent_c, Q_cent_c = kabsch(frame_pos[heavy_idx], comref_heavy_pos)
            aligned_comref = (frame_pos[heavy_idx] - P_cent_c) @ R_c + Q_cent_c
            rmsd_comref = np.sqrt(np.mean((aligned_comref - comref_heavy_pos) ** 2))

            time_ps = frame_idx * 0.05

            rows.append({
                "system": sys_name,
                "species": species,
                "isomer": comp["isomer"],
                "handedness": comp["handedness"],
                "conformer": conformer,
                "frame": frame_idx,
                "time_ps": round(time_ps, 2),
                "rmsd_heavy_A": rmsd_heavy,
                "rmsd_eu9_A": rmsd_eu9,
                "rmsd_heavy_comref_A": rmsd_comref,
            })

        elapsed = time.time() - t0
        print(f" {n_frames} frames, {n_heavy} heavy, Eu+9={eu9_align_idx}, "
              f"comref={len(comref_heavy_idx)} heavy, {elapsed:.1f}s")

        del u
        gc.collect()

    df = pd.DataFrame(rows)
    csv_path = out_dir / "rmsd_xtb.csv"
    df.to_csv(csv_path, index=False)
    print(f"  Saved: {csv_path} ({len(df)} rows)")
    return df


# ───────────────────────────────────────────────────────────────
# Task 1b: solv_md RMSD computation
# ───────────────────────────────────────────────────────────────

def compute_solv_md_rmsd(solv_dir: Path, out_dir: Path, sap_to_tsap: dict,
                        com_dir: Path, systems=None):
    """
    Compute per-frame heavy-atom RMSD for solv_md trajectories.

    Three RMSD columns per frame:
    1. rmsd_heavy_A: all-heavy alignment to own frame-0 reference (primary)
    2. rmsd_eu8_A: Eu+8 alignment to own frame-0 reference (legacy)
    3. rmsd_heavy_comref_A: all-heavy alignment to com_md frame-0 common reference

    Reference for columns 1-2 = frame 0 of own trajectory.
    Common reference for column 3 = com_md frame-0 MOL heavy atoms.
    Frame-by-frame processing to avoid OOM.
    """
    if systems is None:
        systems = SYSTEMS

    rows = []
    for sys_name in systems:
        t0 = time.time()
        tpr_path = solv_dir / sys_name / "prod_0.tpr"
        xtc_path = solv_dir / sys_name / "solv_all.xtc"

        if not tpr_path.exists() or not xtc_path.exists():
            print(f"  SKIP {sys_name}: solv_md files not found", file=sys.stderr)
            continue

        comp = parse_system_name(sys_name)
        conformer = comp["conformer"]
        species = comp["species"]

        comref_result = load_comref_heavy(com_dir, sys_name)
        if comref_result is None:
            print(f"  SKIP {sys_name}: comref not available", file=sys.stderr)
            continue
        comref_heavy_pos, comref_heavy_idx = comref_result

        print(f"  solv_md: {sys_name} ...", end="", flush=True)

        u = mda.Universe(str(tpr_path), str(xtc_path))
        mol = u.select_atoms("resname MOL")
        n_mol = mol.n_atoms

        eu8_align_idx = get_eu8_align_idx(conformer, species, sap_to_tsap)

        heavy_idx = get_heavy_indices(mol)
        n_heavy = len(heavy_idx)

        if len(heavy_idx) != len(comref_heavy_idx):
            print(f"  SKIP {sys_name}: heavy count mismatch "
                  f"(solv={len(heavy_idx)}, comref={len(comref_heavy_idx)})",
                  file=sys.stderr)
            del u
            gc.collect()
            continue

        u.trajectory[0]
        ref_mol_pos = mol.positions.copy()
        ref_heavy = ref_mol_pos[heavy_idx]
        ref_eu8_align = ref_mol_pos[eu8_align_idx]

        n_frames = len(u.trajectory)

        for frame_idx in range(n_frames):
            u.trajectory[frame_idx]
            mol_pos = mol.positions.copy()

            R_h, P_cent_h, Q_cent_h = kabsch(mol_pos[heavy_idx], ref_heavy)
            aligned_heavy = (mol_pos[heavy_idx] - P_cent_h) @ R_h + Q_cent_h
            rmsd_heavy = np.sqrt(np.mean((aligned_heavy - ref_heavy) ** 2))

            R_e, P_cent_e, Q_cent_e = kabsch(mol_pos[eu8_align_idx], ref_eu8_align)
            aligned_eu8 = (mol_pos - P_cent_e) @ R_e + Q_cent_e
            rmsd_eu8 = np.sqrt(np.mean((aligned_eu8[heavy_idx] - ref_heavy) ** 2))

            R_c, P_cent_c, Q_cent_c = kabsch(mol_pos[heavy_idx], comref_heavy_pos)
            aligned_comref = (mol_pos[heavy_idx] - P_cent_c) @ R_c + Q_cent_c
            rmsd_comref = np.sqrt(np.mean((aligned_comref - comref_heavy_pos) ** 2))

            rows.append({
                "system": sys_name,
                "species": species,
                "isomer": comp["isomer"],
                "handedness": comp["handedness"],
                "conformer": conformer,
                "frame": frame_idx,
                "rmsd_heavy_A": rmsd_heavy,
                "rmsd_eu8_A": rmsd_eu8,
                "rmsd_heavy_comref_A": rmsd_comref,
            })

        elapsed = time.time() - t0
        print(f" {n_frames} frames, {n_heavy} heavy (primary), "
              f"Eu+8={eu8_align_idx} (legacy), comref={len(comref_heavy_idx)} heavy, "
              f"{elapsed:.1f}s")

        del u
        gc.collect()

    df = pd.DataFrame(rows)
    csv_path = out_dir / "rmsd_solv_md.csv"
    df.to_csv(csv_path, index=False)
    print(f"  Saved: {csv_path} ({len(df)} rows)")
    return df


# ───────────────────────────────────────────────────────────────
# Task 1c: com_md RMSD computation
# ───────────────────────────────────────────────────────────────

def compute_com_md_rmsd(com_dir: Path, out_dir: Path, sap_to_tsap: dict,
                        systems=None):
    """
    Compute per-frame heavy-atom RMSD for com_md trajectories.

    Three RMSD columns per frame:
    1. rmsd_heavy_A: all-heavy alignment to own frame-0 reference (primary)
    2. rmsd_eu8_A: Eu+8 alignment to own frame-0 reference (legacy)
    3. rmsd_heavy_comref_A: all-heavy alignment to com_md frame-0 common reference

    For com_md, the common reference IS the self-reference, so
    rmsd_heavy_comref_A should equal rmsd_heavy_A within float precision.
    Both are computed independently for verifiability.

    Reference for columns 1-2 = frame 0 of own trajectory.
    Common reference for column 3 = com_md frame-0 MOL heavy atoms.
    Frame-by-frame processing to avoid OOM.
    """
    if systems is None:
        systems = SYSTEMS

    rows = []
    for sys_name in systems:
        t0 = time.time()
        pdb_path = com_dir / f"{sys_name}_fp" / "v1.pdb"
        xtc_path = com_dir / f"{sys_name}_fp" / "v1.xtc"

        if not pdb_path.exists() or not xtc_path.exists():
            print(f"  SKIP {sys_name}: com_md files not found", file=sys.stderr)
            continue

        comp = parse_system_name(sys_name)
        conformer = comp["conformer"]
        species = comp["species"]

        comref_result = load_comref_heavy(com_dir, sys_name)
        if comref_result is None:
            print(f"  SKIP {sys_name}: comref not available", file=sys.stderr)
            continue
        comref_heavy_pos, comref_heavy_idx = comref_result

        print(f"  com_md: {sys_name} ...", end="", flush=True)

        u = mda.Universe(str(pdb_path), str(xtc_path))
        mol = u.select_atoms("resname MOL")
        n_mol = mol.n_atoms

        eu8_align_idx = get_eu8_align_idx(conformer, species, sap_to_tsap)

        heavy_idx = get_heavy_indices(mol)
        n_heavy = len(heavy_idx)

        if len(heavy_idx) != len(comref_heavy_idx):
            print(f"  SKIP {sys_name}: heavy count mismatch "
                  f"(com={len(heavy_idx)}, comref={len(comref_heavy_idx)})",
                  file=sys.stderr)
            del u
            gc.collect()
            continue

        u.trajectory[0]
        ref_mol_pos = mol.positions.copy()
        ref_heavy = ref_mol_pos[heavy_idx]
        ref_eu8_align = ref_mol_pos[eu8_align_idx]

        n_frames = len(u.trajectory)

        for frame_idx in range(n_frames):
            u.trajectory[frame_idx]
            mol_pos = mol.positions.copy()

            R_h, P_cent_h, Q_cent_h = kabsch(mol_pos[heavy_idx], ref_heavy)
            aligned_heavy = (mol_pos[heavy_idx] - P_cent_h) @ R_h + Q_cent_h
            rmsd_heavy = np.sqrt(np.mean((aligned_heavy - ref_heavy) ** 2))

            R_e, P_cent_e, Q_cent_e = kabsch(mol_pos[eu8_align_idx], ref_eu8_align)
            aligned_eu8 = (mol_pos - P_cent_e) @ R_e + Q_cent_e
            rmsd_eu8 = np.sqrt(np.mean((aligned_eu8[heavy_idx] - ref_heavy) ** 2))

            R_c, P_cent_c, Q_cent_c = kabsch(mol_pos[heavy_idx], comref_heavy_pos)
            aligned_comref = (mol_pos[heavy_idx] - P_cent_c) @ R_c + Q_cent_c
            rmsd_comref = np.sqrt(np.mean((aligned_comref - comref_heavy_pos) ** 2))

            rows.append({
                "system": sys_name,
                "species": species,
                "isomer": comp["isomer"],
                "handedness": comp["handedness"],
                "conformer": conformer,
                "frame": frame_idx,
                "rmsd_heavy_A": rmsd_heavy,
                "rmsd_eu8_A": rmsd_eu8,
                "rmsd_heavy_comref_A": rmsd_comref,
            })

        elapsed = time.time() - t0
        print(f" {n_frames} frames, {n_heavy} heavy (primary), "
              f"Eu+8={eu8_align_idx} (legacy), comref={len(comref_heavy_idx)} heavy, "
              f"{elapsed:.1f}s")

        del u
        gc.collect()

    df = pd.DataFrame(rows)
    csv_path = out_dir / "rmsd_com_md.csv"
    df.to_csv(csv_path, index=False)
    print(f"  Saved: {csv_path} ({len(df)} rows)")
    return df


# ───────────────────────────────────────────────────────────────
# Summary statistics
# ───────────────────────────────────────────────────────────────

def build_summary(xtb_df, solv_df, com_df, out_dir: Path):
    """Build rmsd_summary.csv with per-system per-method statistics."""
    summary_rows = []

    def _summarize(df, method, rmsd_col):
        for sys_name in df["system"].unique():
            sub = df[df["system"] == sys_name]
            vals = sub[rmsd_col].values
            summary_rows.append({
                "system": sys_name,
                "method": method,
                "n_frames": len(sub),
                "n_heavy": None,
                "rmsd_mean": np.mean(vals),
                "rmsd_std": np.std(vals, ddof=1),
                "rmsd_median": np.median(vals),
                "rmsd_min": np.min(vals),
                "rmsd_max": np.max(vals),
                "rmsd_q10": np.percentile(vals, 10),
                "rmsd_q90": np.percentile(vals, 90),
                "rmsd_range": np.max(vals) - np.min(vals),
            })

    if xtb_df is not None and len(xtb_df) > 0:
        _summarize(xtb_df, "xtb_heavy", "rmsd_heavy_A")
        _summarize(xtb_df, "xtb_eu9", "rmsd_eu9_A")
        if "rmsd_heavy_comref_A" in xtb_df.columns:
            _summarize(xtb_df, "xtb_comref", "rmsd_heavy_comref_A")
    if solv_df is not None and len(solv_df) > 0:
        _summarize(solv_df, "solv_md", "rmsd_heavy_A")
        _summarize(solv_df, "solv_md_eu8", "rmsd_eu8_A")
        if "rmsd_heavy_comref_A" in solv_df.columns:
            _summarize(solv_df, "solv_md_comref", "rmsd_heavy_comref_A")
    if com_df is not None and len(com_df) > 0:
        _summarize(com_df, "com_md", "rmsd_heavy_A")
        _summarize(com_df, "com_md_eu8", "rmsd_eu8_A")
        if "rmsd_heavy_comref_A" in com_df.columns:
            _summarize(com_df, "com_md_comref", "rmsd_heavy_comref_A")

    df_summary = pd.DataFrame(summary_rows)
    csv_path = out_dir / "rmsd_summary.csv"
    df_summary.to_csv(csv_path, index=False)
    print(f"  Saved: {csv_path} ({len(df_summary)} rows)")
    return df_summary


# ───────────────────────────────────────────────────────────────
# Validation checks
# ───────────────────────────────────────────────────────────────

def validate_results(xtb_df, solv_df, com_df, summary_df):
    """Print validation checks and flag any issues."""
    issues = []

    print("\n=== Validation Checks ===")

    # 1. All 16 systems processed per method
    for method_name, df in [("xtb", xtb_df), ("solv_md", solv_df), ("com_md", com_df)]:
        if df is None or len(df) == 0:
            print(f"  [SKIP] {method_name}: no data (not computed)")
            continue
        n_sys = df["system"].nunique()
        status = "PASS" if n_sys == 16 else "FAIL"
        if n_sys != 16:
            issues.append(f"{method_name}: only {n_sys}/16 systems processed")
        print(f"  [{status}] {method_name}: {n_sys}/16 systems")

    # 2. Frame counts
    for method_name, df, expected in [
        ("xtb", xtb_df, 2000),
        ("solv_md", solv_df, 4004),
        ("com_md", com_df, 4004),
    ]:
        if df is None or len(df) == 0:
            continue
        for sys_name in df["system"].unique():
            n = len(df[df["system"] == sys_name])
            ok = (n == expected) or (method_name == "com_md" and sys_name == "me_rrrD_sap" and n >= 3591)
            if not ok:
                issues.append(f"{method_name}/{sys_name}: {n} frames (expected {expected})")
        print(f"  [INFO] {method_name}: frame counts checked")

    # 3. Frame-0 RMSD ~ 0 (within float precision ~1e-6)
    for method_name, df, col in [
        ("xtb_heavy", xtb_df, "rmsd_heavy_A"),
        ("xtb_eu9", xtb_df, "rmsd_eu9_A"),
        ("solv_md", solv_df, "rmsd_heavy_A"),
        ("com_md", com_df, "rmsd_heavy_A"),
    ]:
        if df is None or len(df) == 0:
            continue
        f0 = df[df["frame"] == 0]
        max_f0 = f0[col].max()
        status = "PASS" if max_f0 < 1e-6 else "FAIL"
        if max_f0 >= 1e-6:
            issues.append(f"{method_name}: frame-0 RMSD max = {max_f0:.2e} (> 1e-6)")
        print(f"  [{status}] {method_name}: frame-0 RMSD max = {max_f0:.2e}")

    # Also check legacy Eu+8 columns
    for method_name, df, col in [
        ("solv_md_eu8", solv_df, "rmsd_eu8_A"),
        ("com_md_eu8", com_df, "rmsd_eu8_A"),
    ]:
        if df is None or len(df) == 0:
            continue
        if "rmsd_eu8_A" not in df.columns:
            continue
        f0 = df[df["frame"] == 0]
        max_f0 = f0[col].max()
        status = "PASS" if max_f0 < 1e-6 else "FAIL"
        if max_f0 >= 1e-6:
            issues.append(f"{method_name}: frame-0 RMSD max = {max_f0:.2e} (> 1e-6)")
        print(f"  [{status}] {method_name}: frame-0 RMSD max = {max_f0:.2e}")

    # Common-reference (comref) validation checks
    # com_md frame-0 self-consistency: rmsd_heavy_comref_A at frame 0 should be ~0
    if com_df is not None and len(com_df) > 0 and "rmsd_heavy_comref_A" in com_df.columns:
        f0_comref = com_df[com_df["frame"] == 0]["rmsd_heavy_comref_A"].max()
        status = "PASS" if f0_comref < 1e-6 else "FAIL"
        if f0_comref >= 1e-6:
            issues.append(f"com_md_comref: frame-0 self-consistency max = {f0_comref:.2e} (> 1e-6)")
        print(f"  [{status}] com_md_comref: frame-0 self-consistency max = {f0_comref:.2e}")

    # xtb/solv_md frame-0 comref offset: should be > 0 (structural offset) and < 5.0
    for method_name, df in [("xtb_comref", xtb_df), ("solv_md_comref", solv_df)]:
        if df is None or len(df) == 0:
            continue
        if "rmsd_heavy_comref_A" not in df.columns:
            continue
        f0_comref = df[df["frame"] == 0]["rmsd_heavy_comref_A"]
        if len(f0_comref) == 0:
            continue
        all_positive = (f0_comref > 0).all()
        max_val = f0_comref.max()
        within_bound = max_val < 5.0
        status = "PASS" if (all_positive and within_bound) else "WARN"
        if not all_positive:
            issues.append(f"{method_name}: frame-0 comref has non-positive values")
        if not within_bound:
            issues.append(f"{method_name}: frame-0 comref max = {max_val:.4f} (> 5.0 A)")
        print(f"  [{status}] {method_name}: frame-0 comref offset "
              f"(min={f0_comref.min():.4f}, max={max_val:.4f}) A")

    # comref RMSD range sanity
    for method_name, df, hi in [
        ("xtb_comref", xtb_df, 15.0),
        ("solv_md_comref", solv_df, 30.0),
        ("com_md_comref", com_df, 15.0),
    ]:
        if df is None or len(df) == 0:
            continue
        if "rmsd_heavy_comref_A" not in df.columns:
            continue
        overall_max = df["rmsd_heavy_comref_A"].max()
        status = "PASS" if overall_max <= hi else "WARN"
        if overall_max > hi:
            issues.append(f"{method_name}: comref RMSD max = {overall_max:.4f} (> {hi})")
        print(f"  [{status}] {method_name}: comref RMSD max = {overall_max:.4f} A "
              f"(bound {hi})")

    # com_md cross-column consistency: rmsd_heavy_comref_A ≈ rmsd_heavy_A at all frames
    if com_df is not None and len(com_df) > 0 and "rmsd_heavy_comref_A" in com_df.columns:
        diff = (com_df["rmsd_heavy_comref_A"] - com_df["rmsd_heavy_A"]).abs().max()
        status = "PASS" if diff < 1e-5 else "FAIL"
        if diff >= 1e-5:
            issues.append(f"com_md cross-column: |comref - heavy|_max = {diff:.2e} (> 1e-5)")
        print(f"  [{status}] com_md cross-column: |comref - heavy|_max = {diff:.2e}")

    # 4. RMSD range sanity
    for method_name, df, col, hi in [
        ("xtb_heavy", xtb_df, "rmsd_heavy_A", 15.0),
        ("xtb_eu9", xtb_df, "rmsd_eu9_A", 15.0),
        ("solv_md", solv_df, "rmsd_heavy_A", 30.0),  # solv_md can be large
        ("com_md", com_df, "rmsd_heavy_A", 15.0),
    ]:
        if df is None or len(df) == 0:
            continue
        overall_min = df[col].min()
        overall_max = df[col].max()
        out_of_range = overall_max > hi
        status = "PASS" if not out_of_range else "WARN"
        if out_of_range:
            issues.append(f"{method_name}: RMSD range [{overall_min:.4f}, {overall_max:.4f}] exceeds {hi}")
        print(f"  [{status}] {method_name}: RMSD range [{overall_min:.4f}, {overall_max:.4f}] Å")

    # Legacy Eu+8-aligned range checks
    for method_name, df, col, hi in [
        ("solv_md_eu8", solv_df, "rmsd_eu8_A", 45.0),  # Eu+8 aligned values can be large
        ("com_md_eu8", com_df, "rmsd_eu8_A", 15.0),
    ]:
        if df is None or len(df) == 0:
            continue
        if "rmsd_eu8_A" not in df.columns:
            continue
        overall_min = df[col].min()
        overall_max = df[col].max()
        out_of_range = overall_max > hi
        status = "PASS" if not out_of_range else "WARN"
        if out_of_range:
            issues.append(f"{method_name}: RMSD range [{overall_min:.4f}, {overall_max:.4f}] exceeds {hi}")
        print(f"  [{status}] {method_name}: RMSD range [{overall_min:.4f}, {overall_max:.4f}] Å")

    # 5. No NaN
    for method_name, df in [("xtb", xtb_df), ("solv_md", solv_df), ("com_md", com_df)]:
        if df is None or len(df) == 0:
            continue
        n_nan = df.isnull().sum().sum()
        status = "PASS" if n_nan == 0 else "FAIL"
        if n_nan > 0:
            issues.append(f"{method_name}: {n_nan} NaN values")
        print(f"  [{status}] {method_name}: {n_nan} NaN values")

    # 6. CSV row counts
    for method_name, df, expected_approx in [
        ("xtb", xtb_df, 32000),
        ("solv_md", solv_df, 64064),
        ("com_md", com_df, 63951),
    ]:
        if df is None or len(df) == 0:
            continue
        n = len(df)
        pct = abs(n - expected_approx) / expected_approx * 100
        status = "PASS" if pct < 5 else "WARN"
        print(f"  [{status}] {method_name}: {n} rows (expected ~{expected_approx}, {pct:.1f}% diff)")

    if issues:
        print(f"\n  ISSUES ({len(issues)}):")
        for iss in issues:
            print(f"    - {iss}")
    else:
        print("\n  All checks passed!")

    return issues


# ───────────────────────────────────────────────────────────────
# Main
# ───────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(
        description="Compute per-frame RMSD for 16 chelate systems across 3 MD methods"
    )
    parser.add_argument("--data-dir", type=Path, default=Path("data"),
                        help="Path to xtb data directory (default: data)")
    parser.add_argument("--solv-dir", type=Path, default=Path("solv_md"),
                        help="Path to solv_md directory (default: solv_md)")
    parser.add_argument("--com-dir", type=Path, default=Path("com_md"),
                        help="Path to com_md directory (default: com_md)")
    parser.add_argument("--out-dir", type=Path, default=Path("rmsd_analysis/csv"),
                        help="Output directory for CSV files (default: rmsd_analysis/csv)")
    parser.add_argument("--mappings", type=Path, default=Path("analysis/atom_mappings.json"),
                        help="Path to atom_mappings.json (default: analysis/atom_mappings.json)")
    parser.add_argument("--method", choices=["xtb", "solv_md", "com_md", "all"],
                        default="all", help="Which method to compute (default: all)")
    parser.add_argument("--system", type=str, default=None,
                        help="Compute only this system (default: all 16)")

    args = parser.parse_args()
    args.out_dir.mkdir(parents=True, exist_ok=True)

    systems = [args.system] if args.system else SYSTEMS

    # Load TSAP→SAP mappings
    sap_to_tsap = load_tsap_mappings(args.mappings)
    print(f"Loaded TSAP mappings for species: {list(sap_to_tsap.keys())}")
    # Verify the mapping
    for species, inv in sap_to_tsap.items():
        print(f"  {species}: SAP[54](Eu) -> TSAP[{inv[54]}], "
              f"SAP[7](N) -> TSAP[{inv[7]}], "
              f"SAP[63](cap N) -> TSAP[{inv[63]}]")

    xtb_df = None
    solv_df = None
    com_df = None

    print(f"\nRMSD Computation: method={args.method}, systems={len(systems)}")

    # Sequential execution: xtb -> solv_md -> com_md
    if args.method in ("xtb", "all"):
        print("\n--- xtb (dual alignment: all-heavy + Eu+9, comref common reference) ---")
        xtb_df = compute_xtb_rmsd(args.data_dir, args.out_dir, sap_to_tsap,
                                   args.com_dir, systems)

    if args.method in ("solv_md", "all"):
        print("\n--- solv_md (dual: all-heavy primary + Eu+8 legacy, comref common reference) ---")
        solv_df = compute_solv_md_rmsd(args.solv_dir, args.out_dir, sap_to_tsap,
                                        args.com_dir, systems)

    if args.method in ("com_md", "all"):
        print("\n--- com_md (dual: all-heavy primary + Eu+8 legacy, comref common reference) ---")
        com_df = compute_com_md_rmsd(args.com_dir, args.out_dir, sap_to_tsap, systems)

    # Build summary from any computed methods
    print("\n--- Summary Statistics ---")
    summary_df = build_summary(xtb_df, solv_df, com_df, args.out_dir)

    # Validation
    validate_results(xtb_df, solv_df, com_df, summary_df)

    print("\nDone.")


if __name__ == "__main__":
    main()
