#!/usr/bin/env python3
"""
A5_nmd_generation.py — Generate NMD file from eu8_nochrom PCA eigenvectors.

Re-runs the eu8_nochrom PCA via T3_joint_pca_eu8.py to obtain signed
pca.components_ vectors, then writes a VMD Normal Mode Wizard (NMWiz) format
NMD file with displacement arrows for the top 10 PCs on 56 feature atoms.

Reference structure: data/me_rrrD_sap/xtbopt.xyz (135 atoms, SAP ordering).

Output: analysis/joint_pca_eu8_nochrom_top10.nmd

Usage:
    python scripts/A5_nmd_generation.py
"""

import argparse
import os
import shutil
import sys

import numpy as np

# ───────────────────────────────────────────────────────────────
# Add scripts directory to path for T3 import
# ───────────────────────────────────────────────────────────────
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
PROJECT_ROOT = os.path.dirname(SCRIPT_DIR)
sys.path.insert(0, SCRIPT_DIR)

from T3_joint_pca_eu8 import run_pca, build_parser, backup_if_exists


# ───────────────────────────────────────────────────────────────
# Helper: parse single-frame XYZ
# ───────────────────────────────────────────────────────────────
def read_xyz_single(filepath):
    """
    Parse single-frame XYZ file.

    Parameters
    ----------
    filepath : str
        Path to XYZ file.

    Returns
    -------
    elements : list of str
        Element symbols.
    coords : ndarray, shape (n_atoms, 3)
        Cartesian coordinates in Angstrom.
    """
    with open(filepath) as f:
        lines = f.readlines()
    n_atoms = int(lines[0].strip())
    elements = []
    coords = np.empty((n_atoms, 3), dtype=np.float64)
    for i in range(n_atoms):
        parts = lines[2 + i].split()
        elements.append(parts[0])
        coords[i] = [float(parts[1]), float(parts[2]), float(parts[3])]
    return elements, coords


# ───────────────────────────────────────────────────────────────
# Main
# ───────────────────────────────────────────────────────────────
def main():
    parser = argparse.ArgumentParser(
        description="Generate NMD file from eu8_nochrom PCA eigenvectors")
    parser.add_argument("--ref-xyz", default="data/me_rrrD_sap/xtbopt.xyz",
                        help="Reference XYZ file (135 atoms, SAP ordering)")
    parser.add_argument("--out", default="analysis/joint_pca_eu8_nochrom_top10.nmd",
                        help="Output NMD file path")
    parser.add_argument("--n-modes", type=int, default=10,
                        help="Number of PCA modes to include")
    parser.add_argument("--max-disp", type=float, default=2.0,
                        help="Target max displacement (Å) for visualization scaling")
    args = parser.parse_args()

    print("=" * 65)
    print("A5: NMD File Generation (eu8_nochrom)")
    print("=" * 65)

    # ── Step 1: Re-run eu8_nochrom PCA ────────────────────────
    print("\n[1/5] Re-running eu8_nochrom PCA via T3_joint_pca_eu8.py ...")
    pca_parser = build_parser()
    pca_args = pca_parser.parse_args([
        "--tag", "eu8_nochrom",
        "--exclude-features", "102,103,104,105,107,109,113,114,115,116,118,120"
    ])
    result = run_pca(pca_args)

    pca = result["pca"]
    feature_idx = result["feature_idx"]  # list of 56 SAP atom indices
    ref_elements = result["ref_elements"]  # scaffold elements list

    n_feat = len(feature_idx)
    eigenvalues = pca.explained_variance_  # (10,)
    evr = pca.explained_variance_ratio_
    cumulative = np.cumsum(evr)

    # Reshape components from (10, 168) to (10, 56, 3)
    components_3d = pca.components_.reshape(args.n_modes, n_feat, 3)

    print(f"  PCA model obtained: {args.n_modes} modes, {n_feat} feature atoms")
    print(f"  Eigenvalues (PC1-10): {eigenvalues}")

    # ── Step 2: Read reference structure ───────────────────────
    print(f"\n[2/5] Reading reference structure: {args.ref_xyz} ...")
    if not os.path.exists(args.ref_xyz):
        sys.exit(f"ERROR: Reference XYZ not found: {args.ref_xyz}")

    ref_elements_xyz, ref_coords = read_xyz_single(args.ref_xyz)
    n_atoms = len(ref_elements_xyz)
    print(f"  Reference: {n_atoms} atoms")
    if n_atoms != 135:
        print(f"  [WARN] Expected 135 atoms, got {n_atoms}")

    # ── Step 3: Compute per-atom displacement vectors ──────────
    print(f"\n[3/5] Computing per-atom displacement vectors for {args.n_modes} modes ...")

    # Physical displacement = unit eigenvector * sqrt(eigenvalue)
    # This gives RMS fluctuation amplitude
    feat_set = set(feature_idx)
    feat_to_j = {aidx: j for j, aidx in enumerate(feature_idx)}

    # Build displacement array: (n_modes, n_atoms, 3)
    # Non-feature atoms get zero displacement
    all_displacements = np.zeros((args.n_modes, n_atoms, 3), dtype=np.float64)

    for m in range(args.n_modes):
        for aidx in feature_idx:
            if aidx < n_atoms:
                j = feat_to_j[aidx]
                all_displacements[m, aidx] = (
                    components_3d[m, j] * np.sqrt(eigenvalues[m])
                )

    # ── Step 4: Scale displacements so max ~2 Å ───────────────
    print(f"\n[4/5] Scaling displacements (max target = {args.max_disp:.1f} Å) ...")

    # Find global max displacement magnitude across all modes and atoms
    max_mag = 0.0
    for m in range(args.n_modes):
        mags = np.sqrt(np.sum(all_displacements[m] ** 2, axis=1))
        mode_max = mags.max()
        max_mag = max(max_mag, mode_max)
        print(f"  Mode {m+1}: max displacement (raw) = {mode_max:.4f} Å")

    if max_mag > 0:
        scale = args.max_disp / max_mag
        all_displacements *= scale
        print(f"  Scale factor: {scale:.6f} (max raw = {max_mag:.4f} Å)")

    # Verify scaled max
    for m in range(args.n_modes):
        mags = np.sqrt(np.sum(all_displacements[m] ** 2, axis=1))
        print(f"  Mode {m+1}: max displacement (scaled) = {mags.max():.4f} Å")

    # ── Step 5: Write NMD file ────────────────────────────────
    print(f"\n[5/5] Writing NMD file: {args.out} ...")

    out_dir = os.path.dirname(args.out)
    if out_dir:
        os.makedirs(out_dir, exist_ok=True)

    backup_if_exists(args.out)

    with open(args.out, "w") as f:
        # Header
        f.write("nmwiz_format nmd\n")
        f.write("title Joint PCA Modes (eu8_nochrom) - top 10 PCs\n")

        # Atom names (element symbols from xtbopt.xyz)
        f.write("atomnames " + " ".join(ref_elements_xyz) + "\n")

        # Residue IDs (all 1 - single residue molecule)
        f.write("resids " + " ".join(["1"] * n_atoms) + "\n")

        # Residue names (all LIG)
        f.write("resnames " + " ".join(["LIG"] * n_atoms) + "\n")

        # Coordinates: flattened x,y,z for all atoms
        flat_coords = ref_coords.flatten()  # (405,)
        f.write("coordinates " + " ".join(f"{c:.6f}" for c in flat_coords) + "\n")

        # Modes: each mode with eigenvalue and displacement vectors
        for m in range(args.n_modes):
            ev = eigenvalues[m]
            f.write(f"mode {m+1} {ev:.6f}\n")
            # Displacements: flattened dx,dy,dz for all 135 atoms
            flat_disp = all_displacements[m].flatten()  # (405,)
            f.write(" ".join(f"{d:.8f}" for d in flat_disp) + "\n")

    # ── Summary diagnostics ────────────────────────────────────
    print("\n" + "=" * 65)
    print("NMD GENERATION SUMMARY")
    print("=" * 65)
    print(f"  Total atoms in NMD:       {n_atoms}")
    print(f"  Feature atoms (non-zero): {n_feat}")
    print(f"  Reference structure:      {args.ref_xyz}")
    print(f"  Number of modes:          {args.n_modes}")
    print(f"  Max displacement target:  {args.max_disp:.1f} Å")
    print(f"\n  Eigenvalues and variance:")
    for i in range(args.n_modes):
        print(f"    PC{i+1:2d}: eigenvalue = {eigenvalues[i]:.6f}, "
              f"variance = {evr[i]*100:.1f}%, "
              f"cumulative = {cumulative[i]*100:.1f}%")
    print(f"\n  Per-mode max displacement (scaled):")
    for m in range(args.n_modes):
        mags = np.sqrt(np.sum(all_displacements[m] ** 2, axis=1))
        feat_mags = mags[list(feature_idx)]
        max_feat = feat_mags.max()
        print(f"    Mode {m+1:2d}: max = {max_feat:.4f} Å "
              f"(atom {feature_idx[feat_mags.argmax()]})")
    print(f"\n  Output: {args.out}")
    file_size = os.path.getsize(args.out)
    with open(args.out) as f:
        n_lines = sum(1 for _ in f)
    print(f"  File size: {file_size:,} bytes")
    print(f"  Line count: {n_lines}")
    print("=" * 65)


if __name__ == "__main__":
    main()
