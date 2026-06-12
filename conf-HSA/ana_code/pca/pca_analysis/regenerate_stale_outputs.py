#!/usr/bin/env python3
"""
Regenerate stale pc_modes_top10.nmd and reference_frame.npy from cached intermediates.

FIX-11 changed combined_pca.py to set ref_coords = mean.reshape(natoms, 3).astype(np.float64)
after the alignment loop (instead of using frame 0). The cached mean_coords.npy already
contains the correct mean — only reference_frame.npy and pc_modes_top10.nmd need updating.

This script avoids re-running the full pipeline (alignment + covariance + eigendecomposition).
"""
import sys
from pathlib import Path
import numpy as np

# Import write_nmd_file directly from combined_pca.py to ensure identical output
SCRIPT_DIR = Path(__file__).resolve().parent
sys.path.insert(0, str(SCRIPT_DIR))
from combined_pca import write_nmd_file

OUTDIR = SCRIPT_DIR  # pca_analysis/

def main():
    # Load cached intermediates
    mean = np.load(OUTDIR / "mean_coords.npy")
    eigvals = np.load(OUTDIR / "eigenvalues.npy")
    eigvecs = np.load(OUTDIR / "eigenvectors.npy")

    natoms = mean.shape[0] // 3
    print(f"natoms = {natoms}")
    print(f"eigenvalues shape: {eigvals.shape}")
    print(f"eigenvectors shape: {eigvecs.shape}")

    # Compute ref_coords the FIX-11 way: converged mean structure
    ref_coords = mean.reshape(natoms, 3).astype(np.float64)

    # Overwrite reference_frame.npy
    np.save(OUTDIR / "reference_frame.npy", ref_coords)
    print(f"Saved reference_frame.npy  shape={ref_coords.shape} dtype={ref_coords.dtype}")

    # Regenerate NMD file using the same function as combined_pca.py
    nmd_path = write_nmd_file(eigvals, eigvecs, ref_coords, n_modes=10,
                               outpath=str(OUTDIR / "pc_modes_top10.nmd"))
    print(f"Saved {nmd_path}")

    # Quick sanity checks
    ref_check = np.load(OUTDIR / "reference_frame.npy")
    assert ref_check.shape == (natoms, 3), f"Unexpected shape: {ref_check.shape}"
    assert np.allclose(ref_check, ref_coords), "reference_frame mismatch after save/load"
    print("Sanity checks passed.")

if __name__ == "__main__":
    main()
