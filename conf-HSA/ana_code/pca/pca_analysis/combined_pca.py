import argparse
import os
import sys
import json
import itertools
import logging
from pathlib import Path
from typing import List, Tuple

import numpy as np
np.random.seed(42)  # Reproducibility for stochastic analyses
import MDAnalysis as mda
import scipy.linalg
from scipy.ndimage import gaussian_filter
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.ticker import MaxNLocator

logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")


def build_system_metadata():
    CHIRALITIES = ["me", "phe"]
    STEREOS = ["sss", "rrr"]
    ENANTIOMERS = ["L", "D"]
    LIGANDS = ["sap", "tsap"]
    return [
        {"name": f"{ch}_{st}{en}_{lg}", "chirality": ch, "stereo": st, "enantiomer": en, "ligand": lg}
        for ch, st, en, lg in itertools.product(CHIRALITIES, STEREOS, ENANTIOMERS, LIGANDS)
    ]


def load_systems(input_dir, systems_meta):
    BACKBONE_SEL = "backbone"
    frame_arrays = []
    system_indices = [0]
    natoms = None
    for sys in systems_meta:
        u = mda.Universe(f"{input_dir}/{sys['name']}/fp/v1.pdb", f"{input_dir}/{sys['name']}/fp/v1.xtc")
        bb = u.select_atoms(BACKBONE_SEL)
        if natoms is None:
            natoms = len(bb)
        else:
            assert len(bb) == natoms, f"{sys['name']}: expected {natoms} backbone atoms, got {len(bb)}"
        sys_frames = np.array(
            [bb.positions.astype(np.float32).ravel() for ts in u.trajectory],
            dtype=np.float32
        )
        frame_arrays.append(sys_frames)
        system_indices.append(system_indices[-1] + len(sys_frames))
        sys["n_frames"] = len(sys_frames)
        logging.info(f"Loaded {sys['name']}: {len(sys_frames)} frames")
    coords = np.vstack(frame_arrays)
    del frame_arrays
    system_indices = np.array(system_indices, dtype=np.int64)
    return coords, system_indices, natoms


def align_to_reference(coords, natoms, chunk_size=5000, max_iter=2):
    N, D = coords.shape
    ref = coords[0].reshape(natoms, 3).astype(np.float64)
    ref_centroid = ref.mean(axis=0)
    ref_centered = ref - ref_centroid
    ref_coords = ref.copy()
    aligned = np.empty_like(coords)
    prev_mean = None
    for iteration in range(max_iter):
        for start in range(0, N, chunk_size):
            end = min(start + chunk_size, N)
            chunk = coords[start:end].reshape(-1, natoms, 3)
            centroids = chunk.mean(axis=1, keepdims=True)
            chunk_centered = chunk - centroids
            H = np.einsum('nmk,mj->nkj', chunk_centered, ref_centered)
            U, S, Vh = np.linalg.svd(H)
            R = U @ Vh
            det = np.linalg.det(R)
            sign = np.where(det < 0, -1.0, 1.0)
            Vh[:, -1, :] *= sign[:, np.newaxis]
            R = U @ Vh
            rotated = chunk_centered @ R
            rotated += ref_centroid[np.newaxis, np.newaxis, :]
            aligned[start:end] = rotated.reshape(-1, D)
        mean = np.mean(aligned, axis=0).astype(np.float64)
        if prev_mean is not None:
            rmsd_change = np.sqrt(np.mean((mean - prev_mean) ** 2)) * np.sqrt(3)
            logging.info(f"Alignment iteration {iteration+1}: mean atom-RMSD change = {rmsd_change:.4f} Å")
            if rmsd_change < 0.01:
                logging.info(f"Alignment converged at iteration {iteration+1}")
                break
        else:
            logging.info(f"Alignment iteration 1 complete")
        prev_mean = mean.copy()
        if iteration < max_iter - 1:
            ref_centered = (mean.reshape(natoms, 3) - ref_centroid)
    ref_coords = mean.reshape(natoms, 3).astype(np.float64)
    return aligned, ref_coords, mean


def compute_covariance(aligned, mean, chunk_size=5000):
    N, D = aligned.shape
    cov = np.zeros((D, D), dtype=np.float64)
    for start in range(0, N, chunk_size):
        end = min(start + chunk_size, N)
        centered = (aligned[start:end] - mean).astype(np.float64)
        cov += centered.T @ centered
        del centered
    cov /= (N - 1)
    cov = (cov + cov.T) / 2
    total_variance = np.trace(cov)
    return cov, total_variance


def eigendecompose(cov, top_k=0):
    cov = (cov + cov.T) / 2
    if top_k == 0 or top_k >= cov.shape[0]:
        eigvals, eigvecs = scipy.linalg.eigh(cov)
    else:
        from scipy.sparse.linalg import eigsh
        eigvals, eigvecs = eigsh(cov, k=top_k, which="LM")
    idx = np.argsort(eigvals)[::-1]
    return eigvals[idx], eigvecs[:, idx]


def save_intermediates(outdir, systems_meta, system_indices, aligned, mean,
                       eigvals, eigvecs, ref_coords, total_variance, cov=None):
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    np.save(outdir / "aligned_coords.npy", aligned)
    np.save(outdir / "mean_coords.npy", mean)
    np.save(outdir / "eigenvalues.npy", eigvals)
    np.save(outdir / "eigenvectors.npy", eigvecs)
    np.save(outdir / "total_variance.npy", total_variance)
    np.save(outdir / "system_indices.npy", system_indices)
    np.save(outdir / "reference_frame.npy", ref_coords)
    metadata_json = [
        {"name": s["name"], "chirality": s["chirality"], "stereo": s["stereo"],
         "enantiomer": s["enantiomer"], "ligand": s["ligand"], "n_frames": s.get("n_frames", 0)}
        for s in systems_meta
    ]
    json.dump(metadata_json, (outdir / "system_metadata.json").open("w"), indent=2)
    if cov is not None:
        np.save(outdir / "covariance_matrix.npy", cov)
    cum_var = np.cumsum(eigvals) / total_variance
    logging.info(f"Top 2 PCs explain {cum_var[1]*100:.1f}% variance")
    logging.info(f"Top 10 PCs explain {cum_var[9]*100:.1f}% variance")
    logging.info(f"Eigenvalues computed: {len(eigvals)}, sum: {eigvals.sum():.2f}")
    logging.info(f"Total variance (trace of cov): {total_variance:.2f}")


def write_eigenvalue_csv(eigvals, total_variance, n_modes=10, outpath="eigenvalue_top10.csv"):
    assert total_variance > 0, f"Total variance is {total_variance:.6e} — covariance was degenerate"
    top_eigvals = eigvals[:n_modes]
    frac = top_eigvals / total_variance
    cumulative = np.cumsum(frac)
    with open(outpath, "w") as f:
        f.write("#Mode  Frac     Cumulative  Eigenval\n")
        for i in range(n_modes):
            f.write(f"{i+1:<6}{frac[i]:.5f}  {cumulative[i]:.5f}     {top_eigvals[i]:.5f}\n")
    return outpath


def write_nmd_file(eigvals, eigvecs, ref_coords, n_modes=10, outpath="pc_modes_top10.nmd"):
    assert ref_coords.ndim == 2 and ref_coords.shape[1] == 3, f"ref_coords shape {ref_coords.shape}, expected (N, 3)"
    natoms = ref_coords.shape[0]
    assert eigvecs.shape[0] == 3 * natoms, f"eigvecs rows {eigvecs.shape[0]} != 3*natoms {3*natoms}"
    top_eigvals = np.maximum(eigvals[:n_modes], 0.0)
    ref_flat = ref_coords.flatten()
    with open(outpath, "w") as f:
        f.write("name PC_combined_top10\n")
        f.write("coordinates " + " ".join(f"{c:.6f}" for c in ref_flat) + "\n")
        for i in range(n_modes):
            scale = np.sqrt(top_eigvals[i])
            mode_vec = eigvecs[:, i]  # Already flat: x1,y1,z1,x2,y2,z2,...
            f.write(f"mode {i+1} {scale:.6f} " + " ".join(f"{v:.6f}" for v in mode_vec) + "\n")
    return outpath


STEREO_KEYS = ["L_sap", "L_tsap", "D_sap", "D_tsap"]

COLOR_MAP = {"L_sap": "#E69F00", "L_tsap": "#56B4E9", "D_sap": "#009E73", "D_tsap": "#CC79A7"}

from pca_utils import R_GAS  # IUPAC 2019 value: 0.008314463 kJ/(mol·K)
FES_BINS = 20
FES_EPS = 1e-10
FES_SMOOTH_SIGMA = 1.0

DEFAULT_N_PCS = 3

def _make_pc_pairs(n_pcs):
    return [(i, j) for i in range(1, n_pcs + 1) for j in range(i + 1, n_pcs + 1)]

PC_PAIRS = _make_pc_pairs(DEFAULT_N_PCS)


def _build_groups(system_names):
    return {
        prefix: [s for s in system_names if s.startswith(prefix)]
        for prefix in ["me_sss", "me_rrr", "phe_sss", "phe_rrr"]
    }


def compute_projections(aligned, mean, eigvecs, n_pcs=3):
    centered = aligned - mean
    projections = centered @ eigvecs[:, :n_pcs]
    return projections


def save_projections_csv(projections, system_labels, frame_indices, outpath="projections_all.csv"):
    df = pd.DataFrame({
        "system_label": system_labels,
        "frame_index": frame_indices,
        "PC1": projections[:, 0],
        "PC2": projections[:, 1],
        "PC3": projections[:, 2],
    })
    df.to_csv(outpath, index=False)
    logging.info(f"Saved {outpath}: {len(df)} rows")
    return df


def make_scatter(proj_df, group_name, group_systems, pc_x, pc_y, outdir):
    outdir = Path(outdir)
    df = proj_df[proj_df["system_label"].isin(group_systems)].copy()
    fig, ax = plt.subplots(figsize=(8, 6))
    for stereo_key in STEREO_KEYS:
        stereo_systems = [s for s in group_systems if s.endswith(stereo_key)]
        if not stereo_systems:
            continue
        stereo_df = df[df["system_label"].isin(stereo_systems)]
        ax.scatter(
            stereo_df[f"PC{pc_x}"],
            stereo_df[f"PC{pc_y}"],
            c=COLOR_MAP[stereo_key],
            label=stereo_key,
            s=8,
            alpha=0.4,
            edgecolors="none",
        )
    ax.set_xlabel(f"PC{pc_x}")
    ax.set_ylabel(f"PC{pc_y}")
    legend_elements = [
        Line2D([0], [0], marker="o", color="w", markerfacecolor=COLOR_MAP[sk], markersize=8, label=sk)
        for sk in STEREO_KEYS
    ]
    ax.legend(handles=legend_elements, loc="best", framealpha=0.9, fontsize=9)
    ax.grid(True, linestyle="--", alpha=0.3)
    fig.tight_layout()
    fig.savefig(outdir / f"scatter_{group_name}_{pc_x}{pc_y}.png", dpi=150)
    plt.close(fig)
    df[["system_label", "frame_index", "PC1", "PC2", "PC3"]].to_csv(
        outdir / f"scatter_{group_name}_{pc_x}{pc_y}.csv", index=False
    )


def compute_free_energy(x, y, x_range, y_range, bins=FES_BINS, kT=2.4942,
                        smooth=False, sigma=FES_SMOOTH_SIGMA, eps=FES_EPS):
    """Compute free energy surface. Delegates to pca_utils.compute_fes for
    the core FES computation; retains smoothing as a local extension."""
    from pca_utils import compute_fes
    if smooth:
        # Smoothing path: still use local implementation
        # (compute_fes does not support Gaussian smoothing)
        xedges = np.linspace(x_range[0], x_range[1], bins + 1)
        yedges = np.linspace(y_range[0], y_range[1], bins + 1)
        hist_raw, _, _ = np.histogram2d(x, y, bins=(xedges, yedges))
        hist = gaussian_filter(hist_raw, sigma=sigma)
        total = hist.sum()
        prob = hist / total
        G = -kT * np.log(prob + eps)
        G[prob < 1.0 / (2 * total)] = np.nan  # Less than 0.5 expected counts
        G -= G[np.isfinite(G)].min()
        X, Y = np.meshgrid((xedges[:-1] + xedges[1:]) / 2,
                           (yedges[:-1] + yedges[1:]) / 2)
    else:
        # Standard path: delegate to shared compute_fes
        G, xedges, yedges = compute_fes(x, y, n_bins=bins, kT=kT,
                                         x_range=x_range, y_range=y_range)
        X, Y = np.meshgrid((xedges[:-1] + xedges[1:]) / 2,
                           (yedges[:-1] + yedges[1:]) / 2)
    return X, Y, G


def plot_fes(X, Y, G, group_name, pc_x, pc_y, vmin, vmax, outdir):
    outdir = Path(outdir)
    pair_label = f"{pc_x}{pc_y}"
    fig, ax = plt.subplots(figsize=(8, 6))
    cmap = plt.get_cmap("cubehelix")
    im = ax.pcolormesh(X, Y, G, cmap=cmap, shading="auto", vmin=vmin, vmax=vmax)
    cbar = fig.colorbar(im, ax=ax)
    ax.set_xlabel(f"PC{pc_x}", fontsize=12)
    ax.set_ylabel(f"PC{pc_y}", fontsize=12)
    cbar.set_label("Free Energy (kJ/mol)", fontsize=12)
    ax.yaxis.set_major_locator(MaxNLocator())
    fig.tight_layout()
    png_path = outdir / f"fes_{group_name}_{pair_label}.png"
    fig.savefig(png_path, dpi=150)
    plt.close(fig)
    csv_path = outdir / f"fes_{group_name}_{pair_label}.csv"
    df_csv = pd.DataFrame({"X": X.ravel(), "Y": Y.ravel(), "G": G.ravel()})
    df_csv.to_csv(csv_path, index=False)
    logging.info(f"Saved FES: {png_path.name} + {csv_path.name}")
    return png_path, csv_path


def generate_all_fes(projections, system_indices, systems_meta, outdir,
                     temperature=300.0, smooth=False, fes_bins=FES_BINS):
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    kT = R_GAS * temperature
    logging.info(f"FES: kT = {kT:.4f} kJ/mol (T = {temperature} K)")
    system_names = [s["name"] for s in systems_meta]
    GROUPS = _build_groups(system_names)
    pc_global_min = projections.min(axis=0)
    pc_global_max = projections.max(axis=0)
    fes_results = []
    pc_pairs = _make_pc_pairs(projections.shape[1])
    for group_name, (pc_x, pc_y) in itertools.product(GROUPS, pc_pairs):
        group_systems = GROUPS[group_name]
        mask = np.zeros(projections.shape[0], dtype=bool)
        for i, sys in enumerate(systems_meta):
            if sys["name"] in group_systems:
                start = int(system_indices[i])
                end = int(system_indices[i + 1])
                mask[start:end] = True
        group_proj = projections[mask]
        x = group_proj[:, pc_x - 1]
        y = group_proj[:, pc_y - 1]
        x_range = (pc_global_min[pc_x - 1], pc_global_max[pc_x - 1])
        y_range = (pc_global_min[pc_y - 1], pc_global_max[pc_y - 1])
        X, Y, G = compute_free_energy(x, y, x_range, y_range, bins=fes_bins, kT=kT, smooth=smooth)
        fes_results.append((group_name, pc_x, pc_y, X, Y, G))
        logging.info(f"FES computed: {group_name} PC{pc_x}vsPC{pc_y} "
                      f"({group_proj.shape[0]} frames, G range "
                      f"[{np.nanmin(G):.2f}, {np.nanmax(G):.2f}] kJ/mol)")
    pc_pair_ranges = {}
    for pc_x, pc_y in pc_pairs:
        pair_label = f"{pc_x}{pc_y}"
        group_fes = [G for gn, px, py, _, _, G in fes_results
                     if px == pc_x and py == pc_y]
        vmin = min(np.nanmin(G) for G in group_fes)
        vmax = max(np.nanmax(G) for G in group_fes)
        pc_pair_ranges[pair_label] = (vmin, vmax)
        logging.info(f"FES range for PC{pair_label}: [{vmin:.2f}, {vmax:.2f}] kJ/mol")
    for group_name, pc_x, pc_y, X, Y, G in fes_results:
        pair_label = f"{pc_x}{pc_y}"
        vmin, vmax = pc_pair_ranges[pair_label]
        plot_fes(X, Y, G, group_name, pc_x, pc_y, vmin, vmax, outdir)


def generate_all_plots(projections, system_indices, systems_meta, outdir):
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    system_names = [s["name"] for s in systems_meta]
    n_frames = projections.shape[0]
    system_labels = np.empty(n_frames, dtype=object)
    frame_indices = np.empty(n_frames, dtype=int)
    for i in range(len(system_names)):
        start = int(system_indices[i])
        end = int(system_indices[i + 1])
        system_labels[start:end] = system_names[i]
        frame_indices[start:end] = np.arange(end - start)
    df_proj = save_projections_csv(projections, system_labels, frame_indices, outpath=outdir / "projections_all.csv")
    GROUPS = _build_groups(system_names)
    pc_pairs = _make_pc_pairs(projections.shape[1])
    for group_name, (pc_x, pc_y) in itertools.product(GROUPS, pc_pairs):
        make_scatter(df_proj, group_name, GROUPS[group_name], pc_x, pc_y, outdir)


def parse_args():
    parser = argparse.ArgumentParser(description="Combined PCA of 16 MD systems")
    parser.add_argument("--input-dir",
                        default=os.environ.get('PCA_MD_DIR', '../com_md/'))
    parser.add_argument("--output-dir",
                        default=os.environ.get('PCA_OUTPUT_DIR', 'pca_analysis/'))
    parser.add_argument("--top-k-eigs", type=int, default=100,
                        help="Number of top eigenvectors to compute (0=all)")
    parser.add_argument("--align-chunk", type=int, default=5000,
                        help="Frames per chunk during alignment")
    parser.add_argument("--cov-chunk", type=int, default=5000,
                        help="Frames per chunk during covariance")
    parser.add_argument("--save-cov", action="store_true", default=False,
                        help="Save full covariance matrix (390 MB); skip to save disk")
    parser.add_argument("--force", action="store_true",
                        help="Overwrite existing output files")
    parser.add_argument("--temperature", type=float, default=300.0,
                        help="Temperature in K for free energy calculation (default: 300)")
    parser.add_argument("--smooth", action="store_true", default=False,
                        help="Apply Gaussian smoothing (sigma=1.0 bin) to FES histograms")
    parser.add_argument("--fes-bins", type=int, default=FES_BINS,
                        help=f"Bins per dimension for FES (default: {FES_BINS})")
    parser.add_argument("--n-align-iter", type=int, default=2,
                        help="Number of alignment iterations (1=single-pass, 2+=iterative, default: 2)")
    parser.add_argument("--n-modes", type=int, default=10,
                        help="Number of modes for eigenvalue CSV and NMD output (default: 10)")
    parser.add_argument("--n-pcs", type=int, default=3,
                        help="Number of PCs for projections, scatter plots, and FES (default: 3)")
    return parser.parse_args()


def main():
    args = parse_args()
    outdir = Path(args.output_dir)
    outdir.mkdir(parents=True, exist_ok=True)
    if not args.force:
        existing = list(outdir.glob("*.npy"))
        if existing:
            logging.error(
                f"Output directory contains {len(existing)} .npy files "
                f"(e.g. {existing[0].name}). Use --force to overwrite."
            )
            sys.exit(1)
    systems_meta = build_system_metadata()
    coords, system_indices, natoms = load_systems(args.input_dir, systems_meta)
    aligned, ref_coords, mean = align_to_reference(
        coords, natoms, chunk_size=args.align_chunk, max_iter=args.n_align_iter
    )
    del coords
    cov, total_variance = compute_covariance(aligned, mean, chunk_size=args.cov_chunk)
    eigvals, eigvecs = eigendecompose(cov, top_k=args.top_k_eigs)
    if not args.save_cov:
        del cov
    try:
        save_intermediates(outdir, systems_meta, system_indices, aligned, mean, eigvals, eigvecs, ref_coords, total_variance, cov if args.save_cov else None)
        write_eigenvalue_csv(eigvals, total_variance, n_modes=args.n_modes,
                             outpath=outdir / f"eigenvalue_top{args.n_modes}.csv")
        write_nmd_file(eigvals, eigvecs, ref_coords, n_modes=args.n_modes,
                       outpath=outdir / f"pc_modes_top{args.n_modes}.nmd")
    except IOError as e:
        logging.error(f"Failed to save results: {e}")
        raise
    projections = compute_projections(aligned, mean, eigvecs, n_pcs=args.n_pcs)
    del aligned
    generate_all_plots(projections, system_indices, systems_meta, outdir)
    generate_all_fes(projections, system_indices, systems_meta, outdir,
                     temperature=args.temperature, smooth=args.smooth,
                     fes_bins=args.fes_bins)
    logging.info("All tasks completed.")


if __name__ == "__main__":
    main()
