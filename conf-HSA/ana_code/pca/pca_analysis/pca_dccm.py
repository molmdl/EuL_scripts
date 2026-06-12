import json
import logging
import os
import sys
from pathlib import Path

import numpy as np
np.random.seed(42)  # Reproducibility for stochastic analyses
import MDAnalysis as mda
from pca_utils import compute_dccm, DCCM_BINDING_SITE_START, DCCM_BINDING_SITE_END, DCCM_BINDING_SITE_RESID_OFFSET
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle

logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")

N_BACKBONE = 2327
N_CA = 582
N_COORDS_PER_CA = 3
CA_FLAT_DIM = N_CA * N_COORDS_PER_CA
# Binding site constants imported from pca_utils; aliased for backward compat
BINDING_SITE_START = DCCM_BINDING_SITE_START   # 376
BINDING_SITE_END = DCCM_BINDING_SITE_END       # 487
RESID_OFFSET = DCCM_BINDING_SITE_RESID_OFFSET  # 3
GROUP_PREFIXES = ["me_sss", "me_rrr", "phe_sss", "phe_rrr"]


def build_ca_flat_indices():
    ca_atom_indices = np.arange(1, N_BACKBONE, 4)
    ca_flat_indices = np.sort(np.concatenate(
        [3 * ca_atom_indices + k for k in range(N_COORDS_PER_CA)]
    ))
    return ca_flat_indices


def build_groups(systems_meta):
    groups = {}
    for prefix in GROUP_PREFIXES:
        groups[prefix] = [i for i, s in enumerate(systems_meta) if s["name"].startswith(prefix)]
    return groups


def get_frame_ranges(system_indices, sys_idx_list):
    return [(int(system_indices[i]), int(system_indices[i + 1])) for i in sys_idx_list]


def extract_ca_coords(aligned_mmap, frame_ranges, ca_flat_indices):
    total_frames = sum(end - start for start, end in frame_ranges)
    ca_coords = np.empty((total_frames, len(ca_flat_indices)), dtype=np.float32)
    offset = 0
    for start, end in frame_ranges:
        n = end - start
        ca_coords[offset:offset + n] = aligned_mmap[start:end, ca_flat_indices]
        offset += n
    return ca_coords


# compute_dccm() is now imported from pca_utils.py (was duplicated here
# and in pca_sap_tsap_figures.py). The pca_utils version uses the same
# algorithm but accepts n_ca_atoms and n_coords_per_ca as parameters.


def compute_dccm_from_eigvecs(eigvecs_path, eigvals_path, ca_flat_indices, n_atoms=N_CA):
    eigvecs = np.load(eigvecs_path, mmap_mode="r")
    eigvals = np.load(eigvals_path)
    n_eigs = eigvals.shape[0]
    if n_eigs < 3 * n_atoms:
        raise ValueError(
            f"Cannot compute reliable global DCCM from partial eigendecomposition "
            f"({n_eigs} eigenvectors < 3*{n_atoms} = {3*n_atoms}). "
            f"Use --top-k-eigs 0 for full decomposition, or compute DCCM "
            f"directly from aligned coordinates instead."
        )
    V_ca = eigvecs[ca_flat_indices, :].astype(np.float64)
    cov_ca = (V_ca * eigvals[np.newaxis, :]) @ V_ca.T
    cov_4d = cov_ca.reshape(n_atoms, N_COORDS_PER_CA, n_atoms, N_COORDS_PER_CA)
    numerator = np.einsum("imjm->ij", cov_4d)
    var_i = numerator.diagonal().copy()
    denom = np.sqrt(np.outer(var_i, var_i))
    denom[denom == 0] = 1.0
    return numerator / denom


def save_dccm_csv(dccm, outpath, resid_offset=RESID_OFFSET):
    resids = np.arange(dccm.shape[0]) + resid_offset
    header = ",".join(str(r) for r in resids)
    np.savetxt(outpath, dccm, delimiter=",", header=header, comments="", fmt="%.6f")


def plot_dccm_overview(dccm_dict, outpath, binding_start=BINDING_SITE_START,
                       binding_end=BINDING_SITE_END, resid_offset=RESID_OFFSET):
    fig, axes = plt.subplots(2, 2, figsize=(11.69, 8.27))
    axes = axes.ravel()
    n = N_CA
    tick_positions = np.linspace(0, n - 1, 6, dtype=int)
    tick_labels = [str(p + resid_offset) for p in tick_positions]
    im = None
    for idx, group_name in enumerate(GROUP_PREFIXES):
        ax = axes[idx]
        dccm = dccm_dict[group_name]
        im = ax.imshow(dccm, cmap="bwr", vmin=-1.0, vmax=1.0,
                       interpolation="nearest", origin="lower", aspect="equal")
        rect = Rectangle((binding_start, binding_start),
                          binding_end - binding_start,
                          binding_end - binding_start,
                          linewidth=1.5, edgecolor="black", facecolor="none",
                          linestyle="--")
        ax.add_patch(rect)
        ax.set_title(f"{chr(65 + idx)} DCCM {group_name}", fontsize=11)
        ax.set_xticks(tick_positions)
        ax.set_xticklabels(tick_labels, fontsize=8)
        ax.set_yticks(tick_positions)
        ax.set_yticklabels(tick_labels, fontsize=8)
        ax.set_xlabel("Residue ID", fontsize=9)
        ax.set_ylabel("Residue ID", fontsize=9)
    fig.subplots_adjust(right=0.88, wspace=0.35, hspace=0.35)
    cbar_ax = fig.add_axes([0.90, 0.15, 0.02, 0.7])
    fig.colorbar(im, cax=cbar_ax, label="Cross-correlation")
    fig.savefig(outpath, dpi=300, bbox_inches="tight")
    plt.close(fig)
    logging.info(f"Saved {outpath}")


def plot_dccm_differences(dccm_dict, dccm_global, outpath,
                          binding_start=BINDING_SITE_START,
                          binding_end=BINDING_SITE_END,
                          resid_offset=RESID_OFFSET):
    diff_pairs = [
        ("me_sss", "me_rrr", "Stereo effect (me): sss \u2212 rrr"),
        ("phe_sss", "phe_rrr", "Stereo effect (phe): sss \u2212 rrr"),
        ("me_sss", "phe_sss", "Chirality effect (sss): me \u2212 phe"),
        ("me_rrr", "phe_rrr", "Chirality effect (rrr): me \u2212 phe"),
    ]
    fig, axes = plt.subplots(2, 2, figsize=(11.69, 8.27))
    axes = axes.ravel()
    n = N_CA
    tick_positions = np.linspace(0, n - 1, 6, dtype=int)
    tick_labels = [str(p + resid_offset) for p in tick_positions]
    im = None
    for idx, (g1, g2, title) in enumerate(diff_pairs):
        ax = axes[idx]
        diff = dccm_dict[g1] - dccm_dict[g2]
        im = ax.imshow(diff, cmap="bwr", vmin=-0.5, vmax=0.5,
                       interpolation="nearest", origin="lower", aspect="equal")
        rect = Rectangle((binding_start, binding_start),
                          binding_end - binding_start,
                          binding_end - binding_start,
                          linewidth=1.5, edgecolor="black", facecolor="none",
                          linestyle="--")
        ax.add_patch(rect)
        ax.set_title(f"{chr(65 + idx)} {title}", fontsize=10)
        ax.set_xticks(tick_positions)
        ax.set_xticklabels(tick_labels, fontsize=8)
        ax.set_yticks(tick_positions)
        ax.set_yticklabels(tick_labels, fontsize=8)
        ax.set_xlabel("Residue ID", fontsize=9)
        ax.set_ylabel("Residue ID", fontsize=9)
    fig.subplots_adjust(right=0.88, wspace=0.35, hspace=0.35)
    cbar_ax = fig.add_axes([0.90, 0.15, 0.02, 0.7])
    fig.colorbar(im, cax=cbar_ax, label="\u0394Cross-correlation")
    fig.savefig(outpath, dpi=300, bbox_inches="tight")
    plt.close(fig)
    logging.info(f"Saved {outpath}")


def compute_correlation_stats(dccm_dict, dccm_global, outpath):
    rows = []
    group_names = list(dccm_dict.keys())
    for i, g1 in enumerate(group_names):
        for j, g2 in enumerate(group_names):
            if i >= j:
                continue
            v1 = dccm_dict[g1][np.triu_indices(N_CA, k=1)]
            v2 = dccm_dict[g2][np.triu_indices(N_CA, k=1)]
            r = np.corrcoef(v1, v2)[0, 1]
            rmse = np.sqrt(np.mean((v1 - v2) ** 2))
            rows.append({"pair": f"{g1} vs {g2}", "pearson_r": round(r, 6),
                          "rmse": round(rmse, 6)})
    for g in group_names:
        v1 = dccm_dict[g][np.triu_indices(N_CA, k=1)]
        v2 = dccm_global[np.triu_indices(N_CA, k=1)]
        r = np.corrcoef(v1, v2)[0, 1]
        rmse = np.sqrt(np.mean((v1 - v2) ** 2))
        rows.append({"pair": f"{g} vs global", "pearson_r": round(r, 6),
                      "rmse": round(rmse, 6)})
    with open(outpath, "w") as f:
        f.write("pair,pearson_r,rmse\n")
        for row in rows:
            f.write(f"{row['pair']},{row['pearson_r']},{row['rmse']}\n")
    logging.info(f"Saved {outpath}: {len(rows)} pairs")
    return rows


def parse_args():
    import argparse
    parser = argparse.ArgumentParser(description="Compute DCCM from PCA results")
    parser.add_argument("--force", action="store_true",
                        help="Overwrite existing DCCM output files")
    parser.add_argument("--output-dir", type=str, default=None,
                        help="Output directory (default: same as script)")
    return parser.parse_args()


def main():
    args = parse_args()
    data_dir = Path(args.output_dir) if args.output_dir else Path(__file__).parent

    # Prevent silent overwrite of DCCM data (RQ-REP1)
    stats_file = data_dir / "dccm_correlation_stats.csv"
    if stats_file.exists() and not args.force:
        logging.error(
            f"DCCM output files exist in {data_dir} "
            f"(e.g. {stats_file.name}). Use --force to overwrite."
        )
        sys.exit(1)

    ca_flat_indices = build_ca_flat_indices()
    logging.info(f"C\u03b1 flat coordinate indices: {len(ca_flat_indices)}")

    # Validate CA extraction: positional method vs name-based method
    md_dir = os.environ.get('PCA_MD_DIR', '../com_md/')
    pdb_path = Path(md_dir) / "me_rrrD_sap/fp/v1.pdb"
    if not pdb_path.exists():
        logging.warning(f"PDB for CA validation not found: {pdb_path}. Skipping name-based CA validation.")
    else:
        u = mda.Universe(pdb_path)
        bb = u.select_atoms('backbone')
        ca_by_name = np.where(bb.names == 'CA')[0]
        ca_by_arange = np.arange(1, N_BACKBONE, 4)
        if not np.array_equal(ca_by_name, ca_by_arange):
            logging.error("CA extraction mismatch — positional method invalid for this protein. Use name-based extraction instead.")
            sys.exit(1)
        else:
            logging.info("CA extraction validated: positional and name-based methods agree (%d CA atoms)", len(ca_by_name))

    aligned_mmap = np.load(data_dir / "aligned_coords.npy", mmap_mode="r")
    system_indices = np.load(data_dir / "system_indices.npy")
    with open(data_dir / "system_metadata.json") as f:
        systems_meta = json.load(f)

    groups = build_groups(systems_meta)
    logging.info(f"Groups: { {k: len(v) for k, v in groups.items()} }")

    dccm_dict = {}
    for group_name in GROUP_PREFIXES:
        sys_idx_list = groups[group_name]
        frame_ranges = get_frame_ranges(system_indices, sys_idx_list)
        total_frames = sum(end - start for start, end in frame_ranges)
        logging.info(f"Computing DCCM for {group_name}: {total_frames} frames "
                      f"from {len(sys_idx_list)} systems")
        ca_coords = extract_ca_coords(aligned_mmap, frame_ranges, ca_flat_indices)
        dccm = compute_dccm(ca_coords)
        dccm_dict[group_name] = dccm
        np.save(data_dir / f"dccm_{group_name}.npy", dccm)
        save_dccm_csv(dccm, data_dir / f"dccm_{group_name}.csv")
        off_diag = dccm[np.triu_indices(N_CA, k=1)]
        logging.info(f"  {group_name}: range=[{dccm.min():.3f}, {dccm.max():.3f}], "
                      f"mean|corr|={np.mean(np.abs(off_diag)):.4f}")

    logging.info("Computing global DCCM from eigenvectors")
    dccm_global = compute_dccm_from_eigvecs(
        data_dir / "eigenvectors.npy", data_dir / "eigenvalues.npy", ca_flat_indices
    )
    np.save(data_dir / "dccm_global.npy", dccm_global)
    save_dccm_csv(dccm_global, data_dir / "dccm_global.csv")
    off_diag_g = dccm_global[np.triu_indices(N_CA, k=1)]
    logging.info(f"  global: range=[{dccm_global.min():.3f}, {dccm_global.max():.3f}], "
                  f"mean|corr|={np.mean(np.abs(off_diag_g)):.4f}")

    plot_dccm_overview(dccm_dict, data_dir / "figure_dccm_overview.png")
    plot_dccm_differences(dccm_dict, dccm_global, data_dir / "figure_dccm_differences.png")

    compute_correlation_stats(dccm_dict, dccm_global, data_dir / "dccm_correlation_stats.csv")

    logging.info("All DCCM tasks completed.")


if __name__ == "__main__":
    main()
