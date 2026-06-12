import warnings
# Suppress only specific noisy warnings, not all (M6 fix)
warnings.filterwarnings('ignore', category=FutureWarning)
warnings.filterwarnings('ignore', category=UserWarning, module='matplotlib')

import argparse
import os
import numpy as np
np.random.seed(42)  # Reproducibility for stochastic analyses
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from pathlib import Path
import MDAnalysis as mda


def build_residue_map(pdb_path):
    u = mda.Universe(pdb_path)
    bb = u.select_atoms('backbone')
    assert len(bb) == 2327, f"Expected 2327 backbone atoms, got {len(bb)}"
    unique_resindices = np.unique(bb.resindices)
    resindex_to_seq = {ri: i for i, ri in enumerate(unique_resindices)}
    residue_ids = np.array([resindex_to_seq[ri] for ri in bb.resindices], dtype=int)
    resids = bb.residues.resids
    resnames = bb.residues.resnames
    atoms_per_res = np.bincount(residue_ids)
    atom_indices = np.arange(len(bb))
    return residue_ids, resids, resnames, atoms_per_res, atom_indices


def compute_residue_contribution(eigvecs, residue_ids, atoms_per_res, n_pcs=5, eigvals=None):
    n_residues = len(atoms_per_res)
    vecs = eigvecs[:, :n_pcs].copy().reshape(-1, 3, n_pcs)
    atom_contribs = np.sum(vecs ** 2, axis=1)
    res_contrib = np.zeros((n_residues, n_pcs))
    for pc in range(n_pcs):
        np.add.at(res_contrib[:, pc], residue_ids, atom_contribs[:, pc])
    frac = res_contrib / res_contrib.sum(axis=0, keepdims=True)
    if eigvals is not None:
        absolute = frac * eigvals[:n_pcs]
    else:
        absolute = None
    return frac, absolute


def save_scree_data(eigvals, total_variance, outpath):
    fraction = eigvals / total_variance
    cumulative = np.cumsum(fraction)
    df = pd.DataFrame({
        'pc_number': np.arange(1, len(eigvals) + 1),
        'eigenvalue_angstrom2': eigvals,
        'fraction': fraction,
        'cumulative_fraction': cumulative
    })
    df.to_csv(outpath, index=False)


def plot_scree(eigvals, cumfrac, outpath, n_show=50):
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 4))
    pc_numbers = np.arange(1, n_show + 1)

    ax1.plot(pc_numbers, eigvals[:n_show], 'o-', markersize=3, linewidth=1)
    ax1.set_xlabel('PC number')
    ax1.set_ylabel('Eigenvalue (Å²)')
    ax1.axvline(x=5, color='gray', linestyle='--', linewidth=0.8)
    ax1.annotate('Elbow', xy=(5, eigvals[4]), xytext=(8, eigvals[4] * 0.8),
                 arrowprops=dict(arrowstyle='->', color='gray'), fontsize=8, color='gray')

    ax2.plot(pc_numbers, cumfrac[:n_show], 'o-', markersize=3, linewidth=1)
    ax2.set_xlabel('PC number')
    ax2.set_ylabel('Cumulative variance fraction')
    for val, label in [(0.5, '50%'), (0.75, '75%'), (0.9, '90%')]:
        ax2.axhline(y=val, color='gray', linestyle='--', linewidth=0.8)
        ax2.text(n_show * 0.98, val + 0.02, label, ha='right', fontsize=8, color='gray')

    crossing_50 = int(np.searchsorted(cumfrac, 0.5)) + 1
    crossing_75 = int(np.searchsorted(cumfrac, 0.75)) + 1
    crossing_90 = int(np.searchsorted(cumfrac, 0.9)) + 1
    if crossing_50 <= n_show:
        ax2.axvline(x=crossing_50, color='gray', linestyle=':', linewidth=0.6)
        ax2.annotate(f'{crossing_50} PCs: {cumfrac[crossing_50-1]*100:.1f}%',
                     xy=(crossing_50, 0.5), fontsize=7, color='gray')
    if crossing_75 <= n_show:
        ax2.axvline(x=crossing_75, color='gray', linestyle=':', linewidth=0.6)
        ax2.annotate(f'{crossing_75} PCs: {cumfrac[crossing_75-1]*100:.1f}%',
                     xy=(crossing_75, 0.75), fontsize=7, color='gray')
    if crossing_90 > n_show:
        ax2.annotate(f'{crossing_90} PCs: {cumfrac[crossing_90-1]*100:.1f}%',
                     xy=(n_show * 0.95, 0.9), fontsize=7, color='gray', ha='right')

    fig.tight_layout()
    fig.savefig(outpath, dpi=150)
    plt.close(fig)


def plot_residue_contribution(res_df, outpath):
    fig, ax = plt.subplots(figsize=(12, 4))
    colors = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd']
    frac_cols = [c for c in res_df.columns if c.endswith('_frac') and c.startswith('PC')]
    for i, col in enumerate(frac_cols):
        ax.plot(res_df['resid'], res_df[col], color=colors[i % len(colors)],
                label=col.replace('_frac', ''), linewidth=1)
    baseline = 1.0 / len(res_df)
    ax.axhline(y=baseline, color='gray', linestyle='--', linewidth=0.8, label='Baseline (1/N)')
    top10_idx = res_df['PC1_frac'].nlargest(10).index
    for idx in top10_idx:
        row = res_df.loc[idx]
        ax.annotate(f"{row['resname']}{int(row['resid'])}",
                    xy=(row['resid'], row['PC1_frac']),
                    fontsize=6, rotation=90, va='bottom', ha='center')
    ax.set_xlabel('Residue ID')
    ax.set_ylabel('Fractional contribution')
    ax.legend(loc='upper right', fontsize=8)
    fig.tight_layout()
    fig.savefig(outpath, dpi=150)
    plt.close(fig)


def plot_variance_decomposition(eigvals, cumfrac, res_df, outpath, n_show=50):
    fig = plt.figure(figsize=(8.27, 11.69))
    gs = gridspec.GridSpec(2, 2, figure=fig, hspace=0.35, wspace=0.30)
    ax_scree = fig.add_subplot(gs[0, 0])
    ax_cumul = fig.add_subplot(gs[0, 1])
    ax_full = fig.add_subplot(gs[1, 0])
    ax_zoom = fig.add_subplot(gs[1, 1])

    pc_numbers = np.arange(1, n_show + 1)

    ax_scree.plot(pc_numbers, eigvals[:n_show], 'o-', markersize=3, linewidth=1)
    ax_scree.set_xlabel('PC number')
    ax_scree.set_ylabel('Eigenvalue (Å²)')
    ax_scree.axvline(x=5, color='gray', linestyle='--', linewidth=0.8)
    ax_scree.text(0.02, 0.98, 'A', transform=ax_scree.transAxes, va='top',
                  fontweight='bold', fontsize=11)

    ax_cumul.plot(pc_numbers, cumfrac[:n_show], 'o-', markersize=3, linewidth=1)
    ax_cumul.set_xlabel('PC number')
    ax_cumul.set_ylabel('Cumulative variance fraction')
    for val, label in [(0.5, '50%'), (0.75, '75%'), (0.9, '90%')]:
        ax_cumul.axhline(y=val, color='gray', linestyle='--', linewidth=0.8)
        ax_cumul.text(n_show * 0.98, val + 0.02, label, ha='right', fontsize=8, color='gray')
    ax_cumul.text(0.02, 0.98, 'B', transform=ax_cumul.transAxes, va='top',
                  fontweight='bold', fontsize=11)

    colors = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd']
    frac_cols = [c for c in res_df.columns if c.endswith('_frac') and c.startswith('PC')]
    baseline = 1.0 / len(res_df)

    for i, col in enumerate(frac_cols):
        ax_full.plot(res_df['resid'], res_df[col], color=colors[i % len(colors)],
                     label=col.replace('_frac', ''), linewidth=1)
    ax_full.axhline(y=baseline, color='gray', linestyle='--', linewidth=0.8)
    ax_full.set_xlabel('Residue ID')
    ax_full.set_ylabel('Fractional contribution')
    ax_full.legend(loc='upper right', fontsize=7)
    ax_full.text(0.02, 0.98, 'C', transform=ax_full.transAxes, va='top',
                 fontweight='bold', fontsize=11)

    mask = (res_df['resid'] >= 370) & (res_df['resid'] <= 500)
    zoom_df = res_df[mask].copy()
    for i, col in enumerate(frac_cols):
        lw = 2.0 if 'PC3' in col else 1.2
        ax_zoom.plot(zoom_df['resid'], zoom_df[col], color=colors[i % len(colors)],
                     label=col.replace('_frac', ''), linewidth=lw)
    ax_zoom.axhline(y=baseline, color='gray', linestyle='--', linewidth=0.8)
    key_residues = [379, 382, 383, 384, 386, 387, 390, 391, 414, 485, 489, 490]
    for kr in key_residues:
        ax_zoom.axvline(x=kr, color='lightgray', alpha=0.5, linewidth=0.8)
    pc3_col = 'PC3_frac' if 'PC3_frac' in zoom_df.columns else None
    if pc3_col:
        for kr in key_residues:
            kr_rows = zoom_df[zoom_df['resid'] == kr]
            if not kr_rows.empty and kr_rows[pc3_col].values[0] > baseline:
                ax_zoom.annotate(
                    f"{kr_rows['resname'].values[0]}{int(kr)}",
                    xy=(kr, kr_rows[pc3_col].values[0]),
                    fontsize=6, rotation=90, va='bottom', ha='center', color='#2ca02c')
    # Use the standardized binding site range from pca_utils (S3 fix)
    from pca_utils import BINDING_SITE_RESID_RANGE
    bs_mask = (res_df['resid'] >= BINDING_SITE_RESID_RANGE[0]) & (res_df['resid'] <= BINDING_SITE_RESID_RANGE[1])
    binding_frac = res_df.loc[bs_mask, 'PC1_frac'].sum()
    expected_frac = bs_mask.sum() / len(res_df)
    ax_zoom.text(0.98, 0.02,
                 f"Binding site: {binding_frac*100:.1f}% of PC1\n({expected_frac*100:.1f}% expected if uniform)\nNote: y-axis scale differs from C",
                 transform=ax_zoom.transAxes, ha='right', va='bottom', fontsize=7,
                 bbox=dict(boxstyle='round,pad=0.3', facecolor='wheat', alpha=0.5))
    ax_zoom.set_xlabel('Residue ID')
    ax_zoom.set_ylabel('Fractional contribution')
    ax_zoom.text(0.02, 0.98, 'D', transform=ax_zoom.transAxes, va='top',
                 fontweight='bold', fontsize=11)

    fig.tight_layout()
    fig.savefig(outpath, dpi=300)
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--input-dir', default='pca_analysis/')
    parser.add_argument('--pdb-path',
                        default=os.environ.get('PCA_PDB_PATH', '../com_md/me_sssL_sap/fp/v1.pdb'))
    parser.add_argument('--n-pcs', type=int, default=5)
    parser.add_argument('--n-show', type=int, default=50)
    args = parser.parse_args()

    outdir = Path(args.input_dir)

    eigvecs = np.load(outdir / "eigenvectors.npy", mmap_mode='r')
    eigvals = np.load(outdir / "eigenvalues.npy")
    total_var = np.load(outdir / "total_variance.npy")

    residue_ids, resids, resnames, atoms_per_res, atom_indices = build_residue_map(args.pdb_path)

    np.savez(outdir / "residue_map.npz",
             atom_indices=atom_indices,
             residue_ids=residue_ids,
             residue_names=np.array(resnames, dtype='U4'),
             n_residues=np.array(len(resids)),
             n_atoms_per_res=atoms_per_res)
    pd.DataFrame({'resid': resids, 'resname': resnames}).to_csv(
        outdir / "residue_info.csv", index=False)

    frac, absolute = compute_residue_contribution(
        eigvecs, residue_ids, atoms_per_res, n_pcs=args.n_pcs, eigvals=eigvals)
    del eigvecs

    res_df = pd.DataFrame()
    res_df['resid'] = resids
    res_df['resname'] = resnames
    for i in range(args.n_pcs):
        res_df[f'PC{i+1}_frac'] = frac[:, i]
        res_df[f'PC{i+1}_abs_angstrom2'] = absolute[:, i]
    res_df.to_csv(outdir / "residue_contribution.csv", index=False)

    save_scree_data(eigvals, total_var, outdir / "scree_data.csv")

    fraction = eigvals / total_var
    cumfrac = np.cumsum(fraction)

    plot_scree(eigvals, cumfrac, outdir / "scree_plot.png", n_show=args.n_show)
    plot_residue_contribution(res_df, outdir / "residue_contribution.png")
    plot_variance_decomposition(eigvals, cumfrac, res_df,
                                outdir / "figure_variance_decomposition.png", n_show=args.n_show)


if __name__ == '__main__':
    main()
