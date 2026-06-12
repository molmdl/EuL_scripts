#!/usr/bin/env python3
"""
Block averaging for DCCM convergence analysis.

Computes DCCM from non-overlapping trajectory blocks and quantifies
inter-block variability. Targets phe_sssD_tsap (non-convergent) but
runs on all 4 key phe systems for comparison.

Usage:
    python pca_block_average_dccm.py
    python pca_block_average_dccm.py --systems phe_sssD_tsap
    python pca_block_average_dccm.py --block-sizes 500
    python pca_block_average_dccm.py --input-dir pca_analysis/ --output-dir pca_analysis/block_averaging_data/
"""

import sys
import os
import json
import argparse
import logging
import time

import numpy as np
np.random.seed(42)  # Reproducibility
import pandas as pd
from scipy import stats

import matplotlib
matplotlib.use('Agg')  # Headless backend
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle

# Import from existing pipeline
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from pca_utils import compute_dccm, DCCM_BINDING_SITE_START, DCCM_BINDING_SITE_END
from pca_dccm import build_ca_flat_indices

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------
N_BACKBONE = 2327
N_CA = 582
N_COORDS_PER_CA = 3
CA_FLAT_DIM = N_CA * N_COORDS_PER_CA
BINDING_SITE_START = DCCM_BINDING_SITE_START   # 376
BINDING_SITE_END = DCCM_BINDING_SITE_END       # 487

logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")


# ---------------------------------------------------------------------------
# 1. Core computation functions
# ---------------------------------------------------------------------------

def compute_block_dccms(ca_coords, block_size, n_ca_atoms=N_CA):
    """Compute DCCM for each non-overlapping block of frames.

    Parameters
    ----------
    ca_coords : ndarray, shape (n_frames, 1746)
        CA coordinates for the entire system.
    block_size : int
        Number of frames per block.
    n_ca_atoms : int
        Number of CA atoms (default 582).

    Returns
    -------
    dccm_blocks : list of ndarray, each shape (582, 582)
        DCCM computed from each block.
    n_blocks : int
        Number of complete blocks used.
    """
    n_frames = ca_coords.shape[0]
    n_blocks = n_frames // block_size
    n_used = n_blocks * block_size
    n_discarded = n_frames - n_used
    if n_discarded > 0:
        logging.info(f"    Discarding {n_discarded} trailing frames "
                      f"({n_frames} mod {block_size} = {n_discarded})")

    dccm_blocks = []
    for k in range(n_blocks):
        start = k * block_size
        end = start + block_size
        block_coords = ca_coords[start:end].astype(np.float64)
        dccm_k = compute_dccm(block_coords, n_ca_atoms=n_ca_atoms)
        dccm_blocks.append(dccm_k)
        logging.info(f"    Block {k+1}/{n_blocks}: frames {start}-{end-1}")

    return dccm_blocks, n_blocks


def block_average_dccm(dccm_blocks):
    """Element-wise mean and std of block DCCMs.

    Parameters
    ----------
    dccm_blocks : list of ndarray, each shape (582, 582)
        DCCM computed from each block.

    Returns
    -------
    dccm_avg : ndarray, shape (582, 582)
        Element-wise mean across blocks.
    dccm_std : ndarray, shape (582, 582)
        Element-wise sample standard deviation across blocks.
    """
    stacked = np.stack(dccm_blocks, axis=0)  # (K, 582, 582)
    dccm_avg = stacked.mean(axis=0)
    dccm_std = stacked.std(axis=0, ddof=1)   # Sample std
    return dccm_avg, dccm_std


def compare_dccm_pair(dccm1, dccm2, n_ca=N_CA):
    """Compare two DCCM matrices via Pearson r and RMSE on upper triangle.

    Parameters
    ----------
    dccm1, dccm2 : ndarray, shape (n_ca, n_ca)
        Two DCCM matrices to compare.
    n_ca : int
        Number of CA atoms.

    Returns
    -------
    dict with keys 'pearson_r', 'rmse'
    """
    triu_idx = np.triu_indices(n_ca, k=1)
    v1 = dccm1[triu_idx]
    v2 = dccm2[triu_idx]
    r, p = stats.pearsonr(v1, v2)
    rmse = np.sqrt(np.mean((v1 - v2) ** 2))
    return {'pearson_r': r, 'rmse': rmse}


def interblock_correlation_matrix(dccm_blocks, n_ca=N_CA):
    """K x K Pearson r matrix between all pairs of block DCCMs.

    Parameters
    ----------
    dccm_blocks : list of ndarray, each shape (n_ca, n_ca)
    n_ca : int

    Returns
    -------
    r_matrix : ndarray, shape (K, K)
        Pairwise Pearson r between block DCCM upper triangles.
    summary : dict
        Keys: mean_r_off_diag, min_r_off_diag, max_r_off_diag
    """
    K = len(dccm_blocks)
    triu_idx = np.triu_indices(n_ca, k=1)
    vectors = [dccm_k[triu_idx] for dccm_k in dccm_blocks]

    r_matrix = np.ones((K, K))
    for i in range(K):
        for j in range(i + 1, K):
            r, _ = stats.pearsonr(vectors[i], vectors[j])
            r_matrix[i, j] = r
            r_matrix[j, i] = r

    # Off-diagonal stats
    off_diag_mask = ~np.eye(K, dtype=bool)
    off_diag = r_matrix[off_diag_mask]
    summary = {
        'mean_r_off_diag': off_diag.mean(),
        'min_r_off_diag': off_diag.min(),
        'max_r_off_diag': off_diag.max(),
    }
    return r_matrix, summary


def flyvbjerg_petersen_analysis(ca_coords, block_sizes, n_ca=N_CA):
    """Compute mean DCCM element variance at each block size.

    If variance plateaus (stops increasing) beyond some block size,
    that block size is sufficient for independent sampling.

    Parameters
    ----------
    ca_coords : ndarray, shape (n_frames, 1746)
        CA coordinates for the entire system.
    block_sizes : list of int
        Block sizes to test.
    n_ca : int

    Returns
    -------
    results : list of dict
        Each dict has: block_size, n_blocks, mean_var, median_var,
        max_var, p90_var, binding_site_mean_var
    """
    triu_idx = np.triu_indices(n_ca, k=1)
    results = []

    for B in block_sizes:
        logging.info(f"    Flyvbjerg-Petersen: block_size={B}")
        dccm_blocks, n_blocks = compute_block_dccms(ca_coords, B)
        stacked = np.stack(dccm_blocks, axis=0)  # (K, 582, 582)
        var_matrix = stacked.var(axis=0, ddof=1)  # Element-wise variance

        # Summary statistics on upper triangle
        triu_var = var_matrix[triu_idx]
        binding_site_var = var_matrix[BINDING_SITE_START:BINDING_SITE_END,
                                      BINDING_SITE_START:BINDING_SITE_END]

        results.append({
            'block_size': B,
            'n_blocks': n_blocks,
            'mean_var': triu_var.mean(),
            'median_var': np.median(triu_var),
            'max_var': triu_var.max(),
            'p90_var': np.percentile(triu_var, 90),
            'binding_site_mean_var': binding_site_var.mean(),
        })

        # Free memory
        del dccm_blocks, stacked, var_matrix

    return results


# ---------------------------------------------------------------------------
# 2. Plotting functions
# ---------------------------------------------------------------------------

def plot_block_dccms_overview(dccm_blocks, system_name, block_size, outpath):
    """Multi-panel figure of per-block DCCM heatmaps.

    Parameters
    ----------
    dccm_blocks : list of ndarray, each shape (582, 582)
    system_name : str
    block_size : int
    outpath : str
        Output PNG path.
    """
    K = len(dccm_blocks)
    ncols = min(4, K)
    nrows = (K + ncols - 1) // ncols

    fig, axes = plt.subplots(nrows, ncols, figsize=(4 * ncols, 3.5 * nrows))
    if K == 1:
        axes = np.array([[axes]])
    elif nrows == 1:
        axes = axes.reshape(1, -1)
    axes = axes.ravel()

    tick_positions = np.linspace(0, N_CA - 1, 6, dtype=int)
    tick_labels = [str(p + 3) for p in tick_positions]

    for k in range(K):
        ax = axes[k]
        im = ax.imshow(dccm_blocks[k], cmap='bwr', vmin=-1.0, vmax=1.0,
                        interpolation='nearest', origin='lower', aspect='equal')
        rect = Rectangle((BINDING_SITE_START, BINDING_SITE_START),
                         BINDING_SITE_END - BINDING_SITE_START,
                         BINDING_SITE_END - BINDING_SITE_START,
                         linewidth=1.5, edgecolor='black', facecolor='none',
                         linestyle='--')
        ax.add_patch(rect)
        ax.set_title(f'Block {k+1}', fontsize=9)
        ax.set_xticks(tick_positions)
        ax.set_xticklabels(tick_labels, fontsize=6)
        ax.set_yticks(tick_positions)
        ax.set_yticklabels(tick_labels, fontsize=6)

    # Hide unused subplots
    for k in range(K, len(axes)):
        axes[k].set_visible(False)

    fig.suptitle(f'Per-block DCCMs: {system_name} (B={block_size})', fontsize=12)
    fig.subplots_adjust(right=0.88, wspace=0.35, hspace=0.4)
    cbar_ax = fig.add_axes([0.90, 0.15, 0.02, 0.7])
    fig.colorbar(im, cax=cbar_ax, label='Cross-correlation')
    fig.savefig(outpath, dpi=300, bbox_inches='tight')
    plt.close(fig)
    logging.info(f"  Saved {outpath}")


def plot_block_avg_vs_full(dccm_avg, dccm_full, dccm_std, system_name, outpath):
    """Three-panel figure: block-averaged DCCM, full-sample DCCM, inter-block std.

    Parameters
    ----------
    dccm_avg : ndarray, shape (582, 582)
        Block-averaged DCCM.
    dccm_full : ndarray, shape (582, 582)
        Full-sample DCCM.
    dccm_std : ndarray, shape (582, 582)
        Inter-block standard deviation.
    system_name : str
    outpath : str
    """
    fig, axes = plt.subplots(1, 3, figsize=(16, 5))

    tick_positions = np.linspace(0, N_CA - 1, 6, dtype=int)
    tick_labels = [str(p + 3) for p in tick_positions]

    # Panel (a): Block-averaged DCCM
    im0 = axes[0].imshow(dccm_avg, cmap='bwr', vmin=-1.0, vmax=1.0,
                          interpolation='nearest', origin='lower', aspect='equal')
    rect0 = Rectangle((BINDING_SITE_START, BINDING_SITE_START),
                       BINDING_SITE_END - BINDING_SITE_START,
                       BINDING_SITE_END - BINDING_SITE_START,
                       linewidth=1.5, edgecolor='black', facecolor='none',
                       linestyle='--')
    axes[0].add_patch(rect0)
    axes[0].set_title('A Block-averaged DCCM', fontsize=11)
    axes[0].set_xticks(tick_positions)
    axes[0].set_xticklabels(tick_labels, fontsize=7)
    axes[0].set_yticks(tick_positions)
    axes[0].set_yticklabels(tick_labels, fontsize=7)
    fig.colorbar(im0, ax=axes[0], fraction=0.046, pad=0.04)

    # Panel (b): Full-sample DCCM
    im1 = axes[1].imshow(dccm_full, cmap='bwr', vmin=-1.0, vmax=1.0,
                          interpolation='nearest', origin='lower', aspect='equal')
    rect1 = Rectangle((BINDING_SITE_START, BINDING_SITE_START),
                       BINDING_SITE_END - BINDING_SITE_START,
                       BINDING_SITE_END - BINDING_SITE_START,
                       linewidth=1.5, edgecolor='black', facecolor='none',
                       linestyle='--')
    axes[1].add_patch(rect1)
    axes[1].set_title('B Full-sample DCCM', fontsize=11)
    axes[1].set_xticks(tick_positions)
    axes[1].set_xticklabels(tick_labels, fontsize=7)
    axes[1].set_yticks(tick_positions)
    axes[1].set_yticklabels(tick_labels, fontsize=7)
    fig.colorbar(im1, ax=axes[1], fraction=0.046, pad=0.04)

    # Panel (c): Inter-block std
    vmax_std = max(np.percentile(dccm_std, 95), 0.05)
    im2 = axes[2].imshow(dccm_std, cmap='viridis', vmin=0, vmax=vmax_std,
                          interpolation='nearest', origin='lower', aspect='equal')
    rect2 = Rectangle((BINDING_SITE_START, BINDING_SITE_START),
                       BINDING_SITE_END - BINDING_SITE_START,
                       BINDING_SITE_END - BINDING_SITE_START,
                       linewidth=1.5, edgecolor='white', facecolor='none',
                       linestyle='--')
    axes[2].add_patch(rect2)
    axes[2].set_title(f'C Inter-block std (max={dccm_std.max():.3f})', fontsize=11)
    axes[2].set_xticks(tick_positions)
    axes[2].set_xticklabels(tick_labels, fontsize=7)
    axes[2].set_yticks(tick_positions)
    axes[2].set_yticklabels(tick_labels, fontsize=7)
    fig.colorbar(im2, ax=axes[2], fraction=0.046, pad=0.04, label='Std')

    fig.suptitle(f'Block-averaged vs Full-sample DCCM: {system_name}', fontsize=13)
    fig.tight_layout(rect=[0, 0, 1, 0.95])
    fig.savefig(outpath, dpi=300, bbox_inches='tight')
    plt.close(fig)
    logging.info(f"  Saved {outpath}")


def plot_variance_plateau(fp_results_all_systems, outpath):
    """Plot variance vs block size for all systems (Flyvbjerg-Petersen).

    Parameters
    ----------
    fp_results_all_systems : dict
        system_name -> list of result dicts from flyvbjerg_petersen_analysis.
    outpath : str
    """
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 5))

    colors = {'phe_sssD_sap': '#2ca02c', 'phe_sssD_tsap': '#d62728',
              'phe_rrrL_sap': '#1f77b4', 'phe_rrrL_tsap': '#ff7f0e'}
    markers = {'phe_sssD_sap': 'o', 'phe_sssD_tsap': 's',
               'phe_rrrL_sap': '^', 'phe_rrrL_tsap': 'D'}

    for sys_name, fp_results in fp_results_all_systems.items():
        color = colors.get(sys_name, '#333333')
        marker = markers.get(sys_name, 'o')
        block_sizes = [r['block_size'] for r in fp_results]
        mean_vars = [r['mean_var'] for r in fp_results]

        ax1.plot(block_sizes, mean_vars, '-', color=color, marker=marker,
                 label=sys_name, markersize=6, linewidth=1.5)
        # 1/block_size axis
        inv_bs = [1.0 / b for b in block_sizes]
        ax2.plot(inv_bs, mean_vars, '-', color=color, marker=marker,
                 label=sys_name, markersize=6, linewidth=1.5)

    ax1.set_xlabel('Block size (frames)', fontsize=11)
    ax1.set_ylabel('Mean element variance', fontsize=11)
    ax1.set_title('A Variance vs Block Size', fontsize=12)
    ax1.legend(fontsize=9)
    ax1.grid(True, alpha=0.3)

    ax2.set_xlabel('1 / Block Size', fontsize=11)
    ax2.set_ylabel('Mean element variance', fontsize=11)
    ax2.set_title('B Flyvbjerg-Petersen Plot', fontsize=12)
    ax2.legend(fontsize=9)
    ax2.grid(True, alpha=0.3)

    # Highlight plateau region: find where variance change < 10%
    # between consecutive block sizes
    for sys_name, fp_results in fp_results_all_systems.items():
        mean_vars = [r['mean_var'] for r in fp_results]
        block_sizes = [r['block_size'] for r in fp_results]
        for i in range(1, len(mean_vars)):
            pct_change = abs(mean_vars[i] - mean_vars[i-1]) / mean_vars[i-1]
            if pct_change < 0.10:
                ax1.axvline(x=block_sizes[i], color=colors.get(sys_name, '#333'),
                           linestyle=':', alpha=0.5)
                break

    fig.suptitle('Flyvbjerg-Petersen Variance Plateau Test', fontsize=13)
    fig.tight_layout(rect=[0, 0, 1, 0.95])
    fig.savefig(outpath, dpi=300, bbox_inches='tight')
    plt.close(fig)
    logging.info(f"  Saved {outpath}")


def plot_interblock_correlation(r_matrix, system_name, block_size, outpath):
    """K x K heatmap of inter-block Pearson r.

    Parameters
    ----------
    r_matrix : ndarray, shape (K, K)
    system_name : str
    block_size : int
    outpath : str
    """
    K = r_matrix.shape[0]
    fig, ax = plt.subplots(figsize=(max(5, K + 1), max(4, K)))

    im = ax.imshow(r_matrix, cmap='YlOrRd', vmin=0.5, vmax=1.0,
                    interpolation='nearest', origin='lower', aspect='equal')

    # Annotate cells with r values
    for i in range(K):
        for j in range(K):
            ax.text(j, i, f'{r_matrix[i, j]:.2f}',
                    ha='center', va='center', fontsize=7,
                    color='white' if r_matrix[i, j] > 0.85 else 'black')

    ax.set_xticks(range(K))
    ax.set_xticklabels([f'B{k+1}' for k in range(K)], fontsize=9)
    ax.set_yticks(range(K))
    ax.set_yticklabels([f'B{k+1}' for k in range(K)], fontsize=9)
    ax.set_title(f'Inter-block Pearson r: {system_name} (B={block_size})', fontsize=12)
    fig.colorbar(im, ax=ax, label='Pearson r')
    fig.tight_layout()
    fig.savefig(outpath, dpi=300, bbox_inches='tight')
    plt.close(fig)
    logging.info(f"  Saved {outpath}")


# ---------------------------------------------------------------------------
# 3. Save and report functions
# ---------------------------------------------------------------------------

def save_results(system_name, block_size, n_blocks, dccm_blocks,
                 dccm_avg, dccm_std, comparison, r_matrix, ib_summary,
                 fp_results, output_dir):
    """Save all outputs for a single system at a single block size.

    Parameters
    ----------
    system_name : str
    block_size : int
    n_blocks : int
    dccm_blocks : list of ndarray
    dccm_avg, dccm_std : ndarray
    comparison : dict from compare_dccm_pair
    r_matrix : ndarray
    ib_summary : dict
    fp_results : list of dict (Flyvbjerg-Petersen for all block sizes)
    output_dir : str
        Base output directory.
    """
    # Per-block DCCMs (stacked)
    stacked_blocks = np.stack(dccm_blocks, axis=0)
    np.save(os.path.join(output_dir, 'per_block_dccm',
                         f'block_dccm_{system_name}_{block_size}.npy'),
            stacked_blocks)

    # Block-averaged DCCM and std
    np.save(os.path.join(output_dir, 'block_averaged_dccm',
                         f'block_avg_dccm_{system_name}_{block_size}.npy'),
            dccm_avg)
    np.save(os.path.join(output_dir, 'block_averaged_dccm',
                         f'block_std_dccm_{system_name}_{block_size}.npy'),
            dccm_std)

    # Inter-block correlation CSV
    K = r_matrix.shape[0]
    r_df = pd.DataFrame(r_matrix,
                        index=[f'block_{k+1}' for k in range(K)],
                        columns=[f'block_{k+1}' for k in range(K)])
    r_df.to_csv(os.path.join(output_dir, 'convergence_metrics',
                             f'interblock_r_{system_name}_{block_size}.csv'))

    # Flyvbjerg-Petersen variance plateau CSV
    fp_df = pd.DataFrame(fp_results)
    fp_df.to_csv(os.path.join(output_dir, 'convergence_metrics',
                              f'variance_plateau_{system_name}.csv'),
                 index=False)

    # Collect per-block-size stats row
    triu_idx = np.triu_indices(N_CA, k=1)
    triu_std = dccm_std[triu_idx]
    bs_std = dccm_std[BINDING_SITE_START:BINDING_SITE_END,
                     BINDING_SITE_START:BINDING_SITE_END]
    frac_high_var = np.sum(triu_std > 0.1) / len(triu_std)

    # Get Flyvbjerg-Petersen stats for this block size
    fp_row = [r for r in fp_results if r['block_size'] == block_size]
    fp_mean_var = fp_row[0]['mean_var'] if fp_row else np.nan
    fp_median_var = fp_row[0]['median_var'] if fp_row else np.nan

    row = {
        'system': system_name,
        'block_size': block_size,
        'n_blocks': n_blocks,
        'pearson_r_block_avg_vs_full': comparison['pearson_r'],
        'rmse_block_avg_vs_full': comparison['rmse'],
        'mean_interblock_r': ib_summary['mean_r_off_diag'],
        'min_interblock_r': ib_summary['min_r_off_diag'],
        'max_interblock_r': ib_summary['max_r_off_diag'],
        'mean_interblock_std': triu_std.mean(),
        'max_interblock_std': triu_std.max(),
        'binding_site_mean_std': bs_std.mean(),
        'fraction_high_variance': frac_high_var,
        'mean_var_flyvbjerg': fp_mean_var,
        'median_var_flyvbjerg': fp_median_var,
    }

    return row


def generate_comparison_report(input_dir, output_dir, summary_csv_path):
    """Generate convergence comparison Markdown report.

    Reads existing dccm_convergence_analysis.csv and autocorrelation_analysis.csv
    and compares with block averaging results.

    Parameters
    ----------
    input_dir : str
        Directory with existing analysis CSVs.
    output_dir : str
        Block averaging output directory.
    summary_csv_path : str
        Path to block_averaging_summary.csv.
    """
    # Load data
    conv_path = os.path.join(input_dir, 'dccm_convergence_analysis.csv')
    auto_path = os.path.join(input_dir, 'autocorrelation_analysis.csv')

    conv_df = pd.read_csv(conv_path) if os.path.exists(conv_path) else None
    auto_df = pd.read_csv(auto_path) if os.path.exists(auto_path) else None
    summary_df = pd.read_csv(summary_csv_path)

    # Focus on block_size=500 for the primary comparison
    primary_df = summary_df[summary_df['block_size'] == 500]

    lines = []
    lines.append("# Block Averaging Comparison Report")
    lines.append("")
    lines.append("## Executive Summary")
    lines.append("")

    # Determine phe_sssD_tsap status
    tsap_row = primary_df[primary_df['system'] == 'phe_sssD_tsap']
    if len(tsap_row) > 0:
        tsap_r = tsap_row['mean_interblock_r'].values[0]
        tsap_min_r = tsap_row['min_interblock_r'].values[0]
        tsap_avg_vs_full = tsap_row['pearson_r_block_avg_vs_full'].values[0]

        if tsap_r > 0.9:
            exec_summary = (
                "Block averaging confirms phe_sssD_tsap DCCM is **reliable** despite slow decorrelation. "
                f"Inter-block mean r = {tsap_r:.3f} indicates consistent DCCM patterns across blocks. "
                "The non-monotonic convergence observed in cumulative-frame analysis was noise, "
                "not structural — no metastable states detected."
            )
        elif tsap_min_r < 0.7:
            exec_summary = (
                "Block averaging reveals **metastable states** in phe_sssD_tsap. "
                f"Inter-block min r = {tsap_min_r:.3f} indicates outlier blocks with divergent DCCM patterns. "
                "The block-averaged DCCM should be reported as the primary result, "
                "and the DCCM should be qualified in the manuscript."
            )
        elif tsap_r < 0.8:
            exec_summary = (
                "Block averaging indicates phe_sssD_tsap DCCM is **not converged**. "
                f"Inter-block mean r = {tsap_r:.3f} shows substantial inconsistency across blocks. "
                "Conclusions involving this system should be strongly qualified, "
                "and extended simulation is recommended."
            )
        else:
            exec_summary = (
                "Block averaging shows phe_sssD_tsap DCCM has **moderate reliability**. "
                f"Inter-block mean r = {tsap_r:.3f}, min r = {tsap_min_r:.3f}. "
                "The DCCM captures the dominant correlation pattern but has some inter-block variability. "
                "Recommend reporting block-averaged DCCM with uncertainty estimates."
            )
    else:
        exec_summary = "Insufficient data for phe_sssD_tsap block averaging assessment."

    lines.append(exec_summary)
    lines.append("")

    # Per-system results
    lines.append("## Per-System Results")
    lines.append("")

    for sys_name in ['phe_sssD_sap', 'phe_sssD_tsap', 'phe_rrrL_sap', 'phe_rrrL_tsap']:
        sys_primary = primary_df[primary_df['system'] == sys_name]
        if len(sys_primary) == 0:
            lines.append(f"### {sys_name}")
            lines.append("- No data available")
            lines.append("")
            continue

        row = sys_primary.iloc[0]

        # Get existing convergence data
        sys_conv = conv_df[conv_df['system'] == sys_name] if conv_df is not None else None
        conv_3k_r = "N/A"
        if sys_conv is not None and len(sys_conv) > 0:
            r3k = sys_conv[sys_conv['n_frames_sub'] == 3000]
            if len(r3k) > 0:
                conv_3k_r = f"{r3k['pearson_r'].values[0]:.4f}"

        # Get autocorrelation data
        sys_auto = auto_df[auto_df['system'] == sys_name] if auto_df is not None else None
        neff_str = "N/A"
        if sys_auto is not None and len(sys_auto) > 0:
            mean_neff = sys_auto['n_eff'].mean()
            mean_tau = sys_auto['tau_autocorr'].mean()
            neff_str = f"{mean_neff:.1f} (mean tau={mean_tau:.0f})"

        # Decision tree assessment
        mean_r = row['mean_interblock_r']
        min_r = row['min_interblock_r']
        avg_vs_full_r = row['pearson_r_block_avg_vs_full']

        if mean_r > 0.9:
            assessment = "Reliable"
            recommendation = "Full-sample DCCM is trustworthy. No qualification needed."
        elif min_r < 0.7:
            assessment = "Metastable states detected"
            recommendation = "Report block-averaged DCCM with uncertainty. Qualify in manuscript."
        elif mean_r < 0.8:
            assessment = "Not converged"
            recommendation = "Strongly qualify conclusions. Recommend extended simulation."
        else:
            assessment = "Moderately reliable"
            recommendation = "Report block-averaged DCCM with inter-block std as uncertainty bars."

        lines.append(f"### {sys_name}")
        lines.append(f"- Full-sample DCCM convergence (3k frames): r = {conv_3k_r}")
        lines.append(f"- Block-averaged DCCM vs full-sample: r = {avg_vs_full_r:.4f}")
        lines.append(f"- Per-block r range: {min_r:.4f} – {row['max_interblock_r']:.4f}")
        lines.append(f"- Mean inter-block r: {mean_r:.4f}")
        lines.append(f"- Mean inter-block std: {row['mean_interblock_std']:.4f}")
        lines.append(f"- Binding site mean std: {row['binding_site_mean_std']:.4f}")
        lines.append(f"- Fraction high-variance entries (std > 0.1): {row['fraction_high_variance']:.4f}")
        lines.append(f"- Autocorrelation N_eff: {neff_str}")
        lines.append(f"- **Assessment:** {assessment}")
        lines.append(f"- **Recommendation:** {recommendation}")
        lines.append("")

    # Flyvbjerg-Petersen variance plateau test
    lines.append("## Flyvbjerg-Petersen Variance Plateau Test")
    lines.append("")

    # Load variance plateau CSVs for each system
    plateau_info = []
    for sys_name in ['phe_sssD_sap', 'phe_sssD_tsap', 'phe_rrrL_sap', 'phe_rrrL_tsap']:
        vp_path = os.path.join(output_dir, 'convergence_metrics',
                               f'variance_plateau_{sys_name}.csv')
        if os.path.exists(vp_path):
            vp_df = pd.read_csv(vp_path)
            # Find plateau: first block size where variance change < 10%
            plateau_bs = None
            for i in range(1, len(vp_df)):
                pct_change = abs(vp_df.iloc[i]['mean_var'] - vp_df.iloc[i-1]['mean_var']) / vp_df.iloc[i-1]['mean_var']
                if pct_change < 0.10:
                    plateau_bs = int(vp_df.iloc[i]['block_size'])
                    break
            plateau_info.append({
                'system': sys_name,
                'plateau_block_size': plateau_bs if plateau_bs else 'Not reached (B<=1000)',
                'var_at_500': vp_df[vp_df['block_size'] == 500]['mean_var'].values[0] if len(vp_df[vp_df['block_size'] == 500]) > 0 else 'N/A',
                'var_at_1000': vp_df[vp_df['block_size'] == 1000]['mean_var'].values[0] if len(vp_df[vp_df['block_size'] == 1000]) > 0 else 'N/A',
            })

    if plateau_info:
        plateau_df = pd.DataFrame(plateau_info)
        # Manual markdown table (no tabulate dependency)
        cols = plateau_df.columns.tolist()
        lines.append('| ' + ' | '.join(str(c) for c in cols) + ' |')
        lines.append('| ' + ' | '.join(['---'] * len(cols)) + ' |')
        for _, row in plateau_df.iterrows():
            lines.append('| ' + ' | '.join(str(row[c]) for c in cols) + ' |')
    else:
        lines.append("Variance plateau data not available.")
    lines.append("")

    # Convergence improvement
    lines.append("## Convergence Improvement")
    lines.append("")
    lines.append("| System | Cumulative r (3k) | Block-avg r (B=500) vs full | Improvement |")
    lines.append("|--------|-------------------|------------------------------|-------------|")

    for sys_name in ['phe_sssD_sap', 'phe_sssD_tsap', 'phe_rrrL_sap', 'phe_rrrL_tsap']:
        sys_primary = primary_df[primary_df['system'] == sys_name]
        if len(sys_primary) == 0:
            lines.append(f"| {sys_name} | N/A | N/A | N/A |")
            continue

        avg_r = sys_primary.iloc[0]['pearson_r_block_avg_vs_full']

        sys_conv = conv_df[conv_df['system'] == sys_name] if conv_df is not None else None
        cum_r = "N/A"
        if sys_conv is not None and len(sys_conv) > 0:
            r3k = sys_conv[sys_conv['n_frames_sub'] == 3000]
            if len(r3k) > 0:
                cum_r = f"{r3k['pearson_r'].values[0]:.4f}"

        improvement = "N/A"
        if conv_df is not None and len(r3k) > 0:
            improvement = f"{avg_r - r3k['pearson_r'].values[0]:+.4f}"
        lines.append(f"| {sys_name} | {cum_r} | {avg_r:.4f} | {improvement} |")

    lines.append("")
    lines.append("Note: Block-avg r measures how well the block-averaged DCCM reproduces")
    lines.append("the full-sample DCCM. Higher r means the block average is a good")
    lines.append("representation of the full-sample DCCM, suggesting consistent")
    lines.append("correlation patterns across blocks.")
    lines.append("")

    # Recommendations
    lines.append("## Recommendations")
    lines.append("")

    tsap_row = primary_df[primary_df['system'] == 'phe_sssD_tsap']
    if len(tsap_row) > 0:
        tsap_mean_r = tsap_row['mean_interblock_r'].values[0]
        tsap_min_r = tsap_row['min_interblock_r'].values[0]
        tsap_bs_std = tsap_row['binding_site_mean_std'].values[0]

        if tsap_mean_r > 0.9:
            lines.append("For **phe_sssD_tsap**: The DCCM is reliable despite non-monotonic convergence.")
            lines.append("Block averaging shows high inter-block consistency (mean r > 0.9).")
            lines.append("The non-monotonic convergence in the cumulative-frame test was noise,")
            lines.append("not a sign of metastable states. No qualification needed in the manuscript.")
        elif tsap_min_r < 0.7:
            lines.append("For **phe_sssD_tsap**: Metastable states detected (min inter-block r < 0.7).")
            lines.append("The DCCM pattern varies significantly between blocks, indicating")
            lines.append("the trajectory visits distinct conformational states.")
            lines.append("Recommend reporting the block-averaged DCCM with inter-block std,")
            lines.append("and qualifying the DCCM interpretation in the manuscript methods.")
        elif tsap_mean_r < 0.8:
            lines.append("For **phe_sssD_tsap**: DCCM is not converged (mean inter-block r < 0.8).")
            lines.append("The correlation patterns differ substantially between blocks.")
            lines.append("Strongly qualify all DCCM-based conclusions for this system.")
            lines.append("Recommend extended simulation (at least 2x current length).")
        else:
            lines.append(f"For **phe_sssD_tsap**: Moderate DCCM reliability (mean inter-block r = {tsap_mean_r:.3f}).")
            lines.append("The dominant correlation pattern is consistent across blocks,")
            lines.append("but some inter-block variability exists.")
            lines.append("Report block-averaged DCCM with inter-block std as uncertainty bars.")
            lines.append(f"Binding site region mean std = {tsap_bs_std:.4f}.")

        lines.append("")

        # Cross-system comparison
        for other in ['phe_sssD_sap', 'phe_rrrL_sap', 'phe_rrrL_tsap']:
            other_row = primary_df[primary_df['system'] == other]
            if len(other_row) > 0:
                other_mean_r = other_row.iloc[0]['mean_interblock_r']
                other_std = other_row.iloc[0]['mean_interblock_std']
                lines.append(f"- {other}: inter-block mean r = {other_mean_r:.4f}, "
                             f"mean std = {other_std:.4f}")
    else:
        lines.append("Insufficient data for phe_sssD_tsap block averaging assessment.")

    lines.append("")

    # Write report
    report_path = os.path.join(output_dir, 'summary', 'convergence_comparison.md')
    with open(report_path, 'w') as f:
        f.write('\n'.join(lines))
    logging.info(f"  Saved comparison report: {report_path}")

    return report_path


# ---------------------------------------------------------------------------
# 4. Main function
# ---------------------------------------------------------------------------

def parse_args():
    parser = argparse.ArgumentParser(
        description='Block averaging for DCCM convergence analysis')
    parser.add_argument('--systems', nargs='+',
                        default=['phe_sssD_sap', 'phe_sssD_tsap',
                                 'phe_rrrL_sap', 'phe_rrrL_tsap'],
                        help='System names to analyze')
    parser.add_argument('--block-sizes', nargs='+', type=int,
                        default=[400, 500, 572, 667, 800, 1000],
                        help='Block sizes for Flyvbjerg-Petersen test')
    parser.add_argument('--input-dir', default='pca_analysis/',
                        help='Directory with aligned_coords.npy, system_indices.npy, system_metadata.json')
    parser.add_argument('--output-dir', default='pca_analysis/block_averaging_data/',
                        help='Output directory for all results')
    parser.add_argument('--force', action='store_true',
                        help='Overwrite existing output directory contents')
    parser.add_argument('--report-only', action='store_true',
                        help='Only generate comparison report (skip computation)')
    return parser.parse_args()


def main():
    args = parse_args()

    input_dir = args.input_dir
    output_dir = args.output_dir

    logging.info("=" * 60)
    logging.info("Block Averaging for DCCM Convergence Analysis")
    logging.info("=" * 60)
    logging.info(f"Input dir:  {input_dir}")
    logging.info(f"Output dir: {output_dir}")
    logging.info(f"Systems:    {args.systems}")
    logging.info(f"Block sizes: {args.block_sizes}")

    # Create output directory structure
    for subdir in ['per_block_dccm', 'block_averaged_dccm',
                   'convergence_metrics', 'figures', 'summary']:
        os.makedirs(os.path.join(output_dir, subdir), exist_ok=True)

    # Check for existing output (unless --force)
    if not args.force and not args.report_only:
        existing_files = []
        for subdir in ['per_block_dccm', 'block_averaged_dccm',
                       'convergence_metrics', 'figures', 'summary']:
            subpath = os.path.join(output_dir, subdir)
            if os.path.exists(subpath):
                for f in os.listdir(subpath):
                    if f.endswith(('.npy', '.csv', '.png', '.md')):
                        existing_files.append(os.path.join(subpath, f))
        if existing_files:
            logging.error(
                f"Output directory contains {len(existing_files)} existing files. "
                f"Use --force to overwrite."
            )
            sys.exit(1)

    # If report-only mode, just generate the report
    summary_csv_path = os.path.join(output_dir, 'summary', 'block_averaging_summary.csv')
    if args.report_only:
        if os.path.exists(summary_csv_path):
            generate_comparison_report(input_dir, output_dir, summary_csv_path)
            logging.info("Report-only mode: comparison report generated.")
            return
        else:
            logging.error("Report-only mode requires existing summary CSV. Run full analysis first.")
            sys.exit(1)

    # Load shared data
    t_start = time.time()
    logging.info("Loading aligned coordinates and system data...")

    aligned_path = os.path.join(input_dir, 'aligned_coords.npy')
    aligned_mmap = np.load(aligned_path, mmap_mode='r')
    system_indices = np.load(os.path.join(input_dir, 'system_indices.npy'))
    with open(os.path.join(input_dir, 'system_metadata.json')) as f:
        systems_meta = json.load(f)

    # Build system name -> index mapping
    sys_name_to_idx = {s['name']: i for i, s in enumerate(systems_meta)}

    # Build CA flat indices
    ca_flat_indices = build_ca_flat_indices()
    logging.info(f"CA flat coordinate indices: {len(ca_flat_indices)}")

    # Primary block size for figures
    primary_block_size = 500

    # Collect all stats rows
    all_stats_rows = []
    fp_results_all_systems = {}

    # Process each system
    for sys_name in args.systems:
        logging.info(f"\n{'='*50}")
        logging.info(f"Processing system: {sys_name}")
        logging.info(f"{'='*50}")

        if sys_name not in sys_name_to_idx:
            logging.warning(f"System {sys_name} not found in metadata. Skipping.")
            continue

        idx = sys_name_to_idx[sys_name]
        start = int(system_indices[idx])
        end = int(system_indices[idx + 1])
        n_frames = end - start

        logging.info(f"  Frame range: {start}-{end} ({n_frames} frames)")

        # Extract CA coordinates for this system
        ca_coords = aligned_mmap[start:end, ca_flat_indices].astype(np.float64)

        # Load or compute full-sample DCCM
        dccm_full_path = os.path.join(input_dir, f'dccm_{sys_name}.npy')
        if os.path.exists(dccm_full_path):
            logging.info(f"  Loading full-sample DCCM from {dccm_full_path}")
            dccm_full = np.load(dccm_full_path)
        else:
            logging.info(f"  Computing full-sample DCCM from all {n_frames} frames...")
            dccm_full = compute_dccm(ca_coords, n_ca_atoms=N_CA)

        # Flyvbjerg-Petersen analysis (all block sizes)
        logging.info(f"  Running Flyvbjerg-Petersen analysis for {sys_name}...")
        fp_results = flyvbjerg_petersen_analysis(ca_coords, args.block_sizes)
        fp_results_all_systems[sys_name] = fp_results

        # Process each block size
        for block_size in args.block_sizes:
            logging.info(f"  Block size: {block_size}")

            # Compute per-block DCCMs
            dccm_blocks, n_blocks = compute_block_dccms(ca_coords, block_size)

            # Block-averaged DCCM
            dccm_avg, dccm_std = block_average_dccm(dccm_blocks)

            # Compare block-averaged DCCM with full-sample DCCM
            comparison = compare_dccm_pair(dccm_avg, dccm_full)
            logging.info(f"    Block-avg vs full: r={comparison['pearson_r']:.4f}, "
                         f"RMSE={comparison['rmse']:.4f}")

            # Inter-block correlation matrix
            r_matrix, ib_summary = interblock_correlation_matrix(dccm_blocks)
            logging.info(f"    Inter-block r: mean={ib_summary['mean_r_off_diag']:.4f}, "
                         f"min={ib_summary['min_r_off_diag']:.4f}")

            # Save results
            row = save_results(sys_name, block_size, n_blocks, dccm_blocks,
                              dccm_avg, dccm_std, comparison, r_matrix,
                              ib_summary, fp_results, output_dir)
            all_stats_rows.append(row)

            # Generate figures for primary block size only
            if block_size == primary_block_size:
                fig_dir = os.path.join(output_dir, 'figures')

                # Per-block DCCM overview
                plot_block_dccms_overview(
                    dccm_blocks, sys_name, block_size,
                    os.path.join(fig_dir, f'fig_block_dccms_{sys_name}_{block_size}.png'))

                # Block-averaged vs full
                plot_block_avg_vs_full(
                    dccm_avg, dccm_full, dccm_std, sys_name,
                    os.path.join(fig_dir, f'fig_block_avg_vs_full_{sys_name}.png'))

                # Inter-block correlation
                plot_interblock_correlation(
                    r_matrix, sys_name, block_size,
                    os.path.join(fig_dir, f'fig_interblock_r_{sys_name}_{block_size}.png'))

            # Free memory for this block size
            del dccm_blocks

        # Save per-system block stats CSV
        sys_rows = [r for r in all_stats_rows if r['system'] == sys_name]
        sys_stats_df = pd.DataFrame(sys_rows)
        sys_stats_df.to_csv(os.path.join(output_dir, 'convergence_metrics',
                                          f'block_dccm_stats_{sys_name}.csv'),
                             index=False)
        logging.info(f"  Saved block stats for {sys_name}")

        # Free system memory
        del ca_coords, dccm_full

    # Cross-system variance plateau figure
    if fp_results_all_systems:
        plot_variance_plateau(
            fp_results_all_systems,
            os.path.join(output_dir, 'figures', 'fig_variance_plateau_all_systems.png'))

    # Save cross-system summary CSV
    if all_stats_rows:
        summary_df = pd.DataFrame(all_stats_rows)
        summary_df.to_csv(summary_csv_path, index=False)
        logging.info(f"Saved summary CSV: {summary_csv_path}")

    # Generate comparison report
    if os.path.exists(summary_csv_path):
        generate_comparison_report(input_dir, output_dir, summary_csv_path)

    elapsed = time.time() - t_start
    logging.info(f"\n{'='*60}")
    logging.info(f"Block averaging analysis complete ({elapsed:.0f}s)")
    logging.info(f"{'='*60}")

    # Print summary table
    if all_stats_rows:
        summary_df = pd.DataFrame(all_stats_rows)
        primary = summary_df[summary_df['block_size'] == primary_block_size]
        print("\n" + "=" * 80)
        print(f"BLOCK AVERAGING SUMMARY (block_size={primary_block_size})")
        print("=" * 80)
        for _, row in primary.iterrows():
            print(f"\n  {row['system']}:")
            print(f"    Block-avg vs full r:    {row['pearson_r_block_avg_vs_full']:.4f}")
            print(f"    RMSE:                   {row['rmse_block_avg_vs_full']:.4f}")
            print(f"    Mean inter-block r:     {row['mean_interblock_r']:.4f}")
            print(f"    Min inter-block r:      {row['min_interblock_r']:.4f}")
            print(f"    Mean inter-block std:   {row['mean_interblock_std']:.4f}")
            print(f"    Binding site mean std:  {row['binding_site_mean_std']:.4f}")
            print(f"    Fraction high-var:     {row['fraction_high_variance']:.4f}")


if __name__ == '__main__':
    main()
