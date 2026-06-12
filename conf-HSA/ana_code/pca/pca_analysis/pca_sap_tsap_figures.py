import argparse
import json
import os
import sys

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
import matplotlib.lines as mlines
import matplotlib.patches as mpatches
from matplotlib.patches import Rectangle
import numpy as np
np.random.seed(42)  # Reproducibility for stochastic analyses
import pandas as pd

FONT = {
    'tick': 11,
    'axis': 13,
    'panel': 14,
    'panel_weight': 'bold',
    'legend': 10,
    'annotation': 10,
    'annotation_bold': 12,
    'colorbar': 10,
    'heatmap_cell': 9,
    'heatmap_tick': 9,
    'subtitle': 12,
}

GROUP_COLORS = {
    'me_sss': '#E69F00',
    'me_rrr': '#56B4E9',
    'phe_sss': '#009E73',
    'phe_rrr': '#CC79A7',
}

SAP_COLORS = {
    'phe_sss': '#009E73',
    'phe_rrr': '#CC79A7',
}

TSAP_COLORS = {
    'phe_sss': '#006D5B',
    'phe_rrr': '#7B2D8E',
}

GROUP_LABELS = {
    'me_sss': 'Me-SSS',
    'me_rrr': 'Me-RRR',
    'phe_sss': 'Phe-SSS',
    'phe_rrr': 'Phe-RRR',
}

GROUP_ORDER = ['me_sss', 'me_rrr', 'phe_sss', 'phe_rrr']

RC_PARAMS = {
    'font.size': FONT['tick'],
    'font.family': 'sans-serif',
    'axes.labelsize': FONT['axis'],
    'axes.titlesize': FONT['axis'],
    'xtick.labelsize': FONT['tick'],
    'ytick.labelsize': FONT['tick'],
    'legend.fontsize': FONT['legend'],
    'figure.dpi': 300,
    'savefig.dpi': 300,
    'axes.linewidth': 0.8,
    'lines.linewidth': 1.0,
}

DPI = 300
# kB / R_GAS and kT are imported from pca_utils (via compute_fes default kT=2.4942)

KEY_PHE_SYSTEMS = ['phe_sssD_sap', 'phe_sssD_tsap', 'phe_rrrL_sap', 'phe_rrrL_tsap']

ALL_8_PAIRS = [
    ('phe_sssD_sap', 'phe_sssD_tsap'),
    ('phe_sssL_sap', 'phe_sssL_tsap'),
    ('phe_rrrD_sap', 'phe_rrrD_tsap'),
    ('phe_rrrL_sap', 'phe_rrrL_tsap'),
    ('me_sssD_sap', 'me_sssD_tsap'),
    ('me_sssL_sap', 'me_sssL_tsap'),
    ('me_rrrD_sap', 'me_rrrD_tsap'),
    ('me_rrrL_sap', 'me_rrrL_tsap'),
]

N_BACKBONE = 2327
N_CA = 582
N_COORDS_PER_CA = 3
# Binding site and DCCM constants imported from shared module (S3/M4 fix)
from pca_utils import (compute_dccm, DCCM_BINDING_SITE_START, DCCM_BINDING_SITE_END,
                        DCCM_BINDING_SITE_RESID_OFFSET, compute_fes, R_GAS)
BINDING_SITE_START = DCCM_BINDING_SITE_START   # 376
BINDING_SITE_END = DCCM_BINDING_SITE_END       # 487
RESID_OFFSET = DCCM_BINDING_SITE_RESID_OFFSET  # 3

DISPLAY_NAMES = {
    'phe_sssD_sap': 'Phe-SSS-D sap',
    'phe_sssD_tsap': 'Phe-SSS-D tsap',
    'phe_rrrL_sap': 'Phe-RRR-L sap',
    'phe_rrrL_tsap': 'Phe-RRR-L tsap',
}


def system_display_name(label):
    parts = label.split('_')
    lig = parts[0].capitalize()
    stereo_enant = parts[1]
    stereo = stereo_enant[:-1].upper()
    enant = stereo_enant[-1].upper()
    bind = parts[2]
    return f"{lig}-{stereo}-{enant} {bind}"


def system_to_group(system_label):
    parts = system_label.split('_')
    return f"{parts[0]}_{parts[1][:3]}"


def add_panel_label(ax, label, x=-0.08, y=1.10):
    ax.text(x, y, f"{label.upper()}", transform=ax.transAxes,
            fontsize=FONT['panel'], fontweight=FONT['panel_weight'], va='bottom')


def save_figure(fig, outpath_base):
    fig.savefig(outpath_base + '.png', dpi=DPI, bbox_inches='tight')
    fig.savefig(outpath_base + '.svg', format='svg', bbox_inches='tight')
    plt.close(fig)


# compute_fes() is now imported from pca_utils.py (was duplicated here
# and in pca_si_figures.py). The pca_utils version has the same behavior
# (shift to minimum = 0, NaN for empty bins) with an explicit kT parameter.


def plot_phe_scatter(projections_path, scree_path, outpath_base):
    df = pd.read_csv(projections_path)
    scree = pd.read_csv(scree_path)
    pc1_frac = scree.iloc[0]['fraction']
    pc2_frac = scree.iloc[1]['fraction']
    df = df[df['system_label'].isin(KEY_PHE_SYSTEMS)]
    fig, ax = plt.subplots(figsize=(6, 5))
    scatter_specs = [
        ('phe_sssD_sap', SAP_COLORS['phe_sss']),
        ('phe_sssD_tsap', TSAP_COLORS['phe_sss']),
        ('phe_rrrL_sap', SAP_COLORS['phe_rrr']),
        ('phe_rrrL_tsap', TSAP_COLORS['phe_rrr']),
    ]
    for sys_label, color in scatter_specs:
        sub = df[df['system_label'] == sys_label]
        ax.scatter(sub['PC1'], sub['PC2'],
                   c=color, marker='o', s=8, alpha=0.5, rasterized=True,
                   label='_nolegend_')
    legend_elements = []
    for group in ['phe_sss', 'phe_rrr']:
        label_base = 'Phe-SSS-D' if group == 'phe_sss' else 'Phe-RRR-L'
        legend_elements.append(mpatches.Patch(facecolor=SAP_COLORS[group],
                                              label=f'{label_base} sap'))
        legend_elements.append(mpatches.Patch(facecolor=TSAP_COLORS[group],
                                              label=f'{label_base} tsap'))
    ax.legend(handles=legend_elements, loc='best', framealpha=0.9, fontsize=FONT['legend'])
    ax.set_xlabel(f'PC1 ({pc1_frac * 100:.1f}%)')
    ax.set_ylabel(f'PC2 ({pc2_frac * 100:.1f}%)')
    add_panel_label(ax, 'a')
    save_figure(fig, outpath_base)
    print(f"Saved {outpath_base}.png/.svg")


def pc_label(pc_num, var_fracs):
    """Format a PC axis label with variance percentage, e.g. 'PC1 (26.5%)'."""
    return f"PC{pc_num} ({var_fracs[pc_num] * 100:.1f}%)"


def plot_fes_4panel(projections_path, scree_path, outpath_base):
    """3×4 grid: 3 PC-pair rows × 4 system columns for FES of key Phe systems."""
    df = pd.read_csv(projections_path)
    scree = pd.read_csv(scree_path)
    var_fracs = {i + 1: scree.iloc[i]['fraction'] for i in range(3)}

    systems = KEY_PHE_SYSTEMS
    PC_PAIRS = [(1, 2), (1, 3), (2, 3)]  # rows: top, middle, bottom
    n_rows = len(PC_PAIRS)
    n_cols = len(systems)

    T = 298.15
    kT = R_GAS * T  # 2.478957 kJ/mol

    key_sub = df[df["system_label"].isin(systems)]
    pc_ranges = {}
    for pc in [1, 2, 3]:
        vals = key_sub[f"PC{pc}"].values
        pc_ranges[pc] = (vals.min(), vals.max())

    # Compute FES for all 12 panels: (pc_pair, system) → (Gm, xe, ye)
    fes_data = {}  # key = (row_idx, col_idx)
    for row_idx, (pc_x, pc_y) in enumerate(PC_PAIRS):
        for col_idx, sys_label in enumerate(systems):
            sub = df[df['system_label'] == sys_label]
            G, xe, ye = compute_fes(sub[f'PC{pc_x}'].values, sub[f'PC{pc_y}'].values,
                                     n_bins=20, kT=kT,
                                     x_range=pc_ranges[pc_x], y_range=pc_ranges[pc_y])
            Gm = np.ma.masked_invalid(G)
            fes_data[(row_idx, col_idx)] = (Gm, xe, ye)

    # Per-row vmax: one per PC pair, computed across all 4 systems in that row
    row_vmax = {}
    for row_idx in range(n_rows):
        row_vmax[row_idx] = max(fes_data[(row_idx, c)][0].max() for c in range(n_cols))

    global_vmax = max(row_vmax.values())

    fig = plt.figure(figsize=(14, 8))
    gs = GridSpec(n_rows, n_cols, hspace=0.38, wspace=0.30,
                  left=0.07, right=0.88, top=0.95, bottom=0.07)

    row_images = {}  # track last image per row for colorbar

    for row_idx, (pc_x, pc_y) in enumerate(PC_PAIRS):
        for col_idx, sys_label in enumerate(systems):
            ax = fig.add_subplot(gs[row_idx, col_idx])
            Gm, xe, ye = fes_data[(row_idx, col_idx)]

            im = ax.pcolormesh(xe, ye, Gm, cmap='cubehelix', shading='auto',
                               vmin=0, vmax=global_vmax)
            row_images[row_idx] = im  # store for colorbar

            try:
                ax.contour(0.5 * (xe[:-1] + xe[1:]), 0.5 * (ye[:-1] + ye[1:]),
                           Gm, levels=8, colors='white', linewidths=0.3, alpha=0.5)
            except ValueError:
                pass

            name = DISPLAY_NAMES.get(sys_label, system_display_name(sys_label))
            ax.set_title(name, fontsize=FONT['subtitle'])
            ax.set_xlim(pc_ranges[pc_x])
            ax.set_ylim(pc_ranges[pc_y])

            # Axis labels: every subplot gets x-label; y-label only on left column
            ax.set_xlabel(pc_label(pc_x, var_fracs), fontsize=FONT['axis'])
            if col_idx == 0:
                ax.set_ylabel(pc_label(pc_y, var_fracs), fontsize=FONT['axis'])

    # Single colorbar for entire figure
    all_axes = [fig.axes[i] for i in range(n_rows * n_cols)]
    cbar = fig.colorbar(row_images[0], ax=all_axes, shrink=0.6,
                        label='\u0394G (kJ/mol)', pad=0.02)
    cbar.ax.tick_params(labelsize=FONT['tick'])
    cbar.set_label('\u0394G (kJ/mol)', fontsize=FONT['axis'])

    save_figure(fig, outpath_base)
    print(f"Saved {outpath_base}.png/.svg")


def plot_s6_all_fes_pairs(projections_path, outpath_base):
    df = pd.read_csv(projections_path)
    fig = plt.figure(figsize=(6, 8))
    gs = GridSpec(8, 2, hspace=0.60, wspace=0.3,
                  left=0.10, right=0.88, top=0.96, bottom=0.04)
    panel_idx = 0
    for row_idx, (sys_sap, sys_tsap) in enumerate(ALL_8_PAIRS):
        for col_idx, sys_label in enumerate([sys_sap, sys_tsap]):
            ax = fig.add_subplot(gs[row_idx, col_idx])
            sub = df[df['system_label'] == sys_label]
            G, xe, ye = compute_fes(sub['PC1'].values, sub['PC2'].values)
            Gm = np.ma.masked_invalid(G)
            sub_other = df[df['system_label'] == ([sys_sap, sys_tsap][1 - col_idx])]
            G_other, _, _ = compute_fes(sub_other['PC1'].values, sub_other['PC2'].values)
            Gm_other = np.ma.masked_invalid(G_other)
            shared_vmax = max(Gm.max(), Gm_other.max())
            im = ax.pcolormesh(xe, ye, Gm, cmap='cubehelix', shading='auto',
                                vmin=0, vmax=shared_vmax)
            if col_idx == 1:
                fig.colorbar(im, ax=ax, shrink=0.7)
            title = system_display_name(sys_label)
            ax.set_title(title, fontsize=FONT['subtitle'])
            panel_label = chr(65 + panel_idx)
            add_panel_label(ax, panel_label, x=-0.15, y=1.10)
            if row_idx == 7:
                ax.set_xlabel('PC1 (\u00c5)')
            else:
                ax.set_xticklabels([])
            if col_idx == 0:
                ax.set_ylabel('PC2 (\u00c5)')
            else:
                ax.set_yticklabels([])
            panel_idx += 1
    save_figure(fig, outpath_base)
    print(f"Saved {outpath_base}.png/.svg")


def plot_s7_binding_bar(mmpbsa_path, contacts_path, outpath_base):
    df_mmpbsa = pd.read_csv(mmpbsa_path)
    df_contacts = pd.read_csv(contacts_path)
    df_mmpbsa = df_mmpbsa[df_mmpbsa['system'].isin(KEY_PHE_SYSTEMS)]
    df_contacts = df_contacts[df_contacts['system'].isin(KEY_PHE_SYSTEMS)]
    pair_keys = ['phe_sssD', 'phe_rrrL']
    sap_systems = [f"{p}_sap" for p in pair_keys]
    tsap_systems = [f"{p}_tsap" for p in pair_keys]
    fig = plt.figure(figsize=(6, 4))
    gs = GridSpec(1, 2, wspace=0.45)
    ax1 = fig.add_subplot(gs[0, 0])
    ax2 = fig.add_subplot(gs[0, 1])
    x = np.arange(len(pair_keys))
    bar_width = 0.35
    offset = bar_width / 2
    colors_sap = [SAP_COLORS['phe_sss'], SAP_COLORS['phe_rrr']]
    colors_tsap = [TSAP_COLORS['phe_sss'], TSAP_COLORS['phe_rrr']]
    dg_sap = []
    dg_tsap = []
    dg_sem_sap = []
    dg_sem_tsap = []
    for s, t in zip(sap_systems, tsap_systems):
        row_s = df_mmpbsa[df_mmpbsa['system'] == s].iloc[0]
        row_t = df_mmpbsa[df_mmpbsa['system'] == t].iloc[0]
        dg_sap.append(row_s['TOTAL_corrected'])
        dg_tsap.append(row_t['TOTAL_corrected'])
        dg_sem_sap.append(row_s['SEM'])
        dg_sem_tsap.append(row_t['SEM'])
    bars_sap = ax1.bar(x - offset, dg_sap, bar_width, color=colors_sap,
                       edgecolor='black', linewidth=0.5, label='SAP')
    bars_tsap = ax1.bar(x + offset, dg_tsap, bar_width, color=colors_tsap,
                        edgecolor='black', linewidth=0.5, label='TSAP')
    ax1.errorbar(x - offset, dg_sap, yerr=dg_sem_sap, fmt='none', ecolor='black', capsize=2, linewidth=0.5)
    ax1.errorbar(x + offset, dg_tsap, yerr=dg_sem_tsap, fmt='none', ecolor='black', capsize=2, linewidth=0.5)
    for bar, val in zip(bars_sap, dg_sap):
        ax1.text(bar.get_x() + bar.get_width() / 2, val,
                 f'{val:.1f}', ha='center', va='bottom', fontsize=FONT['annotation'])
    for bar, val in zip(bars_tsap, dg_tsap):
        ax1.text(bar.get_x() + bar.get_width() / 2, val,
                 f'{val:.1f}', ha='center', va='bottom', fontsize=FONT['annotation'])
    ax1.set_xticks(x)
    ax1.set_xticklabels(['Phe-SSS-D', 'Phe-RRR-L'])
    ax1.set_ylabel('\u0394G_corrected (kJ/mol)')
    legend_elements_s7 = []
    for group, lbl in [('phe_sss', 'Phe-SSS-D'), ('phe_rrr', 'Phe-RRR-L')]:
        legend_elements_s7.append(mpatches.Patch(facecolor=SAP_COLORS[group],
                                                  edgecolor='black', linewidth=0.5,
                                                  label=f'{lbl} sap'))
        legend_elements_s7.append(mpatches.Patch(facecolor=TSAP_COLORS[group],
                                                  edgecolor='black', linewidth=0.5,
                                                  label=f'{lbl} tsap'))
    ax1.legend(handles=legend_elements_s7, loc='best', fontsize=FONT['legend'])
    add_panel_label(ax1, 'a')
    ct_sap = []
    ct_tsap = []
    ct_std_sap = []
    ct_std_tsap = []
    for s, t in zip(sap_systems, tsap_systems):
        row_s = df_contacts[df_contacts['system'] == s].iloc[0]
        row_t = df_contacts[df_contacts['system'] == t].iloc[0]
        ct_sap.append(row_s['mean_total_contacts'])
        ct_tsap.append(row_t['mean_total_contacts'])
        ct_std_sap.append(row_s['std_total_contacts'])
        ct_std_tsap.append(row_t['std_total_contacts'])
    bars_sap2 = ax2.bar(x - offset, ct_sap, bar_width, color=colors_sap,
                        edgecolor='black', linewidth=0.5, label='SAP')
    bars_tsap2 = ax2.bar(x + offset, ct_tsap, bar_width, color=colors_tsap,
                         edgecolor='black', linewidth=0.5, label='TSAP')
    ax2.errorbar(x - offset, ct_sap, yerr=ct_std_sap, fmt='none', ecolor='black', capsize=2, linewidth=0.5)
    ax2.errorbar(x + offset, ct_tsap, yerr=ct_std_tsap, fmt='none', ecolor='black', capsize=2, linewidth=0.5)
    for bar, val in zip(bars_sap2, ct_sap):
        ax2.text(bar.get_x() + bar.get_width() / 2, val,
                 f'{val:.0f}', ha='center', va='bottom', fontsize=FONT['annotation'])
    for bar, val in zip(bars_tsap2, ct_tsap):
        ax2.text(bar.get_x() + bar.get_width() / 2, val,
                 f'{val:.0f}', ha='center', va='bottom', fontsize=FONT['annotation'])
    ax2.set_xticks(x)
    ax2.set_xticklabels(['Phe-SSS-D', 'Phe-RRR-L'])
    ax2.set_ylabel('Mean total contacts')
    ax2.legend(handles=legend_elements_s7, loc='best', fontsize=FONT['legend'])
    add_panel_label(ax2, 'b')
    save_figure(fig, outpath_base)
    print(f"Saved {outpath_base}.png/.svg")


def build_ca_flat_indices():
    ca_atom_indices = np.arange(1, N_BACKBONE, 4)
    # Validate CA extraction: expect N_CA indices matching the constant
    assert len(ca_atom_indices) == N_CA, (
        f"CA extraction validation failed: expected {N_CA} CA atoms, "
        f"got {len(ca_atom_indices)} from np.arange(1, {N_BACKBONE}, 4)"
    )
    ca_flat_indices = np.sort(np.concatenate(
        [3 * ca_atom_indices + k for k in range(N_COORDS_PER_CA)]
    ))
    return ca_flat_indices


# compute_dccm() is now imported from pca_utils.py (M4 fix: was duplicated
# here and in pca_dccm.py). The pca_utils version accepts n_ca_atoms and
# n_coords_per_ca as parameters but uses the same algorithm.


def plot_s8_phe_dccm(aligned_path, system_indices_path, metadata_path, outpath_base):
    ca_flat_indices = build_ca_flat_indices()
    aligned_mmap = np.load(aligned_path, mmap_mode='r')
    system_indices = np.load(system_indices_path)
    with open(metadata_path) as f:
        systems_meta = json.load(f)
    name_to_idx = {s['name']: i for i, s in enumerate(systems_meta)}
    out_dir = os.path.dirname(outpath_base)
    dccm_dict = {}
    for sys_label in KEY_PHE_SYSTEMS:
        idx = name_to_idx[sys_label]
        start = int(system_indices[idx])
        end = int(system_indices[idx + 1])
        ca_coords = aligned_mmap[start:end, ca_flat_indices]
        dccm = compute_dccm(ca_coords)
        dccm_dict[sys_label] = dccm
        npy_path = os.path.join(out_dir, f'dccm_{sys_label}.npy')
        np.save(npy_path, dccm)
        print(f"Saved {npy_path}")
    fig = plt.figure(figsize=(6, 6))
    gs = GridSpec(2, 2, hspace=0.35, wspace=0.3)
    tick_positions = np.linspace(0, N_CA - 1, 6, dtype=int)
    tick_labels = [str(p + RESID_OFFSET) for p in tick_positions]
    im = None
    axes = []
    for i, sys_label in enumerate(KEY_PHE_SYSTEMS):
        ax = fig.add_subplot(gs[i // 2, i % 2])
        axes.append(ax)
        dccm = dccm_dict[sys_label]
        im = ax.imshow(dccm, cmap='RdBu_r', vmin=-1, vmax=1,
                       interpolation='nearest', origin='lower', aspect='equal')
        rect = Rectangle((BINDING_SITE_START, BINDING_SITE_START),
                          BINDING_SITE_END - BINDING_SITE_START,
                          BINDING_SITE_END - BINDING_SITE_START,
                          linewidth=1.5, edgecolor='#F0E442', facecolor='none',
                          linestyle='--')
        ax.add_patch(rect)
        ax.set_xticks(tick_positions)
        ax.set_xticklabels(tick_labels, fontsize=FONT['heatmap_tick'])
        ax.set_yticks(tick_positions)
        ax.set_yticklabels(tick_labels, fontsize=FONT['heatmap_tick'])
        if i // 2 == 1:
            ax.set_xlabel('Residue ID')
        else:
            ax.set_xticklabels([])
        if i % 2 == 0:
            ax.set_ylabel('Residue ID')
        else:
            ax.set_yticklabels([])
        title = system_display_name(sys_label)
        ax.set_title(title, fontsize=FONT['subtitle'])
        add_panel_label(ax, chr(65 + i), x=-0.08, y=1.10)
    fig.colorbar(im, ax=axes, shrink=0.6, label='Cross-correlation')
    save_figure(fig, outpath_base)
    print(f"Saved {outpath_base}.png/.svg")


def plot_s9_dccm_difference(dccm_dir, outpath_base):
    dccm_sssD_sap = np.load(os.path.join(dccm_dir, 'dccm_phe_sssD_sap.npy'))
    dccm_sssD_tsap = np.load(os.path.join(dccm_dir, 'dccm_phe_sssD_tsap.npy'))
    dccm_rrrL_sap = np.load(os.path.join(dccm_dir, 'dccm_phe_rrrL_sap.npy'))
    dccm_rrrL_tsap = np.load(os.path.join(dccm_dir, 'dccm_phe_rrrL_tsap.npy'))
    diff_sssD = dccm_sssD_sap - dccm_sssD_tsap
    diff_rrrL = dccm_rrrL_sap - dccm_rrrL_tsap
    vlim = max(np.abs(diff_sssD).max(), np.abs(diff_rrrL).max())
    fig = plt.figure(figsize=(6, 3))
    gs = GridSpec(1, 2, wspace=0.35)
    ax1 = fig.add_subplot(gs[0, 0])
    ax2 = fig.add_subplot(gs[0, 1])
    tick_positions = np.linspace(0, N_CA - 1, 6, dtype=int)
    tick_labels = [str(p + RESID_OFFSET) for p in tick_positions]
    im1 = ax1.imshow(diff_sssD, cmap='RdBu_r', vmin=-vlim, vmax=vlim,
                      interpolation='nearest', origin='lower', aspect='equal')
    rect1 = Rectangle((BINDING_SITE_START, BINDING_SITE_START),
                        BINDING_SITE_END - BINDING_SITE_START,
                        BINDING_SITE_END - BINDING_SITE_START,
                        linewidth=1.5, edgecolor='#F0E442', facecolor='none',
                        linestyle='--')
    ax1.add_patch(rect1)
    ax1.set_xticks(tick_positions)
    ax1.set_xticklabels(tick_labels, fontsize=FONT['heatmap_tick'])
    ax1.set_yticks(tick_positions)
    ax1.set_yticklabels(tick_labels, fontsize=FONT['heatmap_tick'])
    ax1.set_xlabel('Residue ID')
    ax1.set_ylabel('Residue ID')
    ax1.set_title('Phe-SSS-D (SAP \u2212 TSAP)', fontsize=FONT['subtitle'])
    add_panel_label(ax1, 'a', x=-0.08, y=1.10)
    im2 = ax2.imshow(diff_rrrL, cmap='RdBu_r', vmin=-vlim, vmax=vlim,
                      interpolation='nearest', origin='lower', aspect='equal')
    rect2 = Rectangle((BINDING_SITE_START, BINDING_SITE_START),
                        BINDING_SITE_END - BINDING_SITE_START,
                        BINDING_SITE_END - BINDING_SITE_START,
                        linewidth=1.5, edgecolor='#F0E442', facecolor='none',
                        linestyle='--')
    ax2.add_patch(rect2)
    ax2.set_xticks(tick_positions)
    ax2.set_xticklabels(tick_labels, fontsize=FONT['heatmap_tick'])
    ax2.set_yticks(tick_positions)
    ax2.set_yticklabels(tick_labels, fontsize=FONT['heatmap_tick'])
    ax2.set_xlabel('Residue ID')
    ax2.set_ylabel('Residue ID')
    ax2.set_title('Phe-RRR-L (SAP \u2212 TSAP)', fontsize=FONT['subtitle'])
    add_panel_label(ax2, 'b', x=-0.08, y=1.10)
    fig.colorbar(im1, ax=[ax1, ax2], shrink=0.8, label='\u0394Cross-correlation')
    save_figure(fig, outpath_base)
    print(f"Saved {outpath_base}.png/.svg")


def main():
    matplotlib.rcParams.update(RC_PARAMS)
    parser = argparse.ArgumentParser()
    parser.add_argument('--input-dir', default='pca_analysis/')
    parser.add_argument('--output-dir', default=None)
    args = parser.parse_args()
    in_dir = args.input_dir
    out_dir = args.output_dir if args.output_dir else in_dir

    proj_path = os.path.join(in_dir, 'projections_all.csv')
    scree_path = os.path.join(in_dir, 'scree_data.csv')
    mmpbsa_path = os.path.join(in_dir, 'mmpbsa_energies_all_systems.csv')
    contacts_path = os.path.join(in_dir, 'contacts_summary_all_systems.csv')
    aligned_path = os.path.join(in_dir, 'aligned_coords.npy')
    sysidx_path = os.path.join(in_dir, 'system_indices.npy')
    meta_path = os.path.join(in_dir, 'system_metadata.json')

    plot_phe_scatter(proj_path, scree_path, os.path.join(out_dir, 'fig_C_phe_sap_tsap_scatter'))
    plot_fes_4panel(proj_path, scree_path, os.path.join(out_dir, 'fig_D_fes_phe_key_systems'))
    plot_s6_all_fes_pairs(proj_path, os.path.join(out_dir, 'si_fig_S6_sap_tsap_fes'))
    plot_s7_binding_bar(mmpbsa_path, contacts_path,
                        os.path.join(out_dir, 'si_fig_S7_phe_binding_bar'))
    plot_s8_phe_dccm(aligned_path, sysidx_path, meta_path,
                     os.path.join(out_dir, 'si_fig_S8_phe_dccm'))
    plot_s9_dccm_difference(out_dir, os.path.join(out_dir, 'si_fig_S9_dccm_difference'))

    print("All 6 SAP/TSAP figures completed.")


if __name__ == '__main__':
    main()
