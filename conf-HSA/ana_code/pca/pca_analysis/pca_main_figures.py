import argparse
import os

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
import matplotlib.patches as mpatches
import numpy as np
np.random.seed(42)  # Reproducibility for stochastic analyses
import pandas as pd
from scipy import stats

FONT = {
    'tick': 9,
    'axis': 11,
    'panel': 12,
    'panel_weight': 'bold',
    'legend': 8,
    'annotation': 8,
    'annotation_bold': 10,
    'colorbar': 8,
    'heatmap_cell': 7,
    'heatmap_tick': 7,
    'subtitle': 10,
}

SAP_COLORS = {
    'me_sss': '#E69F00',
    'me_rrr': '#56B4E9',
    'phe_sss': '#009E73',
    'phe_rrr': '#CC79A7',
}
TSAP_COLORS = {
    'me_sss': '#D55E00',
    'me_rrr': '#0072B2',
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

PC_COLORS = ['#E69F00', '#56B4E9', '#009E73', '#CC79A7', '#D55E00']

PC_WIDTHS = [2.0, 1.5, 1.5, 1.5, 1.5]


def system_to_group(system_label):
    parts = system_label.split('_')
    return f"{parts[0]}_{parts[1][:3]}"


def system_to_ligand(system_label):
    return system_label.split('_')[-1]


def add_panel_label(ax, label, x=-0.08, y=1.10):
    ax.text(x, y, f"{label.upper()}", transform=ax.transAxes,
            fontsize=FONT['panel'], fontweight=FONT['panel_weight'], va='bottom')


def save_figure(fig, outpath_base):
    fig.savefig(outpath_base + '.png', dpi=DPI, bbox_inches='tight')
    fig.savefig(outpath_base + '.svg', format='svg', bbox_inches='tight')
    plt.close(fig)


def plot_regression_ci(ax, x, y):
    r_val, p_val = stats.pearsonr(x, y)
    m, b = np.polyfit(x, y, 1)
    x_pred = np.linspace(x.min(), x.max(), 100)
    y_pred = m * x_pred + b
    n = len(x)
    x_mean = x.mean()
    ss_x = np.sum((x - x_mean) ** 2)
    se_y = np.sqrt(np.sum((y - (m * x + b)) ** 2) / (n - 2))
    se_pred = se_y * np.sqrt(1 / n + (x_pred - x_mean) ** 2 / ss_x)
    t_crit = stats.t.ppf(0.975, n - 2)
    ax.fill_between(x_pred, y_pred - t_crit * se_pred,
                     y_pred + t_crit * se_pred, alpha=0.2, color='gray', zorder=1)
    ax.plot(x_pred, y_pred, 'k-', lw=1, alpha=0.7, zorder=2)
    return r_val, p_val


def plot_scatter_overlay(projections_path, metadata_path, scree_path, outpath_base):
    df = pd.read_csv(projections_path)
    scree = pd.read_csv(scree_path)
    pc1_frac = scree.iloc[0]['fraction']
    pc2_frac = scree.iloc[1]['fraction']
    df['group'] = df['system_label'].apply(system_to_group)
    df['ligand'] = df['system_label'].apply(system_to_ligand)
    fig = plt.figure(figsize=(6, 5))
    ax = fig.add_subplot(111)
    for group in GROUP_ORDER:
        grp_mask = df['group'] == group
        grp_df = df[grp_mask]
        for ligand in ['sap', 'tsap']:
            sub = grp_df[grp_df['ligand'] == ligand]
            if len(sub) == 0:
                continue
            color = SAP_COLORS[group] if ligand == 'sap' else TSAP_COLORS[group]
            ax.scatter(sub['PC1'], sub['PC2'],
                       c=color, marker='o',
                       s=4, alpha=0.3, rasterized=True, zorder=2)
    ax.set_xlabel(f"PC1 projection ({pc1_frac * 100:.1f}%) / \u00c5", fontsize=FONT['axis'])
    ax.set_ylabel(f"PC2 projection ({pc2_frac * 100:.1f}%) / \u00c5", fontsize=FONT['axis'])
    ax.tick_params(labelsize=FONT['tick'])
    legend_elements = []
    for group in GROUP_ORDER:
        legend_elements.append(mpatches.Patch(facecolor=SAP_COLORS[group],
                                              edgecolor='black', linewidth=0.5,
                                              label=f'{GROUP_LABELS[group]} sap'))
        legend_elements.append(mpatches.Patch(facecolor=TSAP_COLORS[group],
                                              edgecolor='black', linewidth=0.5,
                                              label=f'{GROUP_LABELS[group]} tsap'))
    ax.legend(handles=legend_elements, loc='upper right', fontsize=6,
              ncol=2, framealpha=0.9)
    add_panel_label(ax, 'a')
    fig.tight_layout()
    save_figure(fig, outpath_base)


def plot_residue_zoom(contrib_path, ax=None, outpath_base=None):
    df = pd.read_csv(contrib_path)
    df_zoom = df[(df['resid'] >= 370) & (df['resid'] <= 500)].copy()
    standalone = ax is None
    if standalone:
        fig = plt.figure(figsize=(6, 3.5))
        ax = fig.add_subplot(111)
    for i in range(5):
        pc_col = f'PC{i + 1}_frac'
        ax.plot(df_zoom['resid'], df_zoom[pc_col],
                color=PC_COLORS[i], label=f'PC{i + 1}',
                lw=PC_WIDTHS[i])
    ax.axvspan(377, 490, alpha=0.1, color='coral', zorder=0)
    binding_mask_full = (df['resid'] >= 377) & (df['resid'] <= 490)
    bs_frac = df.loc[binding_mask_full, 'PC1_frac'].sum()
    expected_uniform = binding_mask_full.sum() / len(df)
    ax.text(0.98, 0.98,
            f"Binding site: {bs_frac * 100:.1f}% of PC1\n"
            f"({expected_uniform * 100:.1f}% expected if uniform)",
            transform=ax.transAxes, fontsize=FONT['annotation_bold'],
            fontweight='bold', ha='right', va='top',
            bbox=dict(boxstyle='round,pad=0.3', fc='lightyellow', alpha=0.9))
    peak_labels = df_zoom.nlargest(4, 'PC1_frac').sort_values('resid')
    offsets_y = [0.0003, 0.0006, 0.0003, 0.0006]
    offsets_x = [0, 0, 3, 3]
    for idx, (_, row) in enumerate(peak_labels.iterrows()):
        ax.annotate(f"{row['resname']}{int(row['resid'])}",
                    xy=(row['resid'], row['PC1_frac']),
                    xytext=(row['resid'] + offsets_x[idx % 4],
                            row['PC1_frac'] + offsets_y[idx % 4]),
                    fontsize=FONT['annotation'], ha='center',
                    arrowprops=dict(arrowstyle='-', color='gray', lw=0.5))
    ax.set_xlabel("Residue number", fontsize=FONT['axis'])
    ax.set_ylabel("Fractional contribution", fontsize=FONT['axis'])
    ax.set_xlim(370, 500)
    ax.set_xticks(range(370, 501, 20))
    ax.tick_params(labelsize=FONT['tick'])
    ax.legend(fontsize=FONT['legend'], loc='upper left')
    add_panel_label(ax, 'a')
    if standalone:
        fig.tight_layout()
        if outpath_base:
            save_figure(fig, outpath_base)


def plot_dg_vs_pc1(mmpbsa_path, ax=None, outpath_base=None):
    df = pd.read_csv(mmpbsa_path)
    df['group'] = df['system'].apply(system_to_group)
    df['ligand'] = df['system'].apply(system_to_ligand)
    standalone = ax is None
    if standalone:
        fig = plt.figure(figsize=(6, 5))
        ax = fig.add_subplot(111)
    x_all = df['mean_PC1'].values
    y_all = df['TOTAL_corrected'].values
    for group in GROUP_ORDER:
        grp_df = df[df['group'] == group]
        for ligand in ['sap', 'tsap']:
            sub = grp_df[grp_df['ligand'] == ligand]
            if len(sub) == 0:
                continue
            color = SAP_COLORS[group] if ligand == 'sap' else TSAP_COLORS[group]
            ax.scatter(sub['mean_PC1'], sub['TOTAL_corrected'],
                       c=color, marker='o',
                       s=55, edgecolors='black', linewidths=0.5,
                       zorder=3)
    r_val, p_val = plot_regression_ci(ax, x_all, y_all)
    # Compute 95% CI via Fisher z-transform for annotation
    n = len(x_all)
    z = np.arctanh(r_val)
    se = 1 / np.sqrt(n - 3)
    z_lo = z - 1.96 * se
    z_hi = z + 1.96 * se
    r_lo, r_hi = np.tanh(z_lo), np.tanh(z_hi)
    ax.text(0.05, 0.95,
            f"Pearson r = {r_val:.2f}\n95% CI: [{r_lo:.2f}, {r_hi:.2f}]\np = {p_val:.2f}",
            transform=ax.transAxes, fontsize=FONT['annotation_bold'],
            fontweight='bold', va='top',
            bbox=dict(boxstyle='round,pad=0.3', fc='wheat', alpha=0.5))
    ax.set_xlabel("System-mean PC1 projection (\u00c5)", fontsize=FONT['axis'])
    ax.set_ylabel("\u0394G$_{corrected}$ (kJ mol$^{-1}$)", fontsize=FONT['axis'])
    ax.tick_params(labelsize=FONT['tick'])
    legend_elements = []
    for group in GROUP_ORDER:
        legend_elements.append(mpatches.Patch(facecolor=SAP_COLORS[group],
                                              edgecolor='black', linewidth=0.5,
                                              label=f'{GROUP_LABELS[group]} sap'))
        legend_elements.append(mpatches.Patch(facecolor=TSAP_COLORS[group],
                                              edgecolor='black', linewidth=0.5,
                                              label=f'{GROUP_LABELS[group]} tsap'))
    ax.legend(handles=legend_elements, fontsize=6,
              loc='lower right', framealpha=0.9, ncol=2)
    add_panel_label(ax, 'b')
    if standalone:
        fig.tight_layout()
        if outpath_base:
            save_figure(fig, outpath_base)


def plot_fig_b(contrib_path, mmpbsa_path, outpath_base):
    fig = plt.figure(figsize=(6, 7))
    gs = GridSpec(2, 1, height_ratios=[1, 1], hspace=0.35,
                  left=0.14, right=0.96, top=0.95, bottom=0.06)
    ax_a = fig.add_subplot(gs[0])
    plot_residue_zoom(contrib_path, ax=ax_a)
    ax_b = fig.add_subplot(gs[1])
    plot_dg_vs_pc1(mmpbsa_path, ax=ax_b)
    save_figure(fig, outpath_base)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--input-dir', default='pca_analysis/')
    parser.add_argument('--output-dir', default=None)
    args = parser.parse_args()
    input_dir = args.input_dir.rstrip('/')
    output_dir = args.output_dir.rstrip('/') if args.output_dir else input_dir
    plt.rcParams.update(RC_PARAMS)
    projections_path = os.path.join(input_dir, 'projections_all.csv')
    metadata_path = os.path.join(input_dir, 'system_metadata.json')
    scree_path = os.path.join(input_dir, 'scree_data.csv')
    contrib_path = os.path.join(input_dir, 'residue_contribution.csv')
    mmpbsa_path = os.path.join(input_dir, 'mmpbsa_energies_all_systems.csv')
    fig_a_base = os.path.join(output_dir, 'fig_A_pca_scatter')
    plot_scatter_overlay(projections_path, metadata_path, scree_path, fig_a_base)
    print(f"Saved: {fig_a_base}.png + .svg")
    fig_b_base = os.path.join(output_dir, 'fig_B_pca_binding')
    plot_fig_b(contrib_path, mmpbsa_path, fig_b_base)
    print(f"Saved: {fig_b_base}.png + .svg")


if __name__ == '__main__':
    main()
