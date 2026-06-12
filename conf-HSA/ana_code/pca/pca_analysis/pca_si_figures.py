import argparse
import json
import os
import sys
import warnings

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from matplotlib.patches import Rectangle
import matplotlib.lines as mlines
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

GROUP_COLORS = {
    'me_sss': '#E69F00',
    'me_rrr': '#56B4E9',
    'phe_sss': '#009E73',
    'phe_rrr': '#CC79A7',
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


def r_ci_fisher(r, n, alpha=0.05):
    z = np.arctanh(r)
    se = 1 / np.sqrt(n - 3)
    z_lo = z - stats.norm.ppf(1 - alpha / 2) * se
    z_hi = z + stats.norm.ppf(1 - alpha / 2) * se
    return np.tanh(z_lo), np.tanh(z_hi)


from pca_utils import R_GAS as kB, compute_fes  # Shared from pca_utils.py
from pca_utils import BINDING_SITE_RESID_RANGE, DCCM_BINDING_SITE_START, DCCM_BINDING_SITE_END
T = 300
kT = kB * T


def plot_si_s1_methodology(scree_path, projections_path, metadata_path, outpath_base):
    scree = pd.read_csv(scree_path)
    projections = pd.read_csv(projections_path)
    projections['group'] = projections['system_label'].apply(system_to_group)
    projections['ligand'] = projections['system_label'].apply(system_to_ligand)

    fig = plt.figure(figsize=(6, 8))
    gs = GridSpec(3, 2, height_ratios=[1.2, 1, 1], hspace=0.40, wspace=0.30,
                  left=0.10, right=0.88, top=0.95, bottom=0.05)

    scree50 = scree.head(50)
    ax = fig.add_subplot(gs[0, :])
    ax.bar(scree50['pc_number'], scree50['eigenvalue_angstrom2'],
           color='steelblue', width=0.8, alpha=0.7)
    ax.set_ylabel("Eigenvalue (\u00c5\u00b2)", fontsize=FONT['axis'])
    ax.set_xlabel("PC number", fontsize=FONT['axis'])
    ax.set_xlim(0, 51)
    ax.tick_params(labelsize=FONT['tick'])

    ax2 = ax.twinx()
    ax2.plot(scree50['pc_number'], scree50['cumulative_fraction'],
             color='coral', lw=1.5, zorder=3)
    ax2.set_ylabel("Cumulative fraction", fontsize=FONT['axis'])
    ax2.tick_params(labelsize=FONT['tick'])
    for yval in [0.5, 0.75, 0.9]:
        ax2.axhline(yval, color='gray', ls='--', lw=0.5, alpha=0.6)

    top3_cum = scree.iloc[2]['cumulative_fraction']
    ax.text(0.98, 0.95, f"Top 3 PCs = {top3_cum * 100:.1f}%",
            transform=ax.transAxes, fontsize=FONT['annotation'],
            ha='right', va='top')
    add_panel_label(ax, 'a')

    positions = [(1, 0), (1, 1), (2, 0), (2, 1)]
    for i, (group, (r, c)) in enumerate(zip(GROUP_ORDER, positions)):
        ax_g = fig.add_subplot(gs[r, c])
        grp_df = projections[projections['group'] == group]
        sap_df = grp_df[grp_df['ligand'] == 'sap']
        tsap_df = grp_df[grp_df['ligand'] == 'tsap']
        ax_g.scatter(sap_df['PC1'], sap_df['PC2'],
                     c=SAP_COLORS[group], s=2, alpha=0.2,
                     rasterized=True, marker='o')
        ax_g.scatter(tsap_df['PC1'], tsap_df['PC2'],
                     c=TSAP_COLORS[group], s=2, alpha=0.2,
                     rasterized=True, marker='o')
        ax_g.set_title(GROUP_LABELS[group], fontsize=FONT['subtitle'])
        ax_g.set_xlabel("PC1 projection (Å)", fontsize=FONT['axis'])
        ax_g.set_ylabel("PC2 projection (Å)", fontsize=FONT['axis'])
        ax_g.tick_params(labelsize=FONT['tick'])
        add_panel_label(ax_g, 'bcde'[i])

    save_figure(fig, outpath_base)


def plot_si_s2_fes(projections_path, outpath_base):
    projections = pd.read_csv(projections_path)
    projections['group'] = projections['system_label'].apply(system_to_group)

    fig = plt.figure(figsize=(6, 6))
    gs = GridSpec(2, 2, hspace=0.35, wspace=0.35,
                  left=0.10, right=0.88, top=0.93, bottom=0.07)

    for i, group in enumerate(GROUP_ORDER):
        ax = fig.add_subplot(gs[i // 2, i % 2])
        grp_df = projections[projections['group'] == group]
        G, xedges, yedges = compute_fes(grp_df['PC1'].values, grp_df['PC2'].values)
        xcenters = (xedges[:-1] + xedges[1:]) / 2
        ycenters = (yedges[:-1] + yedges[1:]) / 2

        G_masked = np.ma.masked_invalid(G)
        im = ax.pcolormesh(xedges, yedges, G_masked, cmap='cubehelix', shading='auto')
        try:
            ax.contour(xcenters, ycenters, G_masked.filled(0), levels=8,
                       colors='white', linewidths=0.3, alpha=0.5)
        except ValueError:
            pass
        cb = plt.colorbar(im, ax=ax, shrink=0.8)
        cb.set_label('G (kJ/mol)', fontsize=FONT['colorbar'])
        cb.ax.tick_params(labelsize=FONT['colorbar'])
        ax.set_title(GROUP_LABELS[group], fontsize=FONT['subtitle'])
        ax.set_xlabel("PC1 projection (Å)", fontsize=FONT['axis'])
        ax.set_ylabel("PC2 projection (Å)", fontsize=FONT['axis'])
        ax.tick_params(labelsize=FONT['tick'])
        add_panel_label(ax, 'abcd'[i])

    save_figure(fig, outpath_base)


def plot_si_s3_residue(contrib_path, scree_path, outpath_base):
    df = pd.read_csv(contrib_path)
    scree = pd.read_csv(scree_path)
    pc_fracs = [f"{scree.iloc[i]['fraction'] * 100:.1f}" for i in range(5)]

    fig = plt.figure(figsize=(6, 8))
    gs = GridSpec(3, 1, height_ratios=[1.5, 1, 0.7], hspace=0.35,
                  left=0.12, right=0.96, top=0.95, bottom=0.05)

    ax_a = fig.add_subplot(gs[0])
    for i in range(5):
        ax_a.plot(df['resid'], df[f'PC{i+1}_frac'],
                  color=PC_COLORS[i],
                  label=f'PC{i+1} ({pc_fracs[i]}%)', lw=1)
    ax_a.axvspan(BINDING_SITE_RESID_RANGE[0], BINDING_SITE_RESID_RANGE[1], alpha=0.08, color='gray', zorder=0)
    top10 = df.nlargest(10, 'PC1_frac')
    for _, row in top10.iterrows():
        ax_a.annotate(f"{row['resname']}{int(row['resid'])}",
                      xy=(row['resid'], row['PC1_frac']),
                      fontsize=FONT['annotation'], rotation=90,
                      ha='center', va='bottom')
    ax_a.legend(fontsize=FONT['legend'], loc='upper right', ncol=2)
    ax_a.set_xlabel("Residue number", fontsize=FONT['axis'])
    ax_a.set_ylabel("Fractional contribution", fontsize=FONT['axis'])
    ax_a.set_xticks(range(0, 601, 50))
    ax_a.tick_params(labelsize=FONT['tick'])
    add_panel_label(ax_a, 'a')

    ax_b = fig.add_subplot(gs[1])
    df_zoom = df[(df['resid'] >= 370) & (df['resid'] <= 500)].copy()
    for i in range(5):
        ax_b.plot(df_zoom['resid'], df_zoom[f'PC{i+1}_frac'],
                  color=PC_COLORS[i], label=f'PC{i+1}', lw=1.2)
    ax_b.axvspan(BINDING_SITE_RESID_RANGE[0], BINDING_SITE_RESID_RANGE[1], alpha=0.1, color='coral', zorder=0)
    binding_mask = (df['resid'] >= BINDING_SITE_RESID_RANGE[0]) & (df['resid'] <= BINDING_SITE_RESID_RANGE[1])
    bs_frac = df.loc[binding_mask, 'PC1_frac'].sum()
    ax_b.text(0.98, 0.98, f"Binding site: {bs_frac*100:.1f}% of PC1",
              transform=ax_b.transAxes, fontsize=FONT['annotation_bold'],
              fontweight='bold', ha='right', va='top',
              bbox=dict(boxstyle='round,pad=0.3', fc='lightyellow', alpha=0.9))
    top_pc1 = df_zoom.nlargest(4, 'PC1_frac')
    for _, row in top_pc1.iterrows():
        ax_b.annotate(f"{row['resname']}{int(row['resid'])}",
                      xy=(row['resid'], row['PC1_frac']),
                      xytext=(row['resid'], row['PC1_frac'] + 0.0003),
                      fontsize=FONT['annotation'],
                      arrowprops=dict(arrowstyle='-', lw=0.5),
                      ha='center', va='bottom')
    ax_b.set_xlabel("Residue number", fontsize=FONT['axis'])
    ax_b.set_ylabel("Fractional contribution", fontsize=FONT['axis'])
    ax_b.set_xlim(370, 500)
    ax_b.set_xticks(range(370, 501, 20))
    ax_b.tick_params(labelsize=FONT['tick'])
    ax_b.legend(fontsize=FONT['legend'], loc='upper left')
    add_panel_label(ax_b, 'b')

    ax_c = fig.add_subplot(gs[2])
    binding_mask = (df['resid'] >= BINDING_SITE_RESID_RANGE[0]) & (df['resid'] <= BINDING_SITE_RESID_RANGE[1])
    bs_frac = df.loc[binding_mask, 'PC1_frac'].sum()
    ax_c.pie([bs_frac, 1 - bs_frac],
             labels=['Binding site', 'Scaffold'],
             colors=['coral', 'lightgray'],
             autopct='%1.1f%%',
             textprops={'fontsize': FONT['legend']},
             startangle=90)
    add_panel_label(ax_c, 'c', x=-0.1, y=1.05)

    save_figure(fig, outpath_base)


def plot_si_s4_dccm(input_dir, stats_path, outpath_base):
    stats_df = pd.read_csv(stats_path)

    fig = plt.figure(figsize=(6, 7))
    gs = GridSpec(3, 2, height_ratios=[1, 1, 0.6], hspace=0.35, wspace=0.25,
                  left=0.08, right=0.80, top=0.95, bottom=0.06)

    last_im = None
    for i, group in enumerate(GROUP_ORDER):
        ax = fig.add_subplot(gs[i // 2, i % 2])
        dccm = np.load(os.path.join(input_dir, f'dccm_{group}.npy'))
        last_im = ax.imshow(dccm, cmap='RdBu_r', vmin=-1, vmax=1,
                            aspect='auto', origin='lower')
        rect = Rectangle((DCCM_BINDING_SITE_START, DCCM_BINDING_SITE_START),
                         DCCM_BINDING_SITE_END - DCCM_BINDING_SITE_START,
                         DCCM_BINDING_SITE_END - DCCM_BINDING_SITE_START,
                         fill=False, ec='yellow', ls='--', lw=0.8)
        ax.add_patch(rect)
        ax.set_title(GROUP_LABELS[group], fontsize=FONT['subtitle'])
        tick_pos = [0, 100, 200, 300, 400, 500]
        ax.set_xticks(tick_pos)
        ax.set_yticks(tick_pos)
        ax.tick_params(labelsize=FONT['heatmap_tick'])
        add_panel_label(ax, 'abcd'[i])

    cax = fig.add_axes([0.81, 0.38, 0.012, 0.55])
    cb = fig.colorbar(last_im, cax=cax)
    cb.set_label('Correlation', fontsize=FONT['colorbar'])
    cb.ax.tick_params(labelsize=FONT['colorbar'])

    ax_e = fig.add_subplot(gs[2, :])
    all_pairs = stats_df.copy()
    labels = all_pairs['pair'].values
    r_vals = all_pairs['pearson_r'].values
    n_inter = 6
    colors = ['steelblue'] * n_inter + ['coral'] * (len(labels) - n_inter)
    y_pos = np.arange(len(labels))
    ax_e.barh(y_pos, r_vals, color=colors, height=0.6)
    ax_e.set_yticks(y_pos)
    ax_e.set_yticklabels(labels, fontsize=FONT['tick'])
    ax_e.set_xlabel("Pearson r", fontsize=FONT['axis'])
    ax_e.axvline(0.8, ls='--', color='gray', lw=0.5)
    ax_e.axvline(0.9, ls='--', color='gray', lw=0.5)
    ax_e.set_xlim(0.7, 1.0)
    ax_e.tick_params(labelsize=FONT['tick'])
    inter_rs = stats_df[~stats_df['pair'].str.contains('global')]['pearson_r']
    global_rs = stats_df[stats_df['pair'].str.contains('global')]['pearson_r']
    ax_e.text(0.98, 0.02,
              f"Inter-group: r = {inter_rs.min():.2f}\u2013{inter_rs.max():.2f}\nvs global: r = {global_rs.min():.2f}\u2013{global_rs.max():.2f}",
              transform=ax_e.transAxes, fontsize=FONT['annotation'],
              ha='right', va='bottom')
    add_panel_label(ax_e, 'e')

    save_figure(fig, outpath_base)


def plot_si_s5a_binding_scatters(mmpbsa_path, contacts_path, outpath_base):
    mmpbsa = pd.read_csv(mmpbsa_path)
    contacts = pd.read_csv(contacts_path)
    mmpbsa['group'] = mmpbsa['system'].apply(system_to_group)
    mmpbsa['ligand'] = mmpbsa['system'].apply(system_to_ligand)
    contacts['group'] = contacts['system'].apply(system_to_group)
    contacts['ligand'] = contacts['system'].apply(system_to_ligand)

    fig = plt.figure(figsize=(6, 8))
    gs = GridSpec(3, 2, hspace=0.35, wspace=0.30,
                  left=0.14, right=0.96, top=0.95, bottom=0.05)

    panels = [
        (0, 0, mmpbsa, 'mean_PC1', 'TOTAL_corrected', '\u0394G (kJ mol\u207b\u00b9)'),
        (0, 1, mmpbsa, 'mean_PC2', 'TOTAL_corrected', '\u0394G (kJ mol\u207b\u00b9)'),
        (1, 0, mmpbsa, 'mean_PC3', 'TOTAL_corrected', '\u0394G (kJ mol\u207b\u00b9)'),
        (1, 1, contacts, 'mean_PC1', 'mean_total_contacts', 'Contacts'),
        (2, 0, contacts, 'mean_PC2', 'mean_total_contacts', 'Contacts'),
        (2, 1, contacts, 'mean_PC3', 'mean_total_contacts', 'Contacts'),
    ]

    _panel_r, _panel_p, _panel_ax = [], [], []
    for idx, (row, col, data, xcol, ycol, ylabel) in enumerate(panels):
        ax = fig.add_subplot(gs[row, col])
        x_all = data[xcol].values
        y_all = data[ycol].values
        for group in GROUP_ORDER:
            grp_df = data[data['group'] == group]
            for ligand in ['sap', 'tsap']:
                sub = grp_df[grp_df['ligand'] == ligand]
                if len(sub) == 0:
                    continue
                ax.scatter(sub[xcol], sub[ycol],
                           c=SAP_COLORS[group] if ligand == 'sap' else TSAP_COLORS[group],
                           marker='o',
                           s=40, edgecolors='black', linewidths=0.3,
                           zorder=3)
        r_val, p_val = plot_regression_ci(ax, x_all, y_all)
        # Store for FDR-BH correction after all panels computed
        _panel_r.append(r_val)
        _panel_p.append(p_val)
        _panel_ax.append(ax)
        pc_num = xcol.replace('mean_', '')
        ax.set_xlabel(pc_num, fontsize=FONT['axis'])
        ax.set_ylabel(ylabel, fontsize=FONT['axis'])
        ax.tick_params(labelsize=8)
        add_panel_label(ax, 'abcdef'[idx])

    # Apply FDR-BH multiple testing correction to 6 scatter panels (S5 fix)
    try:
        from statsmodels.stats.multitest import multipletests
        _, p_adj, _, _ = multipletests(_panel_p, method='fdr_bh')
        for ax, r_val, p_raw, p_corrected in zip(_panel_ax, _panel_r, _panel_p, p_adj):
            sig_mark = '*' if p_corrected < 0.05 else ''
            ax.text(0.05, 0.92,
                    f"r = {r_val:.2f}, p = {p_raw:.2f}\np_adj = {p_corrected:.3f}{sig_mark}",
                    transform=ax.transAxes, fontsize=FONT['annotation'])
    except ImportError:
        bonf = 0.05 / len(_panel_p)
        for ax, r_val, p_val in zip(_panel_ax, _panel_r, _panel_p):
            sig_mark = '*' if p_val < bonf else ''
            ax.text(0.05, 0.92,
                    f"r = {r_val:.2f}, p = {p_val:.2f}{sig_mark}\n(Bonferroni α={bonf:.4f})",
                    transform=ax.transAxes, fontsize=FONT['annotation'])

    save_figure(fig, outpath_base)


def plot_si_s5b_forest_heatmap(corr_path, all_ana_dir, projections_path,
                               metadata_path, outpath_base, skip_perframe=False):
    corr_df = pd.read_csv(corr_path)

    fig = plt.figure(figsize=(6, 5))
    gs = GridSpec(1, 2, width_ratios=[1, 1.2], wspace=0.35,
                  left=0.14, right=0.86, top=0.88, bottom=0.10)

    ax_a = fig.add_subplot(gs[0])
    row_labels = []
    for _, row_data in corr_df.iterrows():
        metric = '\u0394G' if row_data['metric'] == 'TOTAL_corrected' else 'Contacts'
        pc = row_data['PC'].replace('mean_', '')
        row_labels.append(f"{metric} vs {pc}")

    r_vals = corr_df['pearson_r'].values
    n = 16
    y_pos = np.arange(len(r_vals))

    for i, (r, label) in enumerate(zip(r_vals, row_labels)):
        lo, hi = r_ci_fisher(r, n)
        ax_a.hlines(y=i, xmin=lo, xmax=hi, color='#4C72B0', linewidth=2, zorder=2)
        ax_a.scatter(r, i, c='#4C72B0', s=50, zorder=3,
                     edgecolors='black', linewidths=0.5)

    ax_a.axvline(0, color='black', lw=0.8, ls='--', zorder=1)
    ax_a.axhline(2.5, color='gray', lw=0.5, ls=':')
    ax_a.set_yticks(range(6))
    ax_a.set_yticklabels(row_labels, fontsize=FONT['tick'])
    ax_a.set_xlabel("Pearson r (95% CI)", fontsize=FONT['axis'])
    ax_a.set_xlim(-0.8, 0.9)
    ax_a.set_title("PC\u2013Binding Correlations\n(N = 16)",
                   fontsize=FONT['subtitle'], fontweight='bold')
    ax_a.tick_params(labelsize=FONT['tick'])
    add_panel_label(ax_a, 'a')

    ax_b = fig.add_subplot(gs[1])
    if not skip_perframe:
        projections = pd.read_csv(projections_path)
        with open(metadata_path) as f:
            metadata = json.load(f)

        system_names = [s['name'] for s in metadata]
        heatmap_data = np.full((16, 3), np.nan)

        for si, sys_name in enumerate(system_names):
            contacts_list = []
            for trial_idx in range(4):
                cpath = os.path.join(all_ana_dir, 'per_ligand', sys_name,
                                     'per_trial', f'mmpbsa_{trial_idx}',
                                     'contacts', 'contacts_total.csv')
                if os.path.exists(cpath):
                    cdf = pd.read_csv(cpath)
                    contacts_list.append(cdf['total_contacts'].values)

            if len(contacts_list) == 0:
                continue

            # Validate per-trial frame alignment (C2 escalation: per-trial count check)
            sys_meta = next((m for m in metadata if m["name"] == sys_name), None)
            n_frames_total = sys_meta["n_frames"] if sys_meta else 0
            n_trials = len(contacts_list)
            frames_per_trial = n_frames_total // n_trials if n_trials > 0 else 0

            for t, contacts_t in enumerate(contacts_list):
                n_contacts = len(contacts_t)
                if n_contacts == 0:
                    warnings.warn(f"Trial {t} for {sys_name}: empty contacts array")
                elif n_contacts != frames_per_trial and frames_per_trial > 0:
                    warnings.warn(
                        f"Trial {t} for {sys_name}: contacts frames ({n_contacts}) "
                        f"!= expected projection frames ({frames_per_trial}). "
                        f"Alignment may be offset."
                    )

            # NOTE: This concatenation assumes that contacts from trials 0-3 are in the
            # same frame order as the projections from the concatenated trajectory.
            # Per-trial frame count validation has been added above, but the
            # concatenation itself still assumes frame-level alignment within
            # each trial. If equilibration frames were removed from contacts but
            # not projections (or vice versa), correlations will be spurious.
            all_contacts = np.concatenate(contacts_list)
            sys_proj = projections[projections['system_label'] == sys_name]
            n_match = min(len(sys_proj), len(all_contacts))

            if n_match < 10:
                continue

            for pci, pc in enumerate(['PC1', 'PC2', 'PC3']):
                pc_vals = sys_proj[pc].values[:n_match]
                c_vals = all_contacts[:n_match]
                if np.std(pc_vals) > 0 and np.std(c_vals) > 0:
                    r, _ = stats.pearsonr(pc_vals, c_vals)
                    heatmap_data[si, pci] = r

        im = ax_b.imshow(heatmap_data, cmap='RdBu_r', vmin=-1, vmax=1,
                         aspect='auto')
        for i in range(16):
            for j in range(3):
                val = heatmap_data[i, j]
                if not np.isnan(val):
                    color = 'white' if abs(val) > 0.5 else 'black'
                    ax_b.text(j, i, f"{val:.2f}", ha='center', va='center',
                              fontsize=FONT['heatmap_cell'], color=color)
        ax_b.set_xticks([0, 1, 2])
        ax_b.set_xticklabels(['PC1', 'PC2', 'PC3'], fontsize=FONT['tick'])
        ax_b.set_yticks(range(16))
        short_names = [n.replace('_', '-') for n in system_names]
        ax_b.set_yticklabels(short_names, fontsize=FONT['heatmap_tick'])
        cb = plt.colorbar(im, ax=ax_b, label='Pearson r', shrink=0.7, pad=0.02)
        cb.ax.tick_params(labelsize=FONT['colorbar'])
        ax_b.set_title("Per-Frame Contacts\u2013PC Correlations",
                       fontsize=FONT['subtitle'])
    else:
        ax_b.text(0.5, 0.5, "Per-frame data\nskipped",
                  transform=ax_b.transAxes, ha='center', va='center',
                  fontsize=FONT['annotation'])

    add_panel_label(ax_b, 'b')

    save_figure(fig, outpath_base)


def write_si_table_t1(corr_path, dccm_stats_path, outpath_base):
    corr_df = pd.read_csv(corr_path)
    dccm_df = pd.read_csv(dccm_stats_path)

    rows = []
    for _, row_data in corr_df.iterrows():
        metric = '\u0394G' if row_data['metric'] == 'TOTAL_corrected' else 'Contacts'
        pc = row_data['PC'].replace('mean_', '')
        lo, hi = r_ci_fisher(row_data['pearson_r'], 16)
        rows.append({
            'Metric': metric,
            'PC': pc,
            'Pearson r': f"{row_data['pearson_r']:.2f}",
            'p-value': f"{row_data['pearson_p']:.3f}",
            '95% CI lower': f"{lo:.2f}",
            '95% CI upper': f"{hi:.2f}",
            'Spearman \u03c1': f"{row_data['spearman_rho']:.2f}",
            'Spearman p': f"{row_data['spearman_p']:.3f}",
            'Partial r': f"{row_data['partial_r']:.2f}",
            'Partial p': f"{row_data['partial_p']:.3f}",
        })

    out_df = pd.DataFrame(rows)
    out_df = out_df.rename(columns={
        '95% CI lower': r'95\% CI lower',
        '95% CI upper': r'95\% CI upper',
    })
    out_df.to_csv(outpath_base + '.csv', index=False)

    with open(outpath_base + '.tex', 'w') as f:
        f.write("\\begin{table}[htbp]\n\\centering\n")
        f.write("\\caption{PCA--Binding Correlation Statistics}\n")
        f.write("\\label{tab:correlations}\n")
        f.write(out_df.to_latex(index=False, escape=False,
                                column_format='llrrrrrrrr'))
        f.write("\n\\end{table}\n")

    dccm_rows = []
    for _, row_data in dccm_df.iterrows():
        dccm_rows.append({
            'Pair': row_data['pair'],
            'Pearson r': f"{row_data['pearson_r']:.3f}",
            'RMSE': f"{row_data['rmse']:.3f}",
        })

    dccm_out = pd.DataFrame(dccm_rows)
    dccm_out['Pair'] = dccm_out['Pair'].str.replace('_', r'\_')
    dccm_tex_path = outpath_base.replace('T1_correlations', 'T1_dccm_stats')
    with open(dccm_tex_path + '.tex', 'w') as f:
        f.write("\\begin{table}[htbp]\n\\centering\n")
        f.write("\\caption{DCCM Inter-Group Correlation Statistics}\n")
        f.write("\\label{tab:dccm_stats}\n")
        f.write(dccm_out.to_latex(index=False, escape=False,
                                   column_format='lrr'))
        f.write("\n\\end{table}\n")


def plot_si_s10_component_scatters(mmpbsa_path, contacts_path, projections_path, outpath_base):
    """SI figure: 2×3 grid of metric–PC scatter plots with regression CIs and p-values.

    Shows TOTAL_corrected and total_contacts vs PC1, PC2, PC3.
    Each panel includes regression line, 95% CI band, and r/p annotation.

    Note: The original plan specified a 4×3 grid (VDWAALS, GGAS, GSOLV, TOTAL_corrected
    × 3 PCs), but component-level energy breakdowns are not available in the aggregated
    MMPBSA data. Adapted to show the two available metrics (ΔG corrected and contacts)
    × 3 PCs, which directly supports the LOO correlation analysis (RQ-RIG3).
    """
    mmpbsa = pd.read_csv(mmpbsa_path)
    contacts = pd.read_csv(contacts_path)
    projections = pd.read_csv(projections_path)

    # System-level means
    sys_proj = projections.groupby('system_label')[['PC1', 'PC2', 'PC3']].mean()
    systems = sorted(sys_proj.index)

    metrics = [
        ('TOTAL_corrected', mmpbsa, 'ΔG corrected (kJ mol⁻¹)'),
        ('mean_total_contacts', contacts, 'Mean total contacts'),
    ]
    pcs = ['PC1', 'PC2', 'PC3']

    fig, axes = plt.subplots(2, 3, figsize=(7.5, 5), dpi=DPI)

    for i, (metric_col, data_df, ylabel) in enumerate(metrics):
        for j, pc in enumerate(pcs):
            ax = axes[i, j]

            # Get data
            metric_vals = data_df.set_index('system').reindex(systems)[metric_col]
            pc_vals = sys_proj.reindex(systems)[pc]
            valid = ~(metric_vals.isna() | pc_vals.isna())
            x = pc_vals[valid].values
            y = metric_vals[valid].values

            if len(x) < 3:
                ax.text(0.5, 0.5, 'Insufficient data', transform=ax.transAxes,
                        ha='center', fontsize=FONT['annotation'])
                continue

            # Scatter with group coloring
            for group in GROUP_ORDER:
                grp_systems = [s for s in systems if system_to_group(s) == group]
                for ligand in ['sap', 'tsap']:
                    sub_systems = [s for s in grp_systems if system_to_ligand(s) == ligand]
                    if not sub_systems:
                        continue
                    sub_x = pc_vals.reindex(sub_systems).dropna().values
                    sub_y = metric_vals.reindex(sub_systems).dropna().values
                    # Only plot where both are valid
                    common_idx = list(set(sub_systems) & set(metric_vals.dropna().index) & set(pc_vals.dropna().index))
                    if len(common_idx) == 0:
                        continue
                    sub_x = pc_vals.reindex(common_idx).values
                    sub_y = metric_vals.reindex(common_idx).values
                    color = SAP_COLORS[group] if ligand == 'sap' else TSAP_COLORS[group]
                    ax.scatter(sub_x, sub_y, s=40, c=color, edgecolors='black',
                               linewidths=0.5, zorder=3)

            # Regression + CI
            r_val, p_val = plot_regression_ci(ax, x, y)

            # Annotation
            sig = '*' if p_val < 0.05 else ''
            ax.text(0.05, 0.95, f"r={r_val:.2f}, p={p_val:.3f}{sig}",
                    transform=ax.transAxes, fontsize=FONT['annotation'],
                    va='top',
                    bbox=dict(boxstyle='round,pad=0.2', fc='wheat', alpha=0.5))

            # Axis labels
            if i == 1:  # Bottom row
                ax.set_xlabel(f'{pc} (Å)', fontsize=FONT['axis'])
            if j == 0:  # Left column
                ax.set_ylabel(ylabel, fontsize=FONT['annotation'])

            ax.tick_params(labelsize=FONT['heatmap_tick'])

    fig.tight_layout(h_pad=1.5, w_pad=1.0)
    save_figure(fig, outpath_base)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--input-dir', default='pca_analysis/')
    parser.add_argument('--all-ana-dir', default='com_md/all_ana/')
    parser.add_argument('--output-dir', default=None)
    parser.add_argument('--skip-perframe', action='store_true', default=False)
    args = parser.parse_args()

    input_dir = args.input_dir.rstrip('/')
    output_dir = args.output_dir.rstrip('/') if args.output_dir else input_dir
    all_ana_dir = args.all_ana_dir.rstrip('/')

    plt.rcParams.update(RC_PARAMS)

    projections_path = os.path.join(input_dir, 'projections_all.csv')
    metadata_path = os.path.join(input_dir, 'system_metadata.json')
    scree_path = os.path.join(input_dir, 'scree_data.csv')
    contrib_path = os.path.join(input_dir, 'residue_contribution.csv')
    mmpbsa_path = os.path.join(input_dir, 'mmpbsa_energies_all_systems.csv')
    contacts_path = os.path.join(input_dir, 'contacts_summary_all_systems.csv')
    corr_path = os.path.join(input_dir, 'binding_correlation_table.csv')
    dccm_stats_path = os.path.join(input_dir, 'dccm_correlation_stats.csv')

    print("Generating SI Fig S1: PCA methodology...")
    plot_si_s1_methodology(scree_path, projections_path, metadata_path,
                           os.path.join(output_dir, 'si_fig_S1_methodology'))

    print("Generating SI Fig S2: FES landscapes...")
    plot_si_s2_fes(projections_path,
                   os.path.join(output_dir, 'si_fig_S2_fes'))

    print("Generating SI Fig S3: Residue detail...")
    plot_si_s3_residue(contrib_path, scree_path,
                       os.path.join(output_dir, 'si_fig_S3_residue'))

    print("Generating SI Fig S4: DCCM...")
    plot_si_s4_dccm(input_dir, dccm_stats_path,
                    os.path.join(output_dir, 'si_fig_S4_dccm'))

    print("Generating SI Fig S5a: Binding scatters...")
    plot_si_s5a_binding_scatters(mmpbsa_path, contacts_path,
                                 os.path.join(output_dir, 'si_fig_S5a_binding_scatters'))

    print("Generating SI Fig S5b: Forest + heatmap...")
    plot_si_s5b_forest_heatmap(corr_path, all_ana_dir, projections_path,
                               metadata_path,
                               os.path.join(output_dir, 'si_fig_S5b_forest_heatmap'),
                               skip_perframe=args.skip_perframe)

    print("Generating SI Table T1...")
    write_si_table_t1(corr_path, dccm_stats_path,
                      os.path.join(output_dir, 'si_table_T1_correlations'))

    # SI Fig S10: Component-level metric-PC scatters (RQ-RIG4)
    print("Generating SI Fig S10: Component-level metric-PC scatters...")
    plot_si_s10_component_scatters(mmpbsa_path, contacts_path, projections_path,
                                  os.path.join(output_dir, 'si_fig_S10_component_scatters'))

    print("All SI figures and tables generated.")


if __name__ == '__main__':
    main()
