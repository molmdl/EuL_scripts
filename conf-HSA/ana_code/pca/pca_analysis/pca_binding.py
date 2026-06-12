import argparse
import json
import os
import sys
import warnings
from pathlib import Path

import numpy as np
np.random.seed(42)  # Reproducibility for stochastic analyses
import pandas as pd
from scipy import stats
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.lines import Line2D

SYSTEMS = [
    "me_sssL_sap", "me_sssL_tsap", "me_sssD_sap", "me_sssD_tsap",
    "me_rrrL_sap", "me_rrrL_tsap", "me_rrrD_sap", "me_rrrD_tsap",
    "phe_sssL_sap", "phe_sssL_tsap", "phe_sssD_sap", "phe_sssD_tsap",
    "phe_rrrL_sap", "phe_rrrL_tsap", "phe_rrrD_sap", "phe_rrrD_tsap",
]

GROUP_COLORS = {
    "me_sss": "#E69F00",
    "me_rrr": "#56B4E9",
    "phe_sss": "#009E73",
    "phe_rrr": "#CC79A7",
}

from pca_utils import BINDING_SITE_CONTACT_RESIDUES as BINDING_SITE_RESIDUES

KCAL_TO_KJ = 4.184
# CORRECTION = 0.33: Empirical scaling factor for MMPBSA TOTAL energies.
# Applied because raw MMPBSA TOTAL systematically overestimates binding free
# energies for this system. The factor accounts for entropy and other
# contributions not captured by the MM/PBSA protocol. See INTERPRETATION_v3.md §9.5.
# If a published reference for this factor is available, it should be cited here.
CORRECTION = 0.33
CONVERSION = KCAL_TO_KJ * CORRECTION


def parse_averaged_energies(path):
    result = {"avg": {}, "sd": {}, "sem": {}}
    energy_terms = ["VDWAALS", "EEL", "GGAS", "GSOLV", "TOTAL"]
    n_valid = 0

    with open(path) as f:
        lines = f.readlines()

    for line in lines:
        parts = line.split()
        if len(parts) >= 6 and parts[0] in energy_terms:
            result["avg"][parts[0]] = float(parts[1])
            result["sd"][parts[0]] = float(parts[2])
            result["sem"][parts[0]] = float(parts[3])

    in_per_trial = False
    in_sd = False
    for line in lines:
        if "Per-trial statistics:" in line:
            in_per_trial = True
            in_sd = False
            continue
        if "Per-trial SD" in line:
            in_sd = True
            continue
        if "=====" in line:
            in_per_trial = False
            in_sd = False
            continue
        if in_per_trial and not in_sd:
            parts = line.split()
            if len(parts) >= 6 and parts[0].startswith("mmpbsa_"):
                try:
                    float(parts[5])
                    n_valid += 1
                except ValueError:
                    pass

    result["n_valid_trials"] = n_valid
    return result


def load_all_mmpbsa(all_ana_dir, systems):
    rows = []
    for sys_name in systems:
        path = Path(all_ana_dir) / sys_name / "mmpbsa" / "averaged_energies.txt"
        if not path.exists():
            continue
        parsed = parse_averaged_energies(str(path))
        total_kcal = parsed["avg"]["TOTAL"]
        rows.append({
            "system": sys_name,
            "TOTAL_kjmol": total_kcal * KCAL_TO_KJ,
            "TOTAL_corrected": total_kcal * CONVERSION,
            "SD": parsed["sd"]["TOTAL"],
            "SEM": parsed["sem"]["TOTAL"],
            "n_trials": parsed["n_valid_trials"],
        })
    return pd.DataFrame(rows)


def load_system_contacts(all_ana_dir, system):
    base = Path(all_ana_dir) / system / "per_trial"
    trial_totals = []
    trial_bs_fracs = []

    for t in range(4):
        total_path = base / f"mmpbsa_{t}" / "contacts" / "contacts_total.csv"
        summary_path = base / f"mmpbsa_{t}" / "contacts" / "contacts_summary.csv"

        if total_path.exists():
            df_t = pd.read_csv(total_path)
            trial_totals.append(df_t["total_contacts"].mean())

        if summary_path.exists():
            df_s = pd.read_csv(summary_path)
            bs_mask = df_s["residue_id"].isin(BINDING_SITE_RESIDUES)
            bs_sum = df_s.loc[bs_mask, "mean_contacts"].sum()
            all_sum = df_s["mean_contacts"].sum()
            if all_sum > 0:
                trial_bs_fracs.append(bs_sum / all_sum)

    return {
        "mean_total_contacts": np.mean(trial_totals) if trial_totals else np.nan,
        "std_total_contacts": np.std(trial_totals, ddof=1) if len(trial_totals) > 1 else 0.0,
        "bindingsite_frac": np.mean(trial_bs_fracs) if trial_bs_fracs else np.nan,
    }


def load_all_contacts(all_ana_dir, systems):
    rows = []
    for sys_name in systems:
        c = load_system_contacts(all_ana_dir, sys_name)
        rows.append({
            "system": sys_name,
            "mean_total_contacts": c["mean_total_contacts"],
            "std_total_contacts": c["std_total_contacts"],
            "bindingsite_frac": c["bindingsite_frac"],
        })
    return pd.DataFrame(rows)


def compute_system_mean_pc(projections_path, system_indices_path, metadata_path):
    proj = pd.read_csv(projections_path)
    pc_means = proj.groupby("system_label")[["PC1", "PC2", "PC3"]].mean().reset_index()
    pc_means.columns = ["system", "mean_PC1", "mean_PC2", "mean_PC3"]
    return pc_means


def compute_partial_corr(x, y, z):
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    z = np.asarray(z, dtype=float)
    Z = np.column_stack([np.ones(len(z)), z])
    beta_x = np.linalg.lstsq(Z, x, rcond=None)[0]
    res_x = x - Z @ beta_x
    beta_y = np.linalg.lstsq(Z, y, rcond=None)[0]
    res_y = y - Z @ beta_y
    r, p = stats.pearsonr(res_x, res_y)
    dfree = len(x) - 3
    if dfree > 0:
        t_stat = r * np.sqrt(dfree / (1 - r**2 + 1e-12))
        p = 2 * stats.t.sf(np.abs(t_stat), dfree)
    else:
        p = np.nan
    return r, p


def compute_correlations(df, x_cols, y_cols, chirality):
    rows = []
    for y_col in y_cols:
        for x_col in x_cols:
            x = df[x_col].values.astype(float)
            y = df[y_col].values.astype(float)
            z = chirality.astype(float)
            mask = ~(np.isnan(x) | np.isnan(y))
            xc, yc, zc = x[mask], y[mask], z[mask]
            pr, pp = stats.pearsonr(xc, yc)
            sr, sp = stats.spearmanr(xc, yc)
            partial_r, partial_p = compute_partial_corr(xc, yc, zc)
            rows.append({
                "metric": y_col,
                "PC": x_col,
                "pearson_r": pr,
                "pearson_p": pp,
                "spearman_rho": sr,
                "spearman_p": sp,
                "partial_r": partial_r,
                "partial_p": partial_p,
            })
    corr_df = pd.DataFrame(rows)

    # Apply FDR-BH multiple testing correction to primary 6 Pearson tests (S5 fix)
    try:
        from statsmodels.stats.multitest import multipletests
        pearson_p_vals = corr_df["pearson_p"].values
        _, pearson_p_adj, _, _ = multipletests(pearson_p_vals, method='fdr_bh')
        corr_df["pearson_p_fdr"] = pearson_p_adj

        # Also correct Spearman and partial p-values
        spearman_p_vals = corr_df["spearman_p"].values
        _, spearman_p_adj, _, _ = multipletests(spearman_p_vals, method='fdr_bh')
        corr_df["spearman_p_fdr"] = spearman_p_adj

        partial_p_vals = corr_df["partial_p"].values
        _, partial_p_adj, _, _ = multipletests(partial_p_vals, method='fdr_bh')
        corr_df["partial_p_fdr"] = partial_p_adj
    except ImportError:
        # Fallback: add Bonferroni threshold note
        n_tests = len(rows)
        corr_df["pearson_p_fdr"] = corr_df["pearson_p"]  # uncorrected if no statsmodels
        warnings.warn(
            f"statsmodels not available — reporting uncorrected p-values. "
            f"Bonferroni threshold for {n_tests} tests: α = {0.05/n_tests:.4f}"
        )

    return corr_df


def _get_group(sys_name):
    parts = sys_name.split("_")
    return f"{parts[0]}_{parts[1][:3]}"


def _get_marker(sys_name):
    return "^" if "tsap" in sys_name else "o"


def correlation_scatter(ax, x, y, color, label, marker="o"):
    ax.scatter(
        x, y, c=color, marker=marker, s=55,
        edgecolors="black", linewidths=0.5, zorder=3, label=label,
    )


def _add_regression(ax, x, y):
    xv = np.asarray(x, dtype=float)
    yv = np.asarray(y, dtype=float)
    valid = ~(np.isnan(xv) | np.isnan(yv))
    if valid.sum() > 2:
        xv, yv = xv[valid], yv[valid]
        r, p = stats.pearsonr(xv, yv)
        m, b = np.polyfit(xv, yv, 1)
        xl = np.linspace(xv.min(), xv.max(), 100)
        ax.plot(xl, m * xl + b, "k-", lw=1, alpha=0.7, zorder=2)
        ax.text(
            0.05, 0.95, f"r={r:.2f}, p={p:.3f}",
            transform=ax.transAxes, fontsize=8, va="top",
            bbox=dict(boxstyle="round,pad=0.3", fc="wheat", alpha=0.5),
        )


def _draw_scatter_panel(ax, df, pc_col, y_col, ylabel_text):
    group_arr = np.array([_get_group(s) for s in df["system"]])
    unique_groups = sorted(set(group_arr))
    for g in unique_groups:
        mask = group_arr == g
        sub = df.iloc[mask]
        mk = _get_marker(sub["system"].iloc[0])
        correlation_scatter(
            ax, sub[pc_col].values, sub[y_col].values,
            GROUP_COLORS[g], g, marker=mk,
        )
    _add_regression(ax, df[pc_col].values, df[y_col].values)
    ax.set_xlabel(pc_col.replace("mean_", ""), fontsize=9)
    ax.set_ylabel(ylabel_text, fontsize=9)


def plot_binding_correlation(df, outpath):
    fig, axes = plt.subplots(2, 3, figsize=(8.27, 11.69), dpi=300)

    pc_cols = ["mean_PC1", "mean_PC2", "mean_PC3"]
    y_cols = ["TOTAL_corrected", "mean_total_contacts"]
    y_labels = [
        r"TOTAL $\Delta$G corrected (kJ/mol)",
        "Mean total contacts",
    ]

    for row, (yc, yl) in enumerate(zip(y_cols, y_labels)):
        for col, pc in enumerate(pc_cols):
            ax = axes[row, col]
            _draw_scatter_panel(ax, df, pc, yc, yl if col == 0 else "")

    legend_groups = sorted(set(_get_group(s) for s in df["system"]))
    handles = [
        Line2D([0], [0], marker="o", color="w", markerfacecolor=GROUP_COLORS[g],
               markersize=7, label=g)
        for g in legend_groups
    ]
    handles.append(Line2D([0], [0], marker="o", color="w", markerfacecolor="gray",
                          markersize=7, label="sap"))
    handles.append(Line2D([0], [0], marker="^", color="w", markerfacecolor="gray",
                          markersize=7, label="tsap"))

    fig.legend(handles=handles, loc="lower center", ncol=6, fontsize=8,
               bbox_to_anchor=(0.5, 0.01))
    fig.tight_layout(rect=[0, 0.06, 1, 1])
    fig.savefig(outpath, dpi=300, bbox_inches="tight")
    plt.close(fig)


def compute_perframe_correlations(projections_path, all_ana_dir, systems, metadata):
    proj = pd.read_csv(projections_path)
    results = {}

    for sys_name in systems:
        sys_meta = next(m for m in metadata if m["name"] == sys_name)
        n_frames = sys_meta["n_frames"]
        sys_proj = proj[proj["system_label"] == sys_name].iloc[:n_frames]

        contacts_list = []
        for t in range(4):
            cpath = Path(all_ana_dir) / sys_name / "per_trial" / f"mmpbsa_{t}" / "contacts" / "contacts_total.csv"
            if cpath.exists():
                df_c = pd.read_csv(cpath)
                contacts_list.append(df_c["total_contacts"].values)

        if not contacts_list:
            results[sys_name] = {"PC1": np.nan, "PC2": np.nan, "PC3": np.nan}
            continue

        # Validate per-trial frame alignment (C2 escalation: per-trial count check)
        n_trials = len(contacts_list)
        frames_per_trial = n_frames // n_trials if n_trials > 0 else 0

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
        contacts_concat = np.concatenate(contacts_list)
        n_match = min(len(contacts_concat), len(sys_proj))

        r_dict = {}
        for pc in ["PC1", "PC2", "PC3"]:
            pc_vals = sys_proj[pc].values[:n_match]
            c_vals = contacts_concat[:n_match]
            valid = ~(np.isnan(pc_vals) | np.isnan(c_vals))
            if valid.sum() > 2:
                r, _ = stats.pearsonr(pc_vals[valid], c_vals[valid])
                r_dict[pc] = r
            else:
                r_dict[pc] = np.nan

        results[sys_name] = r_dict

    return results


def plot_perframe_heatmap(perframe_r, system_names, outpath):
    data = np.array(
        [[perframe_r[s][pc] for pc in ["PC1", "PC2", "PC3"]] for s in system_names]
    )

    fig, ax = plt.subplots(figsize=(4, 8.5), dpi=300)
    im = ax.imshow(data, cmap="RdBu_r", vmin=-1, vmax=1, aspect="auto")

    ax.set_xticks([0, 1, 2])
    ax.set_xticklabels(["PC1", "PC2", "PC3"], fontsize=9)
    ax.set_yticks(range(len(system_names)))
    ax.set_yticklabels(system_names, fontsize=7)

    for i in range(len(system_names)):
        for j in range(3):
            val = data[i, j]
            if not np.isnan(val):
                color = "white" if abs(val) > 0.5 else "black"
                ax.text(j, i, f"{val:.2f}", ha="center", va="center",
                        fontsize=6, color=color)

    plt.colorbar(im, ax=ax, label="Pearson r", shrink=0.5)
    ax.set_title("Per-frame PC\u2013Contacts correlation", fontsize=10)
    fig.tight_layout()
    fig.savefig(outpath, dpi=300, bbox_inches="tight")
    plt.close(fig)


def plot_correlation_bar(corr_df, outpath):
    fig, ax = plt.subplots(figsize=(8, 4), dpi=300)

    metrics = corr_df["metric"].unique()
    pcs = corr_df["PC"].unique()
    x_pos = np.arange(len(pcs))
    width = 0.35

    bar_colors = ["#4C72B0", "#DD8452"]
    for i, metric in enumerate(metrics):
        subset = corr_df[corr_df["metric"] == metric]
        r_vals = [subset[subset["PC"] == pc]["pearson_r"].values[0] for pc in pcs]
        bars = ax.bar(x_pos + i * width, r_vals, width, label=metric,
                      alpha=0.85, color=bar_colors[i % len(bar_colors)])
        for bar, rv in zip(bars, r_vals):
            ax.text(bar.get_x() + bar.get_width() / 2, bar.get_height(),
                    f"{rv:.2f}", ha="center", va="bottom", fontsize=7)

    ax.set_xticks(x_pos + width / 2)
    ax.set_xticklabels(pcs, fontsize=10)
    ax.set_ylabel("Pearson r", fontsize=10)
    ax.set_title("Correlation: PC scores vs binding metrics", fontsize=11)
    ax.axhline(y=0, color="black", linewidth=0.5)
    ax.legend(fontsize=9)

    fig.tight_layout()
    fig.savefig(outpath, dpi=300, bbox_inches="tight")
    plt.close(fig)


def plot_binding_summary(df, corr_df, perframe_r, system_names, outpath):
    fig = plt.figure(figsize=(8.27, 11.69), dpi=300)
    gs = gridspec.GridSpec(3, 3, figure=fig, hspace=0.45, wspace=0.4)

    pc_cols = ["mean_PC1", "mean_PC2", "mean_PC3"]
    y_cols = ["TOTAL_corrected", "mean_total_contacts"]
    y_labels = [
        r"$\Delta$G corrected (kJ/mol)",
        "Mean contacts",
    ]

    for row, (yc, yl) in enumerate(zip(y_cols, y_labels)):
        for col, pc in enumerate(pc_cols):
            ax = fig.add_subplot(gs[row, col])
            _draw_scatter_panel(ax, df, pc, yc, yl if col == 0 else "")

    ax_bar = fig.add_subplot(gs[2, :2])
    metrics = corr_df["metric"].unique()
    pcs = corr_df["PC"].unique()
    x_pos = np.arange(len(pcs))
    width = 0.35
    bar_colors = ["#4C72B0", "#DD8452"]
    for i, metric in enumerate(metrics):
        subset = corr_df[corr_df["metric"] == metric]
        r_vals = [subset[subset["PC"] == pc]["pearson_r"].values[0] for pc in pcs]
        ax_bar.bar(x_pos + i * width, r_vals, width, label=metric,
                   alpha=0.85, color=bar_colors[i % len(bar_colors)])
    ax_bar.set_xticks(x_pos + width / 2)
    ax_bar.set_xticklabels(pcs, fontsize=8)
    ax_bar.set_ylabel("Pearson r", fontsize=8)
    ax_bar.axhline(y=0, color="black", linewidth=0.5)
    ax_bar.legend(fontsize=7)

    ax_heat = fig.add_subplot(gs[2, 2])
    data = np.array(
        [[perframe_r[s][pc] for pc in ["PC1", "PC2", "PC3"]] for s in system_names]
    )
    im = ax_heat.imshow(data, cmap="RdBu_r", vmin=-1, vmax=1, aspect="auto")
    ax_heat.set_xticks([0, 1, 2])
    ax_heat.set_xticklabels(["PC1", "PC2", "PC3"], fontsize=7)
    ax_heat.set_yticks(range(0, len(system_names), 4))
    ax_heat.set_yticklabels(
        [system_names[i] for i in range(0, len(system_names), 4)], fontsize=6
    )
    plt.colorbar(im, ax=ax_heat, shrink=0.6)

    fig.tight_layout()
    fig.savefig(outpath, dpi=300, bbox_inches="tight")
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input-dir", default="pca_analysis/")
    parser.add_argument(
        "--all-ana-dir",
        default=os.environ.get('PCA_ALL_ANA_DIR', '../com_md/all_ana/per_ligand/'),
    )
    parser.add_argument(
        "--com-md-dir",
        default=os.environ.get('PCA_MD_DIR', '../com_md/'),
    )
    args = parser.parse_args()

    input_dir = Path(args.input_dir)
    all_ana_dir = args.all_ana_dir

    systems = SYSTEMS

    metadata_path = input_dir / "system_metadata.json"
    with open(metadata_path) as f:
        metadata = json.load(f)

    projections_path = input_dir / "projections_all.csv"
    system_indices_path = input_dir / "system_indices.npy"

    df_mmpbsa = load_all_mmpbsa(all_ana_dir, systems)
    df_contacts = load_all_contacts(all_ana_dir, systems)
    df_pc = compute_system_mean_pc(
        str(projections_path), str(system_indices_path), str(metadata_path)
    )

    df_merged = df_mmpbsa.merge(df_contacts, on="system").merge(df_pc, on="system")

    chirality = np.array(
        [1.0 if s.startswith("phe") else 0.0 for s in df_merged["system"]]
    )

    corr_df = compute_correlations(
        df_merged,
        x_cols=["mean_PC1", "mean_PC2", "mean_PC3"],
        y_cols=["TOTAL_corrected", "mean_total_contacts"],
        chirality=chirality,
    )

    out_mmpbsa = df_merged[
        ["system", "TOTAL_kjmol", "TOTAL_corrected", "SD", "SEM", "n_trials",
         "mean_PC1", "mean_PC2", "mean_PC3"]
    ]
    out_mmpbsa.to_csv(input_dir / "mmpbsa_energies_all_systems.csv", index=False)

    out_contacts = df_merged[
        ["system", "mean_total_contacts", "std_total_contacts", "bindingsite_frac",
         "mean_PC1", "mean_PC2", "mean_PC3"]
    ]
    out_contacts.to_csv(input_dir / "contacts_summary_all_systems.csv", index=False)

    corr_df.to_csv(input_dir / "binding_correlation_table.csv", index=False)

    perframe_r = compute_perframe_correlations(
        str(projections_path), all_ana_dir, systems, metadata
    )

    plot_binding_correlation(df_merged, str(input_dir / "figure_binding_correlation.png"))
    plot_perframe_heatmap(
        perframe_r, systems, str(input_dir / "figure_binding_perframe_heatmap.png")
    )
    plot_correlation_bar(corr_df, str(input_dir / "figure_correlation_bar.png"))
    plot_binding_summary(
        df_merged, corr_df, perframe_r, systems,
        str(input_dir / "figure_binding_summary.png"),
    )

    print("pca_binding.py complete.")
    print(f"Output directory: {input_dir}")
    for f in sorted(input_dir.glob("mmpbsa_energies*")):
        print(f"  {f.name}")
    for f in sorted(input_dir.glob("contacts_summary*")):
        print(f"  {f.name}")
    for f in sorted(input_dir.glob("binding_correlation*")):
        print(f"  {f.name}")
    for f in sorted(input_dir.glob("figure_binding*")):
        print(f"  {f.name}")
    for f in sorted(input_dir.glob("figure_correlation*")):
        print(f"  {f.name}")


if __name__ == "__main__":
    main()
