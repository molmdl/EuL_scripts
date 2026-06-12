#!/usr/bin/env python3
"""
C3_group_comparison_3way.py — 3-way group comparison metrics and plots
for xtb, solv_md, and com_md in eu8_nochrom PC space.

Computes:
  - Per-system com_md metrics (centroid, std, N_eff)
  - Group comparisons: SAP vs TSAP, D vs L, me vs phe — for each method
  - Cross-method comparisons: xtb_vs_com, solv_vs_com
  - Bootstrap 95% CIs on centroid distances

Plots:
  - plot_SAP_TSAP_3way.png — grouped bar chart (3 methods)
  - plot_D_L_3way.png — grouped bar chart (3 methods, log scale)
  - plot_me_phe_3way.png — grouped bar chart (3 methods)
  - plot_std_comparison_3way.png — std(PC1)/std(PC2) for all 3 methods

Usage:
  python cross_dataset/com_md/scripts/C3_group_comparison_3way.py \
      --input cross_dataset/com_md/analysis/joint_projection_3way_eu8_nochrom.csv \
      --outdir cross_dataset/com_md/analysis
"""

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats


# ---------- Constants ----------
METHODS = ["xtb", "solv", "com_md"]
METHOD_LABELS = {"xtb": "xTB", "solv": "solv_md", "com_md": "com_md"}
METHOD_COLORS = {"xtb": "#1f77b4", "solv": "#2ca02c", "com_md": "#d62728"}


# ---------- Helpers ----------

def compute_neff(series):
    """Compute effective sample size using lag-1 autocorrelation."""
    arr = np.asarray(series, dtype=float)
    n = len(arr)
    if n < 2:
        return n
    rho = np.corrcoef(arr[:-1], arr[1:])[0, 1]
    if np.isnan(rho) or rho <= 0:
        return n
    return max(1, int(n * (1 - rho) / (1 + rho)))


def block_bootstrap_ci(data1, data2, n_bootstrap=1000, block_size=None,
                       rng_seed=42):
    """Block bootstrap for centroid distance confidence interval.

    Returns 95% CI (ci_low, ci_high).
    """
    data1 = np.asarray(data1, dtype=float)
    data2 = np.asarray(data2, dtype=float)
    n1, n2 = len(data1), len(data2)

    if block_size is None:
        block_size = max(1, int(min(n1, n2) / 10))

    rng = np.random.default_rng(rng_seed)
    distances = []

    for _ in range(n_bootstrap):
        idx1 = _block_resample(n1, block_size, rng)
        idx2 = _block_resample(n2, block_size, rng)

        s1 = data1[idx1]
        s2 = data2[idx2]

        c1 = s1.mean(axis=0)
        c2 = s2.mean(axis=0)
        distances.append(np.linalg.norm(c1 - c2))

    distances = np.array(distances)
    ci_low = float(np.percentile(distances, 2.5))
    ci_high = float(np.percentile(distances, 97.5))
    return ci_low, ci_high


def _block_resample(n, block_size, rng):
    """Generate block-resampled indices."""
    indices = []
    while len(indices) < n:
        start = rng.integers(0, n)
        block = np.arange(start, min(start + block_size, n))
        if len(block) < block_size and start + block_size > n:
            wrapped = np.arange(0, block_size - len(block))
            block = np.concatenate([block, wrapped])
        indices.extend(block.tolist())
    return np.array(indices[:n])


# ---------- Per-system com_md metrics ----------

def compute_com_md_per_system_metrics(df):
    """Compute per-system metrics for com_md only.

    Returns DataFrame with one row per system.
    """
    df_com = df[df["method"] == "com_md"].copy()
    systems = sorted(df_com["system"].unique())
    rows = []

    for sys_name in systems:
        sub = df_com[df_com["system"] == sys_name]
        if sub.empty:
            continue

        meta = sub.iloc[0]
        neff1 = compute_neff(sub["PC1"].values)
        neff2 = compute_neff(sub["PC2"].values)

        row = {
            "system": sys_name,
            "species": meta["species"],
            "isomer": meta["isomer"],
            "handedness": meta["handedness"],
            "conformer": meta["conformer"],
            "n_frames": len(sub),
            "Neff_pt1": neff1,
            "Neff_pt2": neff2,
            "centroid_PC1": sub["PC1"].mean(),
            "centroid_PC2": sub["PC2"].mean(),
            "std_PC1": sub["PC1"].std(ddof=1),
            "std_PC2": sub["PC2"].std(ddof=1),
            "span_PC1": sub["PC1"].max() - sub["PC1"].min(),
            "span_PC2": sub["PC2"].max() - sub["PC2"].min(),
        }
        rows.append(row)

    return pd.DataFrame(rows)


# ---------- Grouped metrics ----------

def compute_group_metrics_3way(df):
    """Compute SAP↔TSAP, D↔L, me↔phe comparisons for all 3 methods.

    Also compute cross-method comparisons: xtb_vs_com, solv_vs_com.

    Returns DataFrame with distance, CI, and ratio columns.
    """
    rows = []

    # ---- SAP vs TSAP ----
    groups = df.groupby(["species", "isomer", "handedness"])
    for (sp, iso, hand), grp in groups:
        for method in METHODS:
            sub = grp[grp["method"] == method]
            sap = sub[sub["conformer"] == "sap"]
            tsap = sub[sub["conformer"] == "tsap"]
            if len(sap) == 0 or len(tsap) == 0:
                continue
            sap_xy = sap[["PC1", "PC2"]].values
            tsap_xy = tsap[["PC1", "PC2"]].values
            c_sap = sap_xy.mean(axis=0)
            c_tsap = tsap_xy.mean(axis=0)
            dist = np.linalg.norm(c_sap - c_tsap)
            ci_low, ci_high = block_bootstrap_ci(sap_xy, tsap_xy)
            rows.append({
                "metric_type": "sap_vs_tsap",
                "group_key": f"{sp}_{iso}{hand}",
                "species": sp,
                "isomer": iso,
                "handedness": hand,
                "conformer": "",
                "method": method,
                "distance": dist,
                "ci_low": ci_low,
                "ci_high": ci_high,
            })

    # ---- me vs phe ----
    for isomer in df["isomer"].unique():
        for hand in df["handedness"].unique():
            for conf in df["conformer"].unique():
                for method in METHODS:
                    sub = df[(df["species"].isin(["me", "phe"])) &
                             (df["isomer"] == isomer) &
                             (df["handedness"] == hand) &
                             (df["conformer"] == conf) &
                             (df["method"] == method)]
                    me_sub = sub[sub["species"] == "me"]
                    phe_sub = sub[sub["species"] == "phe"]
                    if len(me_sub) == 0 or len(phe_sub) == 0:
                        continue
                    me_xy = me_sub[["PC1", "PC2"]].values
                    phe_xy = phe_sub[["PC1", "PC2"]].values
                    c_me = me_xy.mean(axis=0)
                    c_phe = phe_xy.mean(axis=0)
                    dist = np.linalg.norm(c_me - c_phe)
                    ci_low, ci_high = block_bootstrap_ci(me_xy, phe_xy)
                    rows.append({
                        "metric_type": "me_vs_phe",
                        "group_key": f"{isomer}{hand}_{conf}",
                        "species": "me_vs_phe",
                        "isomer": isomer,
                        "handedness": hand,
                        "conformer": conf,
                        "method": method,
                        "distance": dist,
                        "ci_low": ci_low,
                        "ci_high": ci_high,
                    })

    # ---- D vs L ----
    for sp in df["species"].unique():
        for iso in df["isomer"].unique():
            for conf in df["conformer"].unique():
                for method in METHODS:
                    sub = df[(df["species"] == sp) &
                             (df["isomer"] == iso) &
                             (df["conformer"] == conf) &
                             (df["method"] == method)]
                    d_sub = sub[sub["handedness"] == "D"]
                    l_sub = sub[sub["handedness"] == "L"]
                    if len(d_sub) == 0 or len(l_sub) == 0:
                        continue
                    d_xy = d_sub[["PC1", "PC2"]].values
                    l_xy = l_sub[["PC1", "PC2"]].values
                    c_d = d_xy.mean(axis=0)
                    c_l = l_xy.mean(axis=0)
                    dist = np.linalg.norm(c_d - c_l)
                    ci_low, ci_high = block_bootstrap_ci(d_xy, l_xy)
                    rows.append({
                        "metric_type": "D_vs_L",
                        "group_key": f"{sp}_{iso}_{conf}",
                        "species": sp,
                        "isomer": iso,
                        "handedness": "",
                        "conformer": conf,
                        "method": method,
                        "distance": dist,
                        "ci_low": ci_low,
                        "ci_high": ci_high,
                    })

    # ---- Cross-method: xtb_vs_com, solv_vs_com (per system) ----
    systems = sorted(df["system"].unique())
    for sys_name in systems:
        for method_pair in [("xtb", "com_md"), ("solv", "com_md")]:
            m1, m2 = method_pair
            sub1 = df[(df["system"] == sys_name) & (df["method"] == m1)]
            sub2 = df[(df["system"] == sys_name) & (df["method"] == m2)]
            if len(sub1) == 0 or len(sub2) == 0:
                continue
            xy1 = sub1[["PC1", "PC2"]].values
            xy2 = sub2[["PC1", "PC2"]].values
            c1 = xy1.mean(axis=0)
            c2 = xy2.mean(axis=0)
            dist = np.linalg.norm(c1 - c2)
            ci_low, ci_high = block_bootstrap_ci(xy1, xy2)
            meta = sub1.iloc[0]
            rows.append({
                "metric_type": f"{m1}_vs_{m2}",
                "group_key": sys_name,
                "species": meta["species"],
                "isomer": meta["isomer"],
                "handedness": meta["handedness"],
                "conformer": meta["conformer"],
                "method": f"{m1}_vs_{m2}",
                "distance": dist,
                "ci_low": ci_low,
                "ci_high": ci_high,
            })

    all_rows = pd.DataFrame(rows)

    # Pivot within-method comparisons (sap_vs_tsap, me_vs_phe, D_vs_L)
    # to get all 3 method distances side-by-side
    within_methods = all_rows[all_rows["metric_type"].isin(
        ["sap_vs_tsap", "me_vs_phe", "D_vs_L"]
    )].copy()

    idx_cols = [c for c in within_methods.columns
                if c not in ("method", "distance", "ci_low", "ci_high")]

    def pivot_metric(df_long, value_col, prefix):
        pivoted = df_long.pivot_table(
            index=[c for c in idx_cols],
            columns="method", values=value_col
        ).reset_index()
        pivoted.columns.name = None
        rename_map = {}
        for m in METHODS:
            if m in pivoted.columns:
                rename_map[m] = f"{prefix}_{m}"
        pivoted = pivoted.rename(columns=rename_map)
        return pivoted

    piv_dist = pivot_metric(within_methods, "distance", "d")
    piv_ci_low = pivot_metric(within_methods, "ci_low", "d_ci_low")
    piv_ci_high = pivot_metric(within_methods, "ci_high", "d_ci_high")

    # Merge
    all_groups = piv_dist
    for piv in [piv_ci_low, piv_ci_high]:
        common_cols = [c for c in piv.columns if c in all_groups.columns]
        new_cols = [c for c in piv.columns if c not in all_groups.columns]
        if new_cols:
            all_groups = all_groups.merge(
                piv[common_cols + new_cols], on=common_cols, how="left"
            )

    # Add ratios (com_md / xtb, com_md / solv)
    for base_method in ["xtb", "solv"]:
        d_col = f"d_{base_method}"
        d_com = "d_com_md"
        if d_col in all_groups.columns and d_com in all_groups.columns:
            ratio_col = f"ratio_com_over_{base_method}"
            all_groups[ratio_col] = all_groups[d_com] / all_groups[d_col].replace(0, np.nan)

    # Also add solv/xtb ratio for backward compatibility
    if "d_xtb" in all_groups.columns and "d_solv" in all_groups.columns:
        all_groups["ratio_solv_over_xtb"] = all_groups["d_solv"] / all_groups["d_xtb"].replace(0, np.nan)

    # Append cross-method rows
    cross_method = all_rows[all_rows["metric_type"].str.contains("_vs_")].copy()
    cross_method = cross_method[~cross_method["metric_type"].isin(
        ["sap_vs_tsap", "me_vs_phe", "D_vs_L"]
    )]

    return all_groups, cross_method


# ---------- Plots ----------

def plot_3way_grouped_bar(group_metrics, metric_type, output_path, log_scale=False):
    """Plot grouped bar chart with 3 methods side-by-side + CI error bars.

    Parameters
    ----------
    group_metrics : DataFrame (pivoted, with d_xtb, d_solv, d_com_md columns)
    metric_type : str — "sap_vs_tsap", "me_vs_phe", or "D_vs_L"
    output_path : Path
    log_scale : bool — use log scale for y-axis
    """
    sub = group_metrics[group_metrics["metric_type"] == metric_type].copy()
    if sub.empty:
        print(f"  No data for {metric_type}, skipping.")
        return

    sub = sub.sort_values("group_key")
    labels = sub["group_key"].values
    n_groups = len(labels)
    x = np.arange(n_groups)
    bar_width = 0.25

    fig, ax1 = plt.subplots(figsize=(max(12, n_groups * 1.5), 6))

    if log_scale:
        ax1.set_yscale("log")
        ylabel = "Centroid distance (PC units, log scale)"
    else:
        ylabel = "Centroid distance (PC units)"

    for i, method in enumerate(METHODS):
        d_col = f"d_{method}"
        ci_lo_col = f"d_ci_low_{method}"
        ci_hi_col = f"d_ci_high_{method}"

        if d_col not in sub.columns:
            continue

        d_vals = sub[d_col].values

        # Error bars
        if ci_lo_col in sub.columns and ci_hi_col in sub.columns:
            ci_lo = sub[ci_lo_col].values
            ci_hi = sub[ci_hi_col].values
            err = np.vstack([
                d_vals - ci_lo,
                ci_hi - d_vals
            ])
        else:
            err = None

        ax1.bar(
            x + (i - 1) * bar_width, d_vals, bar_width,
            label=METHOD_LABELS[method],
            color=METHOD_COLORS[method], alpha=0.85,
            yerr=err, capsize=3, error_kw={"linewidth": 1.0},
        )

    title_map = {
        "sap_vs_tsap": "SAP–TSAP Centroid Distance: 3-Way Comparison (95% CI)",
        "me_vs_phe": "me vs phe Centroid Distance: 3-Way Comparison (95% CI)",
        "D_vs_L": "D vs L Centroid Distance: 3-Way Comparison (95% CI)",
    }

    ax1.set_xticks(x)
    ax1.set_xticklabels(labels, rotation=45, ha="right", fontsize=9)
    ax1.set_ylabel(ylabel, fontsize=11)
    ax1.set_title(title_map.get(metric_type, metric_type), fontsize=13, fontweight="bold")
    ax1.legend(fontsize=10, loc="best")
    ax1.grid(axis="y", alpha=0.3, linestyle="--")

    fig.tight_layout()
    fig.savefig(output_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved: {output_path}")


def plot_std_comparison_3way(df, output_path):
    """Plot std(PC1) and std(PC2) for all 3 methods across all 16 systems.

    Grouped bar chart with 3 methods per system.
    """
    systems = sorted(df["system"].unique())
    n_systems = len(systems)
    x = np.arange(n_systems)
    bar_width = 0.25

    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(16, 10), sharex=True)

    for ax, pc_col, pc_label in [(ax1, "PC1", "PC1"), (ax2, "PC2", "PC2")]:
        for i, method in enumerate(METHODS):
            df_method = df[df["method"] == method]
            stds = []
            for sys_name in systems:
                sub = df_method[df_method["system"] == sys_name]
                if len(sub) > 1:
                    stds.append(sub[pc_col].std(ddof=1))
                else:
                    stds.append(0)
            ax.bar(
                x + (i - 1) * bar_width, stds, bar_width,
                label=METHOD_LABELS[method],
                color=METHOD_COLORS[method], alpha=0.85,
            )

        ax.set_ylabel(f"std({pc_label}) [PC units]", fontsize=11)
        ax.legend(fontsize=9, loc="best")
        ax.grid(axis="y", alpha=0.3, linestyle="--")

    ax2.set_xticks(x)
    ax2.set_xticklabels(systems, rotation=45, ha="right", fontsize=9)
    ax2.set_xlabel("System", fontsize=11)

    fig.suptitle("Conformational Breadth: std(PC1) and std(PC2) across 3 Methods",
                 fontsize=13, fontweight="bold")
    fig.tight_layout()
    fig.savefig(output_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved: {output_path}")


# ---------- Main ----------

def main():
    parser = argparse.ArgumentParser(
        description="C3: 3-way group comparison metrics and plots"
    )
    parser.add_argument(
        "--input", type=Path, required=True,
        help="3-way joint projection CSV",
    )
    parser.add_argument(
        "--outdir", type=Path, required=True,
        help="Output directory for CSVs and PNGs",
    )
    args = parser.parse_args()

    args.outdir.mkdir(parents=True, exist_ok=True)

    # Load data
    print(f"Loading: {args.input}")
    df = pd.read_csv(args.input)
    print(f"  {len(df)} rows, methods: {df['method'].unique().tolist()}")

    # ── Per-system com_md metrics ──
    print("\n[1/5] Computing per-system com_md metrics...")
    com_metrics = compute_com_md_per_system_metrics(df)
    csv_path = args.outdir / "com_md_per_system_metrics.csv"
    com_metrics.to_csv(csv_path, index=False)
    print(f"  Saved: {csv_path} ({len(com_metrics)} rows)")

    # Report summary
    print(f"  com_md avg std(PC1): {com_metrics['std_PC1'].mean():.2f}")
    print(f"  com_md avg std(PC2): {com_metrics['std_PC2'].mean():.2f}")
    print(f"  com_md avg N_eff(pt1): {com_metrics['Neff_pt1'].mean():.0f}")

    # ── Grouped metrics ──
    print("\n[2/5] Computing 3-way group metrics with bootstrap CIs...")
    group_metrics, cross_method_metrics = compute_group_metrics_3way(df)

    csv_path2 = args.outdir / "three_way_group_metrics.csv"
    # Combine within-method and cross-method into one CSV
    combined = pd.concat([group_metrics, cross_method_metrics], ignore_index=True)
    combined.to_csv(csv_path2, index=False)
    print(f"  Saved: {csv_path2} ({len(combined)} rows)")

    # Report key ratios
    for metric_type in ["sap_vs_tsap", "me_vs_phe", "D_vs_L"]:
        sub = group_metrics[group_metrics["metric_type"] == metric_type]
        if len(sub) > 0:
            for ratio_col in ["ratio_com_over_xtb", "ratio_com_over_solv", "ratio_solv_over_xtb"]:
                if ratio_col in sub.columns:
                    avg = sub[ratio_col].mean()
                    print(f"  {metric_type} avg {ratio_col}: {avg:.2f}")

    # ── Plots ──
    print("\n[3/5] Generating SAP/TSAP comparison plot...")
    plot_3way_grouped_bar(
        group_metrics, "sap_vs_tsap",
        args.outdir / "plot_SAP_TSAP_3way.png",
        log_scale=False,
    )

    print("\n[4/5] Generating D/L and me/phe comparison plots...")
    plot_3way_grouped_bar(
        group_metrics, "D_vs_L",
        args.outdir / "plot_D_L_3way.png",
        log_scale=True,
    )
    plot_3way_grouped_bar(
        group_metrics, "me_vs_phe",
        args.outdir / "plot_me_phe_3way.png",
        log_scale=False,
    )

    print("\n[5/5] Generating std comparison plot...")
    plot_std_comparison_3way(df, args.outdir / "plot_std_comparison_3way.png")

    # ── Summary ──
    print("\n" + "=" * 60)
    print("C3 SUMMARY")
    print("=" * 60)
    print(f"  com_md per-system metrics: {csv_path}")
    print(f"  3-way group metrics: {csv_path2}")
    print(f"  Plots:")
    print(f"    {args.outdir / 'plot_SAP_TSAP_3way.png'}")
    print(f"    {args.outdir / 'plot_D_L_3way.png'}")
    print(f"    {args.outdir / 'plot_me_phe_3way.png'}")
    print(f"    {args.outdir / 'plot_std_comparison_3way.png'}")
    print("=" * 60)


if __name__ == "__main__":
    main()
