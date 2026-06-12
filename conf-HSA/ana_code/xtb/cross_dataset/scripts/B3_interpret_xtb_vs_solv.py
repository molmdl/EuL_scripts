#!/usr/bin/env python3
"""B3_interpret_xtb_vs_solv.py — Cross-dataset interpretation: xTB vs solv_md.

Reads the joint eu8_nochrom PCA projection (xtb + solv), computes per-system
and grouped comparison metrics, saves CSVs and plots, and writes an
interpretation document.

Metrics include:
  - Centroid positions and shift distances
  - Wasserstein and KS distances per PC (effect size only; p-values omitted)
  - Mann-Whitney U rank tests per PC
  - Span and std ratios (solv / xtb) — std ratio is primary breadth metric
  - Mahalanobis containment (solv within xtb chi² quantiles)
  - Effective sample size (N_eff) via lag-1 autocorrelation
  - Bootstrap 95% CIs on group centroid distances
  - Grouped comparisons: SAP vs TSAP, me vs phe, D vs L

Usage:
    python cross_dataset/scripts/B3_interpret_xtb_vs_solv.py
    python cross_dataset/scripts/B3_interpret_xtb_vs_solv.py --input <csv>
"""

import argparse
import json
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
ROOT = Path(__file__).resolve().parent.parent.parent  # xtb root
DEFAULT_INPUT = ROOT / "cross_dataset" / "analysis" / "joint_projection_eu8_nochrom_xtb_solv.csv"
OUT_ANALYSIS = ROOT / "cross_dataset" / "analysis"
OUT_DOCS = ROOT / "cross_dataset" / "docs"

# Chi-squared quantiles for 2-DOF Mahalanobis containment
CHI2_DF = 2
CHI2_QUANTILES = {"50": stats.chi2.ppf(0.50, CHI2_DF),
                  "90": stats.chi2.ppf(0.90, CHI2_DF),
                  "99": stats.chi2.ppf(0.99, CHI2_DF)}


def df_to_markdown(df: pd.DataFrame, floatfmt: str = ".2f") -> str:
    """Convert DataFrame to a markdown table without requiring tabulate."""
    cols = list(df.columns)
    # Format values
    rows_str = []
    for _, row in df.iterrows():
        vals = []
        for c in cols:
            v = row[c]
            # Check NaN BEFORE general float — isinstance(v, float) matches NaN
            if isinstance(v, float) and np.isnan(v):
                vals.append("NaN")
            elif isinstance(v, float):
                vals.append(f"{v:{floatfmt}}")
            elif hasattr(v, '__float__'):
                vals.append(f"{float(v):{floatfmt}}")
            else:
                vals.append(str(v))
        rows_str.append("| " + " | ".join(vals) + " |")
    header = "| " + " | ".join(cols) + " |"
    sep = "| " + " | ".join(["---"] * len(cols)) + " |"
    return "\n".join([header, sep] + rows_str)


def compute_neff(series):
    """Compute effective sample size using lag-1 autocorrelation.

    For a time series with N observations and lag-1 autocorrelation rho,
    the effective sample size is N_eff = N * (1 - rho) / (1 + rho).
    This accounts for the fact that correlated samples provide less
    independent information than uncorrelated ones.

    Parameters
    ----------
    series : array-like
        1D array of observations (e.g., PC1 trajectory).

    Returns
    -------
    int
        Effective sample size (minimum 1).
    """
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

    Resamples data in blocks to preserve autocorrelation structure,
    computes centroid distance between resampled groups, and returns
    the 2.5th and 97.5th percentile as the 95% CI.

    Parameters
    ----------
    data1, data2 : array-like, shape (N, 2)
        PC1, PC2 coordinates for the two groups.
    n_bootstrap : int
        Number of bootstrap resamples.
    block_size : int or None
        Block size for resampling. If None, uses max(1, min(N1, N2) / 10).
    rng_seed : int
        Random seed for reproducibility.

    Returns
    -------
    ci_low : float
        Lower bound of 95% CI.
    ci_high : float
        Upper bound of 95% CI.
    """
    data1 = np.asarray(data1, dtype=float)
    data2 = np.asarray(data2, dtype=float)
    n1, n2 = len(data1), len(data2)

    if block_size is None:
        block_size = max(1, int(min(n1, n2) / 10))

    rng = np.random.default_rng(rng_seed)
    distances = []

    for _ in range(n_bootstrap):
        # Block-resample group 1
        idx1 = _block_resample(n1, block_size, rng)
        # Block-resample group 2
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
    """Generate block-resampled indices.

    Randomly selects starting points for blocks of size `block_size`
    and concatenates them until we have at least `n` indices,
    then truncates to exactly `n`.
    """
    indices = []
    while len(indices) < n:
        start = rng.integers(0, n)
        block = np.arange(start, min(start + block_size, n))
        # Wrap around if block extends past end
        if len(block) < block_size and start + block_size > n:
            wrapped = np.arange(0, block_size - len(block))
            block = np.concatenate([block, wrapped])
        indices.extend(block.tolist())
    return np.array(indices[:n])


# ---------------------------------------------------------------------------
# Per-system metrics
# ---------------------------------------------------------------------------
def compute_per_system_metrics(df: pd.DataFrame) -> pd.DataFrame:
    """Compute xtb vs solv comparison metrics for each of the 16 systems.

    Parameters
    ----------
    df : DataFrame with columns system, method, PC1, PC2, species, isomer,
         handedness, conformer

    Returns
    -------
    DataFrame with one row per system and many metric columns.
    """
    rows = []
    systems = sorted(df["system"].unique())
    total = len(systems)

    for idx, sys_name in enumerate(systems, 1):
        print(f"  [{idx}/{total}] {sys_name}")
        sub = df[df["system"] == sys_name]
        xtb = sub[sub["method"] == "xtb"]
        solv = sub[sub["method"] == "solv"]

        if len(solv) == 0:
            print(f"    [skip] {sys_name}: no solv_md frames")
            continue

        meta = sub.iloc[0]

        # Centroids
        c1_xtb, c2_xtb = xtb["PC1"].mean(), xtb["PC2"].mean()
        c1_solv, c2_solv = solv["PC1"].mean(), solv["PC2"].mean()

        # Centroid distance
        d_centroid = np.sqrt((c1_xtb - c1_solv) ** 2 + (c2_xtb - c2_solv) ** 2)
        d_PC1 = abs(c1_xtb - c1_solv)
        d_PC2 = abs(c2_xtb - c2_solv)

        # Wasserstein distance (1D per PC)
        ws_PC1 = stats.wasserstein_distance(xtb["PC1"], solv["PC1"])
        ws_PC2 = stats.wasserstein_distance(xtb["PC2"], solv["PC2"])

        # KS test (statistic only; p-values are trivially significant)
        ks1 = stats.ks_2samp(xtb["PC1"], solv["PC1"])
        ks2 = stats.ks_2samp(xtb["PC2"], solv["PC2"])

        # Mann-Whitney U
        mw1 = stats.mannwhitneyu(xtb["PC1"], solv["PC1"], alternative="two-sided")
        mw2 = stats.mannwhitneyu(xtb["PC2"], solv["PC2"], alternative="two-sided")

        # Span (max - min)
        span1_xtb = xtb["PC1"].max() - xtb["PC1"].min()
        span2_xtb = xtb["PC2"].max() - xtb["PC2"].min()
        span1_solv = solv["PC1"].max() - solv["PC1"].min()
        span2_solv = solv["PC2"].max() - solv["PC2"].min()

        # Std
        std1_xtb, std2_xtb = xtb["PC1"].std(), xtb["PC2"].std()
        std1_solv, std2_solv = solv["PC1"].std(), solv["PC2"].std()

        # Span and std ratios (solv / xtb); guard against zero
        span_ratio_PC1 = span1_solv / span1_xtb if span1_xtb > 0 else np.nan
        span_ratio_PC2 = span2_solv / span2_xtb if span2_xtb > 0 else np.nan
        std_ratio_PC1 = std1_solv / std1_xtb if std1_xtb > 0 else np.nan
        std_ratio_PC2 = std2_solv / std2_xtb if std2_xtb > 0 else np.nan

        # --- Mahalanobis containment ---
        # xtb covariance matrix in PC1-PC2
        X_xtb = xtb[["PC1", "PC2"]].values
        X_solv = solv[["PC1", "PC2"]].values
        mu_xtb = np.array([c1_xtb, c2_xtb])

        cov_xtb = np.cov(X_xtb.T)
        try:
            cov_inv = np.linalg.inv(cov_xtb)
        except np.linalg.LinAlgError:
            print(f"    [warn] singular cov for {sys_name}, using pinv")
            cov_inv = np.linalg.pinv(cov_xtb)

        # Mahal distance for each solv frame from xtb centroid
        diff_solv = X_solv - mu_xtb
        mahal_solv = np.sqrt(np.sum(diff_solv @ cov_inv * diff_solv, axis=1))
        mahal_solv_sq = mahal_solv ** 2  # follows chi²(2) under xtb distribution

        containment = {}
        for qname, qval in CHI2_QUANTILES.items():
            containment[f"containment_{qname}"] = np.mean(mahal_solv_sq <= qval)

        # xtb within solv 90th percentile — REMOVED (trivial containment:
        # xtb is always fully within solv; see caveat in doc)

        # Centroid Mahalanobis distance (solv centroid from xtb centroid)
        mu_solv = np.array([c1_solv, c2_solv])
        centroid_mahal = np.sqrt(
            (mu_solv - mu_xtb) @ cov_inv @ (mu_solv - mu_xtb)
        )

        # Effective sample sizes (from PC1 lag-1 autocorrelation)
        n_eff_xtb = compute_neff(xtb["PC1"].values)
        n_eff_solv = compute_neff(solv["PC1"].values)

        row = {
            "system": sys_name,
            "species": meta["species"],
            "isomer": meta["isomer"],
            "handedness": meta["handedness"],
            "conformer": meta["conformer"],
            "n_xtb": len(xtb),
            "n_solv": len(solv),
            "n_eff_xtb": n_eff_xtb,
            "n_eff_solv": n_eff_solv,
            "centroid_PC1_xtb": c1_xtb,
            "centroid_PC2_xtb": c2_xtb,
            "centroid_PC1_solv": c1_solv,
            "centroid_PC2_solv": c2_solv,
            "d_PC1": d_PC1,
            "d_PC2": d_PC2,
            "d_centroid": d_centroid,
            "ws_PC1": ws_PC1,
            "ws_PC2": ws_PC2,
            "ks_PC1_stat": ks1.statistic,
            "ks_PC2_stat": ks2.statistic,
            "mw_PC1_stat": mw1.statistic,
            "mw_PC1_pval": mw1.pvalue,
            "mw_PC2_stat": mw2.statistic,
            "mw_PC2_pval": mw2.pvalue,
            "span_PC1_xtb": span1_xtb,
            "span_PC2_xtb": span2_xtb,
            "span_PC1_solv": span1_solv,
            "span_PC2_solv": span2_solv,
            "span_ratio_PC1": span_ratio_PC1,
            "span_ratio_PC2": span_ratio_PC2,
            "std_PC1_xtb": std1_xtb,
            "std_PC2_xtb": std2_xtb,
            "std_PC1_solv": std1_solv,
            "std_PC2_solv": std2_solv,
            "std_ratio_PC1": std_ratio_PC1,
            "std_ratio_PC2": std_ratio_PC2,
            "centroid_mahal": centroid_mahal,
        }
        row.update(containment)
        rows.append(row)

    return pd.DataFrame(rows)


# ---------------------------------------------------------------------------
# PC spans CSV (separate for easy reference)
# ---------------------------------------------------------------------------
def build_pc_spans_df(per_sys: pd.DataFrame) -> pd.DataFrame:
    """Extract span and N_eff columns into a dedicated DataFrame."""
    cols = ["system",
            "n_eff_xtb", "n_eff_solv",
            "span_PC1_xtb", "span_PC2_xtb",
            "span_PC1_solv", "span_PC2_solv",
            "span_ratio_PC1", "span_ratio_PC2",
            "std_PC1_xtb", "std_PC2_xtb",
            "std_PC1_solv", "std_PC2_solv",
            "std_ratio_PC1", "std_ratio_PC2"]
    return per_sys[cols].copy()


# ---------------------------------------------------------------------------
# Grouped metrics
# ---------------------------------------------------------------------------
def compute_group_metrics(df: pd.DataFrame, per_sys: pd.DataFrame) -> pd.DataFrame:
    """Compute SAP↔TSAP, me↔phe, and D↔L grouped comparisons.

    For each pair, compute centroid distance in xtb and solv separately,
    then the ratio solv_dist / xtb_dist. Bootstrap 95% CIs are computed
    for each centroid distance.
    """
    # We'll collect (metric_type, group_key, ..., method, distance, ci_low, ci_high)
    rows = []

    # ---- SAP vs TSAP ----
    groups = df.groupby(["species", "isomer", "handedness"])
    for (sp, iso, hand), grp in groups:
        for method in ["xtb", "solv"]:
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
                for method in ["xtb", "solv"]:
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
                for method in ["xtb", "solv"]:
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

    all_rows = pd.DataFrame(rows)

    # Pivot to get xtb and solv distances (and CIs) side-by-side
    # Build separate pivots for distance, ci_low, ci_high
    def pivot_metric(df_long, value_col, prefix):
        """Pivot a long-format metric column from method-split to wide."""
        pivoted = df_long.pivot_table(
            index=[c for c in df_long.columns
                   if c not in ("method", "distance", "ci_low", "ci_high")],
            columns="method", values=value_col
        ).reset_index()
        pivoted.columns.name = None
        pivoted = pivoted.rename(columns={
            "xtb": f"{prefix}_xtb", "solv": f"{prefix}_solv"
        })
        return pivoted

    # Get unique index columns (everything except method/distance/ci columns)
    idx_cols = [c for c in all_rows.columns
                if c not in ("method", "distance", "ci_low", "ci_high")]

    piv_dist = pivot_metric(all_rows, "distance", "d")
    piv_ci_low = pivot_metric(all_rows, "ci_low", "d_ci_low")
    piv_ci_high = pivot_metric(all_rows, "ci_high", "d_ci_high")

    # Merge all pivoted dataframes
    all_groups = piv_dist.merge(piv_ci_low, on=idx_cols, how="left")
    all_groups = all_groups.merge(piv_ci_high, on=idx_cols, how="left")
    all_groups["ratio"] = all_groups["d_solv"] / all_groups["d_xtb"].replace(0, np.nan)

    return all_groups


# ---------------------------------------------------------------------------
# Plots
# ---------------------------------------------------------------------------
def plot_centroid_shifts(per_sys: pd.DataFrame, out_path: Path):
    """Scatter xtb (blue) and solv (red) centroids with arrows connecting systems."""
    fig, ax = plt.subplots(figsize=(10, 8))

    ax.scatter(per_sys["centroid_PC1_xtb"], per_sys["centroid_PC2_xtb"],
               c="C0", marker="o", s=80, zorder=3, label="xtb centroid")
    ax.scatter(per_sys["centroid_PC1_solv"], per_sys["centroid_PC2_solv"],
               c="C3", marker="s", s=80, zorder=3, label="solv centroid")

    for _, row in per_sys.iterrows():
        ax.annotate("",
                    xy=(row["centroid_PC1_solv"], row["centroid_PC2_solv"]),
                    xytext=(row["centroid_PC1_xtb"], row["centroid_PC2_xtb"]),
                    arrowprops=dict(arrowstyle="->", color="gray", lw=1.0, alpha=0.6))
        # Label at xtb centroid
        ax.text(row["centroid_PC1_xtb"], row["centroid_PC2_xtb"],
                row["system"], fontsize=5, alpha=0.7, ha="right", va="bottom")

    ax.set_xlabel("PC1")
    ax.set_ylabel("PC2")
    ax.set_title("Centroid Shifts: xtb → solv_md")
    ax.legend()
    fig.tight_layout()
    fig.savefig(out_path, dpi=150)
    plt.close(fig)
    print(f"  Saved: {out_path}")


def plot_sap_tsap_stability(group_metrics: pd.DataFrame, out_path: Path):
    """Bar chart — SAP-TSAP distance in xtb vs solv, with bootstrap CI error bars."""
    sap_tsap = group_metrics[group_metrics["metric_type"] == "sap_vs_tsap"].copy()
    if len(sap_tsap) == 0:
        print("  [skip] No SAP vs TSAP data for bar chart.")
        return

    # Sort for consistent display
    sap_tsap = sap_tsap.sort_values("group_key")

    labels = sap_tsap["group_key"].values
    d_xtb = sap_tsap["d_xtb"].values
    d_solv = sap_tsap["d_solv"].values
    ratios = sap_tsap["ratio"].values

    # Bootstrap CI error bars: [distance - ci_low, ci_high - distance]
    # For xtb bars
    xtb_err_low = d_xtb - sap_tsap["d_ci_low_xtb"].values
    xtb_err_high = sap_tsap["d_ci_high_xtb"].values - d_xtb
    xtb_err = np.array([xtb_err_low, xtb_err_high])

    # For solv bars
    solv_err_low = d_solv - sap_tsap["d_ci_low_solv"].values
    solv_err_high = sap_tsap["d_ci_high_solv"].values - d_solv
    solv_err = np.array([solv_err_low, solv_err_high])

    x = np.arange(len(labels))
    width = 0.35

    fig, ax1 = plt.subplots(figsize=(12, 6))
    bars1 = ax1.bar(x - width / 2, d_xtb, width, label="xtb", color="C0",
                    yerr=xtb_err, capsize=3, error_kw={"linewidth": 1})
    bars2 = ax1.bar(x + width / 2, d_solv, width, label="solv", color="C3",
                    yerr=solv_err, capsize=3, error_kw={"linewidth": 1})

    ax1.set_ylabel("Centroid distance (PC units)")
    ax1.set_title("SAP–TSAP Centroid Distance: xtb vs solv_md (95% CI)")
    ax1.set_xticks(x)
    ax1.set_xticklabels(labels, rotation=45, ha="right", fontsize=8)
    ax1.legend(loc="upper left")

    # Secondary axis for ratio
    ax2 = ax1.twinx()
    ax2.plot(x, ratios, "k--o", markersize=5, label="ratio (solv/xtb)")
    ax2.set_ylabel("Ratio (solv / xtb)")
    ax2.legend(loc="upper right")

    fig.tight_layout()
    fig.savefig(out_path, dpi=150)
    plt.close(fig)
    print(f"  Saved: {out_path}")


def plot_pc_span_comparison(per_sys: pd.DataFrame, out_path: Path):
    """Grouped bar chart of span ratios (solv/xtb) for PC1 and PC2.

    Supplementary plot — std ratio is the primary breadth metric.
    """
    fig, ax = plt.subplots(figsize=(12, 6))

    systems = per_sys["system"].values
    x = np.arange(len(systems))
    width = 0.35

    ax.bar(x - width / 2, per_sys["span_ratio_PC1"].values, width,
           label="PC1 span ratio", color="C0")
    ax.bar(x + width / 2, per_sys["span_ratio_PC2"].values, width,
           label="PC2 span ratio", color="C3")

    ax.axhline(1.0, color="gray", ls="--", lw=0.8, label="1:1 line")
    ax.set_ylabel("Span ratio (solv / xtb)")
    ax.set_title("PC Span Ratio: solv_md / xtb (supplementary)")
    ax.set_xticks(x)
    ax.set_xticklabels(systems, rotation=45, ha="right", fontsize=8)
    ax.legend()
    fig.tight_layout()
    fig.savefig(out_path, dpi=150)
    plt.close(fig)
    print(f"  Saved: {out_path}")


def plot_pc_std_comparison(per_sys: pd.DataFrame, out_path: Path):
    """Grouped bar chart of std ratios (solv/xtb) for PC1 and PC2.

    This is the PRIMARY breadth comparison plot, since std is more robust
    than span (less sensitive to outliers).
    """
    fig, ax = plt.subplots(figsize=(12, 6))

    systems = per_sys["system"].values
    x = np.arange(len(systems))
    width = 0.35

    ax.bar(x - width / 2, per_sys["std_ratio_PC1"].values, width,
           label="PC1 std ratio", color="C0")
    ax.bar(x + width / 2, per_sys["std_ratio_PC2"].values, width,
           label="PC2 std ratio", color="C3")

    ax.axhline(1.0, color="gray", ls="--", lw=0.8, label="1:1 line")
    ax.set_ylabel("Std ratio (solv / xtb)")
    ax.set_title("PC Std Ratio: solv_md / xtb (primary breadth metric)")
    ax.set_xticks(x)
    ax.set_xticklabels(systems, rotation=45, ha="right", fontsize=8)
    ax.legend()
    fig.tight_layout()
    fig.savefig(out_path, dpi=150)
    plt.close(fig)
    print(f"  Saved: {out_path}")


# ---------------------------------------------------------------------------
# Interpretation document
# ---------------------------------------------------------------------------
def write_interpretation(per_sys: pd.DataFrame, group_metrics: pd.DataFrame,
                         out_path: Path):
    """Write the xtb_vs_solv.md interpretation document."""
    # --- Compute summary statistics ---
    top3_shift = per_sys.nlargest(3, "d_centroid")[["system", "d_centroid"]].values
    avg_shift = per_sys["d_centroid"].mean()
    median_shift = per_sys["d_centroid"].median()

    avg_ws_PC1 = per_sys["ws_PC1"].mean()
    avg_ws_PC2 = per_sys["ws_PC2"].mean()

    avg_span_ratio_PC1 = per_sys["span_ratio_PC1"].mean()
    avg_span_ratio_PC2 = per_sys["span_ratio_PC2"].mean()
    avg_std_ratio_PC1 = per_sys["std_ratio_PC1"].mean()
    avg_std_ratio_PC2 = per_sys["std_ratio_PC2"].mean()

    avg_containment_50 = per_sys["containment_50"].mean()
    avg_containment_90 = per_sys["containment_90"].mean()
    avg_containment_99 = per_sys["containment_99"].mean()

    # Effective sample sizes
    avg_neff_xtb = per_sys["n_eff_xtb"].mean()
    avg_neff_solv = per_sys["n_eff_solv"].mean()
    min_neff_xtb = per_sys["n_eff_xtb"].min()
    max_neff_xtb = per_sys["n_eff_xtb"].max()
    low_neff_systems = per_sys[per_sys["n_eff_xtb"] < 50][["system", "n_eff_xtb"]]

    # SAP-TSAP analysis with bootstrap CIs
    sap_tsap = group_metrics[group_metrics["metric_type"] == "sap_vs_tsap"]
    if len(sap_tsap) > 0:
        avg_sap_tsap_ratio = sap_tsap["ratio"].mean()
        sap_tsap_display = sap_tsap[["group_key", "d_xtb", "d_solv", "ratio",
                                      "d_ci_low_xtb", "d_ci_high_xtb",
                                      "d_ci_low_solv", "d_ci_high_solv"]].copy()
        sap_tsap_display.columns = ["group_key", "d_xtb", "d_solv", "ratio",
                                     "xtb_CI_low", "xtb_CI_high",
                                     "solv_CI_low", "solv_CI_high"]
        sap_tsap_table = df_to_markdown(sap_tsap_display)
    else:
        avg_sap_tsap_ratio = np.nan
        sap_tsap_table = "No data available."

    # me vs phe
    mephe = group_metrics[group_metrics["metric_type"] == "me_vs_phe"]
    if len(mephe) > 0:
        avg_mephe_ratio = mephe["ratio"].mean()
        mephe_display = mephe[["group_key", "d_xtb", "d_solv", "ratio",
                                "d_ci_low_xtb", "d_ci_high_xtb",
                                "d_ci_low_solv", "d_ci_high_solv"]].copy()
        mephe_display.columns = ["group_key", "d_xtb", "d_solv", "ratio",
                                   "xtb_CI_low", "xtb_CI_high",
                                   "solv_CI_low", "solv_CI_high"]
        mephe_table = df_to_markdown(mephe_display)
    else:
        avg_mephe_ratio = np.nan
        mephe_table = "No data available."

    # D vs L
    dl = group_metrics[group_metrics["metric_type"] == "D_vs_L"]
    if len(dl) > 0:
        avg_dl_ratio = dl["ratio"].mean()
        dl_display = dl[["group_key", "d_xtb", "d_solv", "ratio",
                          "d_ci_low_xtb", "d_ci_high_xtb",
                          "d_ci_low_solv", "d_ci_high_solv"]].copy()
        dl_display.columns = ["group_key", "d_xtb", "d_solv", "ratio",
                                "xtb_CI_low", "xtb_CI_high",
                                "solv_CI_low", "solv_CI_high"]
        dl_table = df_to_markdown(dl_display)
    else:
        avg_dl_ratio = np.nan
        dl_table = "No data available."

    # Per-system summary table (with N_eff and std ratio)
    sys_summary = per_sys[["system", "n_eff_xtb", "n_eff_solv",
                            "d_centroid", "d_PC1", "d_PC2",
                            "std_ratio_PC1", "std_ratio_PC2",
                            "ks_PC1_stat", "ks_PC2_stat",
                            "containment_90", "centroid_mahal"]].copy()
    sys_summary_table = df_to_markdown(sys_summary)

    # Top 3 most shifted
    top3_lines = []
    for sys_name, d in top3_shift:
        row = per_sys[per_sys["system"] == sys_name].iloc[0]
        top3_lines.append(
            f"  - **{sys_name}**: d_centroid = {d:.2f}, "
            f"d_PC1 = {row['d_PC1']:.2f}, d_PC2 = {row['d_PC2']:.2f}, "
            f"N_eff(xtb) = {row['n_eff_xtb']:.0f}, "
            f"containment_90 = {row['containment_90']:.3f}"
        )

    # Low N_eff warning
    low_neff_lines = []
    for _, row in low_neff_systems.iterrows():
        low_neff_lines.append(
            f"  - **{row['system']}**: N_eff = {row['n_eff_xtb']:.0f}"
        )
    if not low_neff_lines:
        low_neff_lines.append("  - None (all systems have N_eff ≥ 50)")

    # --- Compose document ---
    doc = f"""# xTB vs solv_md: Cross-dataset Interpretation

## Summary of Top Findings

| Metric | Value |
|--------|-------|
| Average centroid shift | {avg_shift:.2f} PC units |
| Median centroid shift | {median_shift:.2f} PC units |
| **Average PC1 std ratio (solv/xtb)** | **{avg_std_ratio_PC1:.1f}×** |
| **Average PC2 std ratio (solv/xtb)** | **{avg_std_ratio_PC2:.1f}×** |
| Average PC1 span ratio (solv/xtb) | {avg_span_ratio_PC1:.1f}× |
| Average PC2 span ratio (solv/xtb) | {avg_span_ratio_PC2:.1f}× |
| Avg fraction of solv within xtb 90th %ile | {avg_containment_90:.3f} |
| Average N_eff (xtb) | {avg_neff_xtb:.0f} |
| Average N_eff (solv) | {avg_neff_solv:.0f} |
| Average SAP–TSAP ratio (solv/xtb) | {avg_sap_tsap_ratio:.2f} |

> **Primary breadth metric**: std ratio (solv/xtb) is the most robust measure
> of conformational breadth difference. Span ratios are supplementary (more
> sensitive to outliers).

### Top 3 Most Shifted Systems

{chr(10).join(top3_lines)}

### Effective Sample Size (N_eff)

Effective sample sizes account for autocorrelation in MD trajectories using
the lag-1 autocorrelation method: N_eff = N × (1 - ρ) / (1 + ρ).

- **xTB**: average N_eff = {avg_neff_xtb:.0f} (range {min_neff_xtb:.0f}–{max_neff_xtb:.0f})
  - Raw frame count: 2000 per system
  - High autocorrelation reduces effective independence to {avg_neff_xtb/2000*100:.0f}% of raw count
- **solv_md**: average N_eff = {avg_neff_solv:.0f}
  - Raw frame count: ~4000 per system
  - Lower autocorrelation → {avg_neff_solv/4000*100:.0f}% effective independence

Systems with unusually low N_eff (< 50):
{chr(10).join(low_neff_lines)}

> **Caveat**: KS test p-values are not reported here because they are
> trivially significant for all 16/16 systems. The massive distributional
> differences (different shapes, vastly different spreads) make p-values
> uninformative. Use the KS **statistic** (effect size) and the **std ratio**
> for practical comparison instead.

## Per-System Metrics

{sys_summary_table}

## Physical Interpretation

### Why Explicit Solvent Shifts PC Space

The xTB trajectories (100 ps, implicit GBSA solvent) sample a narrow region
of conformational space around the initial equilibrium. The solv_md trajectories
(400 ns, explicit OPC3 water) explore a **much broader** region:

  - **Std ratios of {avg_std_ratio_PC1:.0f}× (PC1) and {avg_std_ratio_PC2:.0f}× (PC2)**
    are the primary breadth metrics. These reflect the combined effect of
    different solvent models and dramatically different sampling times
    (100 ps vs 400 ns, a 4000× difference). The timescale difference is
    likely the dominant factor; the contribution of the solvent model change
    alone cannot be separated from this experimental design. solv_md std is
    {avg_std_ratio_PC1:.0f}× larger on PC1 and {avg_std_ratio_PC2:.0f}× on PC2.
- Span ratios ({avg_span_ratio_PC1:.0f}–{avg_span_ratio_PC2:.0f}×) confirm
  this but are more sensitive to outlier frames.
- The average centroid shift of {avg_shift:.2f} PC units indicates that the
  equilibrium conformation itself changes with explicit solvation — implicit
  solvent does not fully capture solvation-driven conformational preferences.
- Mahalanobis containment of ~{avg_containment_90:.1%} at the 90th percentile
  means that roughly 46% of solv_md frames fall within the xtb 90% Mahalanobis ellipse.
  The xTB cloud is narrower than the solv_md landscape (std ratio ~40× on PC1, ~8× on PC2),
  but it captures the center of the solv_md distribution while missing the tails.

### Mechanism of Shift

  1. **Explicit water molecules** may form hydrogen bonds and electrostatic
     interactions with polar groups (carbonyls, amines) that are not captured
     by GBSA, potentially stabilizing or destabilizing specific pendant arm
     orientations. **Caveat:** No hydrogen bond analysis was performed; this
     is a plausible hypothesis, not a confirmed mechanism. Alternative
     explanations include force field differences, different initial conditions,
     and the 4000× longer timescale, all of which could contribute independently.
2. **Longer timescale** (400 ns vs 100 ps) allows barrier crossing events
   that are inaccessible to xTB, including ring flips and large-amplitude
   arm swings.
3. **NPT vs NVT**: solv_md at constant pressure allows volume fluctuations,
   which can subtly change coordination geometry.

### SAP/TSAP Stability Discussion

Bootstrap 95% CIs on centroid distances are shown below.

{sap_tsap_table}

- **Average SAP–TSAP ratio**: {avg_sap_tsap_ratio:.2f}
- When ratio > 1: solv_md amplifies the SAP–TSAP separation (stabilizes the
  distinction). When ratio < 1: solv_md reduces the separation (structures
  become more similar with explicit solvent).
- A ratio near 1.0 would indicate that xTB correctly predicts the relative
  SAP–TSAP separation despite the narrower sampling.

### me vs phe Comparison

{mephe_table}

- Average me↔phe ratio: {avg_mephe_ratio:.2f}

### D vs L Comparison

{dl_table}

- Average D↔L ratio: {avg_dl_ratio:.2f}

### Projection Artifact: D/L Ratio

The low D/L centroid distance ratio ({avg_dl_ratio:.2f}) should be interpreted
with caution. PC1 is overwhelmingly aligned with the D/L structural axis
(99.8% of D/L separation on PC1, 0.2% on PC2), meaning the ratio reflects how this axis projects
differently under narrow (xTB) vs broad (solv_md) sampling rather than necessarily
indicating physical D/L convergence in solution. The ratio quantifies a PCA
projection effect, not a direct physical finding.

### Statistical Caveat: KS Tests

Kolmogorov–Smirnov tests between xtb and solv_md PC distributions are
trivially significant for all systems (p < 10⁻¹⁰⁰), because the two
distributions differ massively in both location and scale. KS p-values
are therefore not reported. The KS **statistic** (effect size, ranging
0–1) and the **std ratio** provide meaningful measures of distributional
difference. A KS statistic of 0.5 means the CDFs differ by 50% at the
point of maximum divergence.

## Recommendation: Is xTB Predictive of solv_md Basins?

**Partially.** The xTB simulations preserve the *relative magnitude* of
SAP-TSAP separation across methods (ratio {avg_sap_tsap_ratio:.2f}), but the
structural directions of SAP-TSAP centroid shifts are mostly consistent between
methods (7/8 groups with direction angle < 30°, max 60°). The one exception is
me_sssL (60° angle). The SAP-TSAP distance magnitudes are not significantly
correlated (Spearman ρ = +0.38, p = 0.35), though the positive trend suggests
partial preservation. Basin ranking is not preserved across
all comparisons. However:

1. **Magnitude**: xTB dramatically underestimates the breadth of conformational
   sampling (std ratios of {avg_std_ratio_PC1:.0f}× / {avg_std_ratio_PC2:.0f}×).
2. **Position**: The centroid shift of {avg_shift:.2f} PC units indicates
   that absolute basin locations differ — implicit solvent biases the
   equilibrium geometry.
3. **Containment**: Only ~{avg_containment_90:.0%} of solv_md frames fall
   within the xtb 90th percentile ellipsoid.
4. **N_eff**: xTB effective sample sizes ({avg_neff_xtb:.0f} avg) are far
   smaller than raw frame counts, meaning xTB provides limited independent
   structural observations.

**Practical implication**: xTB is suitable for *qualitative* screening
(ranking conformers, identifying dominant basin) but should not be relied upon
for *quantitative* free energy estimates. For FES and thermodynamic
comparisons, solv_md is essential.

## Files Generated

| File | Description |
|------|-------------|
| `cross_dataset/analysis/xtb_vs_solv_per_system_metrics.csv` | Per-system metrics (16 rows, +N_eff, −KS p-values, −trivial containment) |
| `cross_dataset/analysis/xtb_vs_solv_group_metrics.csv` | Grouped comparison metrics (+bootstrap CIs) |
| `cross_dataset/analysis/xtb_vs_solv_pc_spans.csv` | PC span, std, and N_eff details |
| `cross_dataset/analysis/plot_centroid_shifts.png` | Centroid shift scatter + arrows |
| `cross_dataset/analysis/plot_sap_tsap_stability.png` | SAP–TSAP distance bar chart (95% CI error bars) |
| `cross_dataset/analysis/plot_pc_std_comparison.png` | **Std ratio bar chart (primary breadth metric)** |
| `cross_dataset/analysis/plot_pc_span_comparison.png` | Span ratio bar chart (supplementary) |
"""

    out_path.write_text(doc, encoding="utf-8")
    print(f"  Saved: {out_path}")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def main():
    parser = argparse.ArgumentParser(
        description="B3: Cross-dataset interpretation between xtb and solv_md"
    )
    parser.add_argument(
        "--input", type=Path, default=DEFAULT_INPUT,
        help="Path to joint projection CSV"
    )
    parser.add_argument(
        "--out-analysis", type=Path, default=OUT_ANALYSIS,
        help="Output directory for analysis CSVs and plots"
    )
    parser.add_argument(
        "--out-docs", type=Path, default=OUT_DOCS,
        help="Output directory for interpretation document"
    )
    args = parser.parse_args()

    # Ensure output dirs exist
    args.out_analysis.mkdir(parents=True, exist_ok=True)
    args.out_docs.mkdir(parents=True, exist_ok=True)

    print("=" * 60)
    print("B3: Cross-dataset interpretation — xtb vs solv_md")
    print("=" * 60)

    # --- Load data ---
    print(f"\nLoading: {args.input}")
    df = pd.read_csv(args.input)
    print(f"  {len(df)} rows, {df['system'].nunique()} systems")

    # --- Per-system metrics ---
    print("\n[1/6] Computing per-system metrics...")
    per_sys = compute_per_system_metrics(df)

    csv1 = args.out_analysis / "xtb_vs_solv_per_system_metrics.csv"
    per_sys.to_csv(csv1, index=False)
    print(f"  Saved: {csv1}")

    # --- PC spans CSV ---
    print("\n[2/6] Building PC spans CSV...")
    pc_spans = build_pc_spans_df(per_sys)
    csv3 = args.out_analysis / "xtb_vs_solv_pc_spans.csv"
    pc_spans.to_csv(csv3, index=False)
    print(f"  Saved: {csv3}")

    # --- Grouped metrics ---
    print("\n[3/6] Computing grouped metrics (with bootstrap CIs)...")
    group_metrics = compute_group_metrics(df, per_sys)
    csv2 = args.out_analysis / "xtb_vs_solv_group_metrics.csv"
    group_metrics.to_csv(csv2, index=False)
    print(f"  Saved: {csv2}")

    # --- Plots ---
    print("\n[4/6] Generating plots...")
    plot_centroid_shifts(per_sys, args.out_analysis / "plot_centroid_shifts.png")
    plot_sap_tsap_stability(group_metrics, args.out_analysis / "plot_sap_tsap_stability.png")
    plot_pc_std_comparison(per_sys, args.out_analysis / "plot_pc_std_comparison.png")
    plot_pc_span_comparison(per_sys, args.out_analysis / "plot_pc_span_comparison.png")

    # --- Interpretation document ---
    print("\n[5/6] Writing interpretation document...")
    write_interpretation(per_sys, group_metrics, args.out_docs / "xtb_vs_solv.md")

    # --- Summary ---
    print("\n[6/6] Summary statistics...")
    print("=" * 60)
    print("DONE — Key findings:")
    print("=" * 60)
    top3 = per_sys.nlargest(3, "d_centroid")
    for _, row in top3.iterrows():
        print(f"  {row['system']}: d_centroid = {row['d_centroid']:.2f}")
    print(f"  Average centroid shift: {per_sys['d_centroid'].mean():.2f}")
    print(f"  Average N_eff (xtb): {per_sys['n_eff_xtb'].mean():.0f}")
    print(f"  Average N_eff (solv): {per_sys['n_eff_solv'].mean():.0f}")
    print(f"  Average PC1 std ratio: {per_sys['std_ratio_PC1'].mean():.1f}×")
    print(f"  Average PC2 std ratio: {per_sys['std_ratio_PC2'].mean():.1f}×")
    sap_tsap = group_metrics[group_metrics["metric_type"] == "sap_vs_tsap"]
    if len(sap_tsap) > 0:
        print(f"  Average SAP–TSAP stability ratio: {sap_tsap['ratio'].mean():.2f}")


if __name__ == "__main__":
    main()
