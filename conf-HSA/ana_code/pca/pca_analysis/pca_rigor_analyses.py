#!/usr/bin/env python3
"""
Statistical rigor analyses for PCA results.

Runs:
1. PCA convergence analysis (subset overlap — RQ-COMP2/RQ-RIG1)
2. Block-bootstrap eigenvalue confidence intervals (RQ-COMP3)
3. Autocorrelation / effective sample size (RQ-COMP5/RQ-RIG2)
4. Two-way ANOVA for stereo × chirality interaction (RQ-COMP6)
5. Leave-one-out correlation analysis (RQ-RIG3)
6. ΔΔG t-tests for SAP/TSAP pairs (RQ-RIG5)
7. Per-system DCCM convergence test (RQ-REP5)

Usage:
    python pca_rigor_analyses.py --input-dir pca_analysis/ --output-dir pca_analysis/
    python pca_rigor_analyses.py --loo              # Run only LOO analysis
    python pca_rigor_analyses.py --ddg-ttest         # Run only ΔΔG t-tests
    python pca_rigor_analyses.py --dccm-convergence  # Run only DCCM convergence
    python pca_rigor_analyses.py --all               # Run all analyses (1–7)
"""

import argparse
import json
import os
import sys
import time

import numpy as np
import pandas as pd

# Add parent dir for imports
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from pca_utils import R_GAS

np.random.seed(42)  # Reproducibility


# ---------------------------------------------------------------------------
# 1. PCA Convergence Analysis
# ---------------------------------------------------------------------------

def run_convergence_analysis(data_dir, output_dir):
    """PCA convergence: recompute on 25%, 50%, 75%, 100% of frames and measure eigenvector overlap.

    Uses sklearn's randomized SVD for efficient top-k computation on
    large data matrices, avoiding the O(d^3) cost of full eigendecomposition.

    RMSIP (Root Mean Square Inner Product) measures subspace overlap:
        RMSIP = sqrt(1/k * sum_i (u_i^T v_i)^2)
    where u_i and v_i are the i-th eigenvectors from the subset and full PCA.

    RMSIP > 0.9 for top 3 PCs indicates stable eigenvectors.
    """
    from sklearn.decomposition import PCA

    print("  Loading aligned coordinates and reference eigenvectors...")
    aligned = np.load(os.path.join(data_dir, 'aligned_coords.npy'), mmap_mode='r')
    eigvecs_full = np.load(os.path.join(data_dir, 'eigenvectors.npy'), mmap_mode='r')

    n_total = aligned.shape[0]
    n_dims = aligned.shape[1]

    # Load full eigenvalues for reference variance fractions
    eigvals_full = np.load(os.path.join(data_dir, 'eigenvalues.npy'))
    total_var_full = eigvals_full.sum()

    fracs = [0.25, 0.50, 0.75, 1.00]
    results = []

    for frac in fracs:
        n_sub = int(n_total * frac)
        # Sample evenly spaced frames to avoid temporal bias
        indices = np.linspace(0, n_total - 1, n_sub, dtype=int)
        sub_coords = aligned[indices].astype(np.float32)

        # Use sklearn PCA with randomized SVD (efficient for tall-thin matrices)
        pca = PCA(n_components=10, svd_solver='randomized', random_state=42)
        pca.fit(sub_coords)

        # Eigenvectors from PCA: pca.components_ shape (10, n_dims), rows = eigenvectors
        # Full eigenvectors: columns = eigenvectors
        eigvecs_sub = pca.components_.T  # (n_dims, 10)

        # Eigenvalues from PCA: pca.explained_variance_
        eigvals = pca.explained_variance_

        # Eigenvector overlap (RMSIP) with full PCA for top k PCs
        for k in [3, 5, 10]:
            v_sub = eigvecs_sub[:, :k]
            v_full_ref = np.array(eigvecs_full[:, :k])
            overlap_matrix = v_full_ref.T @ v_sub
            rmsip = np.sqrt(np.mean(np.diag(overlap_matrix) ** 2))
            results.append({
                'fraction': frac,
                'n_frames': n_sub,
                'n_pcs': k,
                'rmsip': rmsip,
                'variance_pc': None,
                'var_frac': None,
            })

        # Variance fraction for PC1-3
        # Use explained_variance_ratio_ which normalizes by total variance
        # (consistent with full-data reference fractions)
        var_ratios = pca.explained_variance_ratio_
        for i, label in enumerate(['PC1', 'PC2', 'PC3']):
            results.append({
                'fraction': frac,
                'n_frames': n_sub,
                'n_pcs': 0,
                'rmsip': None,
                'variance_pc': label,
                'var_frac': var_ratios[i],
            })

        # Print progress
        rmsip_k3 = results[-6]['rmsip']
        pc1_varfrac = results[-3]['var_frac']
        print(f"    frac={frac:.2f} ({n_sub} frames): "
              f"RMSIP(top3)={rmsip_k3:.4f}, "
              f"PC1 var={pc1_varfrac:.4f}")

        # Add reference (full-data) variance fractions (only for 100%)
        if frac == 1.0:
            for i, label in enumerate(['PC1', 'PC2', 'PC3']):
                results.append({
                    'fraction': frac,
                    'n_frames': n_total,
                    'n_pcs': -1,  # flag: reference
                    'rmsip': None,
                    'variance_pc': f'{label}_ref',
                    'var_frac': eigvals_full[i] / total_var_full,
                })

    df = pd.DataFrame(results)
    out_path = os.path.join(output_dir, 'convergence_analysis.csv')
    df.to_csv(out_path, index=False)
    print(f"  Convergence analysis saved to {out_path}")
    return df


# ---------------------------------------------------------------------------
# 2. Block-Bootstrap for Eigenvalue CIs
# ---------------------------------------------------------------------------

def run_bootstrap_eigenvalues(data_dir, output_dir, n_bootstrap=500, block_size=200):
    """Block-bootstrap for eigenvalue confidence intervals.

    Resamples blocks of frames (preserving temporal correlation structure)
    and recomputes PCA to estimate uncertainty in eigenvalue fractions.
    Uses sklearn's randomized SVD for efficient computation.
    """
    from sklearn.decomposition import PCA

    print("  Loading aligned coordinates for bootstrap...")
    aligned = np.load(os.path.join(data_dir, 'aligned_coords.npy'), mmap_mode='r')
    n_total = aligned.shape[0]

    # Subsample for speed: every 8th frame ≈ 8k frames
    subsample_idx = np.arange(0, n_total, 8)
    coords = aligned[subsample_idx].astype(np.float32)
    n_sub = len(coords)
    print(f"    Subsampled to {n_sub} frames (every 8th)")

    n_blocks = n_sub // block_size
    boot_fracs = {f'PC{i+1}': [] for i in range(5)}

    # Compute observed fractions on the same subsampled data for fair comparison
    pca_obs = PCA(n_components=5, svd_solver='randomized', random_state=42)
    pca_obs.fit(coords)
    observed_ratios = pca_obs.explained_variance_ratio_
    print(f"    Observed fractions on subsample: "
          f"PC1={observed_ratios[0]:.4f}, PC2={observed_ratios[1]:.4f}")

    rng = np.random.RandomState(42)  # Fixed seed for reproducibility

    t0 = time.time()
    for b in range(n_bootstrap):
        # Sample blocks with replacement
        block_starts = rng.randint(0, n_blocks, size=n_blocks)
        boot_idx = np.concatenate([
            np.arange(s * block_size, min((s + 1) * block_size, n_sub))
            for s in block_starts
        ])
        boot_coords = coords[boot_idx]

        # Use sklearn PCA for efficient top-5 computation
        pca = PCA(n_components=5, svd_solver='randomized', random_state=42)
        pca.fit(boot_coords)

        # Use explained_variance_ratio_ which normalizes by TOTAL variance
        # (not just top-5 sum), ensuring consistent denominator with observed
        ratios = pca.explained_variance_ratio_

        for i in range(5):
            boot_fracs[f'PC{i+1}'].append(ratios[i])

        if (b + 1) % 50 == 0:
            elapsed = time.time() - t0
            rate = (b + 1) / elapsed
            eta = (n_bootstrap - b - 1) / rate
            print(f"    Bootstrap {b+1}/{n_bootstrap} "
                  f"({elapsed:.0f}s elapsed, ETA {eta:.0f}s)")

    # Compute 95% CIs
    # Also include full-data fractions for reference
    eigvals_full = np.load(os.path.join(data_dir, 'eigenvalues.npy'))
    total_full = eigvals_full.sum()

    results = []

    for i in range(5):
        label = f'PC{i+1}'
        observed = observed_ratios[i]
        boot_vals = np.array(boot_fracs[label])
        ci_lo = np.percentile(boot_vals, 2.5)
        ci_hi = np.percentile(boot_vals, 97.5)
        results.append({
            'PC': label,
            'observed_frac_subsampled': observed,
            'full_data_frac': eigvals_full[i] / total_full,
            'CI_2.5': ci_lo,
            'CI_97.5': ci_hi,
            'CI_width': ci_hi - ci_lo,
        })

    df = pd.DataFrame(results)
    out_path = os.path.join(output_dir, 'bootstrap_eigenvalue_ci.csv')
    df.to_csv(out_path, index=False)
    print(f"  Bootstrap eigenvalue CIs saved to {out_path}")
    return df


# ---------------------------------------------------------------------------
# 3. Autocorrelation / Effective Sample Size
# ---------------------------------------------------------------------------

def run_autocorrelation_analysis(data_dir, output_dir, max_lag=500):
    """Compute autocorrelation time and effective sample size for each system's PC trajectories.

    Uses FFT-based autocorrelation (fast) and 1/e threshold for integrated
    autocorrelation time. Effective sample size: N_eff = N / (2*tau + 1).
    """
    print("  Loading projections and system metadata...")
    projections = pd.read_csv(os.path.join(data_dir, 'projections_all.csv'))

    with open(os.path.join(data_dir, 'system_metadata.json')) as f:
        systems_meta = json.load(f)

    results = []
    for sys_info in systems_meta:
        sys_name = sys_info['name']
        sys_proj = projections[projections['system_label'] == sys_name]
        n_frames = len(sys_proj)

        for pc in ['PC1', 'PC2', 'PC3']:
            traj = sys_proj[pc].values

            # Compute autocorrelation using FFT (fast)
            x = traj - traj.mean()
            n = len(x)
            fft_x = np.fft.fft(x, n=2 * n)
            acf = np.fft.ifft(fft_x * np.conj(fft_x))[:n].real
            acf = acf / acf[0]  # Normalize

            # Find first crossing of 1/e (integrated autocorrelation time)
            threshold = 1.0 / np.e
            tau = 1  # Minimum: 1 frame (no correlation)
            for lag in range(1, min(max_lag, n)):
                if acf[lag] < threshold:
                    tau = lag
                    break
            else:
                tau = max_lag  # Didn't cross threshold

            # N_eff = n_frames / (2*tau + 1)
            n_eff = n_frames / (2 * tau + 1)

            results.append({
                'system': sys_name,
                'PC': pc,
                'n_frames': n_frames,
                'tau_autocorr': tau,
                'n_eff': n_eff,
                'n_eff_ratio': n_eff / n_frames,
            })

    df = pd.DataFrame(results)
    out_path = os.path.join(output_dir, 'autocorrelation_analysis.csv')
    df.to_csv(out_path, index=False)
    print(f"  Autocorrelation analysis saved to {out_path}")

    # Print summary
    print("\n  === Effective Sample Size Summary ===")
    for pc in ['PC1', 'PC2', 'PC3']:
        pc_df = df[df['PC'] == pc]
        print(f"    {pc}:")
        print(f"      Mean tau_autocorr: {pc_df['tau_autocorr'].mean():.1f} frames")
        print(f"      Mean N_eff: {pc_df['n_eff'].mean():.0f} "
              f"(of {pc_df['n_frames'].mean():.0f} frames)")
        print(f"      Mean N_eff/N ratio: {pc_df['n_eff_ratio'].mean():.3f}")

    return df


# ---------------------------------------------------------------------------
# 4. Two-Way ANOVA for Stereo × Chirality Interaction
# ---------------------------------------------------------------------------

def run_two_way_anova(data_dir, output_dir):
    """Two-way ANOVA for stereo × chirality interaction on PC1, PC2, PC3.

    Uses system-level means (N=16, 4 per cell) with statsmodels.
    Tests whether the stereo × chirality interaction is significant.
    """
    import statsmodels.api as sm
    from statsmodels.formula.api import ols

    print("  Loading projections and system metadata...")
    projections = pd.read_csv(os.path.join(data_dir, 'projections_all.csv'))

    with open(os.path.join(data_dir, 'system_metadata.json')) as f:
        systems_meta = json.load(f)

    # Build DataFrame with stereo, chirality factors
    # Use system-level means (N=16) for the ANOVA
    sys_means = projections.groupby('system_label')[['PC1', 'PC2', 'PC3']].mean()

    anova_data = []
    for sys_info in systems_meta:
        name = sys_info['name']
        if name in sys_means.index:
            row = {
                'system': name,
                'stereo': sys_info['stereo'],       # sss or rrr
                'chirality': sys_info['chirality'],  # me or phe
            }
            for pc in ['PC1', 'PC2', 'PC3']:
                row[pc] = sys_means.loc[name, pc]
            anova_data.append(row)

    anova_df = pd.DataFrame(anova_data)
    print(f"    ANOVA design: {len(anova_df)} systems, "
          f"stereo={sorted(anova_df['stereo'].unique())}, "
          f"chirality={sorted(anova_df['chirality'].unique())}")

    results = []
    for pc in ['PC1', 'PC2', 'PC3']:
        model = ols(f'{pc} ~ C(stereo) * C(chirality)', data=anova_df).fit()
        anova_table = sm.stats.anova_lm(model, typ=2)

        for source in ['C(stereo)', 'C(chirality)', 'C(stereo):C(chirality)',
                        'Residual']:
            if source in anova_table.index:
                row_result = {
                    'PC': pc,
                    'source': source,
                    'sum_sq': anova_table.loc[source, 'sum_sq'],
                    'df': anova_table.loc[source, 'df'],
                    'F': anova_table.loc[source, 'F'],
                    'p_value': anova_table.loc[source, 'PR(>F)'],
                }
                # Handle NaN for Residual
                if source == 'Residual':
                    row_result['F'] = np.nan
                    row_result['p_value'] = np.nan
                results.append(row_result)

    df = pd.DataFrame(results)
    out_path = os.path.join(output_dir, 'anova_interaction.csv')
    df.to_csv(out_path, index=False)
    print(f"  Two-way ANOVA saved to {out_path}")

    # Print key findings
    print("\n  === Stereo × Chirality Interaction ===")
    interaction_df = df[df['source'] == 'C(stereo):C(chirality)']
    for _, row in interaction_df.iterrows():
        sig = "SIGNIFICANT" if row['p_value'] < 0.05 else "not significant"
        print(f"    {row['PC']}: F={row['F']:.2f}, p={row['p_value']:.4f} ({sig})")

    # Also print main effects
    print("\n  === Main Effects ===")
    for effect in ['C(stereo)', 'C(chirality)']:
        effect_df = df[df['source'] == effect]
        for _, row in effect_df.iterrows():
            sig = "SIGNIFICANT" if row['p_value'] < 0.05 else "not significant"
            print(f"    {row['PC']} ~ {effect}: F={row['F']:.2f}, "
                  f"p={row['p_value']:.4f} ({sig})")

    return df


# ---------------------------------------------------------------------------
# 5. Leave-One-Out Correlation Analysis (RQ-RIG3)
# ---------------------------------------------------------------------------

def run_leave_one_out_correlation(data_dir, output_dir):
    """Systematic leave-one-out for all 6 metric–PC correlation pairs.

    For each pair, removes each of the 16 systems one at a time,
    recomputes Pearson r, and reports the range and sign-flip count.
    This addresses the outlier sensitivity concern (e.g. me_rrrD_tsap
    flips ΔG–PC1 from r=−0.04 to r=+0.22).
    """
    from scipy import stats

    print("  Loading projections, MMPBSA, and contacts data...")
    projections = pd.read_csv(os.path.join(data_dir, 'projections_all.csv'))

    mmpbsa_path = os.path.join(data_dir, 'mmpbsa_energies_all_systems.csv')
    if not os.path.exists(mmpbsa_path):
        print(f"  Warning: {mmpbsa_path} not found. Skipping LOO correlation.")
        return None

    mmpbsa = pd.read_csv(mmpbsa_path)

    contacts_path = os.path.join(data_dir, 'contacts_summary_all_systems.csv')
    contacts = None
    if os.path.exists(contacts_path):
        contacts = pd.read_csv(contacts_path)

    # Build system-level data from projections
    sys_proj = projections.groupby('system_label')[['PC1', 'PC2', 'PC3']].mean()

    # Define the 6 pairs to test
    pairs = [
        ('TOTAL_corrected', 'PC1'), ('TOTAL_corrected', 'PC2'), ('TOTAL_corrected', 'PC3'),
        ('mean_total_contacts', 'PC1'), ('mean_total_contacts', 'PC2'), ('mean_total_contacts', 'PC3'),
    ]

    # Get system-level values
    systems = sorted(sys_proj.index.tolist())

    results = []
    for metric, pc in pairs:
        # Get metric values per system
        if metric == 'TOTAL_corrected':
            metric_vals = mmpbsa.set_index('system').reindex(systems)['TOTAL_corrected']
        elif metric == 'mean_total_contacts' and contacts is not None:
            metric_vals = contacts.set_index('system').reindex(systems)['mean_total_contacts']
        else:
            continue

        pc_vals = sys_proj.reindex(systems)[pc]

        # Full-sample correlation
        valid = ~(metric_vals.isna() | pc_vals.isna())
        r_full, p_full = stats.pearsonr(metric_vals[valid], pc_vals[valid])

        # Leave-one-out
        loo_r = []
        loo_p = []
        valid_indices = np.where(valid)[0]
        for vi in valid_indices:
            # Remove system at valid position vi
            x_loo = np.delete(metric_vals[valid].values, vi)
            y_loo = np.delete(pc_vals[valid].values, vi)
            if len(x_loo) >= 3:  # Need at least 3 points
                r, p = stats.pearsonr(x_loo, y_loo)
                loo_r.append(r)
                loo_p.append(p)
                left_out_sys = systems[vi]
                results.append({
                    'metric': metric,
                    'pc': pc,
                    'left_out': left_out_sys,
                    'r_loo': r,
                    'p_loo': p,
                    'r_full': r_full,
                    'p_full': p_full,
                    'delta_r': r - r_full,
                })

        # Summary statistics
        if len(loo_r) > 0:
            loo_r_arr = np.array(loo_r)
            sign_flips = int(np.sum(np.sign(loo_r_arr) != np.sign(r_full)))
            results.append({
                'metric': metric,
                'pc': pc,
                'left_out': 'SUMMARY',
                'r_loo': np.nan,
                'p_loo': np.nan,
                'r_full': r_full,
                'p_full': p_full,
                'delta_r': np.nan,
                'loo_r_min': loo_r_arr.min(),
                'loo_r_max': loo_r_arr.max(),
                'sign_flips': sign_flips,
                'n_systems': len(loo_r),
            })
            print(f"    {metric} vs {pc}: r={r_full:.3f}, "
                  f"LOO range=[{loo_r_arr.min():.3f}, {loo_r_arr.max():.3f}], "
                  f"sign flips={sign_flips}/{len(loo_r)}")

    df = pd.DataFrame(results)
    out_path = os.path.join(output_dir, 'loo_correlation_analysis.csv')
    df.to_csv(out_path, index=False)
    print(f"  LOO correlation analysis saved to {out_path}")
    return df


# ---------------------------------------------------------------------------
# 6. ΔΔG t-tests for SAP/TSAP Pairs (RQ-RIG5)
# ---------------------------------------------------------------------------

def run_ddg_ttests(data_dir, output_dir):
    """Formal Welch's t-tests for ΔΔG between SAP and TSAP for each pair.

    Since only system-level summary statistics are available (mean, SD, n_trials),
    the t-test is computed from summary statistics rather than raw trial data.

    With N=3–4 per group, power is very low. Even if ΔΔG ≈ 2.3 kJ/mol with
    SEM ≈ 7 kJ/mol, the t-test will not be significant (p > 0.5). This
    reinforces the "indistinguishable" conclusion with formal statistics.
    """
    from scipy import stats

    print("  Loading MMPBSA data...")
    mmpbsa_path = os.path.join(data_dir, 'mmpbsa_energies_all_systems.csv')
    mmpbsa = pd.read_csv(mmpbsa_path)

    # SAP/TSAP pairs to test
    SAP_TSAP_PAIRS = [
        ('phe_sssD_sap', 'phe_sssD_tsap'),
        ('phe_rrrL_sap', 'phe_rrrL_tsap'),
    ]

    results = []

    for sap_sys, tsap_sys in SAP_TSAP_PAIRS:
        sap_row = mmpbsa[mmpbsa['system'] == sap_sys]
        tsap_row = mmpbsa[mmpbsa['system'] == tsap_sys]

        if sap_row.empty or tsap_row.empty:
            print(f"  Warning: {sap_sys} or {tsap_sys} not in MMPBSA data. Skipping.")
            continue

        # Extract summary statistics
        mean_sap = sap_row['TOTAL_corrected'].values[0]
        mean_tsap = tsap_row['TOTAL_corrected'].values[0]
        sd_sap = sap_row['SD'].values[0]
        sd_tsap = tsap_row['SD'].values[0]
        n_sap = int(sap_row['n_trials'].values[0])
        n_tsap = int(tsap_row['n_trials'].values[0])

        if n_sap < 2 or n_tsap < 2:
            print(f"  Warning: {sap_sys} or {tsap_sys} has <2 trials. Skipping.")
            continue

        # Welch's t-test from summary statistics
        ddg = mean_sap - mean_tsap
        se_diff = np.sqrt(sd_sap**2 / n_sap + sd_tsap**2 / n_tsap)
        t_stat = ddg / se_diff

        # Welch-Satterthwaite degrees of freedom
        num = (sd_sap**2 / n_sap + sd_tsap**2 / n_tsap)**2
        denom = (sd_sap**2 / n_sap)**2 / (n_sap - 1) + (sd_tsap**2 / n_tsap)**2 / (n_tsap - 1)
        df_approx = num / denom if denom > 0 else np.nan

        # Two-tailed p-value
        p_value = 2 * stats.t.sf(np.abs(t_stat), df_approx) if not np.isnan(df_approx) else np.nan

        # 95% CI for ΔΔG
        if not np.isnan(df_approx) and df_approx > 0:
            t_crit = stats.t.ppf(0.975, df_approx)
            ci_lo = ddg - t_crit * se_diff
            ci_hi = ddg + t_crit * se_diff
        else:
            ci_lo = np.nan
            ci_hi = np.nan

        results.append({
            'pair': f"{sap_sys} vs {tsap_sys}",
            'sap_system': sap_sys,
            'tsap_system': tsap_sys,
            'n_sap': n_sap,
            'n_tsap': n_tsap,
            'mean_dg_sap': mean_sap,
            'mean_dg_tsap': mean_tsap,
            'ddg': ddg,
            'sd_sap': sd_sap,
            'sd_tsap': sd_tsap,
            'se_diff': se_diff,
            'ci_lo': ci_lo,
            'ci_hi': ci_hi,
            't_stat': t_stat,
            'df_approx': df_approx,
            'p_value': p_value,
            'significant_005': p_value < 0.05 if not np.isnan(p_value) else False,
        })

    df = pd.DataFrame(results)
    out_path = os.path.join(output_dir, 'ddg_ttest_results.csv')
    df.to_csv(out_path, index=False)
    print(f"  ΔΔG t-test results saved to {out_path}")

    for _, row in df.iterrows():
        sig = "SIGNIFICANT" if row['significant_005'] else "not significant"
        print(f"    {row['pair']}: ΔΔG={row['ddg']:.2f} kJ/mol, "
              f"95% CI=[{row['ci_lo']:.2f}, {row['ci_hi']:.2f}], "
              f"t={row['t_stat']:.2f}, p={row['p_value']:.3f} ({sig})")

    return df


# ---------------------------------------------------------------------------
# 7. Per-System DCCM Convergence Test (RQ-REP5)
# ---------------------------------------------------------------------------

def run_dccm_convergence(data_dir, output_dir, key_systems=None):
    """Per-system DCCM convergence: measure correlation of partial DCCMs with full-sample DCCM.

    Compute DCCM from the first 1000, 2000, 3000, 4000 frames of each system
    and measure Pearson correlation with the full-sample DCCM. Reports when
    convergence plateau is reached.
    """
    from scipy import stats
    from pca_utils import compute_dccm

    if key_systems is None:
        key_systems = ['phe_sssD_sap', 'phe_sssD_tsap', 'phe_rrrL_sap', 'phe_rrrL_tsap']

    print("  Loading aligned coordinates and system indices...")
    aligned = np.load(os.path.join(data_dir, 'aligned_coords.npy'), mmap_mode='r')
    system_indices = np.load(os.path.join(data_dir, 'system_indices.npy'))

    with open(os.path.join(data_dir, 'system_metadata.json')) as f:
        systems_meta = json.load(f)

    # Build system name -> index mapping
    sys_name_to_idx = {s['name']: i for i, s in enumerate(systems_meta)}

    # CA flat indices (every 4th atom starting from index 1, all 3 coords)
    N_BACKBONE = 2327
    N_CA = 582
    ca_atom_indices = np.arange(1, N_BACKBONE, 4)
    ca_flat_indices = np.sort(np.concatenate(
        [3 * ca_atom_indices + k for k in range(3)]
    ))

    results = []
    for sys_name in key_systems:
        if sys_name not in sys_name_to_idx:
            print(f"    System {sys_name} not found in metadata. Skipping.")
            continue

        idx = sys_name_to_idx[sys_name]
        start = int(system_indices[idx])
        end = int(system_indices[idx + 1])
        n_frames = end - start

        print(f"    Computing full DCCM for {sys_name} ({n_frames} frames)...")
        # Full DCCM
        ca_full = aligned[start:end, ca_flat_indices].astype(np.float64)
        dccm_full = compute_dccm(ca_full)

        # Upper triangle values (for correlation comparison)
        triu_idx = np.triu_indices(N_CA, k=1)
        full_upper = dccm_full[triu_idx]

        # Partial DCCMs at increasing frame counts
        for n_sub in [1000, 2000, 3000, 4000]:
            if n_sub >= n_frames:
                continue
            sub_start = start
            sub_end = start + n_sub
            ca_sub = aligned[sub_start:sub_end, ca_flat_indices].astype(np.float64)
            dccm_sub = compute_dccm(ca_sub)
            sub_upper = dccm_sub[triu_idx]

            # Pearson r between full and partial DCCM upper triangles
            r, p = stats.pearsonr(full_upper, sub_upper)
            rmse = np.sqrt(np.mean((full_upper - sub_upper)**2))

            results.append({
                'system': sys_name,
                'n_frames_full': n_frames,
                'n_frames_sub': n_sub,
                'fraction': n_sub / n_frames,
                'pearson_r': r,
                'rmse': rmse,
            })
            print(f"      {sys_name}: {n_sub}/{n_frames} frames, "
                  f"r={r:.4f}, RMSE={rmse:.4f}")

        # Free memory
        del ca_full, dccm_full

    df = pd.DataFrame(results)
    out_path = os.path.join(output_dir, 'dccm_convergence_analysis.csv')
    df.to_csv(out_path, index=False)
    print(f"  DCCM convergence saved to {out_path}")

    # Print convergence summary
    if len(df) > 0:
        print("\n  === DCCM Convergence Summary ===")
        for sys_name in key_systems:
            sys_df = df[df['system'] == sys_name]
            if len(sys_df) > 0:
                r_3k = sys_df[sys_df['n_frames_sub'] == 3000]['pearson_r']
                r_4k = sys_df[sys_df['n_frames_sub'] == 4000]['pearson_r']
                if len(r_3k) > 0 and len(r_4k) > 0:
                    conv_status = "converged (r>0.95 at 3000)" if r_3k.values[0] > 0.95 else "not yet converged"
                    print(f"    {sys_name}: {conv_status}")

    return df


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description='Statistical rigor analyses for PCA')
    parser.add_argument('--input-dir', default='pca_analysis/',
                        help='Input directory with PCA data files')
    parser.add_argument('--output-dir', default=None,
                        help='Output directory (default: same as input)')
    parser.add_argument('--convergence', action='store_true',
                        help='Run only PCA convergence analysis')
    parser.add_argument('--bootstrap', action='store_true',
                        help='Run only block-bootstrap eigenvalue CIs')
    parser.add_argument('--autocorr', action='store_true',
                        help='Run only autocorrelation / effective sample size')
    parser.add_argument('--anova', action='store_true',
                        help='Run only two-way ANOVA')
    parser.add_argument('--loo', action='store_true',
                        help='Run only leave-one-out correlation analysis')
    parser.add_argument('--ddg-ttest', action='store_true',
                        help='Run only ΔΔG t-tests')
    parser.add_argument('--dccm-convergence', action='store_true',
                        help='Run only DCCM convergence test')
    parser.add_argument('--all', action='store_true',
                        help='Run all analyses (1–7)')
    args = parser.parse_args()

    data_dir = args.input_dir
    output_dir = args.output_dir or data_dir

    os.makedirs(output_dir, exist_ok=True)

    # Determine which analyses to run
    any_flag = any([args.convergence, args.bootstrap, args.autocorr, args.anova,
                    args.loo, args.ddg_ttest, args.dccm_convergence])
    run_all = args.all or not any_flag

    print("=" * 60)
    print("P2 Statistical Rigor Analyses")
    print("=" * 60)

    t_start = time.time()

    output_files = []

    if run_all or args.convergence:
        print("\n--- 1. PCA Convergence Analysis ---")
        conv_df = run_convergence_analysis(data_dir, output_dir)
        output_files.append('convergence_analysis.csv')

    if run_all or args.bootstrap:
        print("\n--- 2. Block-Bootstrap Eigenvalue CIs ---")
        boot_df = run_bootstrap_eigenvalues(data_dir, output_dir)
        output_files.append('bootstrap_eigenvalue_ci.csv')

    if run_all or args.autocorr:
        print("\n--- 3. Autocorrelation / Effective Sample Size ---")
        acf_df = run_autocorrelation_analysis(data_dir, output_dir)
        output_files.append('autocorrelation_analysis.csv')

    if run_all or args.anova:
        print("\n--- 4. Two-Way ANOVA (Stereo × Chirality) ---")
        anova_df = run_two_way_anova(data_dir, output_dir)
        output_files.append('anova_interaction.csv')

    if run_all or args.loo:
        print("\n--- 5. Leave-One-Out Correlation Analysis ---")
        loo_df = run_leave_one_out_correlation(data_dir, output_dir)
        output_files.append('loo_correlation_analysis.csv')

    if run_all or args.ddg_ttest:
        print("\n--- 6. ΔΔG t-tests (SAP vs TSAP) ---")
        ddg_df = run_ddg_ttests(data_dir, output_dir)
        output_files.append('ddg_ttest_results.csv')

    if run_all or args.dccm_convergence:
        print("\n--- 7. DCCM Convergence Test ---")
        dccm_conv_df = run_dccm_convergence(data_dir, output_dir)
        output_files.append('dccm_convergence_analysis.csv')

    elapsed = time.time() - t_start
    print(f"\n{'=' * 60}")
    print(f"Analyses complete ({elapsed:.0f}s). Output files:")
    for f in output_files:
        fpath = os.path.join(output_dir, f)
        if os.path.exists(fpath):
            n_rows = len(pd.read_csv(fpath))
            print(f"  OK  {fpath} ({n_rows} rows)")
        else:
            print(f"  MISSING  {fpath}")


if __name__ == '__main__':
    main()
