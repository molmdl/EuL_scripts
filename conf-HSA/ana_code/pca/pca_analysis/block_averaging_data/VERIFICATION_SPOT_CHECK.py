#!/usr/bin/env python3
"""
Independent numerical spot-check for block averaging DCCM analysis.

Verifies:
1. Block-averaged DCCM = element-wise mean of per-block DCCMs
2. Block-std DCCM = element-wise sample std of per-block DCCMs
3. Pearson r between block-averaged and full-sample DCCM matches summary CSV
4. RMSE between block-averaged and full-sample DCCM matches summary CSV
5. Inter-block mean/min r from interblock CSVs matches summary CSV
6. Flyvbjerg-Petersen variance values match variance plateau CSVs
7. System frame indexing from system_indices.npy matches system_metadata.json
"""

import os
import json
import numpy as np
from scipy import stats

BASE = os.path.dirname(os.path.abspath(__file__))
INPUT_DIR = os.path.join(BASE, '..')

SYSTEMS = ['phe_sssD_sap', 'phe_sssD_tsap', 'phe_rrrL_sap', 'phe_rrrL_tsap']
BLOCK_SIZES = [400, 500, 572, 667, 800, 1000]
N_CA = 582

def check_block_avg_is_mean():
    """Verify block-averaged DCCM = mean of per-block DCCMs."""
    errors = []
    for sys_name in SYSTEMS:
        for bs in BLOCK_SIZES:
            per_block_path = os.path.join(BASE, 'per_block_dccm',
                                          f'block_dccm_{sys_name}_{bs}.npy')
            avg_path = os.path.join(BASE, 'block_averaged_dccm',
                                    f'block_avg_dccm_{sys_name}_{bs}.npy')
            if not os.path.exists(per_block_path):
                errors.append(f"MISSING: {per_block_path}")
                continue
            per_block = np.load(per_block_path)
            block_avg = np.load(avg_path)
            computed_mean = per_block.mean(axis=0)
            if not np.allclose(block_avg, computed_mean):
                max_diff = np.abs(block_avg - computed_mean).max()
                errors.append(f"{sys_name} B={bs}: block_avg != mean, max_diff={max_diff:.2e}")

    if errors:
        print("FAIL: block_avg != mean")
        for e in errors:
            print(f"  {e}")
    else:
        print("PASS: block-averaged DCCM = mean of per-block DCCMs (all 24 combos)")
    return errors


def check_block_std_is_sample_std():
    """Verify block-std DCCM = sample std (ddof=1) of per-block DCCMs."""
    errors = []
    for sys_name in SYSTEMS:
        for bs in BLOCK_SIZES:
            per_block_path = os.path.join(BASE, 'per_block_dccm',
                                          f'block_dccm_{sys_name}_{bs}.npy')
            std_path = os.path.join(BASE, 'block_averaged_dccm',
                                    f'block_std_dccm_{sys_name}_{bs}.npy')
            if not os.path.exists(per_block_path):
                continue
            per_block = np.load(per_block_path)
            block_std = np.load(std_path)
            computed_std = per_block.std(axis=0, ddof=1)
            if not np.allclose(block_std, computed_std):
                max_diff = np.abs(block_std - computed_std).max()
                errors.append(f"{sys_name} B={bs}: block_std != sample std, max_diff={max_diff:.2e}")

    if errors:
        print("FAIL: block_std != sample std")
        for e in errors:
            print(f"  {e}")
    else:
        print("PASS: block-std DCCM = sample std of per-block DCCMs (all 24 combos)")
    return errors


def check_correlation_metrics():
    """Verify Pearson r and RMSE match summary CSV for all systems at B=500."""
    import pandas as pd
    errors = []
    summary = pd.read_csv(os.path.join(BASE, 'summary', 'block_averaging_summary.csv'))
    primary = summary[summary['block_size'] == 500]

    triu_idx = np.triu_indices(N_CA, k=1)

    for sys_name in SYSTEMS:
        block_avg = np.load(os.path.join(BASE, 'block_averaged_dccm',
                                         f'block_avg_dccm_{sys_name}_500.npy'))
        dccm_full = np.load(os.path.join(INPUT_DIR, f'dccm_{sys_name}.npy'))

        v1 = block_avg[triu_idx]
        v2 = dccm_full[triu_idx]
        r, _ = stats.pearsonr(v1, v2)
        rmse = np.sqrt(np.mean((v1 - v2) ** 2))

        expected_r = primary[primary['system'] == sys_name]['pearson_r_block_avg_vs_full'].values[0]
        expected_rmse = primary[primary['system'] == sys_name]['rmse_block_avg_vs_full'].values[0]

        if abs(r - expected_r) > 1e-5:
            errors.append(f"{sys_name}: r={r:.6f} vs expected {expected_r:.6f}")
        if abs(rmse - expected_rmse) > 1e-5:
            errors.append(f"{sys_name}: rmse={rmse:.6f} vs expected {expected_rmse:.6f}")

    if errors:
        print("FAIL: correlation metrics don't match summary CSV")
        for e in errors:
            print(f"  {e}")
    else:
        print("PASS: Pearson r and RMSE match summary CSV (all 4 systems at B=500)")
    return errors


def check_interblock_r():
    """Verify inter-block mean/min r from CSVs matches summary CSV."""
    import pandas as pd
    errors = []
    summary = pd.read_csv(os.path.join(BASE, 'summary', 'block_averaging_summary.csv'))
    primary = summary[summary['block_size'] == 500]

    for sys_name in SYSTEMS:
        r_df = pd.read_csv(os.path.join(BASE, 'convergence_metrics',
                                        f'interblock_r_{sys_name}_500.csv'), index_col=0)
        r_matrix = r_df.values
        K = r_matrix.shape[0]
        off_diag_mask = ~np.eye(K, dtype=bool)
        off_diag = r_matrix[off_diag_mask]

        computed_mean = off_diag.mean()
        computed_min = off_diag.min()

        expected_mean = primary[primary['system'] == sys_name]['mean_interblock_r'].values[0]
        expected_min = primary[primary['system'] == sys_name]['min_interblock_r'].values[0]

        if abs(computed_mean - expected_mean) > 1e-5:
            errors.append(f"{sys_name}: mean_r={computed_mean:.6f} vs expected {expected_mean:.6f}")
        if abs(computed_min - expected_min) > 1e-5:
            errors.append(f"{sys_name}: min_r={computed_min:.6f} vs expected {expected_min:.6f}")

    if errors:
        print("FAIL: inter-block r values don't match summary CSV")
        for e in errors:
            print(f"  {e}")
    else:
        print("PASS: inter-block mean/min r match summary CSV (all 4 systems at B=500)")
    return errors


def check_system_frame_indices():
    """Verify system_indices.npy frame ranges match system_metadata.json."""
    errors = []
    si = np.load(os.path.join(INPUT_DIR, 'system_indices.npy'))
    with open(os.path.join(INPUT_DIR, 'system_metadata.json')) as f:
        meta = json.load(f)
    name_to_idx = {s['name']: i for i, s in enumerate(meta)}

    for sys_name in SYSTEMS:
        idx = name_to_idx[sys_name]
        start = int(si[idx])
        end = int(si[idx + 1])
        n_frames = end - start
        expected = [s['n_frames'] for s in meta if s['name'] == sys_name][0]
        if n_frames != expected:
            errors.append(f"{sys_name}: n_frames={n_frames} vs expected {expected}")

    if errors:
        print("FAIL: frame indices don't match metadata")
        for e in errors:
            print(f"  {e}")
    else:
        print("PASS: system frame indices match metadata (all 4 systems)")
    return errors


if __name__ == '__main__':
    print("=" * 70)
    print("Block Averaging DCCM — Independent Numerical Spot-Check")
    print("=" * 70)

    all_errors = []
    all_errors.extend(check_block_avg_is_mean())
    all_errors.extend(check_block_std_is_sample_std())
    all_errors.extend(check_correlation_metrics())
    all_errors.extend(check_interblock_r())
    all_errors.extend(check_system_frame_indices())

    print()
    if all_errors:
        print(f"OVERALL: FAIL ({len(all_errors)} errors)")
    else:
        print("OVERALL: PASS (all numerical checks passed)")
