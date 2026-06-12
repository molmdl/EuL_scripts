#!/usr/bin/env python3
"""
T2_qc_energy.py

Batch-parse xTB traj.xyz files, perform QC, extract energies,
convert units, and generate publication-quality plots.

Usage:
    python scripts/T2_qc_energy.py --data-dir data --output-dir analysis
"""

import argparse
import glob
import os
import re
import sys

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import gaussian_kde, pearsonr


# CODATA 2018 conversion factor
HA_TO_KJMOL = 2625.4996394799
TIMESTEP_PS = 0.001  # 1 fs = 0.001 ps per frame

# Color scheme
COLORS = {
    'sap': '#4682B4',
    'tsap': '#BA55D3',
    'me': '#2E8B57',
    'phe': '#CD853F',
    'rrr': '#4169E1',
    'sss': '#DC143C',
    'D': '#1E90FF',
    'L': '#FF6347',
}


def parse_system_name(name):
    """Parse system name into metadata components."""
    parts = name.split('_')
    if len(parts) != 3:
        raise ValueError(f"Unexpected system name format: {name}")
    species, stereo, topo = parts
    isomer = stereo[:3]
    handedness = stereo[3]
    return species, isomer, handedness, topo


def discover_systems(data_dir):
    """Discover all systems with traj.xyz files."""
    pattern = os.path.join(data_dir, '*/traj.xyz')
    paths = sorted(glob.glob(pattern))
    if len(paths) != 16:
        raise ValueError(f"Expected 16 traj.xyz files, found {len(paths)}")

    systems = []
    for path in paths:
        name = os.path.basename(os.path.dirname(path))
        species, isomer, handedness, topology = parse_system_name(name)
        systems.append({
            'name': name,
            'path': path,
            'species': species,
            'isomer': isomer,
            'handedness': handedness,
            'topology': topology,
            'expected_natoms': 135 if species == 'me' else 149,
        })
    return systems


def parse_trajectory(system):
    """Stream-parse a single trajectory XYZ file."""
    path = system['path']
    expected_natoms = system['expected_natoms']
    energy_re = re.compile(r'energy:\s*([+-]?\d+(?:\.\d+)?(?:[eE][+-]?\d+)?)')

    frames = []
    frame_idx = 0

    with open(path, 'r') as f:
        while True:
            line = f.readline()
            if not line:
                break
            natoms = int(line.strip())
            if natoms != expected_natoms:
                raise AssertionError(
                    f"{system['name']}: frame {frame_idx} has {natoms} atoms, expected {expected_natoms}"
                )

            comment = f.readline()
            if not comment:
                raise AssertionError(f"{system['name']}: truncated file at frame {frame_idx}")

            match = energy_re.search(comment)
            if not match:
                raise AssertionError(f"{system['name']}: no energy found at frame {frame_idx}")
            energy_hartree = float(match.group(1))

            if not np.isfinite(energy_hartree):
                raise AssertionError(f"{system['name']}: non-finite energy at frame {frame_idx}")

            # Skip coordinate lines
            for _ in range(natoms):
                coord_line = f.readline()
                if not coord_line:
                    raise AssertionError(f"{system['name']}: truncated coordinates at frame {frame_idx}")

            frames.append({
                'frame': frame_idx,
                'natoms': natoms,
                'energy_hartree': energy_hartree,
            })
            frame_idx += 1

    if len(frames) != 2000:
        raise AssertionError(f"{system['name']}: expected 2000 frames, got {len(frames)}")

    return frames


def build_dataframe(systems):
    """Parse all trajectories and build a combined DataFrame."""
    all_rows = []
    qc_lines = []

    print("Parsing trajectories...")
    for system in systems:
        frames = parse_trajectory(system)
        e_arr = np.array([f['energy_hartree'] for f in frames])

        e_min = e_arr.min()
        e_max = e_arr.max()
        e_first = e_arr[0]
        e_last = e_arr[-1]

        # Drift checks
        frame_indices = np.arange(len(frames))
        r, _ = pearsonr(frame_indices, e_arr)
        abs_r = abs(r)

        max_drift_kj = (e_max - e_min) * HA_TO_KJMOL

        drift_warn = False
        warn_reasons = []
        if abs_r > 0.9:
            drift_warn = True
            warn_reasons.append(f"monotonic drift (|r|={abs_r:.3f})")
        if max_drift_kj > 20.0:
            drift_warn = True
            warn_reasons.append(f"max drift >20 kJ/mol ({max_drift_kj:.2f})")

        if drift_warn:
            print(f"  WARNING: {system['name']} - {', '.join(warn_reasons)}")

        qc_lines.append(
            f"{system['name']:20s}  natoms={system['expected_natoms']:3d}  n_frames={len(frames):4d}  "
            f"E_first={e_first:16.6f}  E_last={e_last:16.6f}  "
            f"E_min={e_min:16.6f}  E_max={e_max:16.6f}  "
            f"max_drift_kjmol={max_drift_kj:8.3f}  drift_warn={drift_warn}"
        )

        for f in frames:
            all_rows.append({
                'system': system['name'],
                'species': system['species'],
                'isomer': system['isomer'],
                'handedness': system['handedness'],
                'topology': system['topology'],
                'frame': f['frame'],
                'time_ps': f['frame'] * TIMESTEP_PS,
                'energy_hartree': f['energy_hartree'],
            })

    df = pd.DataFrame(all_rows)

    # Compute relative energies per system
    df['energy_kjmol'] = df['energy_hartree'] * HA_TO_KJMOL
    df['relative_kjmol'] = df.groupby('system')['energy_kjmol'].transform(lambda x: x - x.min())

    return df, qc_lines


def write_qc_report(qc_lines, output_dir):
    """Append QC report."""
    path = os.path.join(output_dir, 'qc_report.txt')
    with open(path, 'w') as f:
        f.write("# T2 QC Report: xTB Energy Parsing\n")
        f.write("=" * 140 + "\n")
        for line in qc_lines:
            f.write(line + "\n")
        f.write("=" * 140 + "\n\n")
    print(f"QC report appended to {path}")


def write_csv(df, output_dir):
    """Write energies CSV."""
    path = os.path.join(output_dir, 'energies_all.csv')
    df.to_csv(path, index=False)
    print(f"CSV written to {path} ({len(df)} rows)")


def plot_timeseries(df, output_dir):
    """Plot 4x4 grid of energy time series."""
    systems = df['system'].unique()

    fig, axes = plt.subplots(4, 4, figsize=(16, 12), sharex=False, sharey=False)
    axes = axes.flatten()

    for idx, sys_name in enumerate(systems):
        ax = axes[idx]
        sys_data = df[df['system'] == sys_name].sort_values('frame')
        topo = sys_data['topology'].iloc[0]
        color = COLORS.get(topo, 'gray')

        ax.plot(sys_data['time_ps'], sys_data['relative_kjmol'], lw=0.5, color=color)
        ax.axhline(0, color='black', linestyle='--', lw=0.5, alpha=0.5)
        ax.set_title(sys_name, fontsize=8)
        ax.set_xlabel('Time (ps)', fontsize=7)
        ax.set_ylabel('E_rel (kJ/mol)', fontsize=7)
        ax.tick_params(labelsize=6)

    fig.suptitle('xTB MD Energy Drift (relative to minimum, kJ/mol)', fontsize=14, y=1.00)
    plt.tight_layout()

    path = os.path.join(output_dir, 'plot_energy_timeseries.png')
    fig.savefig(path, dpi=150, bbox_inches='tight')
    plt.close(fig)
    print(f"Timeseries plot saved to {path}")


def plot_distributions(df, output_dir):
    """Plot grouped KDE/histogram distributions."""
    fig, axes = plt.subplots(2, 2, figsize=(12, 10))

    groupings = [
        ('species', {'me': COLORS['me'], 'phe': COLORS['phe']}),
        ('isomer', {'rrr': COLORS['rrr'], 'sss': COLORS['sss']}),
        ('handedness', {'D': COLORS['D'], 'L': COLORS['L']}),
        ('topology', {'sap': COLORS['sap'], 'tsap': COLORS['tsap']}),
    ]

    for ax, (group, colors) in zip(axes.flatten(), groupings):
        for val, color in colors.items():
            sub = df[df[group] == val]['relative_kjmol']
            if len(sub) < 2:
                continue
            kde = gaussian_kde(sub.values)
            x_range = np.linspace(sub.min(), sub.max(), 500)
            ax.plot(x_range, kde(x_range), color=color, label=f"{val} (n={len(sub)})", lw=2)

        ax.set_xlabel('Relative Energy (kJ/mol)', fontsize=10)
        ax.set_ylabel('Probability Density', fontsize=10)
        ax.set_title(f'By {group.capitalize()}', fontsize=12)
        ax.legend(fontsize=9)
        ax.grid(True, alpha=0.3)

    fig.suptitle('xTB Relative Energy Distributions', fontsize=14, y=1.00)
    plt.tight_layout()

    path = os.path.join(output_dir, 'plot_energy_distributions.png')
    fig.savefig(path, dpi=150, bbox_inches='tight')
    plt.close(fig)
    print(f"Distributions plot saved to {path}")


def plot_drift_summary(df, output_dir):
    """Plot bar chart of max drift per system."""
    drift = df.groupby('system')['relative_kjmol'].agg(['min', 'max'])
    drift['drift'] = drift['max'] - drift['min']

    # Order: species -> topology -> ...
    def sort_key(name):
        parts = name.split('_')
        s = 0 if parts[0] == 'me' else 1
        t = 0 if parts[2] == 'sap' else 1
        return (s, t, name)

    drift = drift.reindex(sorted(drift.index, key=sort_key))

    fig, ax = plt.subplots(figsize=(10, 6))
    colors = [
        COLORS['sap'] if name.split('_')[-1] == 'sap' else COLORS['tsap']
        for name in drift.index
    ]

    ax.barh(range(len(drift)), drift['drift'].values, color=colors)
    ax.set_yticks(range(len(drift)))
    ax.set_yticklabels(drift.index, fontsize=8)
    ax.set_xlabel('Max Energy Drift (kJ/mol)', fontsize=10)
    ax.set_title('Maximum Energy Drift per System', fontsize=12)
    ax.axvline(20, color='red', linestyle='--', lw=1, label='20 kJ/mol threshold')
    ax.legend()
    ax.grid(True, alpha=0.3, axis='x')

    plt.tight_layout()

    path = os.path.join(output_dir, 'plot_energy_drift_summary.png')
    fig.savefig(path, dpi=150, bbox_inches='tight')
    plt.close(fig)
    print(f"Drift summary plot saved to {path}")


def main():
    parser = argparse.ArgumentParser(description='T2: xTB Energy QC and Plotting')
    parser.add_argument('--data-dir', default='data', help='Path to data directory (default: data)')
    parser.add_argument('--output-dir', default='analysis', help='Path to output directory (default: analysis)')
    args = parser.parse_args()

    data_dir = args.data_dir
    output_dir = args.output_dir

    if not os.path.isdir(data_dir):
        print(f"Error: data directory not found: {data_dir}")
        sys.exit(1)

    os.makedirs(output_dir, exist_ok=True)

    # Set plot style
    try:
        plt.style.use('seaborn-v0_8-whitegrid')
    except Exception:
        pass

    systems = discover_systems(data_dir)
    print(f"Discovered {len(systems)} systems")

    df, qc_lines = build_dataframe(systems)
    write_qc_report(qc_lines, output_dir)
    write_csv(df, output_dir)
    plot_timeseries(df, output_dir)
    plot_distributions(df, output_dir)
    plot_drift_summary(df, output_dir)

    print("\nSummary:")
    print(f"  Total systems: {len(systems)}")
    print(f"  Total frames: {len(df)}")
    print(f"  Energy range (Hartree): {df['energy_hartree'].min():.6f} to {df['energy_hartree'].max():.6f}")
    print(f"  Max relative energy (kJ/mol): {df['relative_kjmol'].max():.3f}")
    print("Done.")


if __name__ == '__main__':
    main()
