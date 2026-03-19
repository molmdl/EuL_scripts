#!/usr/bin/env python3
"""
FES Analysis Script

This script analyzes Free Energy Surface (FES) files generated from PLUMED or similar tools.
It provides functionalities to calculate absolute binding free energy, assess convergence,
and perform block analysis with bootstrapping, mimicking the behavior of the TCL script
from https://raw.githubusercontent.com/limresgrp/FMAP_v1/refs/heads/master/funnel_gui/ffs.tcl.

Key Features:
- Supports 1D and 2D FES files.
- Handles 'Infinity' values by treating them as NaN and skipping in computations.
- Outputs energies in kJ/mol and kcal/mol to 2 decimal places.
- Logs all steps and catches errors for robustness.
- Uses space-separated output files.

Dependencies (from env.yml):
- numpy=1.26.4 for numerical computations.
- matplotlib=3.9.2 for plotting.
- No MDanalysis or scipy needed, as operations are on text files.

Usage:
Run with subcommands: calculate, convergence, block.
Example: python fes_analysis.py calculate --fes fes.dat --xmin 0.01 --xmax 1.5 --wref 6 --rcyl 0.1

For full help: python fes_analysis.py --help
"""

import argparse
import glob
import logging
import os
import re

import matplotlib.pyplot as plt
import numpy as np

# Constants (defined at the beginning)
KT = 2.494339  # Precise kT in kJ/mol at 300 K (R * 300 / 1000)
PI = np.pi
C0 = 0.6020  # Standard concentration factor for 1 M
KCAL_CONV = 4.187  # Conversion factor for kcal/mol (kJ/mol / 4.187)
N_BOOTSTRAP = 1000  # Number of bootstrap iterations
N_BLOCKS = 10  # Number of blocks for analysis
MIN_POINTS_PER_BLOCK = 50  # Minimum points per block to proceed


def calculate_log_k(dg_kj):
    """
    Calculate log10 of binding constant K from binding free energy.

    Parameters:
    dg_kj (float): Binding free energy in kJ/mol.

    Returns:
    float: log10(K) where K is the binding (association) constant in M^-1.

    Notes:
    K = exp(-ΔG / RT) where ΔG is the standard binding free energy.
    log10(K) = -ΔG / (RT * ln(10))
    For a negative ΔG (favorable binding), log10(K) > 0.
    """
    try:
        log_k = -dg_kj / (KT * np.log(10))
        return log_k
    except Exception as e:
        logging.error(f"Error calculating log K: {e}")
        raise

np.seterr(divide='ignore', invalid='ignore')

# Setup logging to stdout with timestamp and level
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

def load_fes(fes_path):
    """
    Load and parse an FES file, mimicking TCL conventions.

    Parameters:
    fes_path (str): Path to the FES .dat file.

    Returns:
    dict: Containing 'dim' (1 or 2), 'dr' (list of bin sizes), 'x' (CV1 array),
          'y' (CV2 array or None), 'fes' (FES values array).

    Raises:
    ValueError: If file is empty, unsupported dim, or insufficient columns.
    """
    try:
        logging.info(f"Loading FES file: {fes_path}")
        mins, maxs, nbins = [], [], []
        data_list = []

        with open(fes_path, 'r') as f:
            for line in f:
                line = line.strip()
                if not line or line.startswith('#!'):
                    # Parse headers
                    if re.search(r'min_\S+', line):
                        val_str = line.split()[3]
                        if val_str == 'pi':
                            val = PI
                        elif val_str == '-pi':
                            val = -PI
                        else:
                            val = float(val_str)
                        mins.append(val)
                    elif re.search(r'max_\S+', line):
                        val_str = line.split()[3]
                        if val_str == 'pi':
                            val = PI
                        elif val_str == '-pi':
                            val = -PI
                        else:
                            val = float(val_str)
                        maxs.append(val)
                    elif re.search(r'nbins_\S+', line):
                        nbins.append(int(line.split()[3]))
                    continue

                # Parse data lines
                fields = line.split()
                if not fields:
                    continue
                row = []
                for v in fields:
                    if v == 'Infinity':
                        row.append(np.nan)
                    else:
                        try:
                            row.append(float(v))
                        except ValueError:
                            row.append(np.nan)
                if len(row) >= 3:  # At least CV(s) + FES
                    data_list.append(row)

        if not data_list:
            raise ValueError("No valid data lines found in FES file")

        data = np.array(data_list)
        dim = len(nbins)
        if dim not in (1, 2):
            raise ValueError(f"Unsupported dimension: {dim} (detected {len(nbins)} nbins sets)")

        dr = [(maxs[i] - mins[i]) / nbins[i] for i in range(dim)]

        cv_x_col = 0
        fes_col = 1 if dim == 1 else 2
        cv_y_col = 1 if dim == 2 else None

        if data.shape[1] <= fes_col:
            raise ValueError(f"Insufficient columns in data (need >{fes_col}, got {data.shape[1]})")

        x = data[:, cv_x_col]
        fes = data[:, fes_col]
        y = data[:, cv_y_col] if dim == 2 else None

        logging.info(f"Detected {dim}D FES, dr={' * '.join(f'{d:.6f}' for d in dr)}, valid points: {np.sum(np.isfinite(fes))}")

        return {
            'dim': dim,
            'dr': dr,
            'x': x,
            'y': y,
            'fes': fes
        }
    except FileNotFoundError:
        logging.error(f"File not found: {fes_path}")
        raise
    except Exception as e:
        logging.error(f"Error loading FES: {e}")
        raise

def compute_dg(fes_info, x_min, x_max, y_min=None, y_max=None, wref=None, rcyl=0.1):
    """
    Compute absolute binding free energy from FES data.

    Parameters:
    fes_info (dict): From load_fes().
    x_min, x_max (float): Interval for CV1 (required).
    y_min, y_max (float): Interval for CV2 (optional for 2D).
    wref (float): Reference point for CV1.
    rcyl (float): Cylinder radius in nm (default 0.1).

    Returns:
    tuple: (dg_kj, dg_kcal) in kJ/mol and kcal/mol.

    Raises:
    ValueError: If invalid regions or computations.
    """
    try:
        dim = fes_info['dim']
        dr = fes_info['dr']
        x = fes_info['x']
        y = fes_info['y']
        fes = fes_info['fes']

        if dim == 1 and (y_min is not None or y_max is not None):
            logging.info("1D FES detected; ignoring y intervals.")

        # Compute reference
        ref_mask = (x >= wref) & (x <= wref + dr[0])
        ref_values = fes[ref_mask & np.isfinite(fes)]
        if len(ref_values) == 0:
            raise ValueError("No valid points in reference region")
        ref = np.mean(ref_values)
        logging.info(f"Reference FES: {ref:.2f} (from {len(ref_values)} points)")

        # Integral mask
        mask = (x >= x_min) & (x <= x_max) & np.isfinite(fes)
        if dim == 2:
            if y_min is not None and y_max is not None:
                mask &= (y >= y_min) & (y <= y_max)
            else:
                logging.info("2D FES: using full y range")

        valid_count = np.sum(mask)
        if valid_count == 0:
            raise ValueError("No valid points in integration region")
        logging.info(f"Integration over {valid_count} valid bins")

        boltzmann = np.exp(-(fes[mask] - ref) / KT)
        bin_volume = dr[0] if dim == 1 else dr[0] * dr[1]
        kb = np.sum(boltzmann) * bin_volume

        if kb <= 0:
            raise ValueError("Computed kb <= 0; check bounds or valid data points")

        dg_kj = -KT * np.log(kb * PI * rcyl**2 * C0)
        dg_kcal = dg_kj / KCAL_CONV

        log_k = calculate_log_k(dg_kj)
        return dg_kj, dg_kcal, log_k
    except Exception as e:
        logging.error(f"Error computing DeltaG: {e}")
        raise

def calculate_mode(args):
    """Calculate mode: Compute DeltaG for a single FES file."""
    try:
        fes_info = load_fes(args.fes)
        dg_kj, dg_kcal, log_k = compute_dg(
            fes_info, args.xmin, args.xmax,
            args.ymin, args.ymax, args.wref, args.rcyl
        )
        print(f"DeltaG: {dg_kj:.2f} kJ/mol ({dg_kcal:.2f} kcal/mol)")
        print(f"log K: {log_k:.2f}")
        out_path = os.path.join(args.output_dir, 'deltaG_single.txt')
        with open(out_path, 'w') as f:
            f.write(f"{dg_kj:.2f} {dg_kcal:.2f} {log_k:.2f}\n")
        logging.info(f"Results saved to {out_path}")
    except Exception as e:
        logging.error(f"Calculate mode failed: {e}")

def convergence_mode(args):
    """Convergence mode: Analyze series of FES files for convergence."""
    try:
        fes_files = sorted(glob.glob(os.path.join(args.output_dir, 'fes_*.dat')))
        if not fes_files:
            raise FileNotFoundError("No fes_*.dat files found in output_dir")

        indices = []
        dgs_kj = []
        for fes in fes_files:
            match = re.search(r'fes_(\d+)\.dat', os.path.basename(fes))
            if match:
                idx = int(match.group(1))
            else:
                logging.warning(f"Could not extract index from {fes}, skipping")
                continue
            indices.append(idx)
            fes_info = load_fes(fes)
            dg_kj, _, _ = compute_dg(fes_info, args.xmin, args.xmax, args.ymin, args.ymax, args.wref, args.rcyl)
            dgs_kj.append(dg_kj)

        # Sort by index if needed (glob sorted, but ensure)
        sorted_pairs = sorted(zip(indices, dgs_kj))
        indices, dgs_kj = zip(*sorted_pairs)

        out_path = os.path.join(args.output_dir, 'deltaG.txt')
        with open(out_path, 'w') as f:
            f.write("#! Num. fes   DeltaG   Mean_otf   Std_error\n")
            tot1 = tot2 = tot3 = 0.0
            mean = std = 0.0
            for i in range(len(indices)):
                f.write(f"{indices[i]}   {dgs_kj[i]:.2f}")
                if indices[i] > args.reject:
                    weight = indices[i] - args.reject
                    tot1 += weight * dgs_kj[i]
                    tot2 += weight
                    mean = tot1 / tot2 if tot2 > 0 else 0.0
                    tot3 += weight * (dgs_kj[i] - mean)**2
                    std = np.sqrt(tot3 / tot2) if tot2 > 0 else 0.0
                    f.write(f"   {mean:.2f}   {std:.2f}\n")
                else:
                    f.write("\n")

        # Plot convergence
        plt.figure()
        plt.plot(indices, dgs_kj, 'o-', label='DeltaG')
        plt.xlabel('Num. fes')
        plt.ylabel('DeltaG (kJ/mol)')
        plt.legend()
        plot_path = os.path.join(args.output_dir, 'deltaG.png')
        plt.savefig(plot_path)
        plt.close()
        logging.info(f"Convergence plot saved to {plot_path}")
        logging.info("Convergence analysis completed")
    except Exception as e:
        logging.error(f"Convergence mode failed: {e}")

def block_mode(args):
    """Block mode: Perform block analysis with bootstrapping on deltaG.txt."""
    try:
        deltaG_path = os.path.join(args.output_dir, 'deltaG.txt')
        indices = []
        dgs = []
        with open(deltaG_path, 'r') as f:
            for line in f:
                if line.startswith('#!'):
                    continue
                parts = line.split()
                if len(parts) >= 2:
                    indices.append(int(parts[0]))
                    dgs.append(float(parts[1]))

        reject_idx = next((i for i, idx in enumerate(indices) if idx > args.reject), len(indices))
        post_dgs = np.array(dgs[reject_idx:])
        post_indices = np.array(indices[reject_idx:])
        n_points = len(post_dgs)

        if n_points < N_BLOCKS * MIN_POINTS_PER_BLOCK:
            raise ValueError(f"Too few points for block bootstrap (need >= {N_BLOCKS * MIN_POINTS_PER_BLOCK}, got {n_points})")

        number_points = n_points // N_BLOCKS

        boot_path = os.path.join(args.output_dir, 'bootstrap.txt')
        with open(boot_path, 'w') as f_boot:
            for block in range(N_BLOCKS):
                averages = []
                for _ in range(N_BOOTSTRAP):
                    offset = n_points * block / N_BLOCKS - 1
                    rand_offsets = np.random.rand(number_points) * number_points
                    rand_idxs = np.round(rand_offsets + offset + args.reject).astype(int)
                    rand_idxs = np.clip(rand_idxs, min(indices), max(indices))
                    # Map to array indices
                    valid_mask = np.isin(rand_idxs, indices)
                    rand_idxs = rand_idxs[valid_mask]
                    if len(rand_idxs) == 0:
                        continue  # Skip empty
                    sampled_dgs = [dgs[indices.index(r)] for r in rand_idxs]
                    weights = rand_idxs - args.reject + 1  # Mimic TCL +1
                    weighted_av = np.average(sampled_dgs, weights=weights)
                    averages.append(weighted_av)
                f_boot.write(' '.join(f"{a:.2f}" for a in averages) + '\n')

        # Histograms
        avs_list = []
        with open(boot_path, 'r') as f:
            for line in f:
                avs_list.append([float(v) for v in line.split() if v])

        fig, axs = plt.subplots(2, 5, figsize=(15, 6))
        axs = axs.flatten()
        for b in range(N_BLOCKS):
            if avs_list[b]:
                axs[b].hist(avs_list[b], bins=21)
                axs[b].set_title(f"Block {b+1} Histogram (kJ/mol)")
        plt.tight_layout()
        histo_path = os.path.join(args.output_dir, 'bootstrap_histo.png')
        plt.savefig(histo_path)
        plt.close()
        logging.info(f"Bootstrap histograms saved to {histo_path}")
        logging.info("Block analysis completed")
    except Exception as e:
        logging.error(f"Block mode failed: {e}")

def main():
    """Main entry point with argument parsing."""
    parser = argparse.ArgumentParser(
        description="FES Analysis Tool: Calculate binding free energy, convergence, and block analysis from PLUMED FES files.",
        epilog="Example: python fes_analysis.py calculate --fes fes.dat --xmin 0.01 --xmax 1.5 --ymin 0 --ymax 0.3 --wref 6 --rcyl 0.1"
    )
    parser.add_argument(
        '--output_dir', default='.', help="Output directory for results (default: current directory)."
    )
    parser.add_argument(
        '--reject', type=int, default=0, help="Number of initial points to reject in convergence/block (default: 0)."
    )

    subparsers = parser.add_subparsers(dest='mode', required=True, help="Available modes")

    # Calculate mode
    calc = subparsers.add_parser('calculate', help="Compute DeltaG for a single FES file.")
    calc.add_argument('--fes', required=True, help="Path to the FES .dat file.")
    calc.add_argument('--xmin', type=float, required=True, help="Minimum value for CV1 interval.")
    calc.add_argument('--xmax', type=float, required=True, help="Maximum value for CV1 interval.")
    calc.add_argument('--ymin', type=float, help="Minimum value for CV2 interval (2D only).")
    calc.add_argument('--ymax', type=float, help="Maximum value for CV2 interval (2D only).")
    calc.add_argument('--wref', type=float, required=True, help="Reference point for CV1.")
    calc.add_argument('--rcyl', type=float, default=0.1, help="Cylinder radius in nm (default: 0.1).")

    # Convergence mode
    conv = subparsers.add_parser('convergence', help="Analyze convergence over multiple FES files.")
    conv.add_argument('--xmin', type=float, required=True, help="Minimum value for CV1 interval.")
    conv.add_argument('--xmax', type=float, required=True, help="Maximum value for CV1 interval.")
    conv.add_argument('--ymin', type=float, help="Minimum value for CV2 interval (2D only).")
    conv.add_argument('--ymax', type=float, help="Maximum value for CV2 interval (2D only).")
    conv.add_argument('--wref', type=float, required=True, help="Reference point for CV1.")
    conv.add_argument('--rcyl', type=float, default=0.1, help="Cylinder radius in nm (default: 0.1).")

    # Block mode
    subparsers.add_parser('block', help="Perform block bootstrap analysis on deltaG.txt.")

    args = parser.parse_args()

    try:
        if args.mode == 'calculate':
            calculate_mode(args)
        elif args.mode == 'convergence':
            convergence_mode(args)
        elif args.mode == 'block':
            block_mode(args)
    except Exception as e:
        logging.error(f"Main execution failed: {e}")
        parser.print_help()

if __name__ == '__main__':
    main()
