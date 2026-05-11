#!/usr/bin/env python3
"""
Example integration of SHAPE 2.1 analysis with MD trajectory processing.

This script demonstrates how to:
1. Extract coordination shell coordinates from MD trajectory
2. Generate SHAPE input files
3. Run SHAPE analysis
4. Parse and visualize results

Usage:
    python shape_trajectory_analysis.py --tpr system.tpr --xtc system.xtc --metal-idx 54 --coord-idx 0 1 2 3 4 5 6 7 63
"""

import argparse
import subprocess
import os
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from pathlib import Path

try:
    import MDAnalysis as mda
except ImportError:
    print("Warning: MDAnalysis not installed. Some functions will not work.")
    mda = None


def extract_coordination_frames(tpr_file, xtc_file, metal_idx, coord_idx, 
                                frame_interval=10, max_frames=None):
    """
    Extract coordination shell coordinates from MD trajectory.
    
    Parameters:
    -----------
    tpr_file : str
        GROMACS topology file
    xtc_file : str
        GROMACS trajectory file
    metal_idx : int
        Index of metal atom (0-based)
    coord_idx : list of int
        Indices of coordinating atoms (0-based)
    frame_interval : int
        Extract every Nth frame (default: 10)
    max_frames : int or None
        Maximum number of frames to extract (None = all)
    
    Returns:
    --------
    list of tuples: [(frame_name, coordinates_array), ...]
        coordinates_array: (9, 3) numpy array centered on metal
    """
    if mda is None:
        raise ImportError("MDAnalysis required for trajectory extraction")
    
    u = mda.Universe(tpr_file, xtc_file)
    frames = []
    
    # Get trajectory slice
    traj = u.trajectory[::frame_interval]
    if max_frames is not None:
        traj = list(traj)[:max_frames]
    
    for i, ts in enumerate(traj):
        # Get positions
        metal_pos = u.atoms[metal_idx].position.copy()
        coord_pos = u.atoms[coord_idx].positions.copy()
        
        # Center on metal
        coord_pos_centered = coord_pos - metal_pos
        
        # Convert nm to Angstroms (SHAPE expects Å)
        coord_pos_centered *= 10.0
        
        frame_name = f"frame_{i:05d}"
        frames.append((frame_name, coord_pos_centered))
    
    print(f"Extracted {len(frames)} frames from trajectory")
    return frames


def generate_shape_input(frames, output_file, title="MD trajectory analysis",
                         ligand_label='O', metal_label='Eu'):
    """
    Generate SHAPE .dat input file from extracted frames.
    
    Parameters:
    -----------
    frames : list of tuples
        [(frame_name, coordinates_array), ...]
    output_file : str
        Output .dat filename
    title : str
        Title for SHAPE input
    ligand_label : str
        Element label for ligands
    metal_label : str
        Element label for metal
    """
    with open(output_file, 'w') as f:
        # Header
        f.write(f"$ {title}\n")
        f.write("! Ligands  Metal\n")
        f.write("     9      10\n")
        f.write("! EP-9 CSAPR-9\n")
        f.write("     1       8\n")
        
        # Write each frame
        for frame_name, coords in frames:
            f.write(f"{frame_name}\n")
            
            # Write ligand coordinates
            for pos in coords:
                f.write(f"  {ligand_label}   {pos[0]:10.4f}  {pos[1]:10.4f}  {pos[2]:10.4f}\n")
            
            # Write metal at origin
            f.write(f"  {metal_label}    0.0000    0.0000    0.0000\n")
    
    print(f"Generated SHAPE input: {output_file}")


def run_shape_analysis(input_file, shape_executable="SHAPE_2.1_linux_64/shape_2.1_linux64"):
    """
    Run SHAPE analysis on input file.
    
    Parameters:
    -----------
    input_file : str
        SHAPE .dat input file
    shape_executable : str
        Path to SHAPE executable
    
    Returns:
    --------
    str: Path to output .tab file
    """
    # Get directory and basename
    input_path = Path(input_file)
    work_dir = input_path.parent
    base_name = input_path.stem
    
    # Run SHAPE
    cmd = [shape_executable, input_file]
    
    print(f"Running SHAPE analysis...")
    print(f"  Command: {' '.join(cmd)}")
    print(f"  Working directory: {work_dir}")
    
    result = subprocess.run(cmd, cwd=work_dir, capture_output=True, text=True)
    
    if result.returncode != 0:
        print(f"SHAPE error: {result.stderr}")
        return None
    
    print(result.stdout)
    
    output_tab = work_dir / f"{base_name}.tab"
    if output_tab.exists():
        print(f"SHAPE output: {output_tab}")
        return str(output_tab)
    else:
        print(f"Error: Output file not found: {output_tab}")
        return None


def parse_shape_tab(filename):
    """
    Parse SHAPE .tab output file.
    
    Parameters:
    -----------
    filename : str
        SHAPE .tab file
    
    Returns:
    --------
    dict: {structure_name: {reference_name: chsm_value, ...}, ...}
    """
    results = {}
    ref_names = []
    
    with open(filename, 'r') as f:
        lines = f.readlines()
    
    # Parse reference geometry names
    for line in lines:
        line = line.strip()
        if line and not line.startswith('-') and not line.startswith('Structure'):
            parts = line.split()
            if len(parts) >= 4:
                # Check if it's a reference line (has a number as second field)
                try:
                    int(parts[1])
                    ref_names.append(parts[0])
                except ValueError:
                    pass
    
    # Parse data
    reading_data = False
    for line in lines:
        if 'Structure [ML' in line:
            reading_data = True
            continue
        
        if reading_data:
            line = line.strip()
            if not line or line.startswith('-'):
                continue
            
            parts = line.split(',')
            if len(parts) >= 2:
                name = parts[0].strip()
                values = []
                for v in parts[1:]:
                    try:
                        values.append(float(v.strip()))
                    except ValueError:
                        pass
                
                if len(values) == len(ref_names):
                    results[name] = dict(zip(ref_names, values))
    
    return results


def plot_chsm_timeseries(results, output_file='chsm_timeseries.png'):
    """
    Plot ChSM values over trajectory frames.
    
    Parameters:
    -----------
    results : dict
        Parsed results from parse_shape_tab()
    output_file : str
        Output figure filename
    """
    # Extract time series
    frames = sorted(results.keys())
    chsm_csapr = [results[f]['CSAPR-9'] for f in frames]
    frame_nums = [int(f.split('_')[1]) for f in frames]
    
    # Create plot
    fig, ax = plt.subplots(figsize=(10, 6))
    
    ax.plot(frame_nums, chsm_csapr, 'b-', linewidth=1.5, label='CSAPR-9 ChSM')
    ax.axhline(y=1.0, color='r', linestyle='--', label='Good match (ChSM=1)')
    ax.axhline(y=3.0, color='orange', linestyle='--', label='Moderate match (ChSM=3)')
    
    ax.set_xlabel('Frame Number', fontsize=12)
    ax.set_ylabel('CSAPR-9 ChSM', fontsize=12)
    ax.set_title('Geometry Evolution: Deviation from Ideal Capped Square Antiprism', fontsize=14)
    ax.legend(loc='best')
    ax.grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig(output_file, dpi=150, bbox_inches='tight')
    print(f"Saved plot: {output_file}")
    
    # Print statistics
    print(f"\nChSM Statistics:")
    print(f"  Mean:   {np.mean(chsm_csapr):.3f}")
    print(f"  Std:    {np.std(chsm_csapr):.3f}")
    print(f"  Min:    {np.min(chsm_csapr):.3f}")
    print(f"  Max:    {np.max(chsm_csapr):.3f}")
    
    return fig, ax


def main():
    """Main function demonstrating complete workflow."""
    parser = argparse.ArgumentParser(description='SHAPE analysis for MD trajectory')
    parser.add_argument('--tpr', help='GROMACS topology file')
    parser.add_argument('--xtc', help='GROMACS trajectory file')
    parser.add_argument('--metal-idx', type=int, help='Metal atom index (0-based)')
    parser.add_argument('--coord-idx', type=int, nargs='+', help='Coordinating atom indices')
    parser.add_argument('--interval', type=int, default=10, help='Frame interval')
    parser.add_argument('--max-frames', type=int, help='Maximum frames to analyze')
    parser.add_argument('--output', default='trajectory_shape.dat', help='Output file')
    parser.add_argument('--shape-exe', default='SHAPE_2.1_linux_64/shape_2.1_linux64',
                       help='SHAPE executable path')
    
    args = parser.parse_args()
    
    # Example with trajectory
    if args.tpr and args.xtc:
        print("=" * 60)
        print("SHAPE 2.1 Analysis for 9-Coordinate Eu Complex")
        print("=" * 60)
        
        # Step 1: Extract frames
        print("\nStep 1: Extracting frames from trajectory...")
        frames = extract_coordination_frames(
            args.tpr, args.xtc, args.metal_idx, args.coord_idx,
            frame_interval=args.interval, max_frames=args.max_frames
        )
        
        # Step 2: Generate SHAPE input
        print("\nStep 2: Generating SHAPE input file...")
        generate_shape_input(frames, args.output)
        
        # Step 3: Run SHAPE
        print("\nStep 3: Running SHAPE analysis...")
        tab_file = run_shape_analysis(args.output, args.shape_exe)
        
        if tab_file:
            # Step 4: Parse results
            print("\nStep 4: Parsing results...")
            results = parse_shape_tab(tab_file)
            print(f"Parsed {len(results)} structures")
            
            # Step 5: Visualize
            print("\nStep 5: Creating visualization...")
            plot_chsm_timeseries(results)
            
            # Save results to CSV
            csv_file = args.output.replace('.dat', '_results.csv')
            with open(csv_file, 'w') as f:
                f.write("Frame,EP-9_ChSM,CSAPR-9_ChSM\n")
                for name in sorted(results.keys()):
                    f.write(f"{name},{results[name]['EP-9']:.6f},{results[name]['CSAPR-9']:.6f}\n")
            print(f"Saved results: {csv_file}")
        
        print("\n" + "=" * 60)
        print("Analysis complete!")
        print("=" * 60)
    
    else:
        print("\nExample usage with test data:")
        print("  See generate_test_example() function below\n")
        generate_test_example()


def generate_test_example():
    """Generate a simple test example with synthetic data."""
    print("Creating test example with synthetic data...")
    
    # Generate synthetic frames
    frames = []
    for i in range(20):
        # Create slightly distorted CSAPR-9 geometry
        coords = np.array([
            [0.0, 0.0, 2.5],      # Cap
            [2.33, 0.0, 0.90],    # Bottom square
            [0.0, 2.33, 0.90],
            [-2.33, 0.0, 0.90],
            [0.0, -2.33, 0.90],
            [1.40, 1.40, -1.52],  # Top square
            [-1.40, 1.40, -1.52],
            [-1.40, -1.40, -1.52],
            [1.40, -1.40, -1.52],
        ])
        
        # Add small random distortion
        coords += np.random.randn(9, 3) * 0.05
        
        frames.append((f"frame_{i:03d}", coords))
    
    # Generate SHAPE input
    output_file = "test_synthetic.dat"
    generate_shape_input(frames, output_file, title="Synthetic test data")
    
    # Run SHAPE
    tab_file = run_shape_analysis(output_file)
    
    if tab_file:
        # Parse and plot
        results = parse_shape_tab(tab_file)
        plot_chsm_timeseries(results, 'test_synthetic_timeseries.png')


if __name__ == '__main__':
    main()
