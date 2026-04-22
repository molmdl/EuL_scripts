# Implementation Plan: Protein-Ligand Trajectory Analysis Script

**Date:** March 22, 2026  
**Domain:** GROMACS molecular dynamics trajectory analysis  
**Confidence:** HIGH

---

## Overview

This plan details the implementation of a Python script for comprehensive protein-ligand trajectory analysis. The script will perform five main analyses:

- **Aim 0:** Load topology/trajectory and align to receptor backbone
- **Aim 1:** Calculate RMSD for receptor backbone and ligand heavy atoms
- **Aim 2:** Calculate per-residue contact numbers between receptor and ligand
- **Aim 3:** Calculate per-residue hydrogen bonds between receptor and ligand
- **Aim 4:** Plot MMPBSA results (energy and decomposition)

The script will be command-line driven, generalized for any protein-ligand system, and use efficient vectorized approaches.

---

## 1. Module Structure

### 1.1 File Organization

```
com_ana_trj/
├── analyze_trajectory.py    # Main entry point
├── modules/
│   ├── __init__.py
│   ├── load_align.py        # Aim 0: Load and align trajectory
│   ├── rmsd.py              # Aim 1: RMSD calculations
│   ├── contacts.py          # Aim 2: Contact analysis
│   ├── hydrogen_bonds.py    # Aim 3: Hydrogen bond analysis
│   └── mmpbsa_plots.py      # Aim 4: MMPBSA plotting
├── utils/
│   ├── __init__.py
│   ├── cli.py               # Command-line argument parsing
│   ├── io.py                # File I/O utilities
│   └── plotting.py          # Shared plotting utilities
└── output/                  # Generated output files
```

### 1.2 Main Script Entry Point

```python
#!/usr/bin/env python3
"""
analyze_trajectory.py - Comprehensive protein-ligand trajectory analysis

Performs RMSD, contact, hydrogen bond, and MMPBSA analysis on GROMACS trajectories.
"""

import argparse
import sys
from pathlib import Path

# Add project root to path
sys.path.insert(0, str(Path(__file__).parent))

from modules.load_align import load_and_align_trajectory
from modules.rmsd import calculate_rmsd
from modules.contacts import calculate_contacts
from modules.hydrogen_bonds import calculate_hydrogen_bonds
from modules.mmpbsa_plots import plot_mmpbsa_results
from utils.cli import parse_arguments
from utils.io import setup_output_directory


def main():
    """Main entry point for trajectory analysis."""
    args = parse_arguments()
    
    # Setup output directory
    output_dir = setup_output_directory(args.output_dir)
    
    # Run analyses based on flags
    if args.aim in ['all', '0']:
        # Load and align (always needed first)
        universe, reference = load_and_align_trajectory(
            topology=args.tpr,
            trajectory=args.xtc,
            index_file=args.ndx,
            align_selection=args.receptor_backbone_sel
        )
        print(f"Loaded and aligned trajectory: {universe.trajectory.n_frames} frames")
    
    if args.aim in ['all', '1']:
        # RMSD analysis
        rmsd_results = calculate_rmsd(
            universe=universe,
            receptor_sel=args.receptor_backbone_sel,
            ligand_sel=args.ligand_sel,
            output_dir=output_dir
        )
        print(f"RMSD analysis complete: {len(rmsd_results)} frames")
    
    if args.aim in ['all', '2']:
        # Contact analysis
        contact_results = calculate_contacts(
            universe=universe,
            receptor_sel=args.receptor_sel,
            ligand_sel=args.ligand_sel,
            cutoff=args.contact_cutoff,
            output_dir=output_dir
        )
        print(f"Contact analysis complete: {len(contact_results)} frames")
    
    if args.aim in ['all', '3']:
        # Hydrogen bond analysis
        hbond_results = calculate_hydrogen_bonds(
            universe=universe,
            receptor_sel=args.receptor_sel,
            ligand_sel=args.ligand_sel,
            distance_cutoff=args.hbond_distance_cutoff,
            angle_cutoff=args.hbond_angle_cutoff,
            output_dir=output_dir
        )
        print(f"Hydrogen bond analysis complete")
    
    if args.aim in ['all', '4']:
        # MMPBSA plotting
        plot_mmpbsa_results(
            mmpbsa_dir=args.mmpbsa_dir,
            output_dir=output_dir
        )
        print(f"MMPBSA plots generated")
    
    print(f"\nAll analyses complete. Output saved to: {output_dir}")


if __name__ == '__main__':
    main()
```

---

## 2. Default Variables Configuration

All default parameters are defined in a configuration section at the top of the main script and in `utils/cli.py`:

```python
# ============================================================================
# DEFAULT CONFIGURATION
# ============================================================================

DEFAULTS = {
    # File paths
    'topology': 'com.tpr',
    'trajectory': 'com_traj.xtc',
    'index_file': 'index.ndx',
    'mmpbsa_dir': '.',
    
    # Selection strings (GROMACS index groups or MDAnalysis selection syntax)
    'receptor_backbone_sel': 'protein and name N CA C O',  # For alignment
    'receptor_sel': 'protein',  # For contacts/hbonds
    'ligand_sel': 'resname MOL',  # Ligand selection
    
    # Analysis parameters
    'contact_cutoff': 4.5,  # Angstroms
    'hbond_distance_cutoff': 3.5,  # Angstroms
    'hbond_angle_cutoff': 150,  # Degrees
    
    # Output settings
    'output_dir': 'output',
    'plot_dpi': 300,
    'plot_format': 'png',
}
```

---

## 3. Command Line Arguments

The script uses `argparse` for command-line argument handling:

```python
# utils/cli.py

import argparse
from pathlib import Path


def parse_arguments():
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(
        description='Comprehensive protein-ligand trajectory analysis',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Run all analyses with default settings
  python analyze_trajectory.py --tpr com.tpr --xtc com_traj.xtc

  # Run only RMSD analysis
  python analyze_trajectory.py --aim 1 --tpr com.tpr --xtc com_traj.xtc

  # Custom output directory and ligand selection
  python analyze_trajectory.py --output ./results --ligand-sel "resname LIG"
        """
    )
    
    # Input files (required)
    parser.add_argument('--tpr', type=str, default='com.tpr',
                        help='Path to TPR topology file (default: com.tpr)')
    parser.add_argument('--xtc', type=str, default='com_traj.xtc',
                        help='Path to XTC trajectory file (default: com_traj.xtc)')
    parser.add_argument('--ndx', type=str, default=None,
                        help='Path to GROMACS index file (default: None)')
    
    # Selection strings
    parser.add_argument('--receptor-backbone-sel', type=str,
                        default='protein and name N CA C O',
                        help='Selection for receptor backbone (default: "protein and name N CA C O")')
    parser.add_argument('--receptor-sel', type=str, default='protein',
                        help='Selection for receptor (default: "protein")')
    parser.add_argument('--ligand-sel', type=str, default='resname MOL',
                        help='Selection for ligand (default: "resname MOL")')
    
    # Analysis parameters
    parser.add_argument('--contact-cutoff', type=float, default=4.5,
                        help='Contact distance cutoff in Angstroms (default: 4.5)')
    parser.add_argument('--hbond-distance-cutoff', type=float, default=3.5,
                        help='Hydrogen bond distance cutoff in Angstroms (default: 3.5)')
    parser.add_argument('--hbond-angle-cutoff', type=float, default=150,
                        help='Hydrogen bond angle cutoff in degrees (default: 150)')
    
    # Output settings
    parser.add_argument('--output-dir', type=str, default='output',
                        help='Output directory (default: output)')
    parser.add_argument('--mmpbsa-dir', type=str, default='.',
                        help='Directory containing MMPBSA output files (default: .)')
    
    # Analysis selection
    parser.add_argument('--aim', type=str, default='all',
                        choices=['all', '0', '1', '2', '3', '4'],
                        help='Analysis to run: 0=load/align, 1=RMSD, 2=contacts, 3=hbonds, 4=mmpbsa (default: all)')
    
    # Plotting options
    parser.add_argument('--plot-dpi', type=int, default=300,
                        help='DPI for output plots (default: 300)')
    parser.add_argument('--plot-format', type=str, default='png',
                        choices=['png', 'pdf', 'svg'],
                        help='Output format for plots (default: png)')
    
    # Performance options
    parser.add_argument('--skip-align', action='store_true',
                        help='Skip alignment (trajectory already aligned)')
    parser.add_argument('--n-jobs', type=int, default=1,
                        help='Number of parallel jobs (default: 1)')
    
    return parser.parse_args()
```

---

## 4. Data Flow

### 4.1 Overall Data Flow

```
┌─────────────────────────────────────────────────────────────────────────────┐
│                            MAIN DATA FLOW                                   │
└─────────────────────────────────────────────────────────────────────────────┘

    INPUT FILES                    PROCESSING                    OUTPUT FILES
    ───────────                    ──────────                    ───────────

┌──────────────┐             ┌──────────────────┐         ┌─────────────────┐
│  com.tpr     │────────────▶│                  │         │ rmsd_receptor.csv
│  (topology)  │             │  MODULE:         │         │ rmsd_ligand.csv │
└──────────────┘             │  load_align.py   │         │ rmsd_plot.png   │
                              │                  │         └─────────────────┘
┌──────────────┐             │  - Load Universe │         ┌─────────────────┐
│ com_traj.xtc │────────────▶│  - Atom Select   │         │ contacts.csv    │
│ (trajectory) │             │  - Align         │────────▶│ contacts_plot.png
└──────────────┘             │                  │         └─────────────────┘
                              └──────────────────┘         ┌─────────────────┐
                                                             │ hbond_timeseries.csv
┌──────────────┐             ┌──────────────────┐         │ hbond_per_residue.csv
│ index.ndx    │────────────▶│  MODULE:          │         │ hbond_plot.png  │
│ (optional)   │             │  rmsd.py          │         └─────────────────┘
└──────────────┘             │  contacts.py      │         ┌─────────────────┐
                              │  hydrogen_bonds.py│         │ energy_plot.png │
┌──────────────┐             │  mmpbsa_plots.py  │         │ decomp_plot.png │
│ MMPBSA files │────────────▶│                  │────────▶│ ...             │
│ (.csv, .dat) │             └──────────────────┘         └─────────────────┘
└──────────────┘
```

### 4.2 Per-Module Data Flow

#### Aim 0: Load and Align

```
Input: TPR file + XTC file (+ optional NDX)
       │
       ▼
┌──────────────────┐
│ MDAnalysis       │
│ Universe()       │────▶ Handle atom count mismatch
└──────────────────┘         │
                             ▼
                    ┌──────────────────┐
                    │ Select atoms     │
                    │ (trajectory size)│
                    └──────────────────┘
                             │
                             ▼
                    ┌──────────────────┐
                    │ fit_rot_trans    │
                    │ (align to ref)   │
                    └──────────────────┘
                             │
                             ▼
Output: Aligned Universe object + reference AtomGroup
```

#### Aim 1: RMSD Calculation

```
Input: Aligned Universe + selections
       │
       ▼
┌──────────────────┐
│ RMSD Analysis    │
│ (MDAnalysis)     │
└──────────────────┘
       │
       ├──▶ Receptor backbone RMSD
       │         │
       │         ▼
       │    ┌──────────────────┐
       │    │ Reference: frame 0│
       │    │ Target: all frames│
       │    └──────────────────┘
       │
       ▼
Ligand heavy atoms RMSD
       │
       ▼
┌──────────────────┐
│ Calculate RMSD   │
│ per frame        │
└──────────────────┘
       │
       ▼
Output: CSV (frame, time, rmsd) + PNG plot
```

#### Aim 2: Contact Analysis

```
Input: Aligned Universe + selections
       │
       ▼
┌──────────────────┐
│ Per-Residue      │
│ Contacts         │
└──────────────────┘
       │
       ▼
┌──────────────────┐
│ For each frame:  │
│  - distance_array│
│  - Group by res  │
│  - Count contacts│
└──────────────────┘
       │
       ▼
Output: CSV (time, res1, res2, ...) + PNG (heatmap/bar)
```

#### Aim 3: Hydrogen Bond Analysis

```
Input: Aligned Universe + selections
       │
       ▼
┌──────────────────┐
│ HydrogenBond     │
│ Analysis         │
│ (MDAnalysis)     │
└──────────────────┘
       │
       ▼
┌──────────────────┐
│ between=[rec,lig]│
│ d_a_cutoff=3.5Å │
│ angle=150°       │
└──────────────────┘
       │
       ▼
┌──────────────────┐
│ Aggregate by     │
│ receptor residue │
└──────────────────┘
       │
       ▼
Output: CSV (timeseries + per-residue) + PNG plot
```

#### Aim 4: MMPBSA Plotting

```
Input: MMPBSA CSV/DAT files
       │
       ▼
┌──────────────────┐
│ Parse CSV files  │
│ (results + decomp)│
└──────────────────┘
       │
       ▼
┌──────────────────┐
│ Generate plots   │
│ - Energy line    │
│ - Summary bar    │
│ - Decomp heatmap │
│ - Decomp bar     │
└──────────────────┘
       │
       ▼
Output: PNG/PDF plots
```

---

## 5. Output File Naming Conventions

All output files follow a consistent naming pattern:

```
output/
├── rmsd/
│   ├── rmsd_receptor.csv          # Receptor backbone RMSD time series
│   ├── rmsd_ligand.csv            # Ligand heavy atom RMSD time series
│   ├── rmsd_combined.csv          # Both RMSD values
│   └── rmsd_plot.png              # Overlaid RMSD plot
│
├── contacts/
│   ├── contacts_timeseries.csv    # Wide format: time, res_1, res_2, ...
│   ├── contacts_summary.csv       # Per-residue statistics
│   ├── contacts_total.csv         # Total contacts per frame
│   ├── contacts_heatmap.png       # Heatmap over time
│   ├── contacts_bar.png           # Bar chart of average contacts
│   └── contacts_timeseries.png    # Time series plot
│
├── hydrogen_bonds/
│   ├── hbond_timeseries.csv       # Total hbonds per frame
│   ├── hbond_per_residue.csv      # Per-residue hbond counts
│   ├── hbond_details.csv          # Detailed hbond info per frame
│   ├── hbond_plot.png             # Combined plot (2 panels)
│   └── hbond_timecourse.png       # Time series of hbond count
│
└── mmpbsa/
    ├── energy_timecourse.png      # Energy components over time
    ├── energy_summary.png         # Bar chart of energy terms
    ├── decomposition_bar.png      # Per-residue decomposition
    ├── decomposition_heatmap.png  # Heatmap of residues over time
    └── (PDF versions for publication)
```

### 5.1 CSV Format Specifications

#### RMSD CSV

```csv
# RMSD Analysis Results
# Receptor: protein and name N CA C O
# Ligand: resname MOL
# Units: Angstroms
frame,time_ps,rmsd_receptor,rmsd_ligand
0,0.0,2.45,3.12
1,100.0,2.51,3.08
...
```

#### Contacts CSV (Wide Format)

```csv
# Contact Analysis Results
# Cutoff: 4.5 Angstrom
# Receptor: protein
# Ligand: resname MOL
time_ps,100,101,102,103,...
0.0,5,2,0,3,...
100.0,4,3,1,2,...
```

#### Hydrogen Bond CSV

```csv
# Hydrogen Bond Analysis Results
# Distance cutoff: 3.5 Angstrom
# Angle cutoff: 150 degrees
frame,time_ps,hbond_count
0,0.0,3
1,100.0,4
...
```

---

## 6. Module Implementation Details

### 6.1 Module: load_align.py (Aim 0)

```python
#!/usr/bin/env python3
"""
load_align.py - Load and align GROMACS trajectory

Functions:
    - load_universe: Load MDAnalysis Universe with proper atom selection
    - align_trajectory: Align trajectory to reference structure
    - get_selections: Get receptor and ligand atom groups
"""

import numpy as np
import MDAnalysis as mda
from MDAnalysis.transformations import fit_rot_trans
from typing import Tuple, Optional


def load_universe(topology: str, trajectory: str, index_file: Optional[str] = None) -> mda.Universe:
    """
    Load MDAnalysis Universe with atom count mismatch handling.
    
    Parameters
    ----------
    topology : str
        Path to TPR file
    trajectory : str
        Path to XTC/TRR file
    index_file : str, optional
        Path to GROMACS index file
    
    Returns
    -------
    universe : MDAnalysis.Universe
        Loaded universe object
    """
    # Load universe
    if index_file and Path(index_file).exists():
        u = mda.Universe(topology, trajectory, index=index_file)
    else:
        u = mda.Universe(topology, trajectory)
    
    # Handle atom count mismatch
    u.trajectory.rewind()
    n_atoms_traj = u.trajectory.n_atoms
    n_atoms_topo = u.atoms.n_atoms
    
    if n_atoms_traj < n_atoms_topo:
        print(f"Warning: Topology has {n_atoms_topo} atoms, "
              f"trajectory has {n_atoms_traj}. Truncating...")
        # Select only atoms present in trajectory
        u.atoms = u.atoms[:n_atoms_traj]
    
    return u


def align_trajectory(universe: mda.Universe, 
                    reference_sel: str = 'protein and name N CA C O') -> Tuple[mda.Universe, mda.AtomGroup]:
    """
    Align trajectory to receptor backbone.
    
    Parameters
    ----------
    universe : MDAnalysis.Universe
        Universe with trajectory
    reference_sel : str
        Selection string for reference atoms
    
    Returns
    -------
    universe : MDAnalysis.Universe
        Aligned universe (transformations applied)
    reference_atoms : AtomGroup
        Reference atom group used for alignment
    """
    # Select reference atoms
    reference_atoms = universe.select_atoms(reference_sel)
    
    print(f"Using {len(reference_atoms)} atoms for alignment")
    print(f"Selection: {reference_sel}")
    
    # Create alignment transformation
    alignment = fit_rot_trans(universe, reference_atoms)
    universe.trajectory.add_transformations(alignment)
    
    return universe, reference_atoms


def get_selections(universe: mda.Universe,
                   receptor_sel: str = 'protein',
                   ligand_sel: str = 'resname MOL') -> Tuple[mda.AtomGroup, mda.AtomGroup]:
    """
    Get receptor and ligand atom groups.
    
    Parameters
    ----------
    universe : MDAnalysis.Universe
        Loaded universe
    receptor_sel : str
        Selection string for receptor
    ligand_sel : str
        Selection string for ligand
    
    Returns
    -------
    receptor : AtomGroup
        Receptor atoms
    ligand : AtomGroup
        Ligand atoms
    """
    receptor = universe.select_atoms(receptor_sel)
    ligand = universe.select_atoms(ligand_sel)
    
    print(f"Receptor: {len(receptor)} atoms ({receptor.n_residues} residues)")
    print(f"Ligand: {len(ligand)} atoms")
    
    return receptor, ligand


def load_and_align_trajectory(topology: str, trajectory: str,
                              index_file: Optional[str] = None,
                              align_selection: str = 'protein and name N CA C O',
                              receptor_sel: str = 'protein',
                              ligand_sel: str = 'resname MOL') -> Tuple[mda.Universe, mda.AtomGroup, mda.AtomGroup, mda.AtomGroup]:
    """
    Complete workflow: load trajectory and align.
    
    Returns
    -------
    universe : Universe
        Aligned universe
    reference_atoms : AtomGroup
        Atoms used for alignment
    receptor : AtomGroup
        Receptor atoms
    ligand : AtomGroup
        Ligand atoms
    """
    # Load universe
    universe = load_universe(topology, trajectory, index_file)
    
    # Align trajectory
    universe, reference_atoms = align_trajectory(universe, align_selection)
    
    # Get selections
    receptor, ligand = get_selections(universe, receptor_sel, ligand_sel)
    
    return universe, reference_atoms, receptor, ligand
```

### 6.2 Module: rmsd.py (Aim 1)

```python
#!/usr/bin/env python3
"""
rmsd.py - RMSD calculations for receptor backbone and ligand heavy atoms

Functions:
    - calculate_receptor_rmsd: RMSD of receptor backbone
    - calculate_ligand_rmsd: RMSD of ligand heavy atoms
    - plot_rmsd_overlay: Create overlaid RMSD plot
"""

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from typing import Tuple, List
from MDAnalysis import Universe, AtomGroup
from MDAnalysis.analysis import align
from pathlib import Path


def calculate_receptor_rmsd(universe: Universe, 
                            selection: str = 'protein and name N CA C O') -> pd.DataFrame:
    """
    Calculate RMSD for receptor backbone atoms.
    
    Parameters
    ----------
    universe : Universe
        Aligned universe with trajectory
    selection : str
        Selection string for receptor backbone
    
    Returns
    -------
    df : DataFrame
        Columns: frame, time_ps, rmsd
    """
    # Select receptor backbone
    backbone = universe.select_atoms(selection)
    
    # Use first frame as reference
    reference = backbone.positions.copy()
    
    results = []
    for ts in universe.trajectory:
        # Calculate RMSD
        rmsd = np.sqrt(np.mean(np.sum((backbone.positions - reference)**2, axis=1)))
        results.append({
            'frame': ts.frame,
            'time_ps': ts.time,
            'rmsd': rmsd
        })
    
    df = pd.DataFrame(results)
    return df


def calculate_ligand_rmsd(universe: Universe,
                          selection: str = 'resname MOL and not name H*') -> pd.DataFrame:
    """
    Calculate RMSD for ligand heavy atoms.
    
    Parameters
    ----------
    universe : Universe
        Aligned universe with trajectory
    selection : str
        Selection string for ligand heavy atoms
    
    Returns
    -------
    df : DataFrame
        Columns: frame, time_ps, rmsd
    """
    # Select ligand heavy atoms (exclude hydrogens)
    ligand = universe.select_atoms(selection)
    
    if len(ligand) == 0:
        print("Warning: No ligand atoms found with selection, trying all ligand atoms")
        ligand = universe.select_atoms('resname MOL')
    
    # Use first frame as reference
    reference = ligand.positions.copy()
    
    results = []
    for ts in universe.trajectory:
        # Calculate RMSD
        rmsd = np.sqrt(np.mean(np.sum((ligand.positions - reference)**2, axis=1)))
        results.append({
            'frame': ts.frame,
            'time_ps': ts.time,
            'rmsd': rmsd
        })
    
    df = pd.DataFrame(results)
    return df


def plot_rmsd_overlay(df_receptor: pd.DataFrame, df_ligand: pd.DataFrame,
                      output_path: str, dpi: int = 300) -> None:
    """
    Create overlaid RMSD plot.
    
    Parameters
    ----------
    df_receptor : DataFrame
        Receptor RMSD data
    df_ligand : DataFrame
        Ligand RMSD data
    output_path : str
        Output file path
    dpi : int
        Plot resolution
    """
    fig, ax = plt.subplots(figsize=(10, 6))
    
    # Plot both RMSD curves
    ax.plot(df_receptor['time_ps'], df_receptor['rmsd'], 
            'b-', linewidth=1.5, label='Receptor Backbone', alpha=0.8)
    ax.plot(df_ligand['time_ps'], df_ligand['rmsd'], 
            'r-', linewidth=1.5, label='Ligand Heavy Atoms', alpha=0.8)
    
    # Add moving averages
    for df, color, label in [(df_receptor, 'blue', 'Receptor (MA)'),
                              (df_ligand, 'red', 'Ligand (MA)')]:
        window = min(10, len(df))
        if window > 1:
            ma = df['rmsd'].rolling(window=window, center=True).mean()
            ax.plot(df['time_ps'], ma, color=color, linewidth=2, 
                    linestyle='--', alpha=0.6, label=label)
    
    # Formatting
    ax.set_xlabel('Time (ps)', fontsize=12)
    ax.set_ylabel('RMSD (Å)', fontsize=12)
    ax.set_title('RMSD: Receptor Backbone vs Ligand Heavy Atoms', fontsize=14)
    ax.legend(loc='best', framealpha=0.9)
    ax.grid(True, alpha=0.3)
    ax.set_xlim(df_receptor['time_ps'].min(), df_receptor['time_ps'].max())
    
    # Add statistics text
    text = (f"Receptor: mean={df_receptor['rmsd'].mean():.2f}Å, "
            f"std={df_receptor['rmsd'].std():.2f}Å\n"
            f"Ligand: mean={df_ligand['rmsd'].mean():.2f}Å, "
            f"std={df_ligand['rmsd'].std():.2f}Å")
    ax.text(0.02, 0.98, text, transform=ax.transAxes, fontsize=9,
            verticalalignment='top', bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))
    
    plt.tight_layout()
    plt.savefig(output_path, dpi=dpi, bbox_inches='tight')
    plt.close()
    
    print(f"Saved RMSD plot to {output_path}")


def calculate_rmsd(universe: Universe,
                   receptor_sel: str = 'protein and name N CA C O',
                   ligand_sel: str = 'resname MOL and not name H*',
                   output_dir: str = 'output/rmsd') -> Tuple[pd.DataFrame, pd.DataFrame]:
    """
    Complete RMSD analysis workflow.
    
    Returns
    -------
    df_receptor : DataFrame
        Receptor RMSD results
    df_ligand : DataFrame
        Ligand RMSD results
    """
    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)
    
    # Calculate RMSD
    print("Calculating receptor backbone RMSD...")
    df_receptor = calculate_receptor_rmsd(universe, receptor_sel)
    
    print("Calculating ligand heavy atom RMSD...")
    df_ligand = calculate_ligand_rmsd(universe, ligand_sel)
    
    # Save CSV files
    df_receptor.to_csv(output_path / 'rmsd_receptor.csv', index=False)
    df_ligand.to_csv(output_path / 'rmsd_ligand.csv', index=False)
    
    # Combined CSV
    df_combined = pd.merge(df_receptor, df_ligand, on=['frame', 'time_ps'], 
                           suffixes=('_receptor', '_ligand'))
    df_combined.to_csv(output_path / 'rmsd_combined.csv', index=False)
    
    # Create plot
    plot_rmsd_overlay(df_receptor, df_ligand, 
                      str(output_path / 'rmsd_plot.png'))
    
    return df_receptor, df_ligand
```

### 6.3 Module: contacts.py (Aim 2)

```python
#!/usr/bin/env python3
"""
contacts.py - Per-residue contact analysis between receptor and ligand

Functions:
    - calculate_per_residue_contacts: Calculate contacts per residue
    - plot_contact_heatmap: Create contact heatmap
    - plot_contact_bar: Create bar chart of average contacts
"""

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from MDAnalysis import Universe
from MDAnalysis.analysis import distances
from typing import Dict, List
from pathlib import Path
from collections import defaultdict


def calculate_per_residue_contacts(universe: Universe,
                                   receptor_sel: str = 'protein',
                                   ligand_sel: str = 'resname MOL',
                                   cutoff: float = 4.5) -> Dict:
    """
    Calculate per-residue contact numbers between receptor and ligand.
    
    Parameters
    ----------
    universe : Universe
        Aligned universe with trajectory
    receptor_sel : str
        Selection for receptor
    ligand_sel : str
        Selection for ligand
    cutoff : float
        Distance cutoff in Angstroms
    
    Returns
    -------
    timeseries : dict
        {residue_id: [count_frame0, count_frame1, ...]}
    times : list
        Time points in ps
    """
    # Get atom groups
    receptor = universe.select_atoms(receptor_sel)
    ligand = universe.select_atoms(ligand_sel)
    residues = receptor.residues
    
    # Build residue atom index mapping
    res_atom_indices = {}
    for res in residues:
        res_atom_indices[res.resid] = res.atoms.indices
    
    # Initialize results
    contact_timeseries = {res.resid: [] for res in residues}
    times = []
    
    # Calculate contacts per frame
    for ts in universe.trajectory:
        times.append(ts.time)
        
        # Distance matrix: receptor atoms x ligand atoms
        dist_mat = distances.distance_array(
            receptor.positions,
            ligand.positions,
            box=universe.dimensions
        )
        
        # Calculate contacts per residue
        for resid, atom_indices in res_atom_indices.items():
            # Distances from this residue's atoms to all ligand atoms
            res_distances = dist_mat[atom_indices, :]
            
            # Count: any atom within cutoff (for boolean contact)
            # Or sum of all atom pairs within cutoff (for contact number)
            n_contacts = (res_distances <= cutoff).sum()
            
            contact_timeseries[resid].append(n_contacts)
    
    return contact_timeseries, times


def save_contacts_csv(contact_timeseries: Dict, times: List[float],
                      output_dir: Path) -> None:
    """
    Save contact time series to CSV files.
    """
    # Wide format CSV
    df = pd.DataFrame(contact_timeseries, index=times)
    df.index.name = 'time_ps'
    df.to_csv(output_dir / 'contacts_timeseries.csv')
    
    # Summary statistics
    summary_data = []
    for resid, contacts in contact_timeseries.items():
        contacts_arr = np.array(contacts)
        summary_data.append({
            'residue_id': resid,
            'mean_contacts': np.mean(contacts_arr),
            'std_contacts': np.std(contacts_arr),
            'max_contacts': np.max(contacts_arr),
            'min_contacts': np.min(contacts_arr),
            'fraction_in_contact': np.mean(contacts_arr > 0)
        })
    
    df_summary = pd.DataFrame(summary_data)
    df_summary = df_summary.sort_values('mean_contacts', ascending=False)
    df_summary.to_csv(output_dir / 'contacts_summary.csv', index=False)
    
    # Total contacts per frame
    df_total = pd.DataFrame({'time_ps': times, 
                             'total_contacts': df.sum(axis=1).values})
    df_total.to_csv(output_dir / 'contacts_total.csv', index=False)


def plot_contact_heatmap(contact_timeseries: Dict, times: List[float],
                         output_path: str, top_n: int = 50, dpi: int = 300) -> None:
    """
    Plot heatmap of contacts over time.
    """
    df = pd.DataFrame(contact_timeseries, index=times)
    
    # Select top N residues by average contacts
    avg_contacts = df.mean().sort_values(ascending=False).head(top_n)
    df_top = df[avg_contacts.index]
    
    fig, ax = plt.subplots(figsize=(12, 8))
    
    im = ax.imshow(df_top.values.T, aspect='auto', origin='lower',
                   cmap='YlOrRd', interpolation='nearest')
    
    ax.set_xlabel('Time (ps)', fontsize=12)
    ax.set_ylabel('Residue ID', fontsize=12)
    ax.set_title(f'Receptor-Ligand Contact Heatmap (Top {top_n} Residues)', fontsize=14)
    
    cbar = plt.colorbar(im, ax=ax)
    cbar.set_label('Number of Contacts', fontsize=12)
    
    # Set x-axis ticks
    n_ticks = min(10, len(times))
    tick_positions = np.linspace(0, len(times)-1, n_ticks).astype(int)
    tick_labels = [f'{times[i]:.0f}' for i in tick_positions]
    ax.set_xticks(tick_positions)
    ax.set_xticklabels(tick_labels)
    
    plt.tight_layout()
    plt.savefig(output_path, dpi=dpi, bbox_inches='tight')
    plt.close()
    
    print(f"Saved contact heatmap to {output_path}")


def plot_contact_bar(contact_timeseries: Dict, output_path: str, 
                     top_n: int = 20, dpi: int = 300) -> None:
    """
    Plot bar chart of average contacts per residue.
    """
    df = pd.DataFrame(contact_timeseries)
    avg_contacts = df.mean().sort_values(ascending=False).head(top_n)
    
    fig, ax = plt.subplots(figsize=(10, 6))
    
    y_pos = np.arange(len(avg_contacts))
    ax.barh(y_pos, avg_contacts.values, color='steelblue', edgecolor='black', alpha=0.8)
    
    ax.set_yticks(y_pos)
    ax.set_yticklabels(avg_contacts.index)
    ax.set_xlabel('Average Contact Number', fontsize=12)
    ax.set_title(f'Top {top_n} Residues by Average Contacts', fontsize=14)
    ax.invert_yaxis()
    ax.grid(True, alpha=0.3, axis='x')
    
    # Add value labels
    for i, v in enumerate(avg_contacts.values):
        ax.text(v + 0.5, i, f'{v:.1f}', va='center', fontsize=9)
    
    plt.tight_layout()
    plt.savefig(output_path, dpi=dpi, bbox_inches='tight')
    plt.close()
    
    print(f"Saved contact bar chart to {output_path}")


def calculate_contacts(universe: Universe,
                       receptor_sel: str = 'protein',
                       ligand_sel: str = 'resname MOL',
                       cutoff: float = 4.5,
                       output_dir: str = 'output/contacts') -> Dict:
    """
    Complete contact analysis workflow.
    """
    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)
    
    print(f"Calculating contacts with {cutoff}Å cutoff...")
    contact_timeseries, times = calculate_per_residue_contacts(
        universe, receptor_sel, ligand_sel, cutoff
    )
    
    # Save CSV files
    save_contacts_csv(contact_timeseries, times, output_path)
    
    # Generate plots
    plot_contact_heatmap(contact_timeseries, times, 
                         str(output_path / 'contacts_heatmap.png'))
    plot_contact_bar(contact_timeseries, 
                     str(output_path / 'contacts_bar.png'))
    
    return contact_timeseries
```

### 6.4 Module: hydrogen_bonds.py (Aim 3)

```python
#!/usr/bin/env python3
"""
hydrogen_bonds.py - Per-residue hydrogen bond analysis

Functions:
    - calculate_hydrogen_bonds: Run hydrogen bond analysis
    - aggregate_per_residue: Aggregate hbonds by residue
    - plot_hydrogen_bonds: Create hbond plots
"""

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from MDAnalysis import Universe
from MDAnalysis.analysis.hydrogenbonds.hbond_analysis import HydrogenBondAnalysis as HBA
from typing import Tuple, Dict, List
from pathlib import Path
from collections import defaultdict


def calculate_hydrogen_bonds(universe: Universe,
                              receptor_sel: str = 'protein',
                              ligand_sel: str = 'resname MOL',
                              distance_cutoff: float = 3.5,
                              angle_cutoff: float = 150) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """
    Calculate hydrogen bonds between receptor and ligand.
    
    Parameters
    ----------
    universe : Universe
        Aligned universe with trajectory
    receptor_sel : str
        Selection for receptor
    ligand_sel : str
        Selection for ligand
    distance_cutoff : float
        Donor-acceptor distance cutoff (Å)
    angle_cutoff : float
        Donor-H-acceptor angle cutoff (degrees)
    
    Returns
    -------
    df_time : DataFrame
        Time series: frame, time_ps, hbond_count
    df_residue : DataFrame
        Per-residue: residue, hbond_count, avg_distance, avg_angle
    """
    print(f"Running hydrogen bond analysis (dist={distance_cutoff}Å, angle={angle_cutoff}°)...")
    
    # Initialize hydrogen bond analysis
    hbonds = HBA(
        universe=universe,
        between=[receptor_sel, ligand_sel],
        d_a_cutoff=distance_cutoff,
        d_h_a_angle_cutoff=angle_cutoff,
        update_selctions=False
    )
    
    # Run analysis
    hbonds.run(verbose=False)
    
    hbond_data = hbonds.results.hbonds
    
    if len(hbond_data) == 0:
        print("Warning: No hydrogen bonds found!")
        return pd.DataFrame(), pd.DataFrame()
    
    # Build atom to residue mapping
    atom_to_residue = {}
    for atom in universe.atoms:
        atom_to_residue[atom.index] = {
            'resid': atom.resid,
            'resname': atom.resname,
            'resnum': atom.resnum
        }
    
    # Aggregate per frame
    frames = np.unique(hbond_data[:, 0])
    time_series = []
    
    for frame in frames:
        frame_hbonds = hbond_data[hbond_data[:, 0] == frame]
        n_hbonds = len(frame_hbonds)
        
        # Get time
        universe.trajectory[int(frame)]
        time_ps = universe.trajectory.time
        
        time_series.append({
            'frame': int(frame),
            'time_ps': time_ps,
            'hbond_count': n_hbonds
        })
    
    df_time = pd.DataFrame(time_series)
    
    # Aggregate per residue
    residue_data = defaultdict(lambda: {
        'count': 0,
        'distances': [],
        'angles': []
    })
    
    for hbond in hbond_data:
        frame, donor_idx, hydrogen_idx, acceptor_idx, distance, angle = hbond
        
        # Get residue info for both donor and acceptor
        for idx in [donor_idx, acceptor_idx]:
            res_info = atom_to_residue.get(idx)
            if res_info and res_info['resname'] != ligand_sel.replace('resname ', ''):
                key = f"{res_info['resname']}{res_info['resnum']}"
                residue_data[key]['count'] += 1
                residue_data[key]['distances'].append(distance)
                residue_data[key]['angles'].append(angle)
    
    # Build per-residue DataFrame
    residue_rows = []
    for key, data in residue_data.items():
        residue_rows.append({
            'residue': key,
            'hbond_count': data['count'],
            'avg_distance_A': np.mean(data['distances']) if data['distances'] else 0,
            'avg_angle_deg': np.mean(data['angles']) if data['angles'] else 0
        })
    
    df_residue = pd.DataFrame(residue_rows)
    df_residue = df_residue.sort_values('hbond_count', ascending=False)
    
    return df_time, df_residue


def plot_hydrogen_bonds(df_time: pd.DataFrame, df_residue: pd.DataFrame,
                        output_path: str, dpi: int = 300) -> None:
    """
    Create combined hydrogen bond plot.
    """
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 8))
    
    # Panel 1: Time series
    ax1.plot(df_time['time_ps'], df_time['hbond_count'], 'b-', linewidth=1, alpha=0.7)
    ax1.fill_between(df_time['time_ps'], 0, df_time['hbond_count'], alpha=0.3)
    
    # Moving average
    window = min(10, len(df_time))
    if window > 1:
        ma = df_time['hbond_count'].rolling(window=window, center=True).mean()
        ax1.plot(df_time['time_ps'], ma, 'r-', linewidth=2, label=f'{window}-frame MA')
        ax1.legend()
    
    ax1.set_xlabel('Time (ps)', fontsize=12)
    ax1.set_ylabel('Hydrogen Bond Count', fontsize=12)
    ax1.set_title('Receptor-Ligand Hydrogen Bonds Over Time', fontsize=14)
    ax1.grid(True, alpha=0.3)
    
    # Statistics text
    mean_hb = df_time['hbond_count'].mean()
    std_hb = df_time['hbond_count'].std()
    ax1.axhline(mean_hb, color='red', linestyle='--', label=f'Mean: {mean_hb:.1f}')
    
    # Panel 2: Per-residue bar chart (top 15)
    top_n = 15
    top_res = df_residue.head(top_n)
    
    y_pos = np.arange(len(top_res))
    colors = plt.cm.viridis(np.linspace(0.2, 0.8, len(top_res)))
    
    ax2.barh(y_pos, top_res['hbond_count'], color=colors, edgecolor='black', alpha=0.8)
    ax2.set_yticks(y_pos)
    ax2.set_yticklabels(top_res['residue'])
    ax2.set_xlabel('Total Hydrogen Bond Count', fontsize=12)
    ax2.set_title(f'Top {top_n} Residues by Hydrogen Bond Frequency', fontsize=14)
    ax2.invert_yaxis()
    ax2.grid(True, alpha=0.3, axis='x')
    
    plt.tight_layout()
    plt.savefig(output_path, dpi=dpi, bbox_inches='tight')
    plt.close()
    
    print(f"Saved hydrogen bond plot to {output_path}")


def calculate_hydrogen_bonds_wrapper(universe: Universe,
                                      receptor_sel: str = 'protein',
                                      ligand_sel: str = 'resname MOL',
                                      distance_cutoff: float = 3.5,
                                      angle_cutoff: float = 150,
                                      output_dir: str = 'output/hydrogen_bonds') -> Tuple[pd.DataFrame, pd.DataFrame]:
    """
    Complete hydrogen bond analysis workflow.
    """
    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)
    
    # Calculate hydrogen bonds
    df_time, df_residue = calculate_hydrogen_bonds(
        universe, receptor_sel, ligand_sel, distance_cutoff, angle_cutoff
    )
    
    # Save CSV files
    df_time.to_csv(output_path / 'hbond_timeseries.csv', index=False)
    df_residue.to_csv(output_path / 'hbond_per_residue.csv', index=False)
    
    # Generate plot
    if len(df_time) > 0:
        plot_hydrogen_bonds(df_time, df_residue, str(output_path / 'hbond_plot.png'))
    
    return df_time, df_residue
```

### 6.5 Module: mmpbsa_plots.py (Aim 4)

```python
#!/usr/bin/env python3
"""
mmpbsa_plots.py - MMPBSA results plotting without GUI

Functions:
    - parse_results_csv: Parse FINAL_RESULTS_MMPBSA.csv
    - parse_decomposition_csv: Parse FINAL_DECOMP_MMPBSA.csv
    - plot_energy_timecourse: Plot energy components over time
    - plot_energy_summary: Plot energy summary bar chart
    - plot_decomposition_bar: Plot per-residue decomposition
    - plot_decomposition_heatmap: Plot decomposition heatmap
"""

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from pathlib import Path
from typing import Dict


# Publication-quality matplotlib settings
plt.rcParams.update({
    'font.family': 'sans-serif',
    'font.sans-serif': ['Arial', 'DejaVu Sans'],
    'font.size': 10,
    'axes.labelsize': 11,
    'axes.titlesize': 12,
    'xtick.labelsize': 9,
    'ytick.labelsize': 9,
    'legend.fontsize': 9,
    'figure.dpi': 300,
    'savefig.dpi': 300,
    'axes.linewidth': 0.8,
    'axes.spines.top': False,
    'axes.spines.right': False,
})

ENERGY_COLORS = {
    'VDWAALS': '#2E86AB',
    'EEL': '#A23B72',
    'EGB': '#F18F01',
    'ESURF': '#C73E1D',
    'GGAS': '#3B1F2B',
    'GSOLV': '#95C623',
    'TOTAL': '#1A1A2E',
    'Delta': '#E63946',
}


def parse_results_csv(csv_path: str) -> Dict[str, pd.DataFrame]:
    """
    Parse FINAL_RESULTS_MMPBSA.csv.
    
    Returns
    -------
    data : dict
        Keys: 'complex', 'receptor', 'ligand', 'delta'
    """
    with open(csv_path, 'r') as f:
        content = f.read()
    
    sections = content.split('\n\n')
    data = {}
    
    for section in sections:
        lines = [l.strip() for l in section.strip().split('\n') if l.strip()]
        if not lines:
            continue
            
        if lines[0] == 'Complex Energy Terms':
            data['complex'] = pd.read_csv(Path(csv_path).parent / '\n'.join(lines[2:]))
        elif lines[0] == 'Receptor Energy Terms':
            data['receptor'] = pd.read_csv(Path(csv_path).parent / '\n'.join(lines[2:]))
        elif lines[0] == 'Ligand Energy Terms':
            data['ligand'] = pd.read_csv(Path(csv_path).parent / '\n'.join(lines[2:]))
        elif lines[0] == 'Delta Energy Terms':
            data['delta'] = pd.read_csv(Path(csv_path).parent / '\n'.join(lines[2:]))
    
    return data


def parse_decomposition_csv(csv_path: str) -> Dict[str, pd.DataFrame]:
    """
    Parse FINAL_DECOMP_MMPBSA.csv.
    """
    with open(csv_path, 'r') as f:
        lines = f.readlines()
    
    data = {'complex': [], 'receptor': [], 'ligand': [], 'delta': []}
    current_section = None
    
    for line in lines:
        line = line.strip()
        
        if line == 'Complex:':
            current_section = 'complex'
        elif line == 'Receptor:':
            current_section = 'receptor'
        elif line == 'Ligand:':
            current_section = 'ligand'
        elif line == 'Delta (Complex - Receptor - Ligand):':
            current_section = 'delta'
        elif line.startswith('Frame #,') or not line:
            continue
        
        parts = line.split(',')
        if len(parts) >= 8 and parts[0].isdigit():
            data[current_section].append({
                'frame': int(parts[0]),
                'residue': parts[1],
                'internal': float(parts[2]),
                'vdw': float(parts[3]),
                'electrostatic': float(parts[4]),
                'polar': float(parts[5]),
                'nonpolar': float(parts[6]),
                'total': float(parts[7])
            })
    
    return {k: pd.DataFrame(v) for k, v in data.items() if v}


def plot_energy_timecourse(df_delta: pd.DataFrame, output_path: str,
                           dpi: int = 300) -> None:
    """
    Plot binding energy components over simulation frames.
    """
    fig, ax = plt.subplots(figsize=(4, 2.5))
    
    frames = df_delta['Frame #']
    terms = ['VDWAALS', 'EEL', 'EGB', 'ESURF', 'TOTAL']
    
    for term in terms:
        if term in df_delta.columns:
            lw = 2 if term == 'TOTAL' else 1.2
            ax.plot(frames, df_delta[term], label=term, linewidth=lw,
                   color=ENERGY_COLORS.get(term, '#333'), alpha=0.9)
    
    # Moving average for TOTAL
    if 'TOTAL' in df_delta.columns:
        ma = df_delta['TOTAL'].rolling(window=5, center=True).mean()
        ax.plot(frames, ma, 'k-', linewidth=2, alpha=0.5, label='5-frame MA')
    
    ax.set_xlabel('Frame')
    ax.set_ylabel('Energy (kcal/mol)')
    ax.set_title('Binding Free Energy Decomposition')
    ax.legend(loc='best', framealpha=0.9, fontsize=8)
    ax.grid(True, alpha=0.3, linestyle=':')
    ax.axhline(y=0, color='gray', linestyle='-', lw=0.5)
    
    plt.tight_layout()
    plt.savefig(output_path, dpi=dpi)
    plt.close()


def plot_energy_summary(df_delta: pd.DataFrame, output_path: str,
                        dpi: int = 300) -> None:
    """
    Plot average energy components as bar chart with error bars.
    """
    terms = ['VDWAALS', 'EEL', 'EGB', 'ESURF', 'GGAS', 'GSOLV', 'TOTAL']
    means = []
    stds = []
    labels = []
    colors_list = []
    
    for term in terms:
        col = term if term in df_delta.columns else f'Δ{term}'
        if col in df_delta.columns:
            means.append(df_delta[col].mean())
            stds.append(df_delta[col].std())
            labels.append(term)
            if term == 'TOTAL':
                colors_list.append(ENERGY_COLORS['Delta'])
            else:
                colors_list.append(ENERGY_COLORS.get(term, '#333'))
    
    fig, ax = plt.subplots(figsize=(4, 3))
    
    x = np.arange(len(labels))
    ax.bar(x, means, yerr=stds, capsize=3, color=colors_list,
           edgecolor='black', lw=0.5, alpha=0.85)
    
    ax.set_xticks(x)
    ax.set_xticklabels(labels, rotation=45, ha='right')
    ax.set_ylabel('Energy (kcal/mol)')
    ax.set_title('Binding Energy Summary')
    ax.axhline(y=0, color='black', lw=0.5)
    
    plt.tight_layout()
    plt.savefig(output_path, dpi=dpi)
    plt.close()


def plot_decomposition_bar(df_delta: pd.DataFrame, output_path: str,
                           top_n: int = 15, dpi: int = 300) -> None:
    """
    Plot per-residue energy contribution to binding.
    """
    if df_delta.empty:
        print("No decomposition data to plot")
        return
    
    # Average per residue
    residue_avg = df_delta.groupby('residue')['total'].mean().sort_values()
    top = residue_avg.head(top_n)
    
    # Shorten labels: R:A:ARG:348 -> R348
    labels = []
    for res in top.index:
        parts = res.split(':')
        if len(parts) >= 4:
            labels.append(f"{parts[0]}{parts[3]}")
        else:
            labels.append(res)
    
    fig, ax = plt.subplots(figsize=(4, 4))
    
    colors = ['#2E86AB' if v < 0 else '#E63946' for v in top.values]
    y = np.arange(len(labels))
    
    ax.barh(y, top.values, color=colors, edgecolor='black', lw=0.5)
    ax.set_yticks(y)
    ax.set_yticklabels(labels)
    ax.set_xlabel('Energy (kcal/mol)')
    ax.set_title('Per-Residue Binding Contribution')
    ax.axvline(x=0, color='black', lw=0.5)
    
    plt.tight_layout()
    plt.savefig(output_path, dpi=dpi)
    plt.close()


def plot_decomposition_heatmap(df_delta: pd.DataFrame, output_path: str,
                               top_n: int = 30, dpi: int = 300) -> None:
    """
    Plot per-residue energy contribution as heatmap over frames.
    """
    if df_delta.empty:
        return
    
    # Pivot: frames x residues
    pivot = df_delta.pivot(index='frame', columns='residue', values='total')
    
    # Select top N by variance or absolute value
    top_cols = pivot.var().sort_values(ascending=False).head(top_n).index
    pivot_top = pivot[top_cols]
    
    # Shorten column names
    short_cols = {}
    for col in pivot_top.columns:
        parts = col.split(':')
        if len(parts) >= 4:
            short_cols[col] = f"{parts[2][0]}{parts[3]}"
        else:
            short_cols[col] = col
    pivot_top = pivot_top.rename(columns=short_cols)
    
    fig, ax = plt.subplots(figsize=(8, 4))
    
    vmax = max(abs(pivot_top.values.min()), abs(pivot_top.values.max()))
    im = ax.imshow(pivot_top.values, aspect='auto', cmap='RdBu_r',
                   vmin=-vmax, vmax=vmax)
    
    ax.set_yticks(range(len(pivot_top.index)))
    ax.set_yticklabels(pivot_top.index)
    ax.set_xticks(range(len(pivot_top.columns)))
    ax.set_xticklabels(pivot_top.columns, rotation=90, fontsize=8)
    
    ax.set_xlabel('Residue')
    ax.set_ylabel('Frame')
    ax.set_title('Per-Residue Energy Contribution (kcal/mol)')
    
    plt.colorbar(im, ax=ax, label='Energy (kcal/mol)', shrink=0.8)
    plt.tight_layout()
    plt.savefig(output_path, dpi=dpi)
    plt.close()


def plot_mmpbsa_results(mmpbsa_dir: str = '.',
                        output_dir: str = 'output/mmpbsa',
                        dpi: int = 300) -> None:
    """
    Complete MMPBSA plotting workflow.
    """
    mmpbsa_path = Path(mmpbsa_dir)
    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)
    
    results_csv = mmpbsa_path / 'FINAL_RESULTS_MMPBSA.csv'
    decomp_csv = mmpbsa_path / 'FINAL_DECOMP_MMPBSA.csv'
    
    # Parse and plot energy results
    if results_csv.exists():
        print(f"Parsing {results_csv}...")
        results = parse_results_csv(str(results_csv))
        
        if 'delta' in results:
            plot_energy_timecourse(results['delta'], 
                                   str(output_path / 'energy_timecourse.png'), dpi)
            plot_energy_summary(results['delta'],
                                str(output_path / 'energy_summary.png'), dpi)
    
    # Parse and plot decomposition
    if decomp_csv.exists():
        print(f"Parsing {decomp_csv}...")
        decomp = parse_decomposition_csv(str(decomp_csv))
        
        if 'delta' in decomp:
            plot_decomposition_bar(decomp['delta'],
                                   str(output_path / 'decomposition_bar.png'), dpi)
            plot_decomposition_heatmap(decomp['delta'],
                                       str(output_path / 'decomposition_heatmap.png'), dpi)
    
    print(f"MMPBSA plots saved to {output_path}")
```

---

## 7. Key Implementation Steps for Each Aim

### 7.1 Aim 0: Load and Align

**Steps:**
1. Create `modules/__init__.py` and `load_align.py`
2. Implement `load_universe()` with atom count mismatch handling
3. Implement `align_trajectory()` using `fit_rot_trans`
4. Implement `get_selections()` for receptor/ligand atoms
5. Export `load_and_align_trajectory()` function

**Key Considerations:**
- Handle cases where trajectory has fewer atoms than topology
- Use index.ndx if available for pre-defined groups
- Default receptor selection: `protein and name N CA C O` (backbone)
- Default ligand selection: `resname MOL` or from index group `[MOL]`

### 7.2 Aim 1: RMSD Calculation

**Steps:**
1. Create `modules/rmsd.py`
2. Implement `calculate_receptor_rmsd()` - backbone RMSD
3. Implement `calculate_ligand_rmsd()` - heavy atom RMSD (exclude H*)
4. Implement `plot_rmsd_overlay()` - overlaid plot with moving averages
5. Export `calculate_rmsd()` with CSV and PNG output

**Key Considerations:**
- Use first frame as reference
- Calculate RMSD per frame: `sqrt(mean(sum((positions - reference)^2, axis=1)))`
- Exclude hydrogen atoms from ligand RMSD
- Include moving average (window=10) in plot

### 7.3 Aim 2: Contact Analysis

**Steps:**
1. Create `modules/contacts.py`
2. Implement `calculate_per_residue_contacts()` with distance_array
3. Implement `save_contacts_csv()` - wide format and summary
4. Implement `plot_contact_heatmap()` - heatmap visualization
5. Implement `plot_contact_bar()` - bar chart of top residues
6. Export `calculate_contacts()` with CSV and PNG output

**Key Considerations:**
- Default cutoff: 4.5Å (standard for all-atom simulations)
- Use `distance_array()` with PBC box
- Output wide format CSV: `time_ps, res_1, res_2, ...`
- Include summary statistics: mean, std, max, min, fraction_in_contact

### 7.4 Aim 3: Hydrogen Bond Analysis

**Steps:**
1. Create `modules/hydrogen_bonds.py`
2. Implement `calculate_hydrogen_bonds()` using HydrogenBondAnalysis
3. Implement per-frame and per-residue aggregation
4. Implement `plot_hydrogen_bonds()` - combined 2-panel plot
5. Export `calculate_hydrogen_bonds_wrapper()` with CSV and PNG output

**Key Considerations:**
- Use `between=[receptor_sel, ligand_sel]` to filter cross-group hbonds
- Default distance cutoff: 3.5Å
- Default angle cutoff: 150°
- Build atom-to-residue mapping for aggregation
- Output time series and per-residue summary

### 7.5 Aim 4: MMPBSA Plotting

**Steps:**
1. Create `modules/mmpbsa_plots.py`
2. Implement `parse_results_csv()` for FINAL_RESULTS_MMPBSA.csv
3. Implement `parse_decomposition_csv()` for FINAL_DECOMP_MMPBSA.csv
4. Implement plotting functions:
   - `plot_energy_timecourse()` - line plot
   - `plot_energy_summary()` - bar chart
   - `plot_decomposition_bar()` - per-residue bar
   - `plot_decomposition_heatmap()` - heatmap
5. Export `plot_mmpbsa_results()` - runs all plots

**Key Considerations:**
- Use CSV parsing (more reliable than DAT)
- Create publication-quality plots (300 DPI, proper labels)
- Parse section markers to identify complex/receptor/ligand/delta
- Use diverging colormap (RdBu_r) for decomposition

---

## 8. Testing and Validation

### 8.1 Unit Tests

Create `tests/test_modules.py`:

```python
import pytest
import numpy as np
import pandas as pd
from pathlib import Path

# Test load_align
def test_load_universe():
    """Test universe loading with atom count mismatch handling."""
    from modules.load_align import load_universe
    u = load_universe('com.tpr', 'com_traj.xtc')
    assert u is not None

# Test RMSD
def test_rmsd_calculation():
    """Test RMSD calculation returns expected columns."""
    from modules.rmsd import calculate_receptor_rmsd
    # ... tests

# Test contacts
def test_contacts_csv_format():
    """Test contacts CSV has correct format."""
    # ... tests
```

### 8.2 Integration Test

Run full analysis pipeline:

```bash
python analyze_trajectory.py \
    --tpr com.tpr \
    --xtc com_traj.xtc \
    --ndx index.ndx \
    --output output \
    --aim all
```

Expected outputs:
- `output/rmsd/rmsd_receptor.csv`
- `output/rmsd/rmsd_ligand.csv`
- `output/rmsd/rmsd_plot.png`
- `output/contacts/contacts_timeseries.csv`
- `output/contacts/contacts_heatmap.png`
- `output/hydrogen_bonds/hbond_timeseries.csv`
- `output/hydrogen_bonds/hbond_plot.png`
- `output/mmpbsa/energy_timecourse.png`
- `output/mmpbsa/decomposition_bar.png`

---

## 9. Dependencies

### Required Libraries

```bash
pip install "mdanalysis>=2.10.0" \
            numpy \
            pandas \
            matplotlib
```

### Version Requirements

| Library | Minimum Version | Purpose |
|---------|-----------------|---------|
| MDAnalysis | 2.10.0 | Trajectory I/O, RMSD, contacts, hbonds |
| NumPy | 1.24.0 | Numerical operations |
| Pandas | 2.0.0 | CSV handling, dataframes |
| Matplotlib | 3.7.0 | Plotting |

---

## 10. Summary

This implementation plan provides a comprehensive Python script for protein-ligand trajectory analysis. The script:

1. **Loads and aligns** GROMACS trajectories using MDAnalysis with proper atom selection
2. **Calculates RMSD** for receptor backbone and ligand heavy atoms with overlaid plots
3. **Computes per-residue contacts** using distance-based cutoff (4.5Å)
4. **Analyzes hydrogen bonds** using geometric criteria (3.5Å, 150°)
5. **Generates MMPBSA plots** from CSV output files without GUI

All analyses are configurable via command-line arguments, use vectorized approaches for efficiency, and produce publication-quality outputs.

---

**Next Steps:**
1. Create module files as specified in the plan
2. Test each module independently
3. Run full pipeline with test data
4. Validate outputs against expected results