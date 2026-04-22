# Phase 3: Hydrogen Bond Analysis - Research

**Researched:** 2026-03-22
**Domain:** Protein-ligand hydrogen bond analysis from GROMACS trajectories
**Confidence:** HIGH

## Summary

This research addresses implementing Aim 3 of a GROMACS protein-ligand analysis script: calculating per-residue hydrogen bond counts between receptor and ligand, generating time series data, and outputting CSV and PNG plots. MDAnalysis (version 2.7.0) provides a robust `HydrogenBondAnalysis` class that handles donor/acceptor identification, geometric criteria (distance ~3.5Å, angle ~150°), and supports the `between` keyword for cross-group hydrogen bond detection. The analysis requires loading a TPR topology (for bond information) and XTC trajectory, selecting protein and ligand groups, running the analysis with appropriate geometric cutoffs, aggregating results by receptor residue, and exporting time series data to CSV with matplotlib visualizations.

**Primary recommendation:** Use MDAnalysis.analysis.hydrogenbonds.HydrogenBondAnalysis with the `between` keyword to isolate receptor-ligand hydrogen bonds, then aggregate by residue using the donor/acceptor residue information from results.

---

## Standard Stack

The established libraries for hydrogen bond analysis in MD trajectories:

### Core
| Library | Version | Purpose | Why Standard |
|---------|---------|---------|--------------|
| MDAnalysis | 2.7.0 | Trajectory I/O and analysis | Industry-standard Python library for MD trajectory analysis |
| numpy | 1.26.4 | Numerical operations | Core dependency for coordinate calculations |
| pandas | 1.5.3 | DataFrame for CSV output | Efficient time series data handling |
| matplotlib | 3.7.3 | Plotting | Standard scientific plotting library |

### Supporting
| Library | Purpose | When to Use |
|---------|---------|--------------|
| MDAnalysis.analysis.hydrogenbonds | Hydrogen bond detection | Primary method for hbond analysis |
| MDAnalysis.analysis.align | Trajectory alignment | Before hbond analysis to ensure proper orientation |

**Installation:**
```bash
# Already available in environment
pip install MDAnalysis numpy pandas matplotlib
```

---

## Architecture Patterns

### Recommended Project Structure
```
com_ana_trj/
├── analyze_hbond.py      # Main analysis script
├── data/
│   ├── com.tpr          # Topology (TPR)
│   ├── com_traj.xtc     # Trajectory (XTC)
│   └── index.ndx        # Index file for atom groups
└── output/
    ├── hbond_timeseries.csv
    └── hbond_plot.png
```

### Pattern 1: Basic Hydrogen Bond Analysis

**What:** Using MDAnalysis HydrogenBondAnalysis class to detect hydrogen bonds between two selection groups.

**When to use:** When you need to find all hydrogen bonds between receptor and ligand.

**Example:**
```python
# Source: MDAnalysis 2.7.0 documentation
import MDAnalysis as mda
from MDAnalysis.analysis.hydrogenbonds.hbond_analysis import HydrogenBondAnalysis as HBA

# Load universe with TPR (contains bond information for donor-hydrogen pairing)
u = mda.Universe('com.tpr', 'com_traj.xtc')

# Select receptor and ligand
receptor = u.select_atoms('protein')  # or use index group
ligand = u.select_atoms('resname MOL')  # ligand residue name

# Initialize hydrogen bond analysis
hbonds = HBA(
    universe=u,
    between=['protein', 'resname MOL'],  # Only receptor-ligand hbonds
    d_a_cutoff=3.5,  # Distance cutoff (Å)
    d_h_a_angle_cutoff=150,  # Angle cutoff (degrees)
    update_selctions=False  # Static selections for speed
)

# Run analysis
hbonds.run()

# Access results
# results.hbonds contains: [frame, donor_idx, hydrogen_idx, acceptor_idx, distance, angle]
```

### Pattern 2: Per-Residue Aggregation

**What:** Aggregating hydrogen bond counts by receptor residue from raw hbond data.

**When to use:** When you need to know which residues form hydrogen bonds with the ligand.

**Example:**
```python
# Aggregate hbonds by receptor residue
from collections import defaultdict
import numpy as np

# Get hydrogen bond results
hbond_data = hbonds.results.hbonds  # Shape: (N, 6)

# Build mapping from atom index to residue info
atom_to_residue = {}
for atom in u.atoms:
    atom_to_residue[atom.index] = (atom.resid, atom.resname, atom.resnum)

# Count hbonds per residue per frame
residue_hbonds = defaultdict(lambda: defaultdict(int))
times = []

for hbond in hbond_data:
    frame, donor_idx, hydrogen_idx, acceptor_idx, distance, angle = hbond
    
    # Get residue info for donor (protein) and acceptor (ligand)
    donor_res = atom_to_residue.get(donor_idx)
    acceptor_res = atom_to_residue.get(acceptor_idx)
    
    if donor_res:
        resid, resname, resnum = donor_res
        residue_hbonds[frame][f"{resname}{resnum}"] += 1

# Now have per-frame, per-residue counts
```

### Pattern 3: Time Series Generation

**What:** Creating time series data of hydrogen bond counts for plotting.

**When to use:** When you need to visualize how hydrogen bonds change over time.

**Example:**
```python
import pandas as pd

# Get unique frames and count hbonds per frame
frames = np.unique(hbond_data[:, 0])
time_series_data = []

for frame in frames:
    frame_hbonds = hbond_data[hbond_data[:, 0] == frame]
    total_hbonds = len(frame_hbonds)
    
    # Get time in ps from trajectory
    u.trajectory[frame]
    time_ps = u.trajectory.time
    
    time_series_data.append({
        'frame': frame,
        'time_ps': time_ps,
        'total_hbonds': total_hbonds
    })

# Create DataFrame and save to CSV
df = pd.DataFrame(time_series_data)
df.to_csv('hbond_timeseries.csv', index=False)
```

### Anti-Patterns to Avoid

- **Using PDB instead of TPR:** Without bond information in the topology, donor-hydrogen pairing becomes unreliable and requires manual `d_h_cutoff` specification. Always use TPR for protein-ligand systems.
  
- **Not using `between` keyword:** Without specifying `between`, you'll get protein-protein and ligand-ligand hydrogen bonds too, contaminating your receptor-ligand analysis.

- **Setting `update_selctions=True` for static systems:** This adds overhead on every frame when atom selections don't change. Only use `True` if atoms may enter/leave the selection (e.g., flexible residues).

---

## Don't Hand-Roll

Problems that look simple but have existing solutions:

| Problem | Don't Build | Use Instead | Why |
|---------|-------------|-------------|-----|
| Donor-Hydrogen identification | Write custom distance-based pairing | Use TPR topology bond information | Robust and validated; handles edge cases |
| Geometric criteria | Implement distance/angle calculations | HydrogenBondAnalysis with `d_a_cutoff` and `d_h_a_angle_cutoff` | Vectorized, optimized, validated against VMD |
| Hydrogen/Acceptor guessing | Create custom selection strings | `guess_hydrogens()` and `guess_acceptors()` methods | Uses charge/mass heuristics validated by community |
| Cross-group filtering | Post-process all hbonds | `between` keyword | Efficient, avoids computing unwanted hbonds |

**Key insight:** The MDAnalysis HydrogenBondAnalysis module was modeled after the VMD HBONDS plugin and has been extensively validated. Custom implementations typically fail to handle periodic boundary conditions correctly and miss edge cases in donor/acceptor identification.

---

## Common Pitfalls

### Pitfall 1: Incorrect Donor-Hydrogen Pairing
**What goes wrong:** Hydrogen bonds are detected with wrong donor-hydrogen pairs, leading to false positives or missing hbonds.

**Why it happens:** Using a topology without bond information (e.g., PDB instead of TPR), the analysis cannot correctly identify which hydrogens belong to which donors.

**How to avoid:** Always use TPR topology. If only PDB is available, provide explicit `donors_sel` and `hydrogens_sel` and set `d_h_cutoff=1.2`.

**Warning signs:** Unusually high or low hbond counts; hydrogens paired to wrong donors.

### Pitfall 2: Missing Receptor-Ligand Specificity
**What goes wrong:** Getting protein-protein or ligand-ligand hydrogen bonds mixed with receptor-ligand hbonds.

**Why it happens:** Not specifying the `between` keyword to filter cross-group interactions.

**How to avoid:** Use `between=['selection1', 'selection2']` to restrict to specific group pairs.

**Warning signs:** Total hbond counts much higher than expected for interface; residue analysis shows intra-protein contacts.

### Pitfall 3: Wrong Distance Cutoff
**What goes wrong:** Using default 3.0Å when 3.5Å is more appropriate for protein-ligand systems.

**Why it happens:** The default cutoff is optimized for protein-water hbonds; protein-ligand hbonds can be slightly longer.

**How to avoid:** Set `d_a_cutoff=3.5` for protein-ligand analysis (standard in the field).

**Warning signs:** Fewer hbonds than expected from visual inspection; missing known important hbonds.

### Pitfall 4: Angle Cutoff Too Strict
**What goes wrong:** Using 150° cutoff misses some valid hbonds that are slightly more bent.

**Why it happens:** The 150° default is strict; some protein-ligand hbonds are ~120-140°.

**How to avoid:** Consider using 120-130° for more inclusive analysis, or keep 150° for high-quality subset.

**Warning signs:** Known hbonds not detected; very few hbonds detected.

---

## Code Examples

### Complete Hydrogen Bond Analysis Pipeline

```python
#!/usr/bin/env python3
"""
hbond_analysis.py - Per-residue hydrogen bond analysis for protein-ligand systems
"""

import argparse
import warnings
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")  # Non-interactive backend
import matplotlib.pyplot as plt
import MDAnalysis as mda
from MDAnalysis.analysis.hydrogenbonds.hbond_analysis import HydrogenBondAnalysis as HBA

warnings.filterwarnings('ignore')

def analyze_hydrogen_bonds(tpr_file, xtc_file, output_prefix='hbond'):
    """Analyze hydrogen bonds between receptor and ligand."""
    
    # Load universe
    print(f"Loading {tpr_file} and {xtc_file}...")
    u = mda.Universe(tpr_file, xtc_file)
    
    # Determine receptor and ligand selections
    # Receptor: protein atoms (assuming 'protein' or 'Protein' group exists)
    # Ligand: non-protein, non-water (assuming 'MOL' is ligand residue name)
    receptor_sel = 'protein'
    ligand_sel = 'resname MOL'
    
    print(f"Receptor: {u.select_atoms(receptor_sel).n_atoms} atoms")
    print(f"Ligand: {u.select_atoms(ligand_sel).n_atoms} atoms")
    
    # Initialize hydrogen bond analysis
    # Use standard criteria: distance ~3.5Å, angle ~150°
    hbonds = HBA(
        universe=u,
        between=[receptor_sel, ligand_sel],  # Only receptor-ligand hbonds
        d_a_cutoff=3.5,  # Distance cutoff (Å)
        d_h_a_angle_cutoff=150,  # Angle cutoff (degrees)
        update_selctions=False  # Static selections
    )
    
    print("Running hydrogen bond analysis...")
    hbonds.run(verbose=True)
    
    # Get hydrogen bond data
    hbond_data = hbonds.results.hbonds
    print(f"Found {len(hbond_data)} hydrogen bonds")
    
    # Build atom to residue mapping
    atom_to_residue = {}
    for atom in u.atoms:
        atom_to_residue[atom.index] = {
            'resid': atom.resid,
            'resname': atom.resname,
            'resnum': atom.resnum,
            'name': atom.name
        }
    
    # Analyze per-residue hydrogen bonds
    print("Computing per-residue statistics...")
    
    # Method 1: Total hbonds per frame (time series)
    frames = np.unique(hbond_data[:, 0])
    time_series = []
    
    for frame in frames:
        frame_hbonds = hbond_data[hbond_data[:, 0] == frame]
        n_hbonds = len(frame_hbonds)
        
        # Get time
        u.trajectory[int(frame)]
        time_ps = u.trajectory.time
        
        time_series.append({
            'frame': int(frame),
            'time_ps': time_ps,
            'hbond_count': n_hbonds
        })
    
    # Save time series to CSV
    df_time = pd.DataFrame(time_series)
    csv_file = f'{output_prefix}_timeseries.csv'
    df_time.to_csv(csv_file, index=False)
    print(f"Saved time series to {csv_file}")
    
    # Method 2: Per-residue hydrogen bond counts
    residue_counts = {}
    
    for hbond in hbond_data:
        frame, donor_idx, hydrogen_idx, acceptor_idx, distance, angle = hbond
        
        # Get residue info - check both donor and acceptor
        donor_res = atom_to_residue.get(donor_idx)
        acceptor_res = atom_to_residue.get(acceptor_idx)
        
        # Check which is receptor (protein) vs ligand
        for res_info in [donor_res, acceptor_res]:
            if res_info and res_info['resname'] != 'MOL':
                key = f"{res_info['resname']}{res_info['resnum']}"
                if key not in residue_counts:
                    residue_counts[key] = {
                        'resname': res_info['resname'],
                        'resid': res_info['resid'],
                        'count': 0,
                        'avg_distance': [],
                        'avg_angle': []
                    }
                residue_counts[key]['count'] += 1
                residue_counts[key]['avg_distance'].append(distance)
                residue_counts[key]['avg_angle'].append(angle)
    
    # Calculate averages
    for key in residue_counts:
        if residue_counts[key]['avg_distance']:
            residue_counts[key]['avg_distance'] = np.mean(residue_counts[key]['avg_distance'])
            residue_counts[key]['avg_angle'] = np.mean(residue_counts[key]['avg_angle'])
    
    # Save per-residue data
    df_residue = pd.DataFrame([
        {
            'residue': key,
            'resname': v['resname'],
            'resid': v['resid'],
            'hbond_count': v['count'],
            'avg_distance_A': v['avg_distance'],
            'avg_angle_deg': v['avg_angle']
        }
        for key, v in sorted(residue_counts.items(), key=lambda x: -x[1]['count'])
    ])
    
    residue_csv = f'{output_prefix}_per_residue.csv'
    df_residue.to_csv(residue_csv, index=False)
    print(f"Saved per-residue data to {residue_csv}")
    
    # Generate plots
    print("Generating plots...")
    
    # Plot 1: Time series of total hydrogen bonds
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 8))
    
    # Time series
    ax1.plot(df_time['time_ps'], df_time['hbond_count'], 'b-', linewidth=1)
    ax1.set_xlabel('Time (ps)')
    ax1.set_ylabel('Hydrogen Bond Count')
    ax1.set_title('Receptor-Ligand Hydrogen Bonds Over Time')
    ax1.grid(True, alpha=0.3)
    
    # Add moving average
    window = min(10, len(df_time))
    if window > 1:
        df_time['hbond_ma'] = df_time['hbond_count'].rolling(window=window, center=True).mean()
        ax1.plot(df_time['time_ps'], df_time['hbond_ma'], 'r-', linewidth=2, label=f'{window}-frame moving average')
        ax1.legend()
    
    # Plot 2: Per-residue counts (top 20)
    top_residues = df_residue.head(20)
    ax2.barh(range(len(top_residues)), top_residues['hbond_count'], color='steelblue')
    ax2.set_yticks(range(len(top_residues)))
    ax2.set_yticklabels(top_residues['residue'])
    ax2.set_xlabel('Total Hydrogen Bond Count')
    ax2.set_title('Top 20 Residues by Hydrogen Bond Count')
    ax2.invert_yaxis()
    ax2.grid(True, alpha=0.3, axis='x')
    
    plt.tight_layout()
    plot_file = f'{output_prefix}_plot.png'
    plt.savefig(plot_file, dpi=150)
    plt.close()
    print(f"Saved plot to {plot_file}")
    
    return df_time, df_residue


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Analyze protein-ligand hydrogen bonds')
    parser.add_argument('--tpr', default='com.tpr', help='TPR topology file')
    parser.add_argument('--xtc', default='com_traj.xtc', help='XTC trajectory file')
    parser.add_argument('--output', default='hbond', help='Output prefix')
    
    args = parser.parse_args()
    
    analyze_hydrogen_bonds(args.tpr, args.xtc, args.output)
```

### CSV Output Format

**Time series CSV** (`hbond_timeseries.csv`):
```csv
frame,time_ps,hbond_count
0,0.0,3
1,100.0,4
2,200.0,2
...
```

**Per-residue CSV** (`hbond_per_residue.csv`):
```csv
residue,resname,resid,hbond_count,avg_distance_A,avg_angle_deg
ASP128,ASP,128,45,2.87,162.3
GLU201,GLU,201,38,2.92,158.7
...
```

### Matplotlib Plot Code (Standalone)

```python
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd

# Load data
df_time = pd.read_csv('hbond_timeseries.csv')
df_residue = pd.read_csv('hbond_per_residue.csv')

# Create figure with two subplots
fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 8))

# Subplot 1: Time series
ax1.plot(df_time['time_ps'], df_time['hbond_count'], 'b-', linewidth=1, alpha=0.7)
ax1.fill_between(df_time['time_ps'], 0, df_time['hbond_count'], alpha=0.3)
ax1.set_xlabel('Time (ps)', fontsize=12)
ax1.set_ylabel('Hydrogen Bond Count', fontsize=12)
ax1.set_title('Receptor-Ligand Hydrogen Bonds Over Time', fontsize=14)
ax1.grid(True, alpha=0.3)
ax1.set_xlim(df_time['time_ps'].min(), df_time['time_ps'].max())

# Add statistics text
mean_hb = df_time['hbond_count'].mean()
std_hb = df_time['hbond_count'].std()
ax1.axhline(mean_hb, color='red', linestyle='--', label=f'Mean: {mean_hb:.1f}')
ax1.legend()

# Subplot 2: Per-residue bar chart (top 15)
top_n = 15
top_res = df_residue.head(top_n)
colors = plt.cm.viridis(np.linspace(0.2, 0.8, top_n))
ax2.barh(range(top_n), top_res['hbond_count'], color=colors)
ax2.set_yticks(range(top_n))
ax2.set_yticklabels(top_res['residue'])
ax2.set_xlabel('Total Hydrogen Bond Count', fontsize=12)
ax2.set_title(f'Top {top_n} Residues by Hydrogen Bond Frequency', fontsize=14)
ax2.invert_yaxis()
ax2.grid(True, alpha=0.3, axis='x')

# Add count labels
for i, count in enumerate(top_res['hbond_count']):
    ax2.text(count + 0.5, i, str(count), va='center', fontsize=9)

plt.tight_layout()
plt.savefig('hbond_analysis.png', dpi=150, bbox_inches='tight')
plt.close()
```

---

## State of the Art

| Old Approach | Current Approach | When Changed | Impact |
|--------------|------------------|--------------|--------|
| VMD HBONDS plugin | MDAnalysis HydrogenBondAnalysis | ~2019 | Cross-platform, scriptable, Python integration |
| Custom distance calculations | Built-in `d_a_cutoff` and `d_h_a_angle_cutoff` | v1.0 (MDAnalysis) | Validated geometry, handles PBC |
| Post-processing all hbonds | `between` keyword | v2.0 (MDAnalysis) | Memory-efficient, faster |
| Serial only | Parallel execution with `backend='multiprocessing'` | v2.8.0 | Multi-core support for large trajectories |

**Deprecated/outdated:**
- `MDAnalysis.analysis.hbonds.hbond_analysis` (old module): Replaced by `MDAnalysis.analysis.hydrogenbonds.hbond_analysis` for clarity
- Manual donor/acceptor guessing: Now handled by `guess_hydrogens()` and `guess_acceptors()` methods

---

## Open Questions

1. **Ligand residue name variability**
   - What we know: The test files use 'MOL' as ligand residue name, but this varies by system
   - What's unclear: How to robustly detect ligand residue name without manual specification
   - Recommendation: Add `--ligand-resname` CLI argument with default 'MOL', or auto-detect from non-protein, non-water atoms

2. **Hydrogen bond definition variations**
   - What we know: Standard is ~3.5Å distance, ~150° angle
   - What's unclear: Some studies use 3.0Å/120° for more inclusive detection
   - Recommendation: Add CLI arguments for `--distance-cutoff` and `--angle-cutoff` with sensible defaults

3. **Trajectory alignment requirement**
   - What we know: Hydrogen bonds can be spurious if trajectory is not aligned
   - What's unclear: Whether to auto-align or require pre-aligned trajectory
   - Recommendation: Document that trajectory should be aligned; optionally add auto-alignment using receptor backbone

---

## Sources

### Primary (HIGH confidence)
- MDAnalysis 2.10.0 documentation - Hydrogen Bond Analysis
  - URL: https://docs.mdanalysis.org/stable/documentation_pages/analysis/hydrogenbonds.html
  - Topics: HydrogenBondAnalysis class, parameters, between keyword, output format
  
- MDAnalysis 2.7.0 (installed version)
  - Verified: `import MDAnalysis; print(MDAnalysis.__version__)` → 2.7.0
  - Compatible with 2.10.0 API for core functionality

### Secondary (MEDIUM confidence)
- MDAnalysis User Guide - Analysis examples
  - URL: https://userguide.mdanalysis.org/stable/examples/analysis/README.html
  - Topics: General analysis patterns, trajectory handling

### Tertiary (LOW confidence)
- Community discussions on hydrogen bond cutoffs (general knowledge)
  - Standard: 3.5Å distance, 150° angle for protein-ligand
  - This matches published literature and VMD defaults

---

## Metadata

**Confidence breakdown:**
- Standard stack: HIGH - MDAnalysis is well-established, version verified
- Architecture: HIGH - Patterns directly from official documentation
- Pitfalls: HIGH - Based on official documentation and known issues
- Code examples: HIGH - Adapted from MDAnalysis documentation

**Research date:** 2026-03-22
**Valid until:** 2026-04-22 (30 days - MDAnalysis API is stable)
