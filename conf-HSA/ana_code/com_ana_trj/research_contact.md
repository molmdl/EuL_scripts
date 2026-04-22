# Phase 2: Per-Residue Contact Analysis - Research

**Researched:** 2026-03-22
**Domain:** Protein-ligand contact analysis using MDAnalysis
**Confidence:** HIGH

## Summary

This research covers implementing Aim 2 of a GROMACS protein-ligand analysis script: calculating distance-based per-residue contact numbers between receptor and ligand, generating time series data, and outputting CSV and PNG plots.

Based on analysis of test files in `~/dparker/dp_xinyi/ana_code/com_ana_trj/`:
- **System:** 172,475 atoms total (protein receptor + ligand + solvent + ions)
- **Receptor:** [Protein] group - atoms 1 to ~172,000
- **Ligand:** [MOL] group - atoms 9122-9154 (33 atoms - small molecule)
- **Available selections:** Protein, C-alpha, Backbone, MOL, Water, ions

**Primary recommendation:** Use MDAnalysis `distance_array()` with `contact_matrix()` for per-residue contact calculation. Standard cutoff is 4.5Å for all-atom simulations. Use vectorized NumPy operations for efficiency.

---

## 1. Definition of Contact

### Distance Cutoff Standards

| Simulation Type | Cutoff Distance | Notes |
|-----------------|-----------------|-------|
| All-atom | **4.5 Å** | Standard for protein-ligand contacts |
| Coarse-grained | 6.0 Å | MARTINI-style simulations |
| Strict contacts | 3.5-4.0 Å | For specific interactions (H-bonds, salt bridges) |
| Extended range | 5.0-6.0 Å | For aromatic stacking, weak interactions |

**Source:** MDAnalysis documentation recommends 4.5Å for all-atom and 6.0Å for coarse-grained simulations (see [contacts_within_cutoff.html](https://userguide.mdanalysis.org/stable/examples/analysis/distances_and_contacts/contacts_within_cutoff.html))

### Contact Definition

A **contact** is defined when any atom from a receptor residue is within the cutoff distance of any atom from the ligand. For per-residue analysis:
- Compute all pairwise distances between receptor residue atoms and ligand atoms
- If minimum distance ≤ cutoff → residue is in contact
- Count per-residue contacts over time

---

## 2. Code Snippets for Per-Residue Contact Calculation

### 2.1 Basic Setup - Loading Trajectory

```python
import MDAnalysis as mda
import numpy as np
import pandas as pd

# Load trajectory with GROMACS files
u = mda.Universe('com.tpr', 'com_traj.xtc')

# Select receptor (protein) and ligand groups
# Method 1: Using index groups
receptor = u.select_atoms('[ Protein ]')      # All protein atoms
ligand = u.select_atoms('[ MOL ]')            # Ligand molecule

# Method 2: Alternative selections
# CA atoms only (more efficient for large proteins)
receptor_ca = u.select_atoms('name CA and protein')
# Ligand heavy atoms
ligand_heavy = u.select_atoms('resname MOL and not name H*')
```

### 2.2 Per-Residue Contact Calculation

```python
# Source: Adapted from MDAnalysis distance_array documentation
# https://userguide.mdanalysis.org/stable/examples/analysis/distances_and_contacts/distances_between_selections.html

def calculate_per_residue_contacts(u, receptor_sel, ligand_sel, cutoff=4.5):
    """
    Calculate per-residue contact numbers between receptor and ligand.
    
    Parameters:
    -----------
    u : Universe
        MDAnalysis Universe with trajectory
    receptor_sel : str
        Selection string for receptor residues
    ligand_sel : str
        Selection string for ligand atoms
    cutoff : float
        Distance cutoff in Angstroms (default 4.5)
    
    Returns:
    --------
    dict : {residue_id: [contacts_frame0, contacts_frame1, ...]}
    """
    # Get receptor residue list
    receptor = u.select_atoms(receptor_sel)
    residues = receptor.residues
    
    # Get ligand atoms
    ligand = u.select_atoms(ligand_sel)
    
    # Pre-allocate results dictionary
    contact_timeseries = {res.resid: [] for res in residues}
    times = []
    
    # Iterate through trajectory
    for ts in u.trajectory:
        times.append(ts.time)  # or ts.frame
        
        # Calculate distance array: receptor atoms vs ligand atoms
        # Shape: (n_receptor_atoms, n_ligand_atoms)
        dist_matrix = mda.analysis.distances.distance_array(
            receptor.positions,
            ligand.positions,
            box=u.dimensions
        )
        
        # Reshape to group by residue
        # For each residue, find if ANY atom is within cutoff of ANY ligand atom
        n_residues = len(residues)
        
        # Get atom indices for each residue
        res_atom_indices = []
        for res in residues:
            res_atoms = res.atoms
            res_atom_indices.append(res_atoms.indices)
        
        # Calculate contacts per residue
        for res_idx, atom_indices in enumerate(res_atom_indices):
            # Get distances from this residue's atoms to all ligand atoms
            res_distances = dist_matrix[atom_indices, :]
            
            # Check if ANY atom is within cutoff
            min_distance = res_distances.min()
            n_contacts = (res_distances <= cutoff).any(axis=0).sum() if res_distances.ndim > 1 else int(res_distances <= cutoff)
            
            contact_timeseries[residues[res_idx].resid].append(n_contacts)
    
    return contact_timeseries, times
```

### 2.3 Optimized Vectorized Version

```python
# Source: Adapted from MDAnalysis contacts module
# https://docs.mdanalysis.org/stable/documentation_pages/analysis/contacts.html

def calculate_contacts_vectorized(u, receptor_sel, ligand_sel, cutoff=4.5):
    """
    Vectorized per-residue contact calculation using NumPy broadcasting.
    Much faster than iterating over residues.
    """
    from MDAnalysis.analysis import distances
    
    # Select atoms
    receptor = u.select_atoms(receptor_sel)
    ligand = u.select_atoms(ligand_sel)
    residues = receptor.residues
    n_residues = len(residues)
    n_ligand_atoms = len(ligand)
    
    # Get residue atom mappings
    res_start_idx = []
    res_end_idx = []
    for res in residues:
        res_start_idx.append(res.atoms[0].ix)
        res_end_idx.append(res.atoms[-1].ix)
    
    # Pre-allocate
    contact_matrix = np.zeros((n_residues, n_ligand_atoms), dtype=bool)
    timeseries = {res.resid: [] for res in residues}
    times = []
    
    for ts in u.trajectory:
        times.append(ts.time)
        
        # Calculate full distance matrix
        dist_mat = distances.distance_array(
            receptor.positions,
            ligand.positions,
            box=u.dimensions
        )
        
        # For each residue, check if ANY atom is within cutoff of ANY ligand atom
        # Vectorized approach
        for i, (start, end) in enumerate(zip(res_start_idx, res_end_idx)):
            # Get distances for atoms in this residue
            res_dists = dist_mat[start:end+1, :]
            
            # Contact if ANY atom of residue is within cutoff of ANY ligand atom
            is_contact = (res_dists <= cutoff).any(axis=0).any()
            n_contacts = (res_dists <= cutoff).sum()
            
            timeseries[residues[i].resid].append(n_contacts)
    
    return timeseries, times
```

### 2.4 Using MDAnalysis Contacts Class (Alternative)

```python
# Source: https://userguide.mdanalysis.org/stable/examples/analysis/distances_and_contacts/contacts_within_cutoff.html

from MDAnalysis.analysis import contacts

# Define groups - using CA for efficiency
sel_receptor = "name CA and protein"
sel_ligand = "resname MOL"

receptor_ca = u.select_atoms(sel_receptor)
ligand_mol = u.select_atoms(sel_ligand)

# Use contacts.Contacts class with radius_cut method
ca = contacts.Contacts(u,
                       select=(sel_receptor, sel_ligand),
                       refgroup=(receptor_ca, ligand_mol),
                       radius=4.5,
                       method='radius_cut').run()

# Results in ca.results.timeseries
# First column: frame/time, Second column: number of contacts
```

---

## 3. Identifying Receptor Residues in Contact with Ligand

### 3.1 Method 1: Selection-Based Identification

```python
# Get list of receptor residues in contact at each frame
def get_contacting_residues(u, receptor_sel, ligand_sel, cutoff=4.5):
    """
    Return list of residue IDs that are in contact with ligand.
    """
    receptor = u.select_atoms(receptor_sel)
    ligand = u.select_atoms(ligand_sel)
    residues = receptor.residues
    
    contacting_resids = []
    
    for ts in u.trajectory:
        dist_mat = mda.analysis.distances.distance_array(
            receptor.positions,
            ligand.positions,
            box=u.dimensions
        )
        
        # Find residues with at least one atom within cutoff
        current_contacts = set()
        for res in residues:
            res_atoms = res.atoms
            atom_indices = res_atoms.indices
            
            # Get minimum distance for this residue
            min_dist = dist_mat[atom_indices, :].min()
            if min_dist <= cutoff:
                current_contacts.add(res.resid)
        
        contacting_resids.append(current_contacts)
    
    return contacting_resids
```

### 3.2 Method 2: Using contact_matrix

```python
# Source: https://docs.mdanalysis.org/stable/documentation_pages/analysis/distances.html

from MDAnalysis.analysis import distances

def get_contacts_with_matrix(u, receptor_sel, ligand_sel, cutoff=4.5):
    """
    Use contact_matrix for boolean contact determination.
    """
    receptor = u.select_atoms(receptor_sel)
    ligand = u.select_atoms(ligand_sel)
    
    contacts_per_frame = []
    
    for ts in u.trajectory:
        # Distance array
        dist_mat = distances.distance_array(
            receptor.positions,
            ligand.positions,
            box=u.dimensions
        )
        
        # Boolean contact matrix
        contact_mat = distances.contact_matrix(dist_mat, cutoff=cutoff)
        
        # Sum contacts
        n_contacts = contact_mat.sum()
        contacts_per_frame.append(n_contacts)
    
    return contacts_per_frame
```

### 3.3 With Residue Information

```python
# Get detailed per-residue contact info with residue names
def detailed_contact_analysis(u, receptor_sel, ligand_sel, cutoff=4.5):
    """
    Return DataFrame with per-residue contact details.
    """
    receptor = u.select_atoms(receptor_sel)
    ligand = u.select_atoms(ligand_sel)
    residues = receptor.residues
    
    results = []
    
    for ts in u.trajectory:
        frame = ts.frame
        time = ts.time
        
        dist_mat = mda.analysis.distances.distance_array(
            receptor.positions,
            ligand.positions,
            box=u.dimensions
        )
        
        for res in residues:
            res_atoms = res.atoms
            atom_indices = res_atoms.indices
            
            res_dists = dist_mat[atom_indices, :]
            min_dist = res_dists.min()
            n_within_cutoff = (res_dists <= cutoff).sum()
            
            results.append({
                'frame': frame,
                'time_ps': time,
                'residue_name': res.resname,
                'residue_id': res.resid,
                'min_distance': min_dist,
                'n_atoms_within_cutoff': n_within_cutoff,
                'in_contact': min_dist <= cutoff
            })
    
    return pd.DataFrame(results)
```

---

## 4. CSV Output Format

### 4.1 Time Series Per Residue (Wide Format)

```python
# Source: Adapted from MDAnalysis contacts example
# https://userguide.mdanalysis.org/stable/examples/analysis/distances_and_contacts/contacts_within_cutoff.html

def save_contacts_csv(contact_timeseries, times, output_file):
    """
    Save contact time series to CSV.
    
    Format: wide format with residue IDs as columns
    """
    df = pd.DataFrame(contact_timeseries, index=times)
    df.index.name = 'time_ps'
    df.to_csv(output_file)
    return df

# Example output:
# time_ps,    100,    101,    102,    103,    104, ...
# 0.0,         3,      0,      2,      0,      1, ...
# 100.0,       2,      1,      3,      0,      0, ...
# 200.0,       1,      0,      2,      1,      1, ...
```

### 4.2 Detailed Long Format

```python
def save_detailed_csv(df, output_file):
    """
    Save detailed contact info in long format.
    """
    df.to_csv(output_file, index=False)
    
# Example output:
# frame,time_ps,residue_name,residue_id,min_distance,n_atoms_within_cutoff,in_contact
# 0,0.0,LYS,100,3.2,4,True
# 0,0.0,GLU,101,6.1,0,False
# 0,0.0,ASP,102,4.0,2,True
```

### 4.3 Summary Statistics Format

```python
def save_summary_csv(contact_timeseries, output_file):
    """
    Save summary statistics per residue.
    """
    df = pd.DataFrame(contact_timeseries)
    
    summary = pd.DataFrame({
        'residue_id': df.columns,
        'mean_contacts': df.mean(),
        'std_contacts': df.std(),
        'max_contacts': df.max(),
        'min_contacts': df.min(),
        'fraction_in_contact': (df > 0).mean()
    })
    
    summary.to_csv(output_file, index=False)
    return summary
```

### 4.4 Recommended CSV Structure

**Primary output (time series):**
```
time_ps,res_1,res_2,res_3,...
0.0,5,2,0,...
100.0,4,3,1,...
```

**Metadata header (optional):**
```
# Contact Analysis Results
# Cutoff: 4.5 Angstrom
# Receptor: Protein CA atoms
# Ligand: MOL
# Frames: 100
```

---

## 5. Matplotlib Code for PNG Plots

### 5.1 Time Series Plot

```python
import matplotlib.pyplot as plt
import matplotlib
matplotlib.use('Agg')  # Non-interactive backend

def plot_contact_timeseries(contact_timeseries, times, output_file):
    """
    Plot contact number time series per residue.
    """
    df = pd.DataFrame(contact_timeseries, index=times)
    
    fig, ax = plt.subplots(figsize=(12, 6))
    
    # Plot each residue
    for col in df.columns:
        ax.plot(df.index, df[col], label=col, alpha=0.7)
    
    ax.set_xlabel('Time (ps)', fontsize=12)
    ax.set_ylabel('Number of Contacts', fontsize=12)
    ax.set_title('Per-Residue Contact Number vs Time', fontsize=14)
    ax.legend(bbox_to_anchor=(1.05, 1), loc='upper left', fontsize=8)
    ax.grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    plt.close()
```

### 5.2 Heatmap of Contacts Over Time

```python
def plot_contact_heatmap(contact_timeseries, times, output_file):
    """
    Plot heatmap of contacts over time.
    """
    df = pd.DataFrame(contact_timeseries, index=times)
    
    fig, ax = plt.subplots(figsize=(14, 8))
    
    im = ax.imshow(df.values.T, aspect='auto', origin='lower',
                   cmap='YlOrRd', interpolation='nearest')
    
    # Labels
    ax.set_xlabel('Time (ps)', fontsize=12)
    ax.set_ylabel('Residue ID', fontsize=12)
    ax.set_title('Residue-Ligand Contact Heatmap', fontsize=14)
    
    # Colorbar
    cbar = plt.colorbar(im, ax=ax)
    cbar.set_label('Number of Contacts', fontsize=12)
    
    # Set tick labels
    ax.set_xticks(np.linspace(0, len(times)-1, 10).astype(int))
    ax.set_xticklabels([f'{t:.0f}' for t in np.linspace(times[0], times[-1], 10)])
    
    plt.tight_layout()
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    plt.close()
```

### 5.3 Stacked Area Plot

```python
def plot_stacked_contacts(contact_timeseries, times, output_file):
    """
    Plot stacked area chart of total contacts.
    """
    df = pd.DataFrame(contact_timeseries, index=times)
    
    fig, ax = plt.subplots(figsize=(12, 6))
    
    ax.stackplot(df.index, df.values.T, labels=df.columns, alpha=0.8)
    
    ax.set_xlabel('Time (ps)', fontsize=12)
    ax.set_ylabel('Number of Contacts', fontsize=12)
    ax.set_title('Stacked Per-Residue Contacts vs Time', fontsize=14)
    ax.legend(bbox_to_anchor=(1.05, 1), loc='upper left', fontsize=8)
    ax.grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    plt.close()
```

### 5.4 Bar Plot of Average Contacts

```python
def plot_average_contacts(contact_timeseries, output_file):
    """
    Plot bar chart of average contacts per residue.
    """
    df = pd.DataFrame(contact_timeseries)
    
    avg_contacts = df.mean().sort_values(ascending=False)
    
    fig, ax = plt.subplots(figsize=(12, 6))
    
    ax.bar(avg_contacts.index, avg_contacts.values, color='steelblue')
    ax.set_xlabel('Residue ID', fontsize=12)
    ax.set_ylabel('Average Contacts', fontsize=12)
    ax.set_title('Average Contact Number per Residue', fontsize=14)
    ax.grid(True, alpha=0.3, axis='y')
    
    # Rotate x labels if many residues
    plt.xticks(rotation=45, ha='right')
    
    plt.tight_layout()
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    plt.close()
```

### 5.5 Combined Multi-Panel Figure

```python
def plot_combined_analysis(contact_timeseries, times, output_file):
    """
    Create multi-panel figure with summary visualizations.
    """
    df = pd.DataFrame(contact_timeseries, index=times)
    
    fig = plt.figure(figsize=(16, 10))
    
    # Panel 1: Total contacts over time
    ax1 = fig.add_subplot(2, 2, 1)
    total_contacts = df.sum(axis=1)
    ax1.plot(total_contacts.index, total_contacts.values, 'b-', linewidth=2)
    ax1.set_xlabel('Time (ps)')
    ax1.set_ylabel('Total Contacts')
    ax1.set_title('Total Receptor-Ligand Contacts')
    ax1.grid(True, alpha=0.3)
    
    # Panel 2: Heatmap
    ax2 = fig.add_subplot(2, 2, 2)
    im = ax2.imshow(df.values.T[:50], aspect='auto', origin='lower',
                    cmap='YlOrRd', interpolation='nearest')
    ax2.set_xlabel('Frame')
    ax2.set_ylabel('Residue')
    plt.colorbar(im, ax=ax2, label='Contacts')
    
    # Panel 3: Average contacts bar
    ax3 = fig.add_subplot(2, 2, 3)
    avg = df.mean().sort_values(ascending=False)[:20]
    ax3.bar(range(len(avg)), avg.values, color='steelblue')
    ax3.set_xticks(range(len(avg)))
    ax3.set_xticklabels(avg.index, rotation=45, ha='right')
    ax3.set_xlabel('Residue')
    ax3.set_ylabel('Average Contacts')
    ax3.set_title('Top 20 Residues by Average Contacts')
    
    # Panel 4: Contact frequency
    ax4 = fig.add_subplot(2, 2, 4)
    contact_freq = (df > 0).mean().sort_values(ascending=False)[:20]
    ax4.bar(range(len(contact_freq)), contact_freq.values, color='coral')
    ax4.set_xticks(range(len(contact_freq)))
    ax4.set_xticklabels(contact_freq.index, rotation=45, ha='right')
    ax4.set_xlabel('Residue')
    ax4.set_ylabel('Fraction of Time in Contact')
    ax4.set_title('Contact Persistence (Top 20)')
    
    plt.tight_layout()
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    plt.close()
```

---

## 6. Vectorized/Efficient Approaches

### 6.1 Key Optimization Strategies

| Strategy | Benefit | Implementation |
|----------|---------|----------------|
| **Use CA atoms only** | Reduce computation ~10x | `select_atoms('name CA')` |
| **Pre-allocate arrays** | Avoid repeated memory allocation | `np.zeros()` before loop |
| **Use distance_array with box** | PBC-aware, fast C implementation | Built-in MDAnalysis |
| **Batch residue processing** | NumPy broadcasting | Reshape and compute in one pass |
| **Sparse contact matrix** | Memory efficient for large systems | `returntype='sparse'` |

### 6.2 Efficient Distance Calculation

```python
# Source: https://docs.mdanalysis.org/stable/documentation_pages/analysis/distances.html

# Pre-allocate result array for repeated calls (faster)
result = np.zeros((n_atoms_A, n_atoms_B), dtype=np.float64)

for ts in u.trajectory:
    distances.distance_array(group_A.positions, group_B.positions,
                            box=u.dimensions, result=result)
```

### 6.3 Parallel Execution

```python
# MDAnalysis 2.8.0+ supports parallel execution
# https://docs.mdanalysis.org/stable/documentation_pages/analysis/contacts.html

ca = contacts.Contacts(u,
                       select=(sel_receptor, sel_ligand),
                       refgroup=(receptor_ca, ligand_mol),
                       radius=4.5,
                       method='radius_cut')

# Run with multiprocessing
ca.run(backend='multiprocessing', n_jobs=4)

# Or with dask (requires mdanalysis[dask])
ca.run(backend='dask')
```

### 6.4 Using Numba for Custom Functions

```python
from numba import jit

@jit(nopython=True)
def fast_contact_check(distances, cutoff):
    """Numba-optimized contact checking."""
    n_atoms = distances.shape[0]
    n_ligand = distances.shape[1]
    
    for i in range(n_atoms):
        for j in range(n_ligand):
            if distances[i, j] <= cutoff:
                return True
    return False
```

### 6.5 Memory-Efficient Streaming

```python
# Process trajectory in chunks to avoid memory issues
def process_in_chunks(u, chunk_size=100):
    """Process trajectory in chunks for memory efficiency."""
    n_frames = len(u.trajectory)
    
    for start in range(0, n_frames, chunk_size):
        end = min(start + chunk_size, n_frames)
        
        # Process frames in chunk
        for ts in u.trajectory[start:end]:
            # Your analysis here
            pass
```

---

## Test File Analysis

Based on the files in `~/dparker/dp_xinyi/ana_code/com_ana_trj/`:

- **com.tpr:** GROMACS TPR file (binary topology) - 6.7 MB
- **com_traj.xtc:** GROMACS XTC trajectory - 730 KB (likely 10-50 frames)
- **index.ndx:** Index file with groups:
  - `[ Protein ]` - receptor atoms (1 to ~172,000)
  - `[ MOL ]` - ligand atoms (9122-9154, 33 atoms)
  - `[ C-alpha ]` - CA atoms for representative distance calculations
  - `[ Backbone ]`, `[ SideChain ]` - protein structure selections
  - `[ Water ]`, `[ Ion ]` - solvent and ions

**System size estimate:** ~172,000 atoms (large protein-ligand complex with water box)

**Recommended selections:**
```python
receptor_sel = "name CA and protein"  # CA atoms only for efficiency
ligand_sel = "resname MOL"             # Ligand residue
```

---

## Sources

### Primary (HIGH confidence)
- MDAnalysis User Guide - Distances and Contacts: https://userguide.mdanalysis.org/stable/examples/analysis/distances_and_contacts/
- MDAnalysis API Documentation - contacts module: https://docs.mdanalysis.org/stable/documentation_pages/analysis/contacts.html
- MDAnalysis API Documentation - distances module: https://docs.mdanalysis.org/stable/documentation_pages/analysis/distances.html

### Secondary (MEDIUM confidence)
- Contact analysis within cutoff tutorial: https://userguide.mdanalysis.org/stable/examples/analysis/distances_and_contacts/contacts_within_cutoff.html
- Native contacts fraction tutorial: https://userguide.mdanalysis.org/stable/examples/analysis/distances_and_contacts/contacts_native_fraction.html

### Tertiary (LOW confidence)
- Community discussions on MDAnalysis forums for edge cases
- Custom implementations may be needed for specific protein-ligand systems

---

## Metadata

**Confidence breakdown:**
- Standard stack: HIGH - MDAnalysis is the standard library for MD trajectory analysis
- Architecture: HIGH - distance_array + contact_matrix is well-documented pattern
- Pitfalls: MEDIUM - PBC handling and large system memory are known issues
- Code examples: HIGH - directly from MDAnalysis documentation with adaptations

**Research date:** 2026-03-22
**Valid until:** 2026-04-22 (MDAnalysis APIs are stable, but check for version-specific changes)