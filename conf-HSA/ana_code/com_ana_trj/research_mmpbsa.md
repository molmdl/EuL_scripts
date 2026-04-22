# Research: gmx_MMPBSA Plotting Without GUI

**Researched:** March 22, 2026
**Domain:** Scientific visualization for molecular dynamics binding free energy analysis
**Confidence:** HIGH

## Summary

This research addresses implementing Aim 4: generating plots from gmx_MMPBSA results without GUI interaction. The gmx_MMPBSA tool produces multiple output files (CSV and DAT formats) containing:

- **Energy decomposition results**: Per-frame and summary statistics for binding free energy components (van der Waals, electrostatics, polar solvation, etc.)
- **Per-residue decomposition**: Detailed residue-by-residue energy contributions to binding

The GUI tool `gmx_mmpbsa_ana` generates line plots, bar plots, heatmaps, and PyMOL visualizations. This research provides code patterns to reproduce these plots programmatically using Python, pandas, and matplotlib.

**Primary recommendation:** Parse CSV files for plotting (more reliable than DAT), use matplotlib with publication-quality settings, and create modular plotting functions for different chart types.

---

## 1. File Formats and Data Structure

### 1.1 FINAL_RESULTS_MMPBSA.csv

Per-frame energy data with three sections:

**Structure:**
```
Line 1: GENERALALIZED BORN:
Line 2: Complex Energy Terms
Line 3: Frame #,BOND,ANGLE,DIHED,VDWAALS,EEL,1-4 VDW,1-4 EEL,EGB,ESURF,GGAS,GSOLV,TOTAL
Line 4-24: Frame data (21 frames in example)

Line 25: [blank]
Line 26: Receptor Energy Terms
Line 27: Frame #,BOND,... (header)
Line 28-48: Receptor frame data

Line 49: [blank]
Line 50: Ligand Energy Terms
Line 51: Frame #,... (header)
Line 52-72: Ligand frame data

Line 73: [blank]
Line 74: Delta Energy Terms (Complex - Receptor - Ligand = Binding)
Line 75: Frame #,... (header)
Line 76-96: Delta frame data
```

**Energy Terms:**
| Term | Description |
|------|-------------|
| BOND | Bond energy |
| ANGLE | Angle energy |
| DIHED | Dihedral energy |
| VDWAALS | van der Waals energy |
| EEL | Electrostatic energy |
| 1-4 VDW | 1-4 van der Waals |
| 1-4 EEL | 1-4 electrostatics |
| EGB | Generalized Born polar solvation |
| ESURF | Non-polar solvation (surface area) |
| GGAS | Gas phase energy (BOND+ANGLE+DIHED+VDWAALS+EEL+1-4 VDW+1-4 EEL) |
| GSOLV | Solvation energy (EGB+ESURF) |
| TOTAL | Total binding free energy (GGAS+GSOLV) |

### 1.2 FINAL_RESULTS_MMPBSA.dat

Summary statistics file with header information:

```
Header lines (1-20):
| Run on [date]
| gmx_MMPBSA Version=1.6.4
| Complex Structure file: com.tpr
| ...
| Calculations performed using 21 complex frames
| All units are reported in kcal/mol
| SD - Sample standard deviation, SEM - Sample standard error of the mean

Energy Component       Average     SD(Prop.)         SD   SEM(Prop.)        SEM
BOND                   1724.41         36.63      36.63         7.99       7.99
...
Delta (Complex - Receptor - Ligand):
Energy Component       Average     SD(Prop.)         SD   SEM(Prop.)        SEM
...
ΔTOTAL                  -41.23          0.15       1.83         0.03       0.40
```

### 1.3 FINAL_DECOMP_MMPBSA.csv

Per-frame per-residue decomposition data:

**Structure:**
```
Line 1: "Generalized Born Decomposition Energies"
Line 2: [blank]
Line 3: Complex:
Line 4: Total Decomposition Contribution (TDC)
Line 5: Frame #,Residue,Internal,van der Waals,Electrostatic,Polar Solvation,Non-Polar Solv.,TOTAL
Line 6-26: 21 residues for Frame 1
... (repeats for each frame)

Then sections for:
- Receptor: Total Decomposition Contribution
- Ligand: Total Decomposition Contribution  
- Delta (Complex - Receptor - Ligand): Total Decomposition Contribution
```

**Residue Format:** `R:A:ARG:348` or `L:B:MOL:583`
- R/L: Receptor/Ligand
- A/B: Chain ID
- ARG/MOL: Residue name
- 348/583: Residue number

### 1.4 FINAL_DECOMP_MMPBSA.dat

Summary statistics for per-residue decomposition with sections:
- Total Energy Decomposition
- Sidechain Energy Decomposition
- Backbone Energy Decomposition

Each contains: Residue, Internal (Avg/Std.Dev/Std.Err), van der Waals, Electrostatic, Polar Solvation, Non-Polar Solv., TOTAL

---

## 2. What Plots gmx_mmpbsa_ana Generates

Based on gmx_MMPBSA documentation, the GUI generates:

### 2.1 Line Plots
- Energy component evolution over simulation time
- Can show moving average (rolling mean)
- Shows raw data with optional indicators

### 2.2 Bar Plots  
- Per-residue contribution to binding
- Energy term breakdown (all terms or grouped)
- ΔG Binding summary
- Shows average with error bars (standard deviation)

### 2.3 Heatmaps
- Per-residue contribution per frame
- Inter-residue pair contributions (for per-wise decomposition)
- Best for visualizing residue-residue interactions

### 2.4 PyMOL Visualization
- 3D structure with per-residue energy coloring
- Multiple color palettes available

---

## 3. Code Snippets: Parsing Result Files

### 3.1 Parsing FINAL_RESULTS_MMPBSA.csv

```python
import pandas as pd
import numpy as np

def parse_energy_csv(csv_path):
    """Parse FINAL_RESULTS_MMPBSA.csv and return DataFrames for each section."""
    
    with open(csv_path, 'r') as f:
        lines = f.readlines()
    
    # Find section boundaries
    sections = {}
    current_section = None
    
    for i, line in enumerate(lines):
        line = line.strip()
        if line == "Complex Energy Terms":
            current_section = "complex"
            continue
        elif line == "Receptor Energy Terms":
            current_section = "receptor"
            continue
        elif line == "Ligand Energy Terms":
            current_section = "ligand"
            continue
        elif line == "Delta Energy Terms":
            current_section = "delta"
            continue
        elif line.startswith("Frame #,") and current_section:
            # This is the header line
            header_line = i
            break
    
    # Re-read to get section indices
    data = {}
    with open(csv_path, 'r') as f:
        content = f.read()
    
    # Split by section markers
    parts = content.split('\n\n')
    
    for part in parts:
        lines = part.strip().split('\n')
        if not lines:
            continue
            
        if lines[0] == "Complex Energy Terms":
            # Parse complex data
            df = pd.read_csv(pd.io.common.StringIO('\n'.join(lines[2:])))
            data['complex'] = df
        elif lines[0] == "Receptor Energy Terms":
            df = pd.read_csv(pd.io.common.StringIO('\n'.join(lines[2:])))
            data['receptor'] = df
        elif lines[0] == "Ligand Energy Terms":
            df = pd.read_csv(pd.io.common.StringIO('\n'.join(lines[2:])))
            data['ligand'] = df
        elif lines[0] == "Delta Energy Terms":
            df = pd.read_csv(pd.io.common.StringIO('\n'.join(lines[2:])))
            data['delta'] = df
    
    return data


def parse_decomposition_csv(csv_path):
    """Parse FINAL_DECOMP_MMPBSA.csv for per-residue decomposition."""
    
    with open(csv_path, 'r') as f:
        lines = f.readlines()
    
    # Find header lines
    data = {'complex': [], 'receptor': [], 'ligand': [], 'delta': []}
    current_section = None
    
    for line in lines:
        line = line.strip()
        if line == "Complex:":
            current_section = 'complex'
            continue
        elif line == "Receptor:":
            current_section = 'receptor'
            continue
        elif line == "Ligand:":
            current_section = 'ligand'
            continue
        elif line == "Delta (Complex - Receptor - Ligand):":
            current_section = 'delta'
            continue
            
        if line.startswith("Frame #,"):
            continue
        if not line:
            continue
            
        # Parse data line
        parts = line.split(',')
        if len(parts) >= 8 and parts[0].isdigit():
            frame = int(parts[0])
            residue = parts[1]
            internal = float(parts[2])
            vdw = float(parts[3])
            ele = float(parts[4])
            polar = float(parts[5])
            nonpolar = float(parts[6])
            total = float(parts[7])
            
            data[current_section].append({
                'frame': frame,
                'residue': residue,
                'internal': internal,
                'vdw': vdw,
                'electrostatic': ele,
                'polar_solvation': polar,
                'nonpolar_solvation': nonpolar,
                'total': total
            })
    
    # Convert to DataFrames
    result = {}
    for key, rows in data.items():
        if rows:
            result[key] = pd.DataFrame(rows)
    
    return result
```

### 3.2 Parsing FINAL_RESULTS_MMPBSA.dat

```python
def parse_summary_dat(dat_path):
    """Parse FINAL_RESULTS_MMPBSA.dat for summary statistics."""
    
    with open(dat_path, 'r') as f:
        lines = f.readlines()
    
    data = {'complex': {}, 'receptor': {}, 'ligand': {}, 'delta': {}}
    current_section = None
    
    for line in lines:
        line = line.strip()
        
        if 'Complex:' in line and 'Energy Component' not in line:
            current_section = 'complex'
            continue
        elif 'Receptor:' in line and 'Energy Component' not in line:
            current_section = 'receptor'
            continue
        elif 'Ligand:' in line and 'Energy Component' not in line:
            current_section = 'ligand'
            continue
        elif 'Delta' in line and 'Energy Component' not in line:
            current_section = 'delta'
            continue
            
        if current_section and 'Energy Component' in line:
            continue
            
        if current_section and line and not line.startswith('|') and not line.startswith('-'):
            parts = line.split()
            if len(parts) >= 7 and parts[0] not in ['Energy', 'SD']:
                try:
                    term = ' '.join(parts[:-6])
                    avg = float(parts[-6])
                    sd = float(parts[-4])
                    sem = float(parts[-1])
                    data[current_section][term] = {'average': avg, 'sd': sd, 'sem': sem}
                except:
                    pass
    
    return data
```

---

## 4. Matplotlib Code for Energy and Decomposition Plots

### 4.1 Publication-Quality Style Settings

```python
import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np

# Publication-quality settings
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
    'savefig.bbox': 'tight',
    'axes.linewidth': 0.8,
    'axes.spines.top': False,
    'axes.spines.right': False,
})

# Color scheme for energy components
ENERGY_COLORS = {
    'VDWAALS': '#2E86AB',      # Blue
    'EEL': '#A23B72',          # Magenta
    'EGB': '#F18F01',          # Orange
    'ESURF': '#C73E1D',        # Red
    'GGAS': '#3B1F2B',         # Dark
    'GSOLV': '#95C623',        # Green
    'TOTAL': '#1A1A2E',        # Dark blue
    'Delta': '#E63946',        # Red accent
}

def setup_figure(width=3.5, height=2.5):
    """Create a publication-quality figure."""
    fig, ax = plt.subplots(figsize=(width, height))
    return fig, ax
```

### 4.2 Line Plot: Energy Components Over Time

```python
def plot_energy_line(df_delta, output_path='energy_line.png'):
    """Plot binding energy components over simulation frames."""
    
    fig, ax = setup_figure(width=4, height=2.5)
    
    frames = df_delta['Frame #']
    
    # Plot main energy terms
    terms_to_plot = ['VDWAALS', 'EEL', 'EGB', 'ESURF', 'TOTAL']
    
    for term in terms_to_plot:
        if term in df_delta.columns:
            color = ENERGY_COLORS.get(term, '#333333')
            linewidth = 2 if term == 'TOTAL' else 1.2
            linestyle = '-' if term == 'TOTAL' else '--'
            
            ax.plot(frames, df_delta[term], 
                   label=term, 
                   color=color, 
                   linewidth=linewidth,
                   linestyle=linestyle,
                   alpha=0.9)
    
    # Calculate and plot moving average
    window = 5
    if 'TOTAL' in df_delta.columns:
        ma = df_delta['TOTAL'].rolling(window=window, center=True).mean()
        ax.plot(frames, ma, 'k-', linewidth=2, alpha=0.5, 
               label=f'{window}-frame MA')
    
    ax.set_xlabel('Frame')
    ax.set_ylabel('Energy (kcal/mol)')
    ax.set_title('Binding Free Energy Decomposition')
    ax.legend(loc='best', framealpha=0.9)
    ax.grid(True, alpha=0.3, linestyle=':')
    
    # Add horizontal line at 0
    ax.axhline(y=0, color='gray', linestyle='-', linewidth=0.5)
    
    plt.tight_layout()
    fig.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.close()
    
    return output_path
```

### 4.3 Bar Plot: Energy Terms Summary

```python
def plot_energy_bar(delta_summary, output_path='energy_bar.png'):
    """Plot average energy components as bar chart with error bars."""
    
    fig, ax = setup_figure(width=4, height=3)
    
    # Extract delta terms
    terms = ['ΔVDWAALS', 'ΔEEL', 'ΔEGB', 'ΔESURF', 'ΔGGAS', 'ΔGSOLV', 'ΔTOTAL']
    
    # Get data from summary (parsed dat file)
    values = []
    errors = []
    labels = []
    colors = []
    
    for term in terms:
        if term in delta_summary:
            values.append(delta_summary[term]['average'])
            errors.append(delta_summary[term]['sd'])
            labels.append(term.replace('Δ', ''))
            
            if term == 'ΔTOTAL':
                colors.append(ENERGY_COLORS['Delta'])
            elif 'VDW' in term:
                colors.append(ENERGY_COLORS['VDWAALS'])
            elif 'EEL' in term:
                colors.append(ENERGY_COLORS['EEL'])
            elif 'EGB' in term:
                colors.append(ENERGY_COLORS['EGB'])
            elif 'ESURF' in term:
                colors.append(ENERGY_COLORS['ESURF'])
            elif 'GGAS' in term:
                colors.append(ENERGY_COLORS['GGAS'])
            elif 'GSOLV' in term:
                colors.append(ENERGY_COLORS['GSOLV'])
    
    x = np.arange(len(labels))
    
    # Create bars
    bars = ax.bar(x, values, yerr=errors, capsize=3, 
                  color=colors, edgecolor='black', linewidth=0.5,
                  alpha=0.85)
    
    # Add value labels on bars
    for i, (val, err) in enumerate(zip(values, errors)):
        y_pos = val + err + 1 if val >= 0 else val - err - 2
        ax.text(i, y_pos, f'{val:.1f}', ha='center', va='bottom' if val >= 0 else 'top',
               fontsize=8)
    
    ax.set_xticks(x)
    ax.set_xticklabels(labels, rotation=45, ha='right')
    ax.set_ylabel('Energy (kcal/mol)')
    ax.set_title('Binding Free Energy Components')
    ax.axhline(y=0, color='black', linestyle='-', linewidth=0.5)
    
    plt.tight_layout()
    fig.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.close()
    
    return output_path
```

### 4.4 Bar Plot: Per-Residue Decomposition

```python
def plot_decomposition_bar(df_delta, top_n=15, output_path='decomp_bar.png'):
    """Plot per-residue energy contribution to binding."""
    
    # Get average per residue
    residue_avg = df_delta.groupby('residue')['total'].mean().sort_values()
    
    # Select top contributing residues (most negative = favorable)
    top_favorable = residue_avg.head(top_n)
    
    fig, ax = setup_figure(width=5, height=4)
    
    # Create shortened labels
    labels = []
    for res in top_favorable.index:
        parts = res.split(':')
        short_label = f"{parts[2][0]}{parts[3]}"  # e.g., R348
        labels.append(short_label)
    
    # Color by contribution
    colors = ['#2E86AB' if v < 0 else '#E63946' for v in top_favorable.values]
    
    y_pos = np.arange(len(labels))
    
    bars = ax.barh(y_pos, top_favorable.values, color=colors, 
                   edgecolor='black', linewidth=0.5, alpha=0.85)
    
    ax.set_yticks(y_pos)
    ax.set_yticklabels(labels)
    ax.set_xlabel('Energy Contribution (kcal/mol)')
    ax.set_title('Per-Residue Binding Energy Contribution')
    
    # Add vertical line at 0
    ax.axvline(x=0, color='black', linestyle='-', linewidth=0.5)
    
    # Add value labels
    for i, v in enumerate(top_favorable.values):
        x_pos = v - 2 if v < 0 else v + 2
        ax.text(x_pos, i, f'{v:.1f}', ha='right' if v < 0 else 'left', 
               va='center', fontsize=8)
    
    plt.tight_layout()
    fig.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.close()
    
    return output_path
```

### 4.5 Heatmap: Residue Contributions Over Time

```python
def plot_decomposition_heatmap(df_delta, output_path='decomp_heatmap.png'):
    """Plot per-residue energy contribution as heatmap over frames."""
    
    import matplotlib.colors as mcolors
    
    # Pivot to create matrix: frames x residues
    pivot_data = df_delta.pivot(index='frame', columns='residue', values='total')
    
    # Shorten column names
    short_cols = {}
    for col in pivot_data.columns:
        parts = col.split(':')
        short_cols[col] = f"{parts[2][0]}{parts[3]}"
    pivot_data = pivot_data.rename(columns=short_cols)
    
    fig, ax = plt.subplots(figsize=(8, 4))
    
    # Create diverging colormap (blue = favorable, red = unfavorable)
    vmax = max(abs(pivot_data.values.min()), abs(pivot_data.values.max()))
    
    im = ax.imshow(pivot_data.values, aspect='auto', cmap='RdBu_r',
                   vmin=-vmax, vmax=vmax)
    
    # Set labels
    ax.set_yticks(range(len(pivot_data.index)))
    ax.set_yticklabels(pivot_data.index)
    ax.set_xticks(range(len(pivot_data.columns)))
    ax.set_xticklabels(pivot_data.columns, rotation=90, fontsize=8)
    
    ax.set_xlabel('Residue')
    ax.set_ylabel('Frame')
    ax.set_title('Per-Residue Energy Contribution (kcal/mol)')
    
    # Add colorbar
    cbar = plt.colorbar(im, ax=ax, label='Energy (kcal/mol)', shrink=0.8)
    
    plt.tight_layout()
    fig.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.close()
    
    return output_path
```

---

## 5. Best Practices for Scientific Visualization

### 5.1 Publication Standards

1. **Figure Dimensions**: Use single-column (3.5"), double-column (7"), or full-page (6") widths
2. **Font Sizes**: Minimum 6pt for axis labels, 8pt for tick labels, 10pt for legends
3. **Line Weights**: 0.5-1pt for data lines, 0.5pt for gridlines
4. **Color Accessibility**: Use colorblind-safe palettes (viridis, plasma, or custom)

### 5.2 Error Representation

- Use standard deviation (SD) for sample variability
- Use standard error of mean (SEM) for confidence intervals
- Always include error bars in bar charts
- Consider showing both raw data points and error bars

### 5.3 Axis Labels and Units

- Always include units (kcal/mol for energies)
- Use SI prefixes appropriately
- Label axes clearly with descriptive terms
- Consider dimensionless quantities for comparisons

### 5.4 Color Schemes

```python
# Scientifically appropriate colormaps
from matplotlib.cm import get_cmap

# Sequential (for magnitude): viridis, plasma, inferno
# Diverging (for signed data): RdBu_r, coolwarm, seismic  
# Categorical: Set2, Tableau10, Classic

# Custom accessible palette
SCIENTIFIC_PALETTE = [
    '#0173B2',  # Blue
    '#DE8F05',  # Orange
    '#029E73',  # Green
    '#CC78BC',  # Purple
    '#CA9161',  # Brown
    '#949494',  # Gray
    '#FBAFE4',  # Pink
    '#56B4E9',  # Light blue
]
```

### 5.5 Layout and Composition

```python
# Constrained layout for proper spacing
fig, axes = plt.subplots(2, 2, constrained_layout=True)

# GridSpec for complex layouts  
from matplotlib.gridspec import GridSpec
gs = GridSpec(2, 2, figure=fig, hspace=0.3, wspace=0.3)
```

### 5.6 Export Formats

| Format | Use Case |
|--------|----------|
| PNG | Drafts, web |
| PDF | Publications (vector) |
| SVG | Editable figures |
| EPS | LaTeX documents |

---

## 6. Complete Example Script

```python
#!/usr/bin/env python3
"""
gmx_MMPBSA Results Plotting Script
Generates publication-quality plots from gmx_MMPBSA output files.
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from pathlib import Path


# ============================================================================
# Configuration
# ============================================================================

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

COLORS = {
    'VDWAALS': '#2E86AB',
    'EEL': '#A23B72', 
    'EGB': '#F18F01',
    'ESURF': '#C73E1D',
    'GGAS': '#3B1F2B',
    'GSOLV': '#95C623',
    'TOTAL': '#1A1A2E',
}


# ============================================================================
# Parsing Functions
# ============================================================================

def parse_results_csv(csv_path):
    """Parse FINAL_RESULTS_MMPBSA.csv."""
    
    with open(csv_path, 'r') as f:
        content = f.read()
    
    sections = content.split('\n\n')
    data = {}
    
    for section in sections:
        lines = [l.strip() for l in section.strip().split('\n') if l.strip()]
        if not lines:
            continue
            
        if lines[0] == 'Complex Energy Terms':
            data['complex'] = _parse_energy_table(lines[2:])
        elif lines[0] == 'Receptor Energy Terms':
            data['receptor'] = _parse_energy_table(lines[2:])
        elif lines[0] == 'Ligand Energy Terms':
            data['ligand'] = _parse_energy_table(lines[2:])
        elif lines[0] == 'Delta Energy Terms':
            data['delta'] = _parse_energy_table(lines[2:])
    
    return data


def _parse_energy_table(lines):
    """Helper to parse energy table from CSV lines."""
    if not lines:
        return None
    df = pd.read_csv(lines[0] + '\n' + '\n'.join(lines[1:]))
    return df


def parse_decomposition_csv(csv_path):
    """Parse FINAL_DECOMP_MMPBSA.csv."""
    
    with open(csv_path, 'r') as f:
        lines = f.readlines()
    
    data = {'complex': [], 'receptor': [], 'ligand': [], 'delta': []}
    current_section = None
    
    for line in lines:
        line = line.strip()
        
        if line == 'Complex:':
            current_section = 'complex'
            continue
        elif line == 'Receptor:':
            current_section = 'receptor'
            continue
        elif line == 'Ligand:':
            current_section = 'ligand'
            continue
        elif line == 'Delta (Complex - Receptor - Ligand):':
            current_section = 'delta'
            continue
            
        if line.startswith('Frame #,') or not line:
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


# ============================================================================
# Plotting Functions
# ============================================================================

def plot_binding_energy_timecourse(results, output_dir='plots'):
    """Plot binding energy components over time."""
    
    Path(output_dir).mkdir(exist_ok=True)
    
    if 'delta' not in results:
        print("No delta data found")
        return
    
    df = results['delta']
    fig, ax = plt.subplots(figsize=(4, 2.5))
    
    frames = df['Frame #']
    
    for term in ['VDWAALS', 'EEL', 'EGB', 'ESURF', 'TOTAL']:
        if term in df.columns:
            lw = 2 if term == 'TOTAL' else 1.2
            ax.plot(frames, df[term], label=term, linewidth=lw,
                   color=COLORS.get(term, '#333'), alpha=0.9)
    
    ax.set_xlabel('Frame')
    ax.set_ylabel('Energy (kcal/mol)')
    ax.set_title('Binding Free Energy Decomposition')
    ax.legend(loc='best', framealpha=0.9)
    ax.grid(True, alpha=0.3, linestyle=':')
    ax.axhline(y=0, color='gray', linestyle='-', lw=0.5)
    
    plt.tight_layout()
    fig.savefig(f'{output_dir}/energy_timecourse.png', dpi=300)
    plt.savefig(f'{output_dir}/energy_timecourse.pdf')
    plt.close()
    
    print(f"Saved: {output_dir}/energy_timecourse.png")


def plot_decomposition_bar(decomp, output_dir='plots', top_n=15):
    """Plot per-residue decomposition."""
    
    Path(output_dir).mkdir(exist_ok=True)
    
    if 'delta' not in decomp:
        print("No delta decomposition data found")
        return
    
    df = decomp['delta']
    avg_by_residue = df.groupby('residue')['total'].mean().sort_values()
    
    # Select top contributors
    top = avg_by_residue.head(top_n)
    
    fig, ax = plt.subplots(figsize=(4, 4))
    
    # Shorten labels
    labels = [r.split(':')[2][0] + r.split(':')[3] for r in top.index]
    colors = ['#2E86AB' if v < 0 else '#E63946' for v in top.values]
    
    y = np.arange(len(labels))
    ax.barh(y, top.values, color=colors, edgecolor='black', lw=0.5)
    
    ax.set_yticks(y)
    ax.set_yticklabels(labels)
    ax.set_xlabel('Energy (kcal/mol)')
    ax.set_title('Per-Residue Binding Contribution')
    ax.axvline(x=0, color='black', lw=0.5)
    
    plt.tight_layout()
    fig.savefig(f'{output_dir}/decomposition_bar.png', dpi=300)
    plt.savefig(f'{output_dir}/decomposition_bar.pdf')
    plt.close()
    
    print(f"Saved: {output_dir}/decomposition_bar.png")


def plot_energy_summary(results, output_dir='plots'):
    """Plot energy summary bar chart."""
    
    Path(output_dir).mkdir(exist_ok=True)
    
    if 'delta' not in results:
        return
    
    df = results['delta']
    
    # Calculate means and stds
    terms = ['VDWAALS', 'EEL', 'EGB', 'ESURF', 'GGAS', 'GSOLV', 'TOTAL']
    means = []
    stds = []
    labels = []
    colors_list = []
    
    for term in terms:
        col = term if term in df.columns else f'Δ{term}'
        if col in df.columns:
            means.append(df[col].mean())
            stds.append(df[col].std())
            labels.append(term)
            if term == 'TOTAL':
                colors_list.append('#E63946')
            else:
                colors_list.append(COLORS.get(term, '#333'))
    
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
    fig.savefig(f'{output_dir}/energy_summary.png', dpi=300)
    plt.savefig(f'{output_dir}/energy_summary.pdf')
    plt.close()
    
    print(f"Saved: {output_dir}/energy_summary.png")


# ============================================================================
# Main
# ============================================================================

def main(work_dir):
    """Main function to generate all plots."""
    
    results_csv = Path(work_dir) / 'FINAL_RESULTS_MMPBSA.csv'
    decomp_csv = Path(work_dir) / 'FINAL_DECOMP_MMPBSA.csv'
    
    print(f"Processing: {work_dir}")
    
    # Parse data
    if results_csv.exists():
        results = parse_results_csv(str(results_csv))
        print(f"  Loaded energy results: {list(results.keys())}")
        
        # Generate plots
        plot_binding_energy_timecourse(results)
        plot_energy_summary(results)
    else:
        print(f"  ERROR: {results_csv} not found")
    
    if decomp_csv.exists():
        decomp = parse_decomposition_csv(str(decomp_csv))
        print(f"  Loaded decomposition: {list(decomp.keys())}")
        
        plot_decomposition_bar(decomp)
    else:
        print(f"  WARNING: {decomp_csv} not found")


if __name__ == '__main__':
    import sys
    
    if len(sys.argv) > 1:
        work_dir = sys.argv[1]
    else:
        work_dir = '.'  # Current directory
    
    main(work_dir)
```

---

## 7. Sources

### Primary (HIGH confidence)
- gmx_MMPBSA documentation: https://valdes-tresanco-ms.github.io/gmx_MMPBSA/dev/analyzer/
- Analysis of provided output files in `~/dparker/dp_xinyi/ana_code/com_ana_trj/`

### Secondary (MEDIUM confidence)
- Matplotlib official documentation: https://matplotlib.org/stable/users/explain/customizing.html
- Scientific visualization best practices: Nature, Science publication guidelines

### Tertiary (LOW confidence)
- Community examples on gmx_MMPBSA GitHub discussions

---

## 8. Open Questions

1. **MMXSA file format**: The `COMPACT_MMXSA_RESULTS.mmxsa` file appears to be binary. Further investigation needed to determine if it contains additional data not in CSV files.

2. **Entropy calculations**: If entropy terms (Interaction Entropy, nmode) are calculated, additional parsing would be needed.

3. **Multiple trajectory protocol**: The script assumes single trajectory. Multiple trajectory files would require different parsing logic.

---

## 9. Metadata

**Confidence breakdown:**
- File formats: HIGH - Direct analysis of output files
- Plot types: HIGH - Based on gmx_MMPBSA documentation
- Code snippets: HIGH - Based on standard pandas/matplotlib patterns
- Best practices: MEDIUM - Industry standards with some interpretation

**Research date:** March 22, 2026
**Valid until:** 6 months (matplotlib API stable, gmx_MMPBSA format stable)
