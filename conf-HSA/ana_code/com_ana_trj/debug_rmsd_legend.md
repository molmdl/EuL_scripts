# RMSD Plot Legend Overlap Issue - Debug Report

## Summary

**Issue:** Two legend boxes overlap in the RMSD plot output, obscuring content.

**Root Cause:** The `plot_rmsd_overlay` function creates two separate annotation boxes that occupy the same screen space.

**File:** `~/dparker/dp_xinyi/ana_code/com_ana_trj/analyze_trajectory_standalone.py`

**Function:** `plot_rmsd_overlay` (lines 257-289)

---

## Root Cause Analysis

### The Two Overlapping Boxes

#### Box 1: Main Legend (Line 276)
```python
ax.legend(loc='best', framealpha=0.9)
```
- **Content:** 4 legend entries
  - 'Receptor Backbone' (solid blue line)
  - 'Ligand Heavy Atoms' (solid red line)
  - 'Receptor (MA)' (dashed blue line - moving average)
  - 'Ligand (MA)' (dashed red line - moving average)
- **Location:** matplotlib's 'best' algorithm chooses optimal placement
- **Problem:** For time series data starting from origin, 'best' often chooses **top-left** corner

#### Box 2: Statistics Text Box (Lines 280-285)
```python
text = (f"Receptor: mean={df_receptor['rmsd'].mean():.2f}A, "
        f"std={df_receptor['rmsd'].std():.2f}A\n"
        f"Ligand: mean={df_ligand['rmsd'].mean():.2f}A, "
        f"std={df_ligand['rmsd'].std():.2f}A")
ax.text(0.02, 0.98, text, transform=ax.transAxes, fontsize=9,
        verticalalignment='top', bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))
```
- **Content:** Mean and standard deviation statistics for receptor and ligand
- **Location:** **Explicitly at top-left** corner (0.02, 0.98 in axes coordinates)
- **Problem:** Fixed position conflicts with legend's auto-placement

### Why They Overlap

1. matplotlib's `loc='best'` algorithm evaluates multiple positions and scores them
2. For typical RMSD plots (lines starting near origin), top-left often scores highest
3. The statistics box is hardcoded to top-left (0.02, 0.98)
4. Both boxes have white backgrounds with rounded corners, making overlap obvious

---

## Proposed Fixes

### Option A: Move Statistics Box to Bottom-Left (Recommended)

This keeps the legend at matplotlib's chosen 'best' location and moves the statistics box to an area that won't interfere.

**Change line 284:**
```python
# Before:
ax.text(0.02, 0.98, text, transform=ax.transAxes, fontsize=9,
        verticalalignment='top', bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))

# After:
ax.text(0.02, 0.02, text, transform=ax.transAxes, fontsize=9,
        verticalalignment='bottom', bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))
```

**Pros:**
- Minimal code change
- Statistics box at bottom won't interfere with legend at top
- Bottom area is typically empty in RMSD plots

**Cons:**
- Statistics box might overlap with x-axis label on small plots

---

### Option B: Explicitly Place Legend at Upper Right

Force the legend to upper-right, keeping statistics box at top-left.

**Change line 276:**
```python
# Before:
ax.legend(loc='best', framealpha=0.9)

# After:
ax.legend(loc='upper right', framealpha=0.9)
```

**Pros:**
- Very simple one-line change
- Upper-right often has white space in RMSD plots

**Cons:**
- Less flexible if data shape changes
- Might not be optimal for all trajectory types

---

### Option C: Combine Statistics into Legend

Remove the separate text box and add statistics as additional legend entries or in the title.

**Change lines 276 and 280-285:**
```python
# Create custom legend with statistics
from matplotlib.patches import Patch
from matplotlib.lines import Line2D

legend_elements = [
    Line2D([0], [0], color='blue', linewidth=1.5, label='Receptor Backbone'),
    Line2D([0], [0], color='red', linewidth=1.5, label='Ligand Heavy Atoms'),
    Line2D([0], [0], color='blue', linewidth=2, linestyle='--', label='Receptor (MA)'),
    Line2D([0], [0], color='red', linewidth=2, linestyle='--', label='Ligand (MA)'),
    Patch(facecolor='white', edgecolor='gray', 
          label=f"Stats: R={df_receptor['rmsd'].mean():.2f}±{df_receptor['rmsd'].std():.2f}A"),
    Patch(facecolor='white', edgecolor='gray', 
          label=f"Stats: L={df_ligand['rmsd'].mean():.2f}±{df_ligand['rmsd'].std():.2f}A"),
]
ax.legend(handles=legend_elements, loc='best', framealpha=0.9)

# Remove the separate ax.text() call
```

**Pros:**
- All information in one unified box
- No overlap possible

**Cons:**
- More complex code
- Legend box becomes taller

---

## Recommended Fix

**Use Option A:** Move the statistics text box to the bottom-left corner.

This is the cleanest solution with minimal code change and maintains visual separation between the legend (lines identification) and statistics (numerical summary).

### Implementation

```python
def plot_rmsd_overlay(df_receptor, df_ligand, output_path, dpi=300):
    fig, ax = plt.subplots(figsize=(10, 6))
    
    ax.plot(df_receptor['time_ps'], df_receptor['rmsd'], 
            'b-', linewidth=1.5, label='Receptor Backbone', alpha=0.8)
    ax.plot(df_ligand['time_ps'], df_ligand['rmsd'], 
            'r-', linewidth=1.5, label='Ligand Heavy Atoms', alpha=0.8)
    
    for df, color, label in [(df_receptor, 'blue', 'Receptor (MA)'),
                              (df_ligand, 'red', 'Ligand (MA)')]:
        window = min(10, len(df))
        if window > 1:
            ma = df['rmsd'].rolling(window=window, center=True).mean()
            ax.plot(df['time_ps'], ma, color=color, linewidth=2, 
                    linestyle='--', alpha=0.6, label=label)
    
    ax.set_xlabel('Time (ps)', fontsize=12)
    ax.set_ylabel('RMSD (A)', fontsize=12)
    ax.set_title('RMSD: Receptor Backbone vs Ligand Heavy Atoms', fontsize=14)
    ax.legend(loc='best', framealpha=0.9)
    ax.grid(True, alpha=0.3)
    ax.set_xlim(df_receptor['time_ps'].min(), df_receptor['time_ps'].max())
    
    text = (f"Receptor: mean={df_receptor['rmsd'].mean():.2f}A, "
            f"std={df_receptor['rmsd'].std():.2f}A\n"
            f"Ligand: mean={df_ligand['rmsd'].mean():.2f}A, "
            f"std={df_ligand['rmsd'].std():.2f}A")
    
    # FIX: Move statistics box to bottom-left to avoid legend overlap
    ax.text(0.02, 0.02, text, transform=ax.transAxes, fontsize=9,
            verticalalignment='bottom', bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))
    
    plt.tight_layout()
    plt.savefig(output_path, dpi=dpi, bbox_inches='tight')
    plt.close()
```

---

## Verification Steps

After applying the fix:

1. Run the RMSD analysis:
   ```bash
   python analyze_trajectory_standalone.py --aim 1 --tpr com.tpr --xtc com_traj.xtc
   ```

2. Check the output plot:
   - `~/dparker/dp_xinyi/ana_code/com_ana_trj/output/rmsd/rmsd_plot.png`

3. Verify:
   - [ ] Legend box is clearly visible in upper area
   - [ ] Statistics box is clearly visible in lower-left area
   - [ ] No overlap between boxes
   - [ ] All plot lines and labels are readable