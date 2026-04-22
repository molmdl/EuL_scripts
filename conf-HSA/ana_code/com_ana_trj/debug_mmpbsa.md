# Debug Report: Empty MMPBSA Energy Plots

## Summary
The MMPBSA energy timecourse and energy summary plots are empty because the `parse_results_csv()` function incorrectly parses the CSV file, skipping the header row.

## Root Cause
**Location:** `analyze_trajectory_standalone.py`, lines 597-604 (function `parse_results_csv`)

**Bug:** The parsing logic uses `lines[2:]` to read CSV data, but this skips the header row which is at `lines[1]`.

### Code Analysis

The CSV file structure for each section is:
```
Delta Energy Terms           <- lines[0] (section title)
Frame #,BOND,ANGLE,...       <- lines[1] (HEADER - needed by pandas)
1,0.0,-0.0,-40.15,...        <- lines[2] (first data row)
2,-0.0,-0.0,-43.59,...       <- lines[3] (second data row)
...
```

The current code does:
```python
data['delta'] = pd.read_csv(StringIO('\n'.join(lines[2:])))
```

This passes data starting from `lines[2]` (first data row) to pandas, but pandas has no header to use for column names. As a result, pandas treats the first data row as column names, creating malformed column names like `'1'`, `'0.0'`, `'-0.0.1'`, etc.

### Evidence

When parsing the Delta section:
- Expected columns: `['Frame #', 'BOND', 'ANGLE', 'DIHED', 'VDWAALS', 'EEL', '1-4 VDW', '1-4 EEL', 'EGB', 'ESURF', 'GGAS', 'GSOLV', 'TOTAL']`
- Actual columns: `['1', '0.0', '-0.0', '-0.0.1', '-40.15', '-4.87', '-0.0.2', '-0.0.3', '11.79', '-4.93', '-45.01', '6.85', '-38.16']`

The plotting functions then search for columns like `'VDWAALS'`, `'EEL'`, `'TOTAL'` which don't exist, resulting in empty plots.

## Proposed Fix

Change all `lines[2:]` to `lines[1:]` in the `parse_results_csv()` function:

```python
def parse_results_csv(csv_path):
    with open(csv_path, 'r') as f:
        content = f.read()
    
    sections = content.split('\n\n')
    data = {}
    
    for section in sections:
        lines = [l.strip() for l in section.strip().split('\n') if l.strip()]
        if not lines:
            continue
            
        if lines[0] == 'Complex Energy Terms':
            data['complex'] = pd.read_csv(StringIO('\n'.join(lines[1:])))  # Changed from lines[2:]
        elif lines[0] == 'Receptor Energy Terms':
            data['receptor'] = pd.read_csv(StringIO('\n'.join(lines[1:])))  # Changed from lines[2:]
        elif lines[0] == 'Ligand Energy Terms':
            data['ligand'] = pd.read_csv(StringIO('\n'.join(lines[1:])))    # Changed from lines[2:]
        elif lines[0] == 'Delta Energy Terms':
            data['delta'] = pd.read_csv(StringIO('\n'.join(lines[1:])))      # Changed from lines[2:]
    
    return data
```

## Additional Issue: Complex Section Not Parsed

The Complex section starts with "GENERALIZED BORN:" instead of "Complex Energy Terms":

```
GENERALIZED BORN:             <- lines[0]
Complex Energy Terms          <- lines[1]
Frame #,BOND,ANGLE,...        <- lines[2]
```

The condition `if lines[0] == 'Complex Energy Terms'` never matches because `lines[0]` is "GENERALIZED BORN:".

**However**, this is not critical for plotting because the plotting functions only use 'delta' data. The fix above addresses the primary issue.

## Status: FIXED and VERIFIED

The fix has been applied and tested. All four MMPBSA plot files are now generated correctly:
- `energy_timecourse.png` (102 KB) - Shows energy decomposition over frames
- `energy_summary.png` (66 KB) - Shows bar chart of average energies  
- `decomposition_bar.png` (102 KB) - Shows per-residue binding contribution
- `decomposition_heatmap.png` (128 KB) - Shows per-residue energy heatmap

## Verification Steps (Completed)

1. ✅ Ran: `python analyze_trajectory_standalone.py --aim 4 --mmpbsa-dir . --output-dir ./test_output`
2. ✅ Checked output/mmpbsa/energy_timecourse.png - shows energy decomposition over frames
3. ✅ Checked output/mmpbsa/energy_summary.png - shows bar chart of average energies
4. ✅ Verified parsing returns correct columns: `['Frame #', 'BOND', 'ANGLE', 'DIHED', 'VDWAALS', 'EEL', '1-4 VDW', '1-4 EEL', 'EGB', 'ESURF', 'GGAS', 'GSOLV', 'TOTAL']`

## Files Involved

- `~/dparker/dp_xinyi/ana_code/com_ana_trj/analyze_trajectory_standalone.py` - **FIXED**
- `~/dparker/dp_xinyi/ana_code/com_ana_trj/FINAL_RESULTS_MMPBSA.csv` - Source data file

## MMXSA Binary Format

The `COMPACT_MMXSA_RESULTS.mmxsa` file is a pickled `SimpleNamespace` object from gmx_MMPBSA. It contains:
- Summary statistics (averages, standard deviations)
- File paths and configuration
- But NOT per-frame energy data needed for timecourse plots

The per-frame data must be read from `FINAL_RESULTS_MMPBSA.csv`.
