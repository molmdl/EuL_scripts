# Research: SASA-Based Rescoring for Ligand Selection

**Researched:** March 26, 2026  
**Domain:** Molecular dynamics workflow integration / Solvent Accessible Surface Area calculation  
**Confidence:** MEDIUM (verified via code analysis and algorithm review)

---

## Summary

This research addresses the integration of an approximate SASA (Solvent Accessible Surface Area) calculation into the ligand selection workflow of `dock2com_1.py`. The current workflow selects the best docking pose based on docking scores (e.g., `minimizedAffinity`, `CNNscore`). The goal is to add SASA-based rescoring as an additional selection criterion, where lower SASA values indicate better ligand burial (more favorable binding).

The `quick_sasa.py` script provides the core algorithm: a fast VDW+Gaussian approximation method that computes relative SASA values without requiring external GROMACS tools. This can be integrated into the existing pose selection workflow to rescore docking poses after the initial selection.

**Primary recommendation:** Implement a two-stage selection process: (1) parse all SDF models and compute SASA for each pose using receptor coordinates, (2) select the best pose using a weighted combination or by filtering based on SASA threshold. This preserves all original functionality while adding SASA-based rescoring capability.

---

## Current dock2com_1.py Functionality

### Core Workflow (6 Steps)

| Step | Function | Description |
|------|----------|-------------|
| 1 | `select_best_across_sdfiles()` | Parse all SDF files, select best model by metric |
| 2 | Receptor auto-detection | Find receptor GRO file from SDF filename pattern |
| 3 | `sdf_pose_to_gro()` | Convert SDF pose to GRO using ITP topology |
| 4 | `combine_coordinates()` | Merge receptor + ligand into complex GRO |
| 5 | `extract_receptor_topology()` | Extract receptor ITP from system TOP |
| 6 | `create_system_topology()` | Generate system topology file |

### Model Selection Logic

The current selection is performed in `select_best_across_sdfiles()`:

```python
def select_best_across_sdfiles(sdf_paths, metric=DEFAULT_METRIC):
    best_model = None
    best_score = None
    best_file = None
    reverse = metric not in LOWER_IS_BETTER  # minimizedAffinity: lower is better
    
    for sdf_path in sdf_paths:
        models = parse_sdf(sdf_path)
        for model in models:
            if metric not in model["scores"]:
                continue
            score = model["scores"][metric]
            # Compare and update best_model...
```

**Key metrics available in SDF files:**
- `minimizedAffinity` (kcal/mol) - docking affinity, lower is better
- `CNNscore` - convolutional neural network pose score, higher is better
- `CNNaffinity` - CNN affinity prediction, lower is better
- `intramol` - intramolecular energy

### File Format Support

| Format | Parser Function | Purpose |
|--------|-----------------|---------|
| SDF | `parse_sdf()` | Docking poses with scores in properties |
| MOL2 | `parse_mol2()` | Template for hydrogen reconstruction |
| ITP | `parse_itp()` | Ligand topology (atom ordering) |
| GRO | `write_gro()` | Coordinate output |
| PDB | `parse_pdb()` | Receptor structure input |

---

## quick_sasa.py SASA Approximation Analysis

### Algorithm Overview

The `calculate_sasa_internal()` function implements a fast VDW+Gaussian approximation:

```
SASA_atom ≈ 4π(r_vdw + r_probe)² × exposure_factor
```

Where:
- `r_vdw` = van der Waals radius (Bondi 1964, Mantina 2009)
- `r_probe` = 1.4 Å (water probe radius, configurable)
- `exposure_factor = exp(-λ × neighbor_sum)`
- `neighbor_sum` = Gaussian-weighted count of nearby atoms

### Gaussian Parameters

| Parameter | Value | Purpose |
|-----------|-------|---------|
| λ (lambda) | 0.5 | Scaling factor for burial calculation |
| σ (sigma) | 1.5 Å | Gaussian width for neighbor density |
| neighbor_cutoff | 10.0 Å | Neighbor search radius |

### VDW Radii Table (Key Elements)

```python
vdw_radii_ang = {
    'H': 1.20, 'C': 1.70, 'N': 1.55, 'O': 1.52, 
    'S': 1.80, 'Cl': 1.75, 'Fe': 2.00, 'Zn': 1.39,
    'Eu': 2.05,  # Europium (for lanthanide binding proteins)
    # ... full table in code
}
default_radius = 1.70  # Carbon fallback
```

### Performance Characteristics

- **Speed:** ~0.7ms/frame for typical protein-ligand complexes
- **Accuracy:** 10-20% vs. probe-based methods (suitable for relative comparisons)
- **Dependencies:** numpy, scipy.spatial.cKDTree, MDAnalysis (for trajectory handling)

### Integration Points

The SASA calculation can work in two modes relevant to this project:

1. **Static structure mode:** Calculate SASA for a single complex (ligand + receptor)
2. **Trajectory mode:** Process multiple frames with time-series output

For docking pose selection, we need **mode 1**: calculate SASA for each ligand pose in the context of its receptor.

---

## Integration Approach for SASA Rescoring

### Data Flow

```
SDF files (hsa*-phe_sssD.sdf)
    │
    ▼
parse_sdf() → models with coordinates + scores
    │
    ▼
For each model:
    - Extract ligand coordinates from SDF
    - Load receptor coordinates from hsa*_ali.gro
    - Calculate ligand SASA using VDW+Gaussian method
    - Store SASA value in model metadata
    │
    ▼
select_best_across_sdfiles(metric="SASA")  # or composite
    │
    ▼
Output: Best pose selected by SASA criterion
```

### Implementation Strategy

#### Option 1: Two-Stage Selection (Recommended)

1. **Stage 1 - Initial filtering:** Select top-N poses by docking score
2. **Stage 2 - SASA rescoring:** Re-rank by SASA (lower = better buried)

**Pros:** Combines docking accuracy with binding site complementarity  
**Cons:** Requires two passes through data

#### Option 2: Composite Score

Create a new metric: `composite_score = α × normalized_affinity + β × normalized_SASA`

Where:
- `normalized_affinity` = (score - min) / (max - min)
- `normalized_SASA` = (sasa - min) / (max - min) or inverse
- α, β are configurable weights

**Pros:** Single selection pass, tunable balance  
**Cons:** Weight selection requires validation

#### Option 3: Post-Selection SASA Filter

1. Select pose by original metric
2. Calculate SASA for selected pose
3. Flag if SASA exceeds threshold (poorly buried)

**Pros:** Minimal code change, provides quality check  
**Cons:** Doesn't use SASA for ranking

### Recommended Implementation

**Option 1 with CLI extension:**

```bash
# Original behavior (unchanged)
python dock2com_1.py -i lig.itp -s hsa*.sdf -r topol.top

# With SASA rescoring
python dock2com_1.py -i lig.itp -s hsa*.sdf -r topol.top --sasa-rescore

# With custom SASA metric
python dock2com_1.py -i lig.itp -s hsa*.sdf -r topol.top --metric sasa
```

---

## Key Considerations and Potential Challenges

### 1. Coordinate System Alignment

**Challenge:** SDF poses use Ångstrom coordinates, receptor GRO uses nm.

**Solution:** Convert during SASA calculation (quick_sasa.py uses Å internally).

```python
# In quick_sasa.py:
group_coords = atom_group.positions  # Angstroms
# ... calculate ...
return atom_sasa_array / 100.0  # Convert Å² to nm²
```

### 2. Atom Type Mapping

**Challenge:** GRO/ITP use GROMACS atom types, SDF uses element symbols.

**Solution:** Use element symbols from `parse_sdf()` (field 31:34) which provides element directly.

```python
# From dock2com_1.py parse_sdf():
el = ln[31:34].strip()
if el == "*":
    el = "EU"  # Handle Europium
atoms.append({"el": el, ...})
```

### 3. Receptor Conformation Matching

**Challenge:** Each SDF file corresponds to a specific receptor conformation (hsa0, hsa1, etc.).

**Solution:** Auto-detect receptor GRO from SDF filename pattern (existing `derive_receptor_gro_from_sdf()`).

```python
# Existing code:
rec_gro, prefix = derive_receptor_gro_from_sdf(best_file, ...)
```

### 4. Performance

**Challenge:** SASA calculation for 32+ poses per file × 10 files = 320+ calculations.

**Solution:** The fast VDW+Gaussian method handles ~0.7ms/frame. For ~300 poses:
- Estimated time: ~200ms (acceptable for workflow)
- Can add optional flag to enable/disable for performance

### 5. Missing Receptor Files

**Challenge:** SASA calculation requires receptor coordinates.

**Solution:** Fail gracefully if receptor GRO not available:

```python
if not os.path.exists(rec_gro):
    print(f"WARNING: Cannot calculate SASA without receptor GRO: {rec_gro}")
    print("  Falling back to original selection metric")
    # Fall back to original selection
```

### 6. Hydrogen Treatment

**Challenge:** Quick SASA uses heavy atoms only by default. Docking poses may have implicit hydrogens.

**Solution:** Use same approach as `quick_sasa.py` (element-based, handles H correctly).

---

## Recommended Implementation Strategy

### Phase 1: Core SASA Calculator (Standalone Function)

```python
def calculate_pose_sasa(ligand_atoms, receptor_atoms, probe_radius=1.4):
    """
    Calculate ligand SASA in context of receptor.
    
    Parameters:
    -----------
    ligand_atoms : list of dict
        [{"el": "C", "x": 1.0, "y": 2.0, "z": 3.0}, ...]
    receptor_atoms : list of dict
        Same format as ligand_atoms
    
    Returns:
    --------
    float : SASA in nm²
    """
    # Use VDW+Gaussian method from quick_sasa.py
    # Return summed SASA for all ligand atoms
```

### Phase 2: Integration into Selection

```python
def select_best_with_sasa(sdf_paths, rec_gro_pattern, metric="sasa"):
    """
    Select best model using SASA as primary metric.
    
    For each SDF pose:
    1. Parse coordinates
    2. Load matching receptor GRO
    3. Calculate SASA
    4. Store in model["scores"]["sasa"]
    
    Return best model by SASA (lower = better)
    """
```

### Phase 3: CLI Extensions

```python
# Add to build_parser():
sel_group.add_argument("--sasa-rescore", action="store_true",
                       help="Rescore poses by SASA after initial selection")
sel_group.add_argument("--sasa-threshold", type=float, default=None,
                       help="Maximum SASA (nm²) to accept pose")
sel_group.add_argument("--composite-metrics", nargs="+",
                       help="Combine multiple metrics for selection")
```

### Phase 4: Output Enhancement

Add SASA to summary output:

```python
print(f"  SASA: {model['scores'].get('sasa', 'N/A'):.4f} nm²")
```

---

## Data Format Reference

### SDF Pose Structure (from parse_sdf)

```python
{
    "model_num": 1,
    "title": "optfreq",
    "atoms": [
        {"idx": 1, "el": "O", "x": 9.5202, "y": -4.4198, "z": -26.2499},
        ...
    ],
    "bonds": [...],
    "scores": {
        "minimizedAffinity": -7.96,
        "CNNscore": 0.1776,
        "CNNaffinity": 7.415,
        "intramol": -3.49
    }
}
```

### Receptor GRO Format (hsa*_ali.gro)

Standard GROMACS GRO format with:
- Box dimensions for periodic boundary conditions
- Atomic coordinates in nm
- Residue names and atom names

### Ligand ITP Format (lig_g.itp)

GROMACS ITP format with:
- `[ atoms ]` section: atom types, charges, masses
- `[ bonds ]` section: connectivity
- Used for atom ordering when converting SDF to GRO

---

## Testing Recommendations

### Unit Tests

1. `test_sasa_calculation()` - Verify SASA for known structure (benzene in water)
2. `test_sasa_with_receptor()` - Ligand SASA in presence/absence of receptor
3. `test_coordinate_conversion()` - Verify Å↔nm conversion

### Integration Tests

1. `test_full_workflow_sasa()` - Run with `--sasa-rescore` flag
2. `test_receptor_auto_detect()` - Verify hsa* matching
3. `test_graceful_degradation()` - No receptor GRO available

### Validation

Compare SASA-based selection vs. affinity-based selection for test cases where binding mode is known.

---

## Sources

### Primary (Code Analysis)
- `dock2com_1.py` (lines 40-162) - SDF parsing and model selection
- `quick_sasa.py` (lines 114-300) - VDW+Gaussian SASA algorithm
- `dock2com_1.sh` - Command invocation pattern

### Secondary (Documentation)
- prompt.md - Project requirements and hints
- Sample data files - SDF, GRO, ITP, MOL2 formats

### Tertiary (Context)
- MDAnalysis documentation - For trajectory handling patterns
- GROMACS SASA calculation - For validation of relative SASA values

---

## Metadata

**Confidence breakdown:**
- Current workflow analysis: HIGH (code fully reviewed)
- SASA algorithm: HIGH (algorithm from quick_sasa.py verified)
- Integration approach: MEDIUM (implementation strategy proposed, needs validation)
- Performance estimates: LOW (no benchmark data available)

**Research date:** March 26, 2026  
**Valid until:** Algorithm is stable; only GROMACS tool changes would affect this research