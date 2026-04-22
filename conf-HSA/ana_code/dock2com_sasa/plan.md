# Implementation Plan: SASA-Based Rescoring for dock2com_1.py

**Created:** March 26, 2026  
**Context:** Add SASA (Solvent Accessible Surface Area) as an additional criterion for selecting optimal docking poses  
**Goal:** Lower SASA indicates better ligand burial within the binding site

---

## Overview

The current `dock2com_1.py` selects the best docking pose based solely on docking scores (e.g., `minimizedAffinity`, `CNNscore`). This plan describes two alternative approaches to integrate SASA-based rescoring, enabling selection of poses that are better buried within the binding site.

### Key Integration Points

| Component | Current State | SASA Enhancement |
|-----------|---------------|------------------|
| Model Selection | `select_best_across_sdfiles()` selects by single metric | Add SASA calculation per pose |
| Receptor Loading | Auto-detected from SDF filename | Use for SASA context |
| CLI Interface | `--metric` for selection metric | Add SASA-related flags |
| Output | Docking scores printed | Include SASA value in summary |

### Core Algorithm

The `quick_sasa.py` provides a fast VDW+Gaussian approximation:

```
SASA_atom = 4π(r_vdw + r_probe)² × exp(-λ × neighbor_sum)
```

- **Speed:** ~0.7ms per pose (suitable for hundreds of poses)
- **Accuracy:** 10-20% vs. probe methods (acceptable for relative comparisons)
- **Output:** SASA in nm² (lower = more buried = better)

---

## Plan A: Two-Stage Selection (Primary Recommendation)

### Overview

A two-pass approach that first filters by docking score, then re-ranks by SASA. This preserves the value of docking scores while adding SASA as a secondary criterion.

**Philosophy:** Docking scores capture binding affinity predictions; SASA captures burial quality. Use both.

### Architecture

```
Stage 1: Docking Score Filter
┌─────────────────────────────────────────────────────┐
│  All SDF poses (N poses across M files)             │
│                     │                               │
│                     ▼                               │
│  Sort by docking metric (e.g., minimizedAffinity)   │
│                     │                               │
│                     ▼                               │
│  Keep top-K candidates (configurable, default=10)   │
└─────────────────────────────────────────────────────┘
                        │
                        ▼
Stage 2: SASA Rescoring
┌─────────────────────────────────────────────────────┐
│  For each candidate pose:                           │
│    - Load matching receptor coordinates             │
│    - Calculate ligand SASA in receptor context      │
│    - Store SASA value                               │
│                                                     │
│  Re-rank candidates by SASA (lower = better)        │
│                                                     │
│  Return best pose by SASA                           │
└─────────────────────────────────────────────────────┘
```

### Step-by-Step Implementation Tasks

#### Task 1: Create SASA Calculation Module

**File:** `dock2com_sasa.py` (new module in same directory)

**Purpose:** Extract and adapt SASA calculation logic from `quick_sasa.py` for use with parsed SDF data

**Dependencies:** None (can run in parallel with Task 2)

**Implementation:**
```python
# Key function to implement:
def calculate_pose_sasa(ligand_atoms, receptor_atoms, probe_radius=1.4):
    """
    Calculate total ligand SASA in context of receptor.
    
    Parameters
    ----------
    ligand_atoms : list of dict
        [{"el": "C", "x": 1.0, "y": 2.0, "z": 3.0}, ...]
        Coordinates in Angstroms (from SDF)
    receptor_atoms : list of dict
        Same format, coordinates from GRO (convert nm -> Å)
    probe_radius : float
        Water probe radius in Angstroms (default 1.4)
    
    Returns
    -------
    float
        Total ligand SASA in nm² (lower = more buried)
    """
```

**Key adaptations needed:**
1. Work with parsed atom lists (not MDAnalysis AtomGroups)
2. Handle coordinate unit conversion (GRO uses nm, SDF uses Å)
3. Extract VDW radii table and Gaussian parameters from `quick_sasa.py`
4. Use scipy.spatial.cKDTree for efficient neighbor search

**VDW Radii Table (copy from quick_sasa.py):**
```python
VDW_RADII_ANG = {
    'H': 1.20, 'C': 1.70, 'N': 1.55, 'O': 1.52, 
    'S': 1.80, 'Cl': 1.75, 'Fe': 2.00, 'Zn': 1.39,
    'Eu': 2.05,  # For lanthanide binding proteins
    # ... full table
}
DEFAULT_VDW_RADIUS = 1.70  # Carbon fallback
```

**Verification:**
- Unit test: Calculate SASA for benzene alone → ~100-120 Å²
- Unit test: Calculate SASA for buried ligand → < 50 Å²
- Integration test: Compare output with `quick_sasa.py` for same structure

---

#### Task 2: Add GRO Parser Function

**File:** `dock2com_1.py`

**Purpose:** Parse receptor GRO file to extract atom coordinates and elements for SASA calculation

**Dependencies:** None (can run in parallel with Task 1)

**Implementation:**
```python
def parse_gro_for_sasa(gro_path):
    """
    Parse GRO file and extract atoms for SASA calculation.
    
    Parameters
    ----------
    gro_path : str
        Path to receptor GRO file
    
    Returns
    -------
    list of dict
        [{"el": "C", "x": 10.0, "y": 20.0, "z": 30.0}, ...]
        Coordinates converted to Angstroms
        Element derived from atom name (not stored in GRO)
    """
```

**Element inference logic:**
```python
def infer_element_from_atom_name(atom_name):
    """
    Infer element from GROMACS atom name.
    
    Examples:
        "CA" -> "C" (alpha carbon)
        "NH1" -> "N"
        "HZ" -> "H"
        "EU" -> "Eu" (Europium)
    """
    # Common patterns:
    # - First 1-2 chars are element
    # - Handle special cases: CA, CB, etc. -> Carbon
    # - Metals: FE -> Fe, ZN -> Zn, etc.
```

**Verification:**
- Unit test: Parse known GRO file, verify atom count
- Unit test: Verify coordinate conversion (nm × 10 = Å)
- Unit test: Verify element inference for common atom names

---

#### Task 3: Implement Two-Stage Selection Function

**File:** `dock2com_1.py`

**Purpose:** Replace or augment `select_best_across_sdfiles()` with SASA rescoring

**Dependencies:** Task 1 and Task 2 must complete

**Implementation:**
```python
def select_best_with_sasa_rescore(
    sdf_paths,
    metric=DEFAULT_METRIC,
    top_k=10,
    rec_gro_pattern=DEFAULT_REC_GRO_PATTERN,
    rec_dir=None,
    probe_radius=1.4,
):
    """
    Two-stage selection: filter by docking score, then rescore by SASA.
    
    Parameters
    ----------
    sdf_paths : list
        List of SDF file paths
    metric : str
        Primary docking metric for Stage 1 filtering
    top_k : int
        Number of top candidates to pass to Stage 2 (default: 10)
    rec_gro_pattern : str
        Pattern for auto-detecting receptor GRO files
    rec_dir : str, optional
        Directory containing receptor files
    probe_radius : float
        Water probe radius for SASA calculation
    
    Returns
    -------
    best_model : dict
        Selected model with SASA in scores
    best_file : str
        Path to SDF file containing best model
    best_model_num : int
        Model number within SDF file
    """
```

**Algorithm:**
```python
# Stage 1: Collect and sort all candidates
candidates = []  # (model, sdf_path, score)
for sdf_path in sdf_paths:
    models = parse_sdf(sdf_path)
    for model in models:
        if metric in model["scores"]:
            score = model["scores"][metric]
            candidates.append((model, sdf_path, score))

# Sort by docking metric
reverse = metric not in LOWER_IS_BETTER
candidates.sort(key=lambda x: x[2], reverse=reverse)

# Keep top-K
top_candidates = candidates[:top_k]

# Stage 2: Calculate SASA for each candidate
for model, sdf_path, _ in top_candidates:
    rec_gro, _ = derive_receptor_gro_from_sdf(sdf_path, rec_gro_pattern, rec_dir)
    if os.path.exists(rec_gro):
        receptor_atoms = parse_gro_for_sasa(rec_gro)
        ligand_atoms = model["atoms"]
        sasa = calculate_pose_sasa(ligand_atoms, receptor_atoms, probe_radius)
        model["scores"]["sasa"] = sasa
    else:
        model["scores"]["sasa"] = float('inf')  # Penalty for missing receptor

# Sort by SASA (lower is better)
top_candidates.sort(key=lambda x: x[0]["scores"]["sasa"])

# Return best
best_model, best_file, _ = top_candidates[0]
return best_model, best_file, best_model["model_num"]
```

**Verification:**
- Integration test: Run on test data, verify two-stage selection works
- Integration test: Verify SASA values are reasonable (0-200 nm² range)
- Test: Handle missing receptor GRO gracefully (warning + fallback)

---

#### Task 4: Add CLI Arguments

**File:** `dock2com_1.py`

**Purpose:** Add new command-line flags for SASA rescoring

**Dependencies:** Task 3 must complete

**Implementation:**
```python
# Add to build_parser() in "Model Selection" group:

sel_group.add_argument(
    "--sasa-rescore",
    action="store_true",
    help="Enable two-stage selection: filter by docking metric, then rescore by SASA"
)

sel_group.add_argument(
    "--sasa-top-k",
    type=int,
    default=10,
    metavar="N",
    help="Number of top candidates to consider in SASA rescoring (default: 10)"
)

sel_group.add_argument(
    "--probe-radius",
    type=float,
    default=1.4,
    metavar="R",
    help="Solvent probe radius in Angstroms for SASA calculation (default: 1.4 for water)"
)
```

**CLI Usage Examples:**
```bash
# Original behavior (unchanged)
python dock2com_1.py -i lig.itp -s hsa*.sdf -r topol.top

# With SASA rescoring (two-stage)
python dock2com_1.py -i lig.itp -s hsa*.sdf -r topol.top --sasa-rescore

# Custom top-K and probe radius
python dock2com_1.py -i lig.itp -s hsa*.sdf -r topol.top --sasa-rescore --sasa-top-k 20 --probe-radius 1.5
```

**Verification:**
- Test: `--help` shows new arguments
- Test: Backward compatibility (no flags = original behavior)
- Test: `--sasa-rescore` triggers two-stage selection

---

#### Task 5: Update Main Workflow

**File:** `dock2com_1.py`

**Purpose:** Integrate SASA rescoring into main() function

**Dependencies:** Tasks 3 and 4 must complete

**Implementation:**
```python
def main():
    # ... existing setup ...
    
    print("=" * 60)
    if config.sasa_rescore:
        print("STEP 1: Two-stage selection (docking + SASA)...")
    else:
        print("STEP 1: Selecting best docking model...")
    print("=" * 60)
    
    if config.sasa_rescore:
        best_model, best_file, best_model_num = select_best_with_sasa_rescore(
            config.sdf_paths,
            metric=config.metric,
            top_k=config.sasa_top_k,
            rec_gro_pattern=config.rec_gro_pattern,
            rec_dir=config.rec_dir,
            probe_radius=config.probe_radius,
        )
    else:
        best_model, best_file, best_model_num = select_best_across_sdfiles(
            config.sdf_paths, metric=config.metric
        )
    
    # ... rest of workflow ...
    
    # Update summary output
    print(f"  Scores: minimizedAffinity={...} CNNscore={...} CNNaffinity={...}")
    if 'sasa' in best_model["scores"]:
        print(f"          SASA={best_model['scores']['sasa']:.4f} nm²")
```

**Verification:**
- End-to-end test: Run full workflow with and without `--sasa-rescore`
- Verify output files are created correctly in both modes

---

#### Task 6: Add Config Class Updates

**File:** `dock2com_1.py`

**Purpose:** Store new CLI arguments in Config class

**Dependencies:** Task 4 must complete

**Implementation:**
```python
class Config:
    def __init__(self, args):
        # ... existing fields ...
        self.sasa_rescore = args.sasa_rescore
        self.sasa_top_k = args.sasa_top_k
        self.probe_radius = args.probe_radius
```

**Verification:**
- Unit test: Config correctly stores new arguments

---

### Dependencies Graph

```
Task 1 (SASA module) ─────┐
                          │
Task 2 (GRO parser) ──────┼──► Task 3 (Two-stage selection) ──► Task 5 (Main workflow)
                          │                                          │
Task 4 (CLI args) ────────┘                                          │
       │                                                              │
       └──────────────────────────────────────────────────────────────┘
                                                                          │
Task 6 (Config) ◄─────────────────────────────────────────────────────────┘
```

**Parallel execution:**
- Tasks 1 and 2 can run in parallel (no dependencies)
- Tasks 4 and 6 can be combined (small, related changes)

**Estimated effort:** 3-4 hours total

---

### Testing Approach

#### Unit Tests

| Test | Purpose | Expected Result |
|------|---------|-----------------|
| `test_calculate_pose_sasa_isolated()` | SASA of isolated ligand | High value (> 100 nm² for typical drug) |
| `test_calculate_pose_sasa_buried()` | SASA of buried ligand | Low value (< 50 nm²) |
| `test_parse_gro_for_sasa()` | GRO parsing | Correct atom count, element inference |
| `test_coordinate_conversion()` | nm to Å conversion | Coordinates × 10 |

#### Integration Tests

| Test | Purpose | Expected Result |
|------|---------|-----------------|
| `test_two_stage_selection()` | Full two-stage selection | Best pose has lowest SASA among top-K |
| `test_missing_receptor_graceful()` | Missing receptor GRO | Warning + penalty SASA, not crash |
| `test_cli_sasa_rescore_flag()` | CLI flag behavior | Two-stage triggered when flag set |

#### End-to-End Test

```bash
# Run on test data
python dock2com_1.py -i lig.itp -s hsa*-dzp.sdf -r topol.top --sasa-rescore

# Verify:
# 1. Output files created (best.gro, com.gro, rec.itp, sys.top)
# 2. Summary shows SASA value
# 3. Selected pose differs from non-SASA selection (if SASA affects ranking)
```

---

### CLI Interface Changes Summary

| Flag | Default | Purpose |
|------|---------|---------|
| `--sasa-rescore` | False | Enable two-stage selection with SASA |
| `--sasa-top-k N` | 10 | Number of candidates for Stage 2 |
| `--probe-radius R` | 1.4 | Solvent probe radius (Å) |

**Backward compatibility:** All new flags are optional. Default behavior unchanged.

---

## Plan B: Single-Pass SASA Metric (Alternative)

### Overview

Add "sasa" as a first-class metric alongside `minimizedAffinity`, `CNNscore`, etc. This approach calculates SASA for all poses in a single pass and allows selecting by SASA directly.

**Philosophy:** Simplicity. One metric, one pass, straightforward CLI.

### Architecture

```
Single Pass SASA Calculation
┌─────────────────────────────────────────────────────┐
│  For each SDF pose:                                 │
│    1. Parse coordinates                              │
│    2. Load matching receptor                         │
│    3. Calculate SASA                                 │
│    4. Store in model["scores"]["sasa"]              │
│                                                     │
│  Select best model by chosen metric                  │
│  (sasa, minimizedAffinity, CNNscore, etc.)          │
└─────────────────────────────────────────────────────┘
```

### Step-by-Step Implementation Tasks

#### Task 1: Create SASA Calculation Module (Same as Plan A Task 1)

See Plan A Task 1. Identical implementation.

---

#### Task 2: Add GRO Parser Function (Same as Plan A Task 2)

See Plan A Task 2. Identical implementation.

---

#### Task 3: Modify select_best_across_sdfiles()

**File:** `dock2com_1.py`

**Purpose:** Extend existing selection function to calculate SASA when needed

**Dependencies:** Tasks 1 and 2 must complete

**Implementation:**
```python
# Add LOWER_IS_BETTER constant
LOWER_IS_BETTER = {"minimizedAffinity", "sasa"}  # Added "sasa"

def select_best_across_sdfiles(
    sdf_paths,
    metric=DEFAULT_METRIC,
    rec_gro_pattern=None,  # New parameter
    rec_dir=None,          # New parameter
    probe_radius=1.4,      # New parameter
):
    """
    Select best model across SDF files, optionally calculating SASA.
    
    If metric == "sasa", calculates SASA for all poses.
    Otherwise, calculates SASA only if --list-models is used with sasa in metrics.
    """
    best_model = None
    best_score = None
    best_file = None
    reverse = metric not in LOWER_IS_BETTER
    
    for sdf_path in sdf_paths:
        print(f"Parsing SDF:        {sdf_path}")
        models = parse_sdf(sdf_path)
        print(f"  {len(models)} models found")
        
        # Auto-detect receptor for this SDF file
        rec_gro, prefix = derive_receptor_gro_from_sdf(sdf_path, rec_gro_pattern, rec_dir)
        receptor_atoms = None
        if os.path.exists(rec_gro):
            receptor_atoms = parse_gro_for_sasa(rec_gro)
        
        for model in models:
            # Calculate SASA if needed
            if metric == "sasa" or (metric not in model["scores"]):
                if receptor_atoms:
                    sasa = calculate_pose_sasa(model["atoms"], receptor_atoms, probe_radius)
                    model["scores"]["sasa"] = sasa
                else:
                    model["scores"]["sasa"] = float('nan')
            
            if metric not in model["scores"]:
                continue
            
            score = model["scores"][metric]
            # ... existing comparison logic ...
    
    return best_model, best_file, best_model["model_num"]
```

**Key differences from Plan A:**
- Calculates SASA on-demand, not in a separate stage
- Uses existing selection logic (minimal changes)
- "sasa" is just another metric option

---

#### Task 4: Update CLI and Main

**File:** `dock2com_1.py`

**Purpose:** Add minimal CLI changes for SASA support

**Dependencies:** Task 3 must complete

**Implementation:**
```python
# Add to build_parser() - minimal changes only
sel_group.add_argument(
    "--probe-radius",
    type=float,
    default=1.4,
    metavar="R",
    help="Solvent probe radius for SASA calculation (default: 1.4)"
)

# Update help text for --metric to include "sasa"
sel_group.add_argument(
    "--metric",
    default=DEFAULT_METRIC,
    metavar="NAME",
    help="Score field for model selection. Options: minimizedAffinity, CNNscore, CNNaffinity, sasa"
)

# In main():
best_model, best_file, best_model_num = select_best_across_sdfiles(
    config.sdf_paths,
    metric=config.metric,
    rec_gro_pattern=config.rec_gro_pattern,
    rec_dir=config.rec_dir,
    probe_radius=config.probe_radius,
)

# Update list_models output to show sasa
if 'sasa' in best_model["scores"]:
    print(f"          SASA={best_model['scores']['sasa']:.4f} nm²")
```

**CLI Usage:**
```bash
# Select by SASA directly
python dock2com_1.py -i lig.itp -s hsa*.sdf -r topol.top --metric sasa

# List all metrics including SASA
python dock2com_1.py -i lig.itp -s hsa*.sdf --list-models
```

---

### Dependencies Graph

```
Task 1 (SASA module) ─────┐
                          │
Task 2 (GRO parser) ──────┼──► Task 3 (Modify selection) ──► Task 4 (CLI/Main)
                          │
```

**Parallel execution:**
- Tasks 1 and 2 can run in parallel

**Estimated effort:** 2-3 hours total (simpler than Plan A)

---

### Testing Approach

#### Unit Tests

Same as Plan A (Tasks 1 and 2 are identical).

#### Integration Tests

| Test | Purpose | Expected Result |
|------|---------|-----------------|
| `test_metric_sasa()` | Select by SASA | Lowest SASA pose selected |
| `test_list_models_with_sasa()` | List models shows SASA | SASA column appears in output |

#### End-to-End Test

```bash
# Select by SASA directly
python dock2com_1.py -i lig.itp -s hsa*-dzp.sdf -r topol.top --metric sasa

# Verify SASA selection worked
```

---

### CLI Interface Changes Summary

| Flag | Default | Change |
|------|---------|--------|
| `--metric` | minimizedAffinity | Now accepts "sasa" |
| `--probe-radius` | 1.4 | New (for SASA calculation) |

**Backward compatibility:** Default behavior unchanged. New metric is opt-in.

---

## Plan C: Composite Score (Third Option - Brief)

### Overview

Create a weighted combination of docking score and SASA:

```
composite_score = α × normalized_affinity + β × normalized_sasa
```

Where:
- `normalized_affinity` = (score - min) / (max - min)
- `normalized_sasa` = (sasa - min) / (max - min)  
- α, β are configurable weights (default: α=0.7, β=0.3)

**Pros:**
- Single selection criterion
- Tunable balance between affinity and burial

**Cons:**
- Weight selection requires validation
- Normalization can be sensitive to outliers
- Less intuitive than two-stage approach

**Recommendation:** Do not implement initially. Consider after Plan A is validated and user feedback suggests composite scoring is desired.

---

## Comparison and Recommendation

### Feature Comparison

| Feature | Plan A (Two-Stage) | Plan B (Single Metric) | Plan C (Composite) |
|---------|-------------------|------------------------|-------------------|
| **Complexity** | Medium | Low | High |
| **Lines of code** | ~200 | ~100 | ~250 |
| **CLI flags** | 3 new | 1 new | 4 new |
| **Backward compatible** | Yes | Yes | Yes |
| **Flexibility** | High (tunable top-K) | Medium | High (tunable weights) |
| **Intuitive** | Yes (filter → rank) | Yes (metric option) | No (weights unclear) |
| **Performance** | ~K × SASA calc | ~N × SASA calc | ~N × SASA calc |
| **Validation needed** | Moderate | Low | High (weights) |

### Decision Matrix

| Criterion | Weight | Plan A | Plan B | Plan C |
|-----------|--------|--------|--------|--------|
| Implementation simplicity | 20% | 3/5 | 5/5 | 2/5 |
| Flexibility for future | 25% | 5/5 | 3/5 | 4/5 |
| Intuitiveness | 20% | 5/5 | 5/5 | 2/5 |
| Performance | 15% | 5/5 | 3/5 | 3/5 |
| Maintenance burden | 10% | 4/5 | 5/5 | 3/5 |
| User control | 10% | 4/5 | 3/5 | 5/5 |
| **Weighted Score** | | **4.2/5** | **3.9/5** | **2.9/5** |

### Recommendation

**Implement Plan A (Two-Stage Selection) as primary approach.**

**Rationale:**
1. **Flexibility:** Two-stage approach allows users to tune the trade-off between docking score quality and burial quality via `--sasa-top-k`
2. **Performance:** Only calculates SASA for top-K candidates (typically 10), not all N poses
3. **Intuitiveness:** Clear mental model—"first filter by docking, then pick most buried"
4. **Extensibility:** Easy to add more stages or criteria later
5. **Validation:** No weight tuning required (unlike Plan C)

**Fallback consideration:**
- If users consistently set `--sasa-top-k` to a large value (approaching N), consider adding Plan B's `--metric sasa` as a shortcut
- Plan C should only be considered after user feedback indicates demand for composite scoring

---

## Implementation Order

### Recommended Sequence

```
Phase 1: Core Infrastructure (Day 1)
├── Task 1: SASA calculation module
├── Task 2: GRO parser function
└── Unit tests for both

Phase 2: Integration (Day 1-2)
├── Task 3: Two-stage selection function
├── Task 4: CLI arguments
├── Task 5: Main workflow updates
├── Task 6: Config class
└── Integration tests

Phase 3: Validation (Day 2)
├── End-to-end testing
├── Documentation updates
└── Example usage in script comments
```

### Risk Mitigation

| Risk | Mitigation |
|------|------------|
| SASA calculation too slow | Profile and optimize; add optional flag to disable |
| Element inference fails | Add fallback to carbon radius with warning |
| Missing receptor GRO | Graceful degradation with warning, use docking score only |
| Coordinate unit mismatch | Explicit conversion in SASA function, unit tests |

---

## Output Files After Implementation

```
dock2com_sasa/
├── dock2com_1.py          # Modified (new functions, CLI changes)
├── dock2com_sasa.py       # NEW (SASA calculation module)
├── quick_sasa.py          # Unchanged (reference implementation)
├── research.md            # Unchanged (research findings)
├── plan.md                # This file
└── tests/
    ├── test_sasa_calculation.py   # NEW
    ├── test_gro_parser.py         # NEW
    └── test_integration.py        # NEW
```

---

## Appendix: Code Skeleton

### dock2com_sasa.py (New Module)

```python
#!/usr/bin/env python3
"""
SASA calculation module for docking pose rescoring.

Provides calculate_pose_sasa() for computing ligand SASA in receptor context.
Uses fast VDW+Gaussian approximation from quick_sasa.py.
"""

import numpy as np
from scipy.spatial import cKDTree

# VDW radii (Bondi 1964, Mantina 2009)
VDW_RADII_ANG = {
    'H': 1.20, 'C': 1.70, 'N': 1.55, 'O': 1.52,
    'S': 1.80, 'Cl': 1.75, 'F': 1.47, 'P': 1.80,
    'Fe': 2.00, 'Zn': 1.39, 'Eu': 2.05, 'Br': 1.85,
    'I': 1.98, 'Mg': 1.73, 'Ca': 2.31, 'Cu': 1.40,
}
DEFAULT_VDW_RADIUS = 1.70  # Carbon

# Gaussian parameters (empirically tuned)
GAUSSIAN_LAMBDA = 0.5
GAUSSIAN_SIGMA = 1.5  # Angstroms
NEIGHBOR_CUTOFF = 10.0  # Angstroms


def calculate_pose_sasa(ligand_atoms, receptor_atoms, probe_radius=1.4):
    """
    Calculate total ligand SASA in context of receptor.
    
    Parameters
    ----------
    ligand_atoms : list of dict
        [{"el": "C", "x": 1.0, "y": 2.0, "z": 3.0}, ...]
        Coordinates in Angstroms
    receptor_atoms : list of dict
        Same format
    probe_radius : float
        Water probe radius (default: 1.4 Å)
    
    Returns
    -------
    float
        Total ligand SASA in nm² (lower = more buried)
    """
    # Combine ligand + receptor coordinates
    all_coords = np.array([
        [a["x"], a["y"], a["z"]] 
        for a in ligand_atoms + receptor_atoms
    ])
    
    # Build KDTree
    tree = cKDTree(all_coords)
    
    # Calculate per-atom SASA for ligand atoms
    ligand_coords = np.array([
        [a["x"], a["y"], a["z"]] 
        for a in ligand_atoms
    ])
    ligand_elements = [a["el"] for a in ligand_atoms]
    
    n_ligand = len(ligand_atoms)
    sasa_per_atom = np.zeros(n_ligand)
    
    for i, (coord, elem) in enumerate(zip(ligand_coords, ligand_elements)):
        r_vdw = VDW_RADII_ANG.get(elem, DEFAULT_VDW_RADIUS)
        r_total = r_vdw + probe_radius
        
        # Find neighbors within cutoff
        neighbors = tree.query_ball_point(coord, NEIGHBOR_CUTOFF)
        
        # Calculate neighbor sum (Gaussian weighted)
        neighbor_sum = 0.0
        for j in neighbors:
            if j == i:
                continue
            d = np.linalg.norm(all_coords[j] - coord)
            neighbor_sum += np.exp(-(d ** 2) / (2 * GAUSSIAN_SIGMA ** 2))
        
        # Exposure factor
        exposure = np.exp(-GAUSSIAN_LAMBDA * neighbor_sum)
        
        # SASA for this atom
        sasa_per_atom[i] = 4 * np.pi * (r_total ** 2) * exposure
    
    # Convert Å² to nm² and sum
    total_sasa_nm2 = np.sum(sasa_per_atom) / 100.0
    
    return total_sasa_nm2
```

### parse_gro_for_sasa() (Add to dock2com_1.py)

```python
def parse_gro_for_sasa(gro_path):
    """
    Parse GRO file and extract atoms for SASA calculation.
    
    Parameters
    ----------
    gro_path : str
        Path to GRO file
    
    Returns
    -------
    list of dict
        [{"el": "C", "x": 10.0, "y": 20.0, "z": 30.0}, ...]
        Coordinates in Angstroms (converted from nm)
    """
    atoms = []
    
    with open(gro_path) as f:
        lines = f.readlines()
    
    # Skip header (line 0) and atom count (line 1)
    # Box is last line
    for line in lines[2:-1]:
        if len(line) < 44:
            continue
        
        atom_name = line[10:15].strip()
        x_nm = float(line[20:28])
        y_nm = float(line[28:36])
        z_nm = float(line[36:44])
        
        # Convert nm to Angstroms
        x_ang = x_nm * 10.0
        y_ang = y_nm * 10.0
        z_ang = z_nm * 10.0
        
        # Infer element from atom name
        element = infer_element_from_atom_name(atom_name)
        
        atoms.append({
            "el": element,
            "x": x_ang,
            "y": y_ang,
            "z": z_ang,
        })
    
    return atoms


def infer_element_from_atom_name(atom_name):
    """
    Infer element from GROMACS atom name.
    
    Parameters
    ----------
    atom_name : str
        Atom name from GRO file (e.g., "CA", "NH1", "HZ")
    
    Returns
    -------
    str
        Element symbol (e.g., "C", "N", "H")
    """
    # Special cases
    special_map = {
        "EU": "Eu",
        "FE": "Fe",
        "ZN": "Zn",
        "MG": "Mg",
        "CA": "C",  # Alpha carbon -> Carbon
        "CB": "C",  # Beta carbon -> Carbon
        "CD": "C",  # Delta carbon -> Carbon
        "CG": "C",  # Gamma carbon -> Carbon
        "CE": "C",  # Epsilon carbon -> Carbon
        "CZ": "C",  # Zeta carbon -> Carbon
        "OH": "O",  # Hydroxyl oxygen
        "NH": "N",  # Amide nitrogen
        "OD": "O",  # Oxygen delta
        "OE": "O",  # Oxygen epsilon
        "OG": "O",  # Oxygen gamma
        "SD": "S",  # Sulfur delta
        "SG": "S",  # Sulfur gamma
    }
    
    upper = atom_name.upper()
    if upper in special_map:
        return special_map[upper]
    
    # General rule: first 1-2 characters
    if len(upper) >= 2 and upper[1].isdigit():
        return upper[0]  # H1, H2, C1, etc.
    
    if len(upper) >= 2:
        two_char = upper[:2]
        if two_char in VDW_RADII_ANG or two_char.capitalize() in VDW_RADII_ANG:
            return two_char.capitalize() if two_char[1].islower() else two_char
        return two_char[0]  # Fall back to first character
    
    return upper[0] if upper else "C"
```

---

## Summary

| Item | Value |
|------|-------|
| **Recommended Plan** | Plan A (Two-Stage Selection) |
| **Primary Benefit** | Combines docking accuracy with burial quality |
| **Estimated Effort** | 3-4 hours |
| **Backward Compatible** | Yes (opt-in via `--sasa-rescore`) |
| **New Files** | `dock2com_sasa.py`, test files |
| **Modified Files** | `dock2com_1.py` |
| **New CLI Flags** | `--sasa-rescore`, `--sasa-top-k`, `--probe-radius` |
