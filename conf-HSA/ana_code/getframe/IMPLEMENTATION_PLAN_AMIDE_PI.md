# Implementation Plan: Amide-π Interaction Detection

## Overview

Add amide-π interaction detection to `detect_pi_interactions.py`, along with statistics output across all frames.

---

## 1. New Data Structures

### 1.1 Side-chain Amide Definitions

```python
SIDECHAIN_AMIDES = {
    'ASN': {
        'atoms': {'N': 'ND2', 'H': ['HD21', 'HD22'], 'O': 'OD1'},
        'description': 'Asparagine side-chain amide'
    },
    'GLN': {
        'atoms': {'N': 'NE2', 'H': ['HE21', 'HE22'], 'O': 'OE1'},
        'description': 'Glutamine side-chain amide'
    }
}
```

### 1.2 Geometry Parameters

```python
AMIDE_PI_DISTANCE_CUTOFF = 4.0      # H or O to ring center (Å)
AMIDE_PI_ANGLE_MIN = 150.0          # N-H...ring angle minimum (degrees)
```

### 1.3 Frame Statistics Structure

```python
frame_stats = {
    'parallel_stacking': [],
    't_shaped_stacking': [],
    'cation_pi': [],
    'amide_pi': [],           # NEW
    'total': []
}
```

---

## 2. New Functions

### 2.1 `detect_sidechain_amides(universe)`

**Purpose**: Find ASN/GLN side-chain amide groups in protein.

**Input**: MDAnalysis Universe

**Output**: 
```python
[
    {
        'resname': 'ASN',
        'resid': 45,
        'N_index': 712,      # ND2 atom index
        'H_indices': [713, 714],  # HD21, HD22 indices
        'O_index': 711,      # OD1 atom index
        'source': 'protein'
    },
    ...
]
```

**Logic**:
1. Select all ASN and GLN residues
2. For each residue, find N, H atoms, and O atoms by name
3. Verify all atoms exist
4. Store atom indices for later use

**Estimated Cost**: O(n_residues) - trivial

---

### 2.2 `detect_amide_pi_nh(amide, ring, positions, geo_ring)`

**Purpose**: Detect NH-π interaction (amide N-H → ring).

**Inputs**:
- `amide`: dict with N_index, H_indices
- `ring`: ring information dict
- `positions`: all atom positions array
- `geo_ring`: ring geometry (center, normal, radius)

**Output**: `(is_interaction: bool, details: dict)`

**Logic**:
1. For each H atom in amide['H_indices']:
   - Get H position
   - Calculate distance from H to ring center
   - If distance ≤ 4.0 Å:
     - Calculate N-H vector and H-to-ring-center vector
     - Calculate angle N-H...ring_center
     - If angle ≥ 150°: interaction found
2. Return details with distance, angle, atom info

**Estimated Cost**: O(n_H_atoms) per ring - very fast

---

### 2.3 `detect_amide_pi_npi(amide, ring, positions, geo_ring)`

**Purpose**: Detect n-π* interaction (amide O lone pair → ring).

**Inputs**: Same as `detect_amide_pi_nh`

**Output**: `(is_interaction: bool, details: dict)`

**Logic**:
1. Get O position from amide['O_index']
2. Calculate distance from O to ring center
3. If distance ≤ 4.0 Å: interaction found
4. Return details with distance, atom info

**Note**: No angle criterion for n-π* per research document.

**Estimated Cost**: O(1) per ring - trivial

---

### 2.4 `calculate_statistics(all_frame_data)`

**Purpose**: Calculate mean and max statistics across all frames.

**Input**: List of frame data dicts from `process_trajectory`

**Output**:
```python
{
    'parallel_stacking_mean': 2.34,
    'parallel_stacking_max': 5,
    't_shaped_stacking_mean': 1.12,
    't_shaped_stacking_max': 3,
    'cation_pi_mean': 0.45,
    'cation_pi_max': 2,
    'amide_pi_mean': 1.23,
    'amide_pi_max': 4,
    'total_mean': 5.14,
    'total_max': 12
}
```

**Logic**:
1. Extract counts for each interaction type across all frames
2. Calculate mean and max for each
3. Calculate total mean and max

**Estimated Cost**: O(n_frames) - trivial

---

## 3. Modified Functions

### 3.1 `analyze_frame()` - Modified

**Changes**:
1. Add parameter: `sidechain_amides`
2. Add detection loop for amide-π interactions:
   - For each amide and each ligand ring:
     - Check NH-π interaction
     - Check n-π* interaction
   - Store in `results['amide_pi']`
3. Update total count to include amide_pi
4. Track residues involved in amide-π

**New output structure**:
```python
results = {
    'parallel_stacking': [...],
    't_shaped_stacking': [...],
    'cation_pi': [...],
    'amide_pi': [           # NEW
        {
            'type': 'NH_pi',
            'distance': 3.5,
            'angle': 155.0,
            'amide': 'ASN45',
            'ring': 'PHE123'
        },
        {
            'type': 'n_pi_star',
            'distance': 3.8,
            'amide': 'GLN67',
            'ring': 'TRP89'
        }
    ],
    'residues_involved': set(),
    'total': 8
}
```

---

### 3.2 `process_trajectory()` - Modified

**Changes**:
1. Call `detect_sidechain_amides()` after detecting charge centers
2. Pass `sidechain_amides` to `analyze_frame()`
3. Collect all frame data for statistics
4. Call `calculate_statistics()` at end
5. Return statistics along with other data

**New return signature**:
```python
return (best_frame_idx, best_results, all_frame_data, 
        protein_rings, ligand_rings, u, statistics)
```

---

### 3.3 `extract_frame_with_dummy_atoms()` - Modified

**Changes**:
1. Add REMARK line for amide-π count
2. Add dummy atoms at ring centers involved in amide-π (if not already present)

---

### 3.4 `main()` - Modified

**Changes**:
1. Update result collection to include statistics
2. Write `pi_interactions_all_frames.csv` (optional)
3. Update main CSV columns to include statistics

---

## 4. CSV Output Changes

### 4.1 Main CSV: `pi_interactions_summary.csv`

**Current columns**:
- system, trial_id, frame_num, global_frame
- parallel_stacking, t_shaped_stacking, cation_pi, total
- residues

**New columns** (best frame only, per system):
```
system, best_frame, 
parallel_stacking, t_shaped_stacking, cation_pi, amide_pi, total,
parallel_stacking_mean, parallel_stacking_max,
t_shaped_stacking_mean, t_shaped_stacking_max,
cation_pi_mean, cation_pi_max,
amide_pi_mean, amide_pi_max,
total_mean, total_max,
residues
```

### 4.2 Optional: `pi_interactions_all_frames.csv`

**Columns**:
```
system, frame_num, global_frame,
parallel_stacking, t_shaped_stacking, cation_pi, amide_pi, total,
residues
```

---

## 5. Implementation Order

1. **Add constants and data structures** (lines 27-51)
   - `SIDECHAIN_AMIDES` dictionary
   - `AMIDE_PI_DISTANCE_CUTOFF`, `AMIDE_PI_ANGLE_MIN`

2. **Implement `detect_sidechain_amides()`** (after line 307)
   - New function, ~30 lines

3. **Implement amide-π detection functions** (after above)
   - `detect_amide_pi_nh()` - ~25 lines
   - `detect_amide_pi_npi()` - ~15 lines

4. **Modify `analyze_frame()`** (lines 461-542)
   - Add amide-π detection loop
   - Update total calculation
   - ~20 additional lines

5. **Implement `calculate_statistics()`** (after analyze_frame)
   - New function, ~25 lines

6. **Modify `process_trajectory()`** (lines 545-606)
   - Add amide detection call
   - Add statistics calculation
   - ~10 additional lines

7. **Modify `extract_frame_with_dummy_atoms()`** (lines 609-676)
   - Add amide-π REMARK
   - ~5 additional lines

8. **Modify `main()`** (lines 697-758)
   - Update CSV handling
   - Add all-frames CSV option
   - ~20 additional lines

---

## 6. Estimated Computational Cost

### Per-frame analysis:

| Operation | Current Cost | Added Cost | Notes |
|-----------|-------------|------------|-------|
| Ring geometry | O(n_rings) | 0 | Unchanged |
| Parallel stacking | O(n_protein_rings × n_ligand_rings) | 0 | Unchanged |
| T-shaped stacking | O(n_protein_rings × n_ligand_rings) | 0 | Unchanged |
| Cation-π | O(n_charges × n_rings) | 0 | Unchanged |
| **Amide-π (NEW)** | - | O(n_amides × n_ligand_rings × 3) | 3 = 2 H + 1 O checks |

**Typical values**:
- n_protein_rings: ~10-20
- n_ligand_rings: ~1-3
- n_charges: ~5-10
- n_amides: ~5-15 (ASN + GLN residues)

**Per-frame added cost**: ~15-45 distance calculations + ~30 angle calculations
**Added cost per 1000-frame trajectory**: Negligible (< 1 second)

### Statistics calculation:
- O(n_frames) once per trajectory
- Negligible

---

## 7. Potential Issues and Mitigations

### 7.1 Hydrogen Atom Detection

**Issue**: Some MD trajectories may not have explicit hydrogens.

**Detection**: Check if H atoms exist in `detect_sidechain_amides()`

**Mitigation**: 
```python
if len(h_atoms) == 0:
    print(f"  Warning: No H atoms found for {resname}{resid} amide")
    # Skip NH-π detection, still do n-π*
```

### 7.2 Protonation State

**Issue**: ASN/GLN may have different protonation states.

**Detection**: Check for expected atom names

**Mitigation**: Accept alternative naming (HD2, HE2 as alternatives to HD21/HD22)

### 7.3 Backbone vs Side-chain

**Issue**: Backbone amides might be accidentally detected.

**Mitigation**: Use specific side-chain atom names (ND2, NE2, OD1, OE1) which don't exist in backbone

### 7.4 Duplicate Detection

**Issue**: Same amide might interact with multiple ring positions on the same aromatic.

**Mitigation**: This is valid - each distinct interaction counts separately

### 7.5 Double-counting

**Issue**: NH-π and n-π* from same amide to same ring should count as 2 interactions.

**Decision**: Correct behavior - they are distinct interaction types

---

## 8. Testing Strategy

### 8.1 Unit Tests

1. `test_detect_sidechain_amides()`:
   - Create mock universe with ASN/GLN
   - Verify correct atom indices returned

2. `test_detect_amide_pi_nh()`:
   - Test with H at 3.5 Å, angle 155° → should detect
   - Test with H at 4.5 Å → should not detect
   - Test with H at 3.5 Å, angle 140° → should not detect

3. `test_detect_amide_pi_npi()`:
   - Test with O at 3.5 Å → should detect
   - Test with O at 4.5 Å → should not detect

### 8.2 Integration Tests

1. Run on single known system
2. Verify amide-π interactions match expected geometry
3. Verify statistics calculated correctly

---

## 9. Summary of Changes

| Component | Lines Changed | New Lines |
|-----------|---------------|-----------|
| Constants | 5 | 8 |
| `detect_sidechain_amides()` | 0 | ~30 |
| `detect_amide_pi_nh()` | 0 | ~25 |
| `detect_amide_pi_npi()` | 0 | ~15 |
| `analyze_frame()` | ~20 | ~20 |
| `calculate_statistics()` | 0 | ~25 |
| `process_trajectory()` | ~10 | ~10 |
| `extract_frame_with_dummy_atoms()` | ~5 | ~5 |
| `main()` | ~20 | ~20 |
| **Total** | ~60 | ~158 |

**Net increase**: ~100-160 lines of code

---

## 10. Timeline Estimate

| Task | Time |
|------|------|
| Add constants and data structures | 5 min |
| Implement `detect_sidechain_amides()` | 15 min |
| Implement amide-π detection functions | 20 min |
| Modify `analyze_frame()` | 15 min |
| Implement `calculate_statistics()` | 10 min |
| Modify `process_trajectory()` | 10 min |
| Modify `extract_frame_with_dummy_atoms()` | 5 min |
| Modify `main()` | 15 min |
| Testing | 20 min |
| **Total** | ~2 hours |

---

## 11. Approval Checklist

- [ ] Side-chain amides only (ASN, GLN) - no backbone
- [ ] NH-π detection with distance ≤ 4.0 Å and angle ≥ 150°
- [ ] n-π* detection with distance ≤ 4.0 Å
- [ ] Statistics: mean and max per interaction type
- [ ] Main CSV: best frame with statistics columns
- [ ] Optional: all-frames CSV
- [ ] Existing functionality preserved

---

**Document Created**: May 12, 2026
**Status**: Ready for Implementation
