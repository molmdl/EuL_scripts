# Pi-Interaction Detection Plan

## Objective
Detect pi-interactions (pi-stacking, T-shaped stacking, cation-pi) in MD trajectories and extract frames with maximum pi-interactions for each system.

---

## Trajectory Structure

**Location**: `trj/[system]/fp/`

| File | Description |
|------|-------------|
| `v1.pdb` | Topology reference (9321 atoms, no solvent) |
| `v1.xtc` | Combined trajectory from all trials |
| `com.tpr` | TPR file (has solvent, incompatible with v1.xtc) |

**Frame Counts**:
- Expected: ~4000 frames per system (4 trials × ~1000 frames each)
- Actual: Varies by system (some trials may be shorter)
- Script handles any frame count automatically

**Frame to Trial Mapping**:
- Trial 0: frames 0-999
- Trial 1: frames 1000-1999
- Trial 2: frames 2000-2999
- Trial 3: frames 3000+
- Script calculates trial_id from global frame number

---

## Detection Methods

### 1. Pi-Stacking (Parallel)

**Geometric Criteria** (from BINANA):
- Distance between ring centers: < 7.5 Å
- Angle between ring normal vectors: < 30° (or > 150°)
- Ring center projection falls within ring disk (radius + 0.75 Å padding)

**Algorithm**:
```python
1. Get ring center and normal for each ring
2. Calculate center-to-center distance
3. Calculate angle between normals
4. Project ring center onto opposite ring plane
5. Check if projection is within ring disk
```

### 2. T-Shaped Stacking

**Geometric Criteria** (from BINANA):
- Distance between ring centers: < 7.5 Å
- Angle between ring normal vectors: 90° ± 30°
- Closest atom-atom distance: < 5.0 Å
- Ring center projection within ring disk

**Algorithm**:
```python
1. Check center distance and angle
2. Calculate minimum atom-atom distance between rings
3. Project center onto opposite ring plane
4. Check projection within ring disk
```

### 3. Cation-Pi

**Geometric Criteria** (from BINANA/PLIP):
- Distance between charge center and ring center: < 6.0 Å
- Charge center projection falls within ring disk

**Charge Centers**:
- **Protein**: ARG (guanidinium), LYS (amine), HIS (imidazole N)
- **Ligand**: Tertiary amines, quaternary ammonium, guanidine

---

## Implementation

### Aromatic Ring Detection

**Protein Rings** (residue-based):
| Residue | Ring Atoms |
|---------|-----------|
| PHE | CG, CD1, CD2, CE1, CE2, CZ |
| TYR | CG, CD1, CD2, CE1, CE2, CZ |
| HIS | CG, ND1, CD2, CE1, NE2 |
| TRP | Ring 1: CG, CD1, NE1, CE2, CD2, CE3 / Ring 2: CD2, CE2, CE3, CZ2, CZ3, CH2 |

**Ligand Rings** (connectivity-based):
- Use MDAnalysis bond information
- Identify 5 and 6-membered rings
- Check planarity (dihedral angles < 15°)

### Charge Center Detection

**Protein Charged Groups**:
| Residue | Charge Center | Atoms |
|---------|--------------|-------|
| ARG | Guanidinium C | CZ |
| LYS | Amine N | NZ |
| HIS | Imidazole | NE2 or ND1 |

**Ligand Charged Groups**:
- Tertiary/quaternary amines (N with 3-4 bonds)
- Guanidine groups
- Use MDAnalysis to identify by connectivity

---

## Workflow

```
1. Find all systems in trj/
2. For each system:
   a. Load v1.pdb + v1.xtc
   b. Detect aromatic rings in protein
   c. Detect aromatic rings in ligand
   d. Detect charge centers
   e. For each frame:
      - Count pi-stacking (protein ring - ligand ring)
      - Count T-shaped stacking
      - Count cation-pi (charge - ring)
   f. Find frame with maximum total interactions
   g. Extract frame with MDAnalysis
   h. Add dummy atoms at detected ring centers
   i. Save as PDB
3. Output summary CSV
```

---

## Output Files

### 1. PDB Files (Primary)

**Naming**: `[system]_best_pi.pdb`

**Contents**:
- Protein atoms
- Ligand atoms
- Dummy atoms at ring centers (resname: XPI, atom name: PI)

**PDB REMARK section**:
```
REMARK PI-INTERACTIONS
REMARK Parallel stacking: [count]
REMARK T-shaped stacking: [count]
REMARK Cation-pi: [count]
REMARK Total: [count]
REMARK Residues involved: [list]
```

### 2. CSV File (Secondary)

**Naming**: `pi_interactions_summary.csv`

**Columns**:
```
system              - System name
trial_id            - Trial ID (mmpbsa_0, mmpbsa_1, etc.)
frame_num           - Frame number within trial
global_frame        - Frame number in merged trajectory
parallel_stacking   - Count of parallel pi-stacking interactions
t_shaped_stacking   - Count of T-shaped stacking interactions
cation_pi           - Count of cation-pi interactions
total               - Total pi-interactions
```

---

## Dependencies

All from current conda environment:
- `numpy` - Vector calculations, geometry
- `scipy` - KDTree for distance searches
- `pandas` - Data handling, CSV output
- `MDAnalysis` - Trajectory loading, atom selection

No additional installations required.

---

## Script Structure

```
detect_pi_interactions.py
├── detect_aromatic_rings_protein()    # Residue-based ring detection
├── detect_aromatic_rings_ligand()     # Connectivity-based ring detection
├── detect_charge_centers()            # Identify charged groups
├── get_ring_geometry()                # Calculate center, normal, radius
├── detect_pi_stacking()               # Parallel stacking detection
├── detect_t_shaped()                  # T-shaped stacking detection
├── detect_cation_pi()                 # Cation-pi detection
├── analyze_frame()                    # Analyze single frame
├── process_trajectory()               # Iterate through all frames
├── extract_best_frame()               # Extract with MDAnalysis
├── add_dummy_atoms()                  # Add ring center markers
└── main()                             - Main workflow
```

---

## Performance Considerations

1. **Ring Detection**: Done once per system (rings don't change)
2. **KDTree**: Use for all distance searches (O(n log n) vs O(n²))
3. **Vectorization**: Batch angle calculations
4. **Frame Count**: No hardcoded limit, handles any trajectory length

---

## Status

- [x] Research completed
- [x] Plan created
- [ ] Waiting for approval
- [ ] Implementation

**Ready for executor spawn upon approval.**
