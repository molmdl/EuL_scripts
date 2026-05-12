# Pi-Interaction Detection Methods in Protein-Ligand Complexes

## Executive Summary

This document provides a comprehensive overview of computational methods for detecting pi-stacking, cation-pi, and other pi-related interactions in protein-ligand complexes, with specific focus on BINANA and PLIP implementations.

---

## 1. BINANA (BINding ANAlyzer)

### 1.1 Overview
- **Repository**: https://github.com/durrantlab/binana
- **License**: Apache 2.0
- **Citation**: Durrant JD, McCammon JA. J Mol Graph Model. 2011 Apr; 29(6): 888-893

### 1.2 Aromatic Ring Detection

BINANA uses a recursive approach to identify aromatic rings:

1. **For non-protein residues**: Identifies all 5 or 6-member rings
   - Checks dihedral angles between adjacent ring atoms
   - Checks dihedral angles between ring atoms and ring substituents
   - Rings are considered aromatic if no dihedral deviates from planarity by >15°

2. **For protein residues**: Uses standardized atom names
   - Phenylalanine (PHE): 1 aromatic ring
   - Tyrosine (TYR): 1 aromatic ring
   - Histidine (HIS): 1 aromatic ring
   - Tryptophan (TRP): 2 aromatic rings

3. **Ring Characterization**:
   - Plane defined by 3 ring atoms (preferably 1st, 3rd, 5th)
   - Center = average of all ring atom coordinates
   - Radius = max distance from center to any ring atom
   - Ring disk = center + plane orientation + radius + padding

### 1.3 Pi-Pi Stacking Detection

**Algorithm**:
1. Compare every ligand aromatic ring with every receptor aromatic ring
2. Check if ring centers are within `PI_PI_INTERACTING_DIST_CUTOFF` (7.5 Å)
3. Calculate angle between normal vectors to ring planes
4. If angle suggests planes are within `PI_STACKING_ANGLE_TOLERANCE` (30°) of parallel:
   - Project each ring atom onto the plane of the opposite ring
   - If any projected point falls within the ring disk (radius + padding), pi-stacking is detected

**Key Implementation Detail**: The algorithm projects ALL ring atoms (not just the center) because pi-stacking interactions are often off-center.

### 1.4 T-Shaped (Edge-Face) Detection

**Algorithm**:
1. Compare ligand and receptor aromatic rings
2. Check if ring centers are within `PI_PI_INTERACTING_DIST_CUTOFF` (7.5 Å)
3. Calculate angle between normal vectors
4. If angle suggests planes are within `T_STACKING_ANGLE_TOLERANCE` (30°) of perpendicular:
   - Verify rings come within `T_STACKING_CLOSEST_DIST_CUTOFF` (5.0 Å) at nearest point
   - Project center points onto opposite ring planes
   - If either projected center falls within the ring disk, T-stacking is detected

### 1.5 Cation-Pi Detection

**Algorithm**:
1. Identify charged functional groups:
   - **Proteins**: Lys (amine), Arg (guanidino), His (ring nitrogens), Glu/Asp (carboxylate)
   - **Ligands**: Metal cations, amines, guanidino groups, carboxylates, phosphates, sulfonates
2. For each positive charge and aromatic ring:
   - Check if distance < `CATION_PI_DIST_CUTOFF` (6.0 Å)
   - Project charge center onto aromatic ring plane
   - If projected point falls within ring disk (radius + padding), cation-pi detected

**Note**: BINANA only considers positive charges (no pi-anion interaction).

### 1.6 BINANA Default Parameters

| Parameter | Value | Description |
|-----------|-------|-------------|
| `PI_PADDING_DIST` | 0.75 Å | Added to ring radius for interaction detection |
| `PI_PI_INTERACTING_DIST_CUTOFF` | 7.5 Å | Ring-center distance cutoff |
| `PI_STACKING_ANGLE_TOLERANCE` | 30.0° | Angle tolerance for parallel orientation |
| `T_STACKING_ANGLE_TOLERANCE` | 30.0° | Angle tolerance for perpendicular orientation |
| `T_STACKING_CLOSEST_DIST_CUTOFF` | 5.0 Å | Closest atom-atom distance for T-stacking |
| `CATION_PI_DIST_CUTOFF` | 6.0 Å | Charge-to-ring-center distance |

---

## 2. PLIP (Protein-Ligand Interaction Profiler)

### 2.1 Overview
- **Repository**: https://github.com/pharmai/plip
- **License**: GPL-2.0
- **Citation**: Salentin et al. NAR 2015; Adasme et al. NAR 2021

### 2.2 Aromatic Ring Detection

PLIP uses OpenBabel for ring perception:
1. Uses SSSR (Smallest Set of Smallest Rings) perception
2. If OpenBabel doesn't report aromaticity:
   - Calculates normals of each ring atom to its neighbors
   - Angle between normal pairs must be < `AROMATIC_PLANARITY` (5°)
   - If planar, ring is considered aromatic

### 2.3 Pi-Stacking Detection

**Algorithm**:
1. Ring centers must be within `PISTACK_DIST_MAX` (5.5 Å)
2. Angle between ring planes must deviate no more than `PISTACK_ANG_DEV` (30°) from:
   - 0° (parallel P-stacking) OR
   - 90° (perpendicular T-stacking)
3. Project each ring center onto the opposite ring plane
4. Offset (distance between ring center and projected point) must be < `PISTACK_OFFSET_MAX` (2.0 Å)

**Classification**:
- Parallel stacking: angle deviation from 0° or 180° within tolerance
- T-shaped: angle deviation from 90° within tolerance

### 2.4 Pi-Cation Detection

**Algorithm**:
1. Distance between positive charge and aromatic ring center < `PICATION_DIST_MAX` (6.0 Å)
2. Offset (ring center to charge projection) < `PISTACK_OFFSET_MAX` (2.0 Å)
3. For tertiary amines in ligands, additional angle criterion applied

**Charged Group Detection**:
- **Proteins**:
  - Positive: Arg, His, Lys side chains
  - Negative: Asp, Glu carboxylates
- **Ligands**:
  - Positive: Quaternary ammonium, tertiary amines, sulfonium, guanidine
  - Negative: Phosphate, sulfonate, sulfonic acid, carboxylate

### 2.5 PLIP Default Parameters

| Parameter | Value | Description |
|-----------|-------|-------------|
| `AROMATIC_PLANARITY` | 5.0° | Max deviation from planarity |
| `PISTACK_DIST_MAX` | 5.5 Å | Max ring-center distance |
| `PISTACK_ANG_DEV` | 30° | Max deviation from optimal angle |
| `PISTACK_OFFSET_MAX` | 2.0 Å | Max offset (≈ benzene radius + 0.5 Å) |
| `PICATION_DIST_MAX` | 6.0 Å | Max charge-to-ring-center distance |
| `BS_DIST` | 7.5 Å | Binding site distance cutoff |

---

## 3. Comparison: BINANA vs PLIP

| Feature | BINANA | PLIP |
|---------|--------|------|
| **Ring detection** | Dihedral angle check (15°) | OpenBabel + planarity (5°) |
| **Pi-pi distance cutoff** | 7.5 Å | 5.5 Å |
| **Angle tolerance** | 30° | 30° |
| **Ring offset check** | All atoms projected | Center projected only |
| **T-stacking closest dist** | 5.0 Å | Not explicit |
| **Cation-pi distance** | 6.0 Å | 6.0 Å |
| **Anion-pi** | Not detected | Not detected |
| **Input formats** | PDBQT, PDB | PDB |

### Key Differences:

1. **BINANA uses larger distance cutoff (7.5 Å)** - based on literature showing pi-pi interactions often off-center
2. **BINANA projects all ring atoms** - better catches off-center stacking
3. **PLIP uses OpenBabel** - more robust aromaticity detection
4. **PLIP has smaller offset tolerance (2.0 Å)** - more stringent

---

## 4. Other Pi-Related Interactions

### 4.1 Anion-Pi Interactions

**Literature Criteria**:
- Distance: 3.0-6.0 Å (anion to ring centroid)
- Anion positioned above ring face
- Energy: 2-5 kcal/mol (weaker than cation-pi)

**Detection in tools**: Neither BINANA nor PLIP detect anion-pi interactions by default.

**Implementation approach**:
```python
# Similar to cation-pi but for negative charges
# Check if negative charge center within ~6Å of ring center
# Project charge onto ring plane
# Check if projection falls within ring disk
```

### 4.2 Sulfur-Pi Interactions

**Literature Criteria**:
- Distance: 3.5-5.5 Å (sulfur to ring center)
- Can involve Cys, Met sulfur atoms
- Weak but biologically relevant

### 4.3 Lone Pair-Pi Interactions

**Literature Criteria**:
- Distance: 3.0-4.0 Å
- Involves lone pairs of O, N interacting with pi systems
- Important in protein-ligand binding

---

## 5. Geometric Criteria from Literature

### 5.1 Pi-Pi Stacking

| Source | Distance (Å) | Angle (°) |
|--------|-------------|-----------|
| McGaughey et al., 1998 | 4.5-7.0 | < 30 (parallel) |
| Hunter & Sanders, 1990 | 3.4-4.0 | 0 (offset stacking) |
| Burley & Petsko, 1985 | 4.5-7.5 | < 30 |
| BINANA | 7.5 (centers) | < 30 |
| PLIP | 5.5 (centers) | < 30 |

### 5.2 T-Shaped (Edge-Face)

| Source | Distance (Å) | Angle (°) |
|--------|-------------|-----------|
| Waters et al., 2004 | ~5.0 | 60-90 |
| BINANA | 5.0 (closest atoms) | 60-90 (±30) |
| PLIP | 5.5 (centers) | ~90 (±30) |

### 5.3 Cation-Pi

| Source | Distance (Å) |
|--------|-------------|
| Gallivan & Dougherty, 1999 | < 6.0 |
| Dunietz et al., 2000 | 3.5-6.0 |
| BINANA | 6.0 |
| PLIP | 6.0 |

---

## 6. Python Libraries and Tools

### 6.1 Direct Use Libraries

| Library | Features | Link |
|---------|----------|------|
| **BINANA** | Pi-stacking, T-stacking, cation-pi | github.com/durrantlab/binana |
| **PLIP** | Pi-stacking, cation-pi, comprehensive | github.com/pharmai/plip |
| **OpenBabel** | Aromaticity detection, ring perception | openbabel.org |
| **RDKit** | Aromaticity, ring systems, geometry | rdkit.org |
| **MDAnalysis** | MD trajectory analysis | mdanalysis.org |
| **PyMOL API** | Visualization, structure manipulation | pymol.org |

### 6.2 Using BINANA as Python Library

```python
from binana.load_ligand_receptor import _get_ligand_receptor_aromatic_dists
from binana.interactions._pi_pi import get_pi_pi
from binana.interactions._cat_pi import get_cation_pi

# After loading ligand and receptor structures
pi_pi_results = get_pi_pi(ligand, receptor)
# Returns: {'counts': {}, 'mols': {}, 'labels': {}}

cation_pi_results = get_cation_pi(ligand, receptor)
# Returns: {'counts': {}, 'mol': Mol, 'labels': []}
```

### 6.3 Using PLIP as Python Library

```python
from plip.structure.preparation import PDBComplex

my_mol = PDBComplex()
my_mol.load_pdb('structure.pdb')
my_mol.analyze()
interactions = my_mol.interaction_sets['LIG:A:100']

# Access pi-stacking
for pistack in interactions.pistacking:
    print(f"Residue: {pistack.resnr}")
    print(f"Type: {pystack.type}")  # 'Parallel' or 'Perpendicular'
    print(f"Distance: {pistack.cent_dist}")
    print(f"Angle: {pistack.angle}")
    print(f"Offset: {pistack.offset}")

# Access pi-cation
for picat in interactions.pi_cation_interactions:
    print(f"Residue: {picat.resnr}")
    print(f"Distance: {picat.dist}")
    print(f"Offset: {picat.offset}")
```

---

## 7. Application to PDB Files and MD Trajectories

### 7.1 Single PDB File Analysis

```python
# Using PLIP
from plip.structure.preparation import PDBComplex

def analyze_pdb(pdb_file):
    mol = PDBComplex()
    mol.load_pdb(pdb_file)
    mol.analyze()
    
    results = {}
    for bs_id, interactions in mol.interaction_sets.items():
        results[bs_id] = {
            'pi_stacking': [(p.resnr, p.type, p.cent_dist) 
                           for p in interactions.pistacking],
            'pi_cation': [(p.resnr, p.dist) 
                         for p in interactions.pi_cation_interactions]
        }
    return results
```

### 7.2 MD Trajectory Analysis

```python
import MDAnalysis as mda
from plip.structure.preparation import PDBComplex
import tempfile
import os

def analyze_trajectory(topology, trajectory, ligand_selection):
    u = mda.Universe(topology, trajectory)
    ligand = u.select_atoms(ligand_selection)
    
    pi_interactions = []
    
    for ts in u.trajectory:
        # Write current frame to temp file
        with tempfile.NamedTemporaryFile(suffix='.pdb', delete=False) as f:
            temp_pdb = f.name
        
        u.atoms.write(temp_pdb)
        
        # Analyze with PLIP
        mol = PDBComplex()
        mol.load_pdb(temp_pdb)
        mol.analyze()
        
        frame_data = {'frame': ts.frame, 'time': ts.time}
        for bs_id, interactions in mol.interaction_sets.items():
            frame_data[bs_id] = {
                'pi_stacking_count': len(interactions.pistacking),
                'pi_cation_count': len(interactions.pi_cation_interactions)
            }
        
        pi_interactions.append(frame_data)
        os.unlink(temp_pdb)
    
    return pi_interactions
```

### 7.3 Custom Implementation

For custom implementations, key algorithms needed:

```python
import numpy as np
from scipy.spatial.distance import cdist

def get_ring_normal(ring_coords):
    """Calculate normal vector to ring plane."""
    # Use first three atoms to define plane
    p1, p2, p3 = ring_coords[:3]
    v1 = p2 - p1
    v2 = p3 - p1
    normal = np.cross(v1, v2)
    return normal / np.linalg.norm(normal)

def angle_between_normals(n1, n2):
    """Calculate angle between two normal vectors."""
    cos_angle = np.dot(n1, n2)
    return np.degrees(np.arccos(np.clip(cos_angle, -1, 1)))

def project_point_to_plane(point, plane_point, normal):
    """Project a point onto a plane."""
    v = point - plane_point
    dist = np.dot(v, normal)
    return point - dist * normal

def detect_pi_stacking(ring1_coords, ring2_coords, 
                       dist_cutoff=7.5, angle_tol=30.0, padding=0.75):
    """Detect pi-pi stacking between two rings."""
    center1 = ring1_coords.mean(axis=0)
    center2 = ring2_coords.mean(axis=0)
    
    # Distance check
    if np.linalg.norm(center1 - center2) > dist_cutoff:
        return False, None
    
    # Angle check
    n1 = get_ring_normal(ring1_coords)
    n2 = get_ring_normal(ring2_coords)
    angle = angle_between_normals(n1, n2)
    
    # Check if parallel (0° or 180°)
    min_angle = min(angle, 180 - angle)
    if min_angle > angle_tol:
        return False, None
    
    # Projection check
    radius1 = max(np.linalg.norm(ring1_coords - center1, axis=1))
    
    proj_center2 = project_point_to_plane(center2, center1, n1)
    offset = np.linalg.norm(proj_center2 - center1)
    
    if offset < radius1 + padding:
        return True, {'distance': np.linalg.norm(center1 - center2),
                      'angle': min_angle, 'offset': offset}
    
    return False, None

def detect_t_stacking(ring1_coords, ring2_coords,
                      dist_cutoff=7.5, angle_tol=30.0, 
                      closest_cutoff=5.0, padding=0.75):
    """Detect T-shaped stacking between two rings."""
    center1 = ring1_coords.mean(axis=0)
    center2 = ring2_coords.mean(axis=0)
    
    # Distance check
    if np.linalg.norm(center1 - center2) > dist_cutoff:
        return False, None
    
    # Angle check (should be ~90°)
    n1 = get_ring_normal(ring1_coords)
    n2 = get_ring_normal(ring2_coords)
    angle = angle_between_normals(n1, n2)
    
    if abs(90 - angle) > angle_tol:
        return False, None
    
    # Check closest atom-atom distance
    distances = cdist(ring1_coords, ring2_coords)
    min_dist = distances.min()
    
    if min_dist > closest_cutoff:
        return False, None
    
    # Projection check
    radius1 = max(np.linalg.norm(ring1_coords - center1, axis=1))
    proj_center2 = project_point_to_plane(center2, center1, n1)
    offset = np.linalg.norm(proj_center2 - center1)
    
    if offset < radius1 + padding:
        return True, {'min_distance': min_dist, 'angle': angle}
    
    return False, None

def detect_cation_pi(charge_center, ring_coords,
                     dist_cutoff=6.0, padding=0.75):
    """Detect cation-pi interaction."""
    ring_center = ring_coords.mean(axis=0)
    ring_normal = get_ring_normal(ring_coords)
    ring_radius = max(np.linalg.norm(ring_coords - ring_center, axis=1))
    
    # Distance check
    distance = np.linalg.norm(charge_center - ring_center)
    if distance > dist_cutoff:
        return False, None
    
    # Projection check
    proj_charge = project_point_to_plane(charge_center, ring_center, ring_normal)
    offset = np.linalg.norm(proj_charge - ring_center)
    
    if offset < ring_radius + padding:
        return True, {'distance': distance, 'offset': offset}
    
    return False, None
```

---

## 8. Recommended Approach for Your Use Case

### For PDB File Analysis:
1. **Use PLIP** for comprehensive, well-validated detection
2. **Use BINANA** if working with PDBQT files from AutoDock

### For MD Trajectory Analysis:
1. Use **MDAnalysis** for trajectory handling
2. Integrate with **PLIP** or **BINANA** for per-frame analysis
3. Consider caching ring/charge identifications for efficiency

### For Custom Implementation:
1. Use **RDKit** or **OpenBabel** for aromaticity detection
2. Implement geometric checks following BINANA/PLIP criteria
3. Validate against known complexes

---

## 9. Key Literature References

1. **Cation-Pi**: Gallivan JP, Dougherty DA. "Cation-pi interactions in structural biology." PNAS 1999, 96(17):9459-9464.

2. **Pi-Pi Stacking**: McGaughey GB et al. "Pi-pi interactions in proteins." J Biol Chem 1998, 273:15458-15463.

3. **BINANA**: Durrant JD, McCammon JA. "BINANA: A novel algorithm for ligand-binding characterization." J Mol Graph Model 2011, 29(6):888-893.

4. **PLIP**: Salentin S et al. "PLIP: fully automated protein-ligand interaction profiler." NAR 2015, 43(W1):W443-W447.

5. **T-Shaped**: Waters ML. "Aromatic interactions in model systems." Curr Opin Chem Biol 2002, 6:736-741.

---

## 10. Summary Table: Detection Criteria

| Interaction | Distance (Å) | Angle (°) | Offset (Å) | Additional |
|-------------|-------------|-----------|------------|------------|
| Pi-Pi Stacking | < 7.5 (BINANA) / 5.5 (PLIP) | < 30 from parallel | < radius + 0.75 | - |
| T-Shaped | < 7.5 centers, < 5.0 closest | < 30 from perpendicular | < radius + 0.75 | Min atom-atom < 5.0 Å |
| Cation-Pi | < 6.0 | N/A | < radius + 0.75 | Charge above ring |
| Anion-Pi | < 6.0 | N/A | < radius + 0.75 | Not in tools |

---

*Document generated for pi-interaction detection research*
*Based on BINANA v2.1 and PLIP v3.0.0*
