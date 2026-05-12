#!/usr/bin/env python3
"""
Detect pi-interactions (pi-stacking, T-shaped stacking, cation-pi) in MD trajectories.
Extract frames with maximum pi-interactions for each system.

Geometric criteria based on BINANA:
- Parallel pi-stacking: ring centers <7.5A, angle <30deg, projection in ring disk
- T-shaped stacking: ring centers <7.5A, angle 90+-30deg, closest atoms <5A
- Cation-pi: charge <6A from ring, projection in ring disk
"""

import os
import sys
import warnings
import numpy as np
import pandas as pd
from pathlib import Path
from collections import defaultdict

warnings.filterwarnings('ignore', category=UserWarning, module='MDAnalysis')

import MDAnalysis as mda
from MDAnalysis.analysis import distances
from scipy.spatial import cKDTree


PROTEIN_AROMATIC_RINGS = {
    'PHE': [['CG', 'CD1', 'CD2', 'CE1', 'CE2', 'CZ']],
    'TYR': [['CG', 'CD1', 'CD2', 'CE1', 'CE2', 'CZ']],
    'HIS': [['CG', 'ND1', 'CD2', 'CE1', 'NE2']],
    'TRP': [
        ['CG', 'CD1', 'NE1', 'CE2', 'CD2', 'CE3'],  
        ['CD2', 'CE2', 'CE3', 'CZ2', 'CZ3', 'CH2']  
    ],
}

PROTEIN_CHARGE_CENTERS = {
    'ARG': {'atom': 'CZ', 'name': 'guanidinium'},
    'LYS': {'atom': 'NZ', 'name': 'amine'},
    'HIS': {'atom': 'NE2', 'name': 'imidazole'},  
}

PARALLEL_DIST_CUTOFF = 7.5  
PARALLEL_ANGLE_CUTOFF = 30.0  
TSHAPED_DIST_CUTOFF = 7.5  
TSHAPED_ANGLE_MIN = 60.0  
TSHAPED_ANGLE_MAX = 120.0  
TSHAPED_ATOM_DIST = 5.0  
CATION_PI_DIST_CUTOFF = 6.0  
RING_PADDING = 0.75  

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

AMIDE_PI_DISTANCE_CUTOFF = 4.0
AMIDE_PI_ANGLE_MIN = 150.0  


def get_ring_geometry(coords):
    """
    Calculate ring center, normal vector, and approximate radius.
    
    Parameters:
    -----------
    coords : np.ndarray
        Ring atom coordinates (N x 3)
    
    Returns:
    --------
    dict with 'center', 'normal', 'radius'
    """
    center = np.mean(coords, axis=0)
    
    centered = coords - center
    
    if len(coords) >= 3:
        v1 = centered[1] - centered[0]
        v2 = centered[2] - centered[0]
        normal = np.cross(v1, v2)
        
        if np.linalg.norm(normal) > 1e-6:
            normal = normal / np.linalg.norm(normal)
        else:
            normal = np.array([0.0, 0.0, 1.0])
    else:
        normal = np.array([0.0, 0.0, 1.0])
    
    radii = np.linalg.norm(centered, axis=1)
    radius = np.mean(radii)
    
    return {'center': center, 'normal': normal, 'radius': radius}


def project_point_to_plane(point, plane_point, plane_normal):
    """Project a point onto a plane defined by a point and normal."""
    v = point - plane_point
    dist = np.dot(v, plane_normal)
    return point - dist * plane_normal, abs(dist)


def angle_between_vectors(v1, v2):
    """Calculate angle in degrees between two vectors."""
    v1_norm = v1 / (np.linalg.norm(v1) + 1e-10)
    v2_norm = v2 / (np.linalg.norm(v2) + 1e-10)
    cos_angle = np.clip(np.dot(v1_norm, v2_norm), -1.0, 1.0)
    angle = np.degrees(np.arccos(abs(cos_angle)))  
    return min(angle, 180.0 - angle)


def detect_aromatic_rings_protein(universe):
    """
    Detect aromatic rings in protein residues.
    
    Parameters:
    -----------
    universe : MDAnalysis.Universe
    
    Returns:
    --------
    list of dict: each dict has 'resname', 'resid', 'atoms', 'atom_indices'
    """
    rings = []
    protein = universe.select_atoms('protein')
    
    for resname, ring_definitions in PROTEIN_AROMATIC_RINGS.items():
        residues = protein.select_atoms(f'resname {resname}').residues
        
        for res in residues:
            for ring_idx, ring_atoms in enumerate(ring_definitions):
                atom_group = res.atoms.select_atoms(' or '.join([f'name {a}' for a in ring_atoms]))
                
                if len(atom_group) == len(ring_atoms):
                    atom_indices = atom_group.indices.tolist()
                    rings.append({
                        'resname': resname,
                        'resid': res.resid,
                        'ring_idx': ring_idx,
                        'atom_names': ring_atoms,
                        'atom_indices': atom_indices,
                        'source': 'protein'
                    })
    
    return rings


def detect_aromatic_rings_ligand(universe):
    """
    Detect aromatic rings in ligand using connectivity-based approach.
    
    Parameters:
    -----------
    universe : MDAnalysis.Universe
    
    Returns:
    --------
    list of dict: each dict has ring info
    """
    rings = []
    
    ligand = universe.select_atoms('not protein and not resname HOH SOL WAT')
    if len(ligand) == 0:
        ligand = universe.select_atoms('segid B')
    
    if len(ligand) == 0:
        return rings
    
    bonds = _infer_bonds_from_distances(ligand)
    
    ring_cycles = _find_rings_from_bonds(bonds, ligand)
    
    for cycle in ring_cycles:
        if len(cycle) in [5, 6]:
            coords = ligand[cycle].positions
            
            if _is_planar(coords):
                rings.append({
                    'resname': ligand[0].resname,
                    'resid': ligand[0].resid,
                    'ring_idx': len(rings),
                    'atom_names': [ligand[i].name for i in cycle],
                    'atom_indices': ligand[cycle].indices.tolist(),
                    'source': 'ligand',
                    'cycle_atoms': cycle
                })
    
    return rings


def _infer_bonds_from_distances(atoms, cutoff=1.8):
    """Infer bonds from atom distances."""
    bonds = defaultdict(set)
    coords = atoms.positions
    n_atoms = len(atoms)
    
    for i in range(n_atoms):
        for j in range(i + 1, n_atoms):
            dist = np.linalg.norm(coords[i] - coords[j])
            if dist < cutoff:
                bonds[i].add(j)
                bonds[j].add(i)
    
    return bonds


def _find_rings_from_bonds(bonds, atoms, max_size=7):
    """Find rings using depth-first search."""
    rings = set()
    n_atoms = len(atoms)
    
    def dfs(start, current, path, visited):
        if len(path) > max_size:
            return
        
        for neighbor in bonds[current]:
            if neighbor == start and len(path) >= 4:
                ring = tuple(sorted(path))
                rings.add(ring)
            elif neighbor not in visited:
                dfs(start, neighbor, path + [neighbor], visited | {neighbor})
    
    for start in range(n_atoms):
        dfs(start, start, [start], {start})
    
    unique_rings = []
    for ring in rings:
        if len(ring) in [5, 6]:
            unique_rings.append(list(ring))
    
    return unique_rings


def _is_planar(coords, threshold=15.0):
    """Check if ring atoms are planar."""
    if len(coords) < 4:
        return True
    
    center = np.mean(coords, axis=0)
    centered = coords - center
    
    v1 = centered[1] - centered[0]
    v2 = centered[2] - centered[0]
    normal = np.cross(v1, v2)
    
    if np.linalg.norm(normal) < 1e-6:
        return False
    
    normal = normal / np.linalg.norm(normal)
    
    for i in range(3, len(coords)):
        v = centered[i]
        dist_to_plane = abs(np.dot(v, normal))
        if dist_to_plane > 0.5:  
            return False
    
    return True


def detect_charge_centers(universe, ligand_rings):
    """
    Detect charged groups for cation-pi interactions.
    
    Parameters:
    -----------
    universe : MDAnalysis.Universe
    ligand_rings : list
        Detected ligand rings (to find ligand charged groups)
    
    Returns:
    --------
    dict with 'protein' and 'ligand' charge centers
    """
    charge_centers = {'protein': [], 'ligand': []}
    
    protein = universe.select_atoms('protein')
    
    for resname, info in PROTEIN_CHARGE_CENTERS.items():
        residues = protein.select_atoms(f'resname {resname}').residues
        
        for res in residues:
            charge_atoms = res.atoms.select_atoms(f'name {info["atom"]}')
            
            if len(charge_atoms) > 0:
                charge_centers['protein'].append({
                    'resname': resname,
                    'resid': res.resid,
                    'atom_name': info['atom'],
                    'atom_index': charge_atoms[0].index,
                    'charge_type': info['name'],
                    'source': 'protein'
                })
    
    ligand = universe.select_atoms('not protein and not resname HOH SOL WAT')
    if len(ligand) == 0:
        ligand = universe.select_atoms('segid B')
    
    if len(ligand) > 0:
        nitrogen_atoms = ligand.select_atoms('name N or type N')
        
        for atom in nitrogen_atoms:
            neighbors = _get_bonded_neighbors(ligand, atom.index - ligand[0].index)
            n_bonds = len(neighbors)
            
            if n_bonds >= 3:  
                charge_centers['ligand'].append({
                    'resname': atom.resname,
                    'resid': atom.resid,
                    'atom_name': atom.name,
                    'atom_index': atom.index,
                    'charge_type': 'amine' if n_bonds == 3 else 'quaternary',
                    'source': 'ligand'
                })
    
    return charge_centers


def _get_bonded_neighbors(atoms, atom_idx, cutoff=1.8):
    """Get bonded neighbors of an atom."""
    coords = atoms.positions
    neighbors = []
    
    for i in range(len(atoms)):
        if i != atom_idx:
            dist = np.linalg.norm(coords[atom_idx] - coords[i])
            if dist < cutoff:
                neighbors.append(i)
    
    return neighbors


def detect_sidechain_amides(universe):
    """
    Detect side-chain amide groups (ASN, GLN) in protein.
    
    Parameters:
    -----------
    universe : MDAnalysis.Universe
    
    Returns:
    --------
    list of dict: each dict has amide group info
    """
    amides = []
    protein = universe.select_atoms('protein')
    
    for resname, info in SIDECHAIN_AMIDES.items():
        residues = protein.select_atoms(f'resname {resname}').residues
        
        for res in residues:
            n_atom = res.atoms.select_atoms(f'name {info["atoms"]["N"]}')
            o_atom = res.atoms.select_atoms(f'name {info["atoms"]["O"]}')
            
            if len(n_atom) == 0 or len(o_atom) == 0:
                continue
            
            h_atoms = []
            for h_name in info['atoms']['H']:
                h_atom = res.atoms.select_atoms(f'name {h_name}')
                if len(h_atom) > 0:
                    h_atoms.append(h_atom[0])
            
            amide_data = {
                'resname': resname,
                'resid': res.resid,
                'N_index': n_atom[0].index,
                'O_index': o_atom[0].index,
                'source': 'protein'
            }
            
            if len(h_atoms) > 0:
                amide_data['H_indices'] = [h.index for h in h_atoms]
            else:
                amide_data['H_indices'] = []
            
            amides.append(amide_data)
    
    return amides


def detect_amide_pi_nh(amide, ring, positions, geo_ring):
    """
    Detect NH-π interaction (amide N-H → ring).
    
    Parameters:
    -----------
    amide : dict
        Amide group information
    ring : dict
        Ring information
    positions : np.ndarray
        All atom positions
    geo_ring : dict
        Ring geometry
    
    Returns:
    --------
    tuple: (is_interaction: bool, details: dict)
    """
    if len(amide.get('H_indices', [])) == 0:
        return False, {}
    
    n_pos = positions[amide['N_index']]
    
    for h_idx in amide['H_indices']:
        h_pos = positions[h_idx]
        
        dist = np.linalg.norm(h_pos - geo_ring['center'])
        
        if dist <= AMIDE_PI_DISTANCE_CUTOFF:
            nh_vector = h_pos - n_pos
            h_ring_vector = geo_ring['center'] - h_pos
            
            if np.linalg.norm(nh_vector) < 1e-6 or np.linalg.norm(h_ring_vector) < 1e-6:
                continue
            
            nh_vector = nh_vector / np.linalg.norm(nh_vector)
            h_ring_vector = h_ring_vector / np.linalg.norm(h_ring_vector)
            
            cos_angle = np.dot(nh_vector, h_ring_vector)
            angle = np.degrees(np.arccos(np.clip(cos_angle, -1.0, 1.0)))
            
            if angle >= AMIDE_PI_ANGLE_MIN:
                return True, {
                    'type': 'NH_pi',
                    'distance': dist,
                    'angle': angle,
                    'amide': f"{amide['resname']}{amide['resid']}",
                    'ring': f"{ring['resname']}{ring['resid']}"
                }
    
    return False, {}


def detect_amide_pi_npi(amide, ring, positions, geo_ring):
    """
    Detect n-π* interaction (amide O lone pair → ring).
    
    Parameters:
    -----------
    amide : dict
        Amide group information
    ring : dict
        Ring information
    positions : np.ndarray
        All atom positions
    geo_ring : dict
        Ring geometry
    
    Returns:
    --------
    tuple: (is_interaction: bool, details: dict)
    """
    o_pos = positions[amide['O_index']]
    
    dist = np.linalg.norm(o_pos - geo_ring['center'])
    
    if dist <= AMIDE_PI_DISTANCE_CUTOFF:
        return True, {
            'type': 'n_pi_star',
            'distance': dist,
            'amide': f"{amide['resname']}{amide['resid']}",
            'ring': f"{ring['resname']}{ring['resid']}"
        }
    
    return False, {}


def detect_parallel_stacking(ring1, ring2, positions, geo1, geo2):
    """
    Detect parallel pi-stacking between two rings.
    
    Parameters:
    -----------
    ring1, ring2 : dict
        Ring information dicts
    positions : np.ndarray
        All atom positions
    geo1, geo2 : dict
        Ring geometry (center, normal, radius)
    
    Returns:
    --------
    tuple: (is_stacking: bool, details: dict)
    """
    center_dist = np.linalg.norm(geo1['center'] - geo2['center'])
    
    if center_dist >= PARALLEL_DIST_CUTOFF:
        return False, {}
    
    angle = angle_between_vectors(geo1['normal'], geo2['normal'])
    
    if angle >= PARALLEL_ANGLE_CUTOFF:
        return False, {}
    
    proj1, _ = project_point_to_plane(geo1['center'], geo2['center'], geo2['normal'])
    dist_to_center = np.linalg.norm(proj1 - geo2['center'])
    
    if dist_to_center > geo2['radius'] + RING_PADDING:
        return False, {}
    
    proj2, _ = project_point_to_plane(geo2['center'], geo1['center'], geo1['normal'])
    dist_to_center2 = np.linalg.norm(proj2 - geo1['center'])
    
    if dist_to_center2 > geo1['radius'] + RING_PADDING:
        return False, {}
    
    return True, {
        'center_dist': center_dist,
        'angle': angle,
        'ring1': f"{ring1['resname']}{ring1['resid']}",
        'ring2': f"{ring2['resname']}{ring2['resid']}"
    }


def detect_t_shaped_stacking(ring1, ring2, positions, geo1, geo2):
    """
    Detect T-shaped stacking between two rings.
    
    Parameters:
    -----------
    ring1, ring2 : dict
        Ring information dicts
    positions : np.ndarray
        All atom positions
    geo1, geo2 : dict
        Ring geometry (center, normal, radius)
    
    Returns:
    --------
    tuple: (is_stacking: bool, details: dict)
    """
    center_dist = np.linalg.norm(geo1['center'] - geo2['center'])
    
    if center_dist >= TSHAPED_DIST_CUTOFF:
        return False, {}
    
    angle = angle_between_vectors(geo1['normal'], geo2['normal'])
    
    if angle < TSHAPED_ANGLE_MIN or angle > TSHAPED_ANGLE_MAX:
        return False, {}
    
    atoms1 = positions[ring1['atom_indices']]
    atoms2 = positions[ring2['atom_indices']]
    
    min_atom_dist = np.min(distances.distance_array(atoms1, atoms2))
    
    if min_atom_dist >= TSHAPED_ATOM_DIST:
        return False, {}
    
    proj1, _ = project_point_to_plane(geo1['center'], geo2['center'], geo2['normal'])
    dist_to_center = np.linalg.norm(proj1 - geo2['center'])
    
    if dist_to_center > geo2['radius'] + RING_PADDING:
        return False, {}
    
    return True, {
        'center_dist': center_dist,
        'angle': angle,
        'min_atom_dist': min_atom_dist,
        'ring1': f"{ring1['resname']}{ring1['resid']}",
        'ring2': f"{ring2['resname']}{ring2['resid']}"
    }


def detect_cation_pi(charge, ring, positions, charge_pos, geo_ring):
    """
    Detect cation-pi interaction.
    
    Parameters:
    -----------
    charge : dict
        Charge center information
    ring : dict
        Ring information
    positions : np.ndarray
        All atom positions
    charge_pos : np.ndarray
        Charge center position
    geo_ring : dict
        Ring geometry
    
    Returns:
    --------
    tuple: (is_interaction: bool, details: dict)
    """
    dist = np.linalg.norm(charge_pos - geo_ring['center'])
    
    if dist >= CATION_PI_DIST_CUTOFF:
        return False, {}
    
    proj, height = project_point_to_plane(charge_pos, geo_ring['center'], geo_ring['normal'])
    dist_to_center = np.linalg.norm(proj - geo_ring['center'])
    
    if dist_to_center > geo_ring['radius'] + RING_PADDING:
        return False, {}
    
    return True, {
        'distance': dist,
        'height': height,
        'charge': f"{charge['resname']}{charge['resid']}",
        'ring': f"{ring['resname']}{ring['resid']}"
    }


def analyze_frame(universe, frame_idx, protein_rings, ligand_rings, charge_centers, sidechain_amides=None):
    """
    Analyze a single frame for pi-interactions.
    
    Returns:
    --------
    dict with interaction counts and details
    """
    universe.trajectory[frame_idx]
    positions = universe.atoms.positions
    
    results = {
        'parallel_stacking': [],
        't_shaped_stacking': [],
        'cation_pi': [],
        'amide_pi': [],
        'residues_involved': set()
    }
    
    protein_geometries = []
    for ring in protein_rings:
        coords = positions[ring['atom_indices']]
        geo = get_ring_geometry(coords)
        protein_geometries.append((ring, geo))
    
    ligand_geometries = []
    for ring in ligand_rings:
        coords = positions[ring['atom_indices']]
        geo = get_ring_geometry(coords)
        ligand_geometries.append((ring, geo))
    
    for p_ring, p_geo in protein_geometries:
        for l_ring, l_geo in ligand_geometries:
            is_parallel, details = detect_parallel_stacking(
                p_ring, l_ring, positions, p_geo, l_geo
            )
            if is_parallel:
                results['parallel_stacking'].append(details)
                results['residues_involved'].add(p_ring['resid'])
                results['residues_involved'].add(l_ring['resid'])
                continue
            
            is_tshaped, details = detect_t_shaped_stacking(
                p_ring, l_ring, positions, p_geo, l_geo
            )
            if is_tshaped:
                results['t_shaped_stacking'].append(details)
                results['residues_involved'].add(p_ring['resid'])
                results['residues_involved'].add(l_ring['resid'])
    
    protein_charges = charge_centers['protein']
    ligand_charges = charge_centers['ligand']
    
    for charge in protein_charges:
        charge_pos = positions[charge['atom_index']]
        for l_ring, l_geo in ligand_geometries:
            is_cation_pi, details = detect_cation_pi(
                charge, l_ring, positions, charge_pos, l_geo
            )
            if is_cation_pi:
                results['cation_pi'].append(details)
                results['residues_involved'].add(charge['resid'])
                results['residues_involved'].add(l_ring['resid'])
    
    for charge in ligand_charges:
        charge_pos = positions[charge['atom_index']]
        for p_ring, p_geo in protein_geometries:
            is_cation_pi, details = detect_cation_pi(
                charge, p_ring, positions, charge_pos, p_geo
            )
            if is_cation_pi:
                results['cation_pi'].append(details)
                results['residues_involved'].add(charge['resid'])
                results['residues_involved'].add(p_ring['resid'])
    
    if sidechain_amides:
        for amide in sidechain_amides:
            for l_ring, l_geo in ligand_geometries:
                is_nh_pi, details = detect_amide_pi_nh(
                    amide, l_ring, positions, l_geo
                )
                if is_nh_pi:
                    results['amide_pi'].append(details)
                    results['residues_involved'].add(amide['resid'])
                    results['residues_involved'].add(l_ring['resid'])
                
                is_n_pi, details2 = detect_amide_pi_npi(
                    amide, l_ring, positions, l_geo
                )
                if is_n_pi:
                    results['amide_pi'].append(details2)
                    results['residues_involved'].add(amide['resid'])
                    results['residues_involved'].add(l_ring['resid'])
    
    results['residues_involved'] = sorted(results['residues_involved'])
    results['total'] = (
        len(results['parallel_stacking']) + 
        len(results['t_shaped_stacking']) + 
        len(results['cation_pi']) +
        len(results['amide_pi'])
    )
    
    return results


def calculate_statistics(all_frame_data):
    """
    Calculate mean and max statistics across all frames.
    
    Parameters:
    -----------
    all_frame_data : list
        List of frame data dicts from process_trajectory
    
    Returns:
    --------
    dict: statistics with mean and max for each interaction type
    """
    if not all_frame_data:
        return {}
    
    interaction_types = ['parallel_stacking', 't_shaped_stacking', 'cation_pi', 'amide_pi', 'total']
    stats = {}
    
    for itype in interaction_types:
        values = [fd.get(itype, 0) for fd in all_frame_data]
        stats[f'{itype}_mean'] = np.mean(values)
        stats[f'{itype}_max'] = np.max(values)
    
    return stats


def process_trajectory(system_name, pdb_file, xtc_file):
    """
    Process entire trajectory for a system.
    
    Returns:
    --------
    tuple: (best_frame_idx, best_results, all_frame_data, protein_rings, ligand_rings, u, statistics)
    """
    print(f"  Loading trajectory...")
    u = mda.Universe(pdb_file, xtc_file)
    
    print(f"  Detecting protein aromatic rings...")
    protein_rings = detect_aromatic_rings_protein(u)
    print(f"    Found {len(protein_rings)} protein aromatic rings")
    
    print(f"  Detecting ligand aromatic rings...")
    ligand_rings = detect_aromatic_rings_ligand(u)
    print(f"    Found {len(ligand_rings)} ligand aromatic rings")
    
    print(f"  Detecting charge centers...")
    charge_centers = detect_charge_centers(u, ligand_rings)
    print(f"    Protein: {len(charge_centers['protein'])} charge centers")
    print(f"    Ligand: {len(charge_centers['ligand'])} charge centers")
    
    print(f"  Detecting side-chain amides...")
    sidechain_amides = detect_sidechain_amides(u)
    print(f"    Found {len(sidechain_amides)} side-chain amides (ASN/GLN)")
    
    n_frames = len(u.trajectory)
    print(f"  Analyzing {n_frames} frames...")
    
    all_frame_data = []
    best_frame_idx = 0
    best_total = 0
    best_results = None
    
    for frame_idx in range(n_frames):
        if frame_idx % 500 == 0:
            print(f"    Frame {frame_idx}/{n_frames}...")
        
        results = analyze_frame(u, frame_idx, protein_rings, ligand_rings, charge_centers, sidechain_amides)
        
        trial_id = f"mmpbsa_{frame_idx // 1000}"
        frame_num = frame_idx % 1000
        
        frame_data = {
            'system': system_name,
            'trial_id': trial_id,
            'frame_num': frame_num,
            'global_frame': frame_idx,
            'parallel_stacking': len(results['parallel_stacking']),
            't_shaped_stacking': len(results['t_shaped_stacking']),
            'cation_pi': len(results['cation_pi']),
            'amide_pi': len(results['amide_pi']),
            'total': results['total'],
            'residues': ','.join(map(str, results['residues_involved'])) if results['residues_involved'] else ''
        }
        all_frame_data.append(frame_data)
        
        if results['total'] > best_total:
            best_total = results['total']
            best_frame_idx = frame_idx
            best_results = results
    
    statistics = calculate_statistics(all_frame_data)
    
    print(f"  Best frame: {best_frame_idx} with {best_total} total pi-interactions")
    
    return best_frame_idx, best_results, all_frame_data, protein_rings, ligand_rings, u, statistics


def extract_frame_with_dummy_atoms(u, frame_idx, protein_rings, ligand_rings, best_results, output_file):
    """
    Extract a frame and add dummy atoms at ring centers.
    Uses MDAnalysis (NOT gmx).
    """
    u.trajectory[frame_idx]
    
    protein_ligand = u.select_atoms('protein or not resname HOH SOL WAT')
    if len(protein_ligand) == 0:
        protein_ligand = u.select_atoms('all and not resname HOH SOL WAT')
    
    pdb_lines = []
    pdb_lines.append("REMARK PI-INTERACTIONS")
    pdb_lines.append(f"REMARK Parallel stacking: {len(best_results['parallel_stacking'])}")
    pdb_lines.append(f"REMARK T-shaped stacking: {len(best_results['t_shaped_stacking'])}")
    pdb_lines.append(f"REMARK Cation-pi: {len(best_results['cation_pi'])}")
    pdb_lines.append(f"REMARK Amide-pi: {len(best_results.get('amide_pi', []))}")
    pdb_lines.append(f"REMARK Total: {best_results['total']}")
    pdb_lines.append(f"REMARK Residues involved: {','.join(map(str, best_results['residues_involved']))}")
    pdb_lines.append(f"REMARK Frame: {frame_idx}")
    
    atom_serial = 1
    for atom in protein_ligand:
        line = f"ATOM  {atom_serial:5d} {atom.name:<4s} {atom.resname:>3s} {atom.segid}{atom.resid:4d}    {atom.position[0]:8.3f}{atom.position[1]:8.3f}{atom.position[2]:8.3f}{atom.occupancy:6.2f}{atom.bfactor:6.2f}          {atom.element:>2s}"
        pdb_lines.append(line)
        atom_serial += 1
    
    dummy_positions = []
    
    involved_ring_indices = set()
    for stack in best_results['parallel_stacking']:
        for ring in protein_rings + ligand_rings:
            ring_id = f"{ring['resname']}{ring['resid']}"
            if ring_id in stack['ring1'] or ring_id in stack['ring2']:
                involved_ring_indices.add(id(ring))
    
    for stack in best_results['t_shaped_stacking']:
        for ring in protein_rings + ligand_rings:
            ring_id = f"{ring['resname']}{ring['resid']}"
            if ring_id in stack['ring1'] or ring_id in stack['ring2']:
                involved_ring_indices.add(id(ring))
    
    for cation_pi in best_results['cation_pi']:
        for ring in protein_rings + ligand_rings:
            ring_id = f"{ring['resname']}{ring['resid']}"
            if ring_id in cation_pi['ring']:
                involved_ring_indices.add(id(ring))
    
    for amide_pi in best_results.get('amide_pi', []):
        for ring in protein_rings + ligand_rings:
            ring_id = f"{ring['resname']}{ring['resid']}"
            if ring_id in amide_pi.get('ring', ''):
                involved_ring_indices.add(id(ring))
    
    for ring in protein_rings + ligand_rings:
        coords = u.atoms.positions[ring['atom_indices']]
        geo = get_ring_geometry(coords)
        dummy_positions.append({
            'position': geo['center'],
            'resname': ring['resname'],
            'resid': ring['resid'],
            'source': ring['source']
        })
    
    for dummy in dummy_positions:
        line = f"HETATM{atom_serial:5d}  PI  XPI B   1    {dummy['position'][0]:8.3f}{dummy['position'][1]:8.3f}{dummy['position'][2]:8.3f}  1.00  0.00"
        pdb_lines.append(line)
        atom_serial += 1
    
    pdb_lines.append("END")
    
    with open(output_file, 'w') as f:
        f.write('\n'.join(pdb_lines) + '\n')
    
    print(f"  Saved {output_file}")


def find_systems(trj_dir):
    """Find all systems with fp/v1.pdb and v1.xtc files."""
    systems = []
    
    if os.path.islink(trj_dir) or os.path.isdir(trj_dir):
        for item in os.listdir(trj_dir):
            system_path = os.path.join(trj_dir, item)
            if os.path.isdir(system_path):
                fp_path = os.path.join(system_path, 'fp')
                pdb_file = os.path.join(fp_path, 'v1.pdb')
                xtc_file = os.path.join(fp_path, 'v1.xtc')
                
                if os.path.exists(pdb_file) and os.path.exists(xtc_file):
                    systems.append(item)
    
    return sorted(systems)


def main():
    """Main function to process all systems."""
    trj_dir = 'trj'
    
    print("Finding systems...")
    systems = find_systems(trj_dir)
    print(f"Found {len(systems)} systems: {', '.join(systems)}")
    
    if len(systems) == 0:
        print("No valid systems found. Exiting.")
        return
    
    all_results = []
    
    for system in systems:
        print(f"\n{'='*60}")
        print(f"Processing {system}...")
        print(f"{'='*60}")
        
        pdb_file = os.path.join(trj_dir, system, 'fp', 'v1.pdb')
        xtc_file = os.path.join(trj_dir, system, 'fp', 'v1.xtc')
        output_file = f"{system}_best_pi.pdb"
        
        if os.path.exists(output_file):
            print(f"  Output file {output_file} already exists. Skipping.")
            continue
        
        try:
            best_frame_idx, best_results, frame_data, protein_rings, ligand_rings, u, statistics = process_trajectory(
                system, pdb_file, xtc_file
            )
            
            extract_frame_with_dummy_atoms(
                u, best_frame_idx, protein_rings, ligand_rings, best_results, output_file
            )
            
            best_frame_data = [fd for fd in frame_data if fd['global_frame'] == best_frame_idx][0]
            
            best_frame_data['best_frame'] = best_frame_data.pop('global_frame')
            best_frame_data.pop('trial_id', None)
            best_frame_data.pop('frame_num', None)
            
            for key, value in statistics.items():
                best_frame_data[key] = value
            
            cols_order = ['system', 'best_frame', 'parallel_stacking', 't_shaped_stacking', 
                          'cation_pi', 'amide_pi', 'total',
                          'parallel_stacking_mean', 'parallel_stacking_max',
                          't_shaped_stacking_mean', 't_shaped_stacking_max',
                          'cation_pi_mean', 'cation_pi_max',
                          'amide_pi_mean', 'amide_pi_max',
                          'total_mean', 'total_max', 'residues']
            best_frame_data = {k: best_frame_data.get(k, '') for k in cols_order}
            
            all_results.append(best_frame_data)
            
        except Exception as e:
            print(f"  Error processing {system}: {e}")
            import traceback
            traceback.print_exc()
            continue
    
    if all_results:
        df = pd.DataFrame(all_results)
        csv_file = 'pi_interactions_summary.csv'
        
        if not os.path.exists(csv_file):
            df.to_csv(csv_file, index=False)
            print(f"\nSummary saved to {csv_file}")
        else:
            existing_df = pd.read_csv(csv_file)
            combined_df = pd.concat([existing_df, df], ignore_index=True)
            combined_df.to_csv(csv_file, index=False)
            print(f"\nSummary appended to {csv_file}")
    
    print(f"\n{'='*60}")
    print("Processing complete!")
    print(f"{'='*60}")


if __name__ == '__main__':
    main()
