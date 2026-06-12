#!/usr/bin/env python3
"""
Detect frames where BOTH Arg410 AND Lys414 have cation-pi interactions with ligand.
Distinguishes between interactions with ANY ring vs 4-fused ring systems.
"""

import os
import sys
import warnings
import argparse
import logging
import numpy as np
import pandas as pd
from pathlib import Path
from collections import defaultdict
from itertools import combinations
from scipy.spatial.distance import cdist

warnings.filterwarnings('ignore', category=UserWarning, module='MDAnalysis')

import MDAnalysis as mda
from MDAnalysis.analysis import distances

logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)

CATION_PI_DIST_CUTOFF = 6.0
RING_PADDING = 0.75


def get_ring_geometry(coords):
    """Calculate ring center, normal vector, and approximate radius."""
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


def detect_aromatic_rings_ligand(universe):
    """Detect aromatic rings in ligand using connectivity-based approach."""
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


def check_planarity(points, threshold=0.8):
    """
    Check if points lie on the same plane using SVD.
    
    Returns True if RMSD from best-fit plane < threshold.
    """
    if len(points) < 3:
        return True
    
    centroid = np.mean(points, axis=0)
    centered = points - centroid
    
    U, S, Vt = np.linalg.svd(centered)
    normal = Vt[-1]
    distances = np.abs(np.dot(centered, normal))
    rmsd = np.sqrt(np.mean(distances**2))
    
    return rmsd < threshold


def find_connected_components(graph):
    """Find all connected components in a graph."""
    visited = set()
    components = []
    
    for start in graph:
        if start in visited:
            continue
        
        component = []
        queue = [start]
        while queue:
            node = queue.pop(0)
            if node in visited:
                continue
            visited.add(node)
            component.append(node)
            queue.extend(graph[node] - visited)
        
        components.append(set(component))
    
    return components


def find_4fused_ring_systems_topology(ligand_rings):
    """
    Topology-based detection of 4-fused ring systems.
    
    Two rings are fused if they share >= 2 atoms (a bond).
    """
    n_rings = len(ligand_rings)
    if n_rings < 4:
        return []
    
    # Build fusion adjacency graph
    fusion_graph = {i: set() for i in range(n_rings)}
    
    for i in range(n_rings):
        atoms_i = set(ligand_rings[i]['atom_indices'])
        for j in range(i + 1, n_rings):
            atoms_j = set(ligand_rings[j]['atom_indices'])
            shared = atoms_i & atoms_j
            
            if len(shared) >= 2:
                fusion_graph[i].add(j)
                fusion_graph[j].add(i)
    
    # Find connected components
    components = find_connected_components(fusion_graph)
    
    # Filter for exactly 4 rings
    fused_systems = []
    system_id = 0
    
    for component in components:
        if len(component) == 4:
            ring_indices = sorted(list(component))
            fused_systems.append({
                'system_id': system_id,
                'ring_indices': ring_indices,
                'num_rings': 4,
                'detection_method': 'topology'
            })
            system_id += 1
    
    return fused_systems


def calculate_frame_score(result):
    """
    Calculate quality score for a frame.
    Higher score = better geometry.
    """
    if not result['arg410_details'] or not result['lys414_details']:
        return float('-inf')
    
    arg410_min_dist = min(d['distance'] for d in result['arg410_details'])
    lys414_min_dist = min(d['distance'] for d in result['lys414_details'])
    
    arg410_height = min(d['height'] for d in result['arg410_details'])
    lys414_height = min(d['height'] for d in result['lys414_details'])
    
    def height_penalty(h):
        if h < 2.0:
            return 1.0
        elif h < 3.0:
            return 0.8
        elif h < 4.0:
            return 0.6
        return 0.4
    
    arg410_score = (1.0 / arg410_min_dist) * height_penalty(arg410_height)
    lys414_score = (1.0 / lys414_min_dist) * height_penalty(lys414_height)
    
    return arg410_score + lys414_score


def check_ring_in_fused_system(ring_idx, fused_systems):
    """
    Check if a specific ring belongs to any 4-fused ring system.
    
    Returns:
    --------
    tuple: (is_in_fused_system: bool, system_id: int or None)
    """
    for system in fused_systems:
        if ring_idx in system['ring_indices']:
            return True, system['system_id']
    return False, None


def is_cation_pi(charge_pos, ring_geo):
    """Check single cation-pi interaction."""
    dist = np.linalg.norm(charge_pos - ring_geo['center'])
    if dist >= CATION_PI_DIST_CUTOFF:
        return False, None
    
    proj, height = project_point_to_plane(
        charge_pos, ring_geo['center'], ring_geo['normal']
    )
    dist_to_center = np.linalg.norm(proj - ring_geo['center'])
    if dist_to_center > ring_geo['radius'] + RING_PADDING:
        return False, None
    
    return True, {'distance': dist, 'height': height}


def check_residue_cation_pi_with_fused(charge_pos, ligand_rings, ligand_ring_geos, fused_systems):
    """
    Check cation-pi interactions with distinction for fused ring systems.
    
    Returns:
    --------
    dict: {
        'has_any_interaction': bool,
        'interacting_rings': list of int,
        'has_fused_interaction': bool,
        'interacting_fused_systems': list of int (system_ids),
        'details': list of dict per interaction
    }
    """
    result = {
        'has_any_interaction': False,
        'interacting_rings': [],
        'has_fused_interaction': False,
        'interacting_fused_systems': [],
        'details': []
    }
    
    for ring_idx, (ring, ring_geo) in enumerate(zip(ligand_rings, ligand_ring_geos)):
        is_pi, details = is_cation_pi(charge_pos, ring_geo)
        
        if is_pi:
            result['has_any_interaction'] = True
            result['interacting_rings'].append(ring_idx)
            
            is_fused, system_id = check_ring_in_fused_system(ring_idx, fused_systems)
            
            interaction_detail = {
                'ring_idx': ring_idx,
                'distance': details['distance'],
                'height': details['height'],
                'is_fused': is_fused,
                'system_id': system_id
            }
            result['details'].append(interaction_detail)
            
            if is_fused:
                result['has_fused_interaction'] = True
                if system_id not in result['interacting_fused_systems']:
                    result['interacting_fused_systems'].append(system_id)
    
    return result


def find_target_residues(universe, target_residues):
    """
    Find specific residues (Arg410, Lys414) and return their charge center atoms.
    
    Parameters:
    -----------
    universe : MDAnalysis.Universe
    target_residues : dict, e.g., {'ARG410': 'CZ', 'LYS414': 'NZ'}
    
    Returns:
    --------
    dict: {residue_name: {'atom_index': int, 'position': np.array, 'resid': int}}
    """
    results = {}
    protein = universe.select_atoms('protein')
    
    for res_key, atom_name in target_residues.items():
        resname = res_key[:3]
        resid = int(res_key[3:])
        
        sel = protein.select_atoms(f'resname {resname} and resid {resid} and name {atom_name}')
        
        if len(sel) == 0:
            logger.warning(f"Could not find {res_key} {atom_name} atom")
            continue
        
        results[res_key] = {
            'atom_index': sel[0].index,
            'resid': resid,
            'resname': resname
        }
    
    return results


def analyze_frame_dual_interaction(universe, frame_idx, target_residues, ligand_rings, fused_systems):
    """
    Analyze single frame for dual Arg410+Lys414 cation-pi interactions.
    
    Returns:
    --------
    dict with interaction results
    """
    universe.trajectory[frame_idx]
    positions = universe.atoms.positions
    
    ligand_ring_geos = []
    for ring in ligand_rings:
        coords = positions[ring['atom_indices']]
        geo = get_ring_geometry(coords)
        ligand_ring_geos.append(geo)
    
    results = {
        'frame': frame_idx,
        'arg410_has_any_pi': False,
        'lys414_has_any_pi': False,
        'both_any_pi': False,
        'arg410_has_fused_pi': False,
        'lys414_has_fused_pi': False,
        'both_fused_pi': False,
        'arg410_details': [],
        'lys414_details': [],
        'arg410_fused_details': [],
        'lys414_fused_details': []
    }
    
    if 'ARG410' in target_residues:
        arg410_pos = positions[target_residues['ARG410']['atom_index']]
        arg410_result = check_residue_cation_pi_with_fused(
            arg410_pos, ligand_rings, ligand_ring_geos, fused_systems
        )
        results['arg410_has_any_pi'] = arg410_result['has_any_interaction']
        results['arg410_has_fused_pi'] = arg410_result['has_fused_interaction']
        results['arg410_details'] = arg410_result['details']
        results['arg410_fused_details'] = [d for d in arg410_result['details'] if d['is_fused']]
    
    if 'LYS414' in target_residues:
        lys414_pos = positions[target_residues['LYS414']['atom_index']]
        lys414_result = check_residue_cation_pi_with_fused(
            lys414_pos, ligand_rings, ligand_ring_geos, fused_systems
        )
        results['lys414_has_any_pi'] = lys414_result['has_any_interaction']
        results['lys414_has_fused_pi'] = lys414_result['has_fused_interaction']
        results['lys414_details'] = lys414_result['details']
        results['lys414_fused_details'] = [d for d in lys414_result['details'] if d['is_fused']]
    
    results['both_any_pi'] = results['arg410_has_any_pi'] and results['lys414_has_any_pi']
    results['both_fused_pi'] = results['arg410_has_fused_pi'] and results['lys414_has_fused_pi']
    
    return results


def process_trajectory_dual(traj_path, system_name):
    """
    Process trajectory and find BEST frames with dual interactions.
    Returns single best frame per condition (any-ring and fused-ring).
    """
    pdb_file = os.path.join(traj_path, 'v1.pdb')
    xtc_file = os.path.join(traj_path, 'v1.xtc')
    
    logger.info(f"Loading trajectory for {system_name}...")
    u = mda.Universe(pdb_file, xtc_file)
    
    logger.info("Detecting ligand aromatic rings...")
    ligand_rings = detect_aromatic_rings_ligand(u)
    logger.info(f"  Found {len(ligand_rings)} ligand aromatic rings")
    
    if len(ligand_rings) == 0:
        logger.warning(f"  No ligand aromatic rings found for {system_name}")
        return None
    
    logger.info("Finding 4-fused ring systems (topology-based)...")
    fused_systems = find_4fused_ring_systems_topology(ligand_rings)
    if fused_systems:
        logger.info(f"  Found {len(fused_systems)} 4-fused ring system(s)")
        for fs in fused_systems:
            logger.info(f"    System {fs['system_id']}: rings {fs['ring_indices']}")
    else:
        logger.info("  No 4-fused ring systems detected")
    
    target_residues = find_target_residues(u, {'ARG410': 'CZ', 'LYS414': 'NZ'})
    
    if 'ARG410' not in target_residues or 'LYS414' not in target_residues:
        logger.warning(f"Missing Arg410 or Lys414 in {system_name}")
        return None
    
    n_frames = len(u.trajectory)
    logger.info(f"Analyzing {n_frames} frames...")
    
    best_any_ring_frame = None
    best_any_ring_score = float('-inf')
    any_ring_count = 0
    
    best_fused_ring_frame = None
    best_fused_ring_score = float('-inf')
    fused_ring_count = 0
    
    for frame_idx in range(n_frames):
        if frame_idx % 500 == 0:
            logger.info(f"  Frame {frame_idx}/{n_frames}...")
        
        result = analyze_frame_dual_interaction(
            u, frame_idx, target_residues, ligand_rings, fused_systems
        )
        
        if result['both_any_pi']:
            any_ring_count += 1
            score = calculate_frame_score(result)
            if score > best_any_ring_score:
                best_any_ring_score = score
                best_any_ring_frame = result
        
        if result['both_fused_pi']:
            fused_ring_count += 1
            score = calculate_frame_score(result)
            if score > best_fused_ring_score:
                best_fused_ring_score = score
                best_fused_ring_frame = result
    
    logger.info(f"Found {any_ring_count} frames with dual any-ring interactions")
    logger.info(f"Found {fused_ring_count} frames with dual 4-fused ring interactions")
    
    ligand_ring_centers = []
    for ring in ligand_rings:
        coords = u.atoms[ring['atom_indices']].positions
        center = np.mean(coords, axis=0)
        ligand_ring_centers.append(center.tolist())
    
    return {
        'best_any_ring_frame': best_any_ring_frame,
        'best_fused_ring_frame': best_fused_ring_frame,
        'any_ring_count': any_ring_count,
        'fused_ring_count': fused_ring_count,
        'fused_systems': fused_systems,
        'total_frames': n_frames,
        'universe': u,
        'ligand_rings': ligand_rings,
        'ligand_ring_centers': ligand_ring_centers
    }


def extract_frames_to_pdb(universe, frame_result, output_path, system_name, 
                          interaction_type='any', fused_systems=None, ligand_ring_centers=None):
    """Extract single best frame to a single-model PDB file."""
    if frame_result is None:
        return
    
    pdb_lines = []
    
    frame_idx = frame_result['frame']
    universe.trajectory[frame_idx]
    
    pdb_lines.append(f"REMARK FRAME: {frame_idx}")
    pdb_lines.append(f"REMARK SYSTEM: {system_name}")
    pdb_lines.append(f"REMARK INTERACTION_TYPE: {interaction_type}")
    
    arg410_rings = [d['ring_idx'] for d in frame_result['arg410_details']]
    lys414_rings = [d['ring_idx'] for d in frame_result['lys414_details']]
    pdb_lines.append(f"REMARK ARG410 interacts with rings: {arg410_rings}")
    pdb_lines.append(f"REMARK LYS414 interacts with rings: {lys414_rings}")
    
    if ligand_ring_centers:
        pdb_lines.append("REMARK RING_CENTERS:")
        for i, center in enumerate(ligand_ring_centers):
            pdb_lines.append(f"REMARK   Ring {i}: ({center[0]:.3f}, {center[1]:.3f}, {center[2]:.3f})")
    
    if fused_systems:
        for fs in fused_systems:
            pdb_lines.append(f"REMARK Fused ring system {fs['system_id']}: rings {fs['ring_indices']}")
            if 'ring_centers' in fs:
                for i, rc in enumerate(fs['ring_centers']):
                    pdb_lines.append(f"REMARK   Fused center {i}: ({rc[0]:.3f}, {rc[1]:.3f}, {rc[2]:.3f})")
    
    arg410_fused = any(d['is_fused'] for d in frame_result['arg410_details'])
    lys414_fused = any(d['is_fused'] for d in frame_result['lys414_details'])
    pdb_lines.append(f"REMARK ARG410 fused interaction: {arg410_fused}")
    pdb_lines.append(f"REMARK LYS414 fused interaction: {lys414_fused}")
    
    arg410_min_dist = min(d['distance'] for d in frame_result['arg410_details'])
    lys414_min_dist = min(d['distance'] for d in frame_result['lys414_details'])
    pdb_lines.append(f"REMARK ARG410 min distance: {arg410_min_dist:.2f}")
    pdb_lines.append(f"REMARK LYS414 min distance: {lys414_min_dist:.2f}")
    pdb_lines.append(f"REMARK FRAME_SCORE: {calculate_frame_score(frame_result):.4f}")
    
    protein_ligand = universe.select_atoms('protein or not resname HOH SOL WAT')
    if len(protein_ligand) == 0:
        protein_ligand = universe.select_atoms('all and not resname HOH SOL WAT')
    
    atom_serial = 1
    for atom in protein_ligand:
        line = f"ATOM  {atom_serial:5d} {atom.name:<4s} {atom.resname:>3s} {atom.segid}{atom.resid:4d}    {atom.position[0]:8.3f}{atom.position[1]:8.3f}{atom.position[2]:8.3f}{atom.occupancy:6.2f}{atom.bfactor:6.2f}          {atom.element:>2s}"
        pdb_lines.append(line)
        atom_serial += 1
    
    pdb_lines.append("END")
    
    with open(output_path, 'w') as f:
        f.write('\n'.join(pdb_lines) + '\n')
    
    logger.info(f"Saved {output_path} (frame {frame_idx})")


def generate_csv_output(all_results, output_csv):
    """Generate CSV with best frames per system."""
    csv_rows = []
    
    for system_name, results in all_results.items():
        total_frames = results['total_frames']
        fused_systems = results['fused_systems']
        
        if results['best_any_ring_frame']:
            frame_result = results['best_any_ring_frame']
            row = {
                'system': system_name,
                'frame': frame_result['frame'],
                'interaction_type': 'any_ring',
                'arg410_ring': ','.join(str(d['ring_idx']) for d in frame_result['arg410_details']),
                'arg410_distance': ','.join(f"{d['distance']:.2f}" for d in frame_result['arg410_details']),
                'arg410_height': ','.join(f"{d['height']:.2f}" for d in frame_result['arg410_details']),
                'arg410_is_fused': any(d['is_fused'] for d in frame_result['arg410_details']),
                'lys414_ring': ','.join(str(d['ring_idx']) for d in frame_result['lys414_details']),
                'lys414_distance': ','.join(f"{d['distance']:.2f}" for d in frame_result['lys414_details']),
                'lys414_height': ','.join(f"{d['height']:.2f}" for d in frame_result['lys414_details']),
                'lys414_is_fused': any(d['is_fused'] for d in frame_result['lys414_details']),
                'fused_system_id': '',
                'frame_score': f"{calculate_frame_score(frame_result):.4f}",
                'total_frames': total_frames,
                'matching_frames': results['any_ring_count'],
                'fraction': results['any_ring_count'] / total_frames if total_frames > 0 else 0
            }
            csv_rows.append(row)
        
        if results['best_fused_ring_frame']:
            frame_result = results['best_fused_ring_frame']
            fused_ids = set()
            for d in frame_result['arg410_fused_details'] + frame_result['lys414_fused_details']:
                if d['system_id'] is not None:
                    fused_ids.add(d['system_id'])
            
            row = {
                'system': system_name,
                'frame': frame_result['frame'],
                'interaction_type': '4fused_ring',
                'arg410_ring': ','.join(str(d['ring_idx']) for d in frame_result['arg410_fused_details']),
                'arg410_distance': ','.join(f"{d['distance']:.2f}" for d in frame_result['arg410_fused_details']),
                'arg410_height': ','.join(f"{d['height']:.2f}" for d in frame_result['arg410_fused_details']),
                'arg410_is_fused': True,
                'lys414_ring': ','.join(str(d['ring_idx']) for d in frame_result['lys414_fused_details']),
                'lys414_distance': ','.join(f"{d['distance']:.2f}" for d in frame_result['lys414_fused_details']),
                'lys414_height': ','.join(f"{d['height']:.2f}" for d in frame_result['lys414_fused_details']),
                'lys414_is_fused': True,
                'fused_system_id': ','.join(str(i) for i in sorted(fused_ids)),
                'frame_score': f"{calculate_frame_score(frame_result):.4f}",
                'total_frames': total_frames,
                'matching_frames': results['fused_ring_count'],
                'fraction': results['fused_ring_count'] / total_frames if total_frames > 0 else 0
            }
            csv_rows.append(row)
    
    if csv_rows:
        df = pd.DataFrame(csv_rows)
        df.to_csv(output_csv, index=False)
        logger.info(f"Saved {len(csv_rows)} results to {output_csv}")


def find_systems(trj_dir):
    """Find all systems with fp/v1.pdb and v1.xtc files."""
    systems = []
    
    if os.path.isdir(trj_dir):
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
    parser = argparse.ArgumentParser(
        description='Detect dual Arg410+Lys414 cation-pi interactions'
    )
    parser.add_argument('--trj-dir', default='trj',
                        help='Directory containing trajectory systems')
    parser.add_argument('--output-dir', default='.',
                        help='Output directory for PDB and CSV files')
    parser.add_argument('--system', help='Process single system by name')
    parser.add_argument('--pdb', help='PDB file for single system')
    parser.add_argument('--xtc', help='XTC file for single system')
    parser.add_argument('--min-fused-rings', type=int, default=4,
                        help='Minimum rings in fused system')
    parser.add_argument('--max-fused-rings', type=int, default=4,
                        help='Maximum rings in fused system')
    parser.add_argument('-v', '--verbose', action='store_true',
                        help='Verbose output')
    
    args = parser.parse_args()
    
    if args.verbose:
        logging.getLogger().setLevel(logging.DEBUG)
    
    os.makedirs(args.output_dir, exist_ok=True)
    
    all_results = {}
    
    if args.system and args.pdb and args.xtc:
        logger.info(f"Processing single system: {args.system}")
        fp_dir = os.path.dirname(args.pdb)
        results = process_trajectory_dual(fp_dir, args.system)
        
        if results:
            all_results[args.system] = results
            
            if results['best_any_ring_frame']:
                output_pdb = os.path.join(args.output_dir, f"{args.system}_dual_pi_any_ring.pdb")
                extract_frames_to_pdb(
                    results['universe'], results['best_any_ring_frame'],
                    output_pdb, args.system, 'any', results['fused_systems'],
                    results['ligand_ring_centers']
                )
            
            if results['best_fused_ring_frame'] and results['fused_systems']:
                output_pdb = os.path.join(args.output_dir, f"{args.system}_dual_pi_4fused.pdb")
                extract_frames_to_pdb(
                    results['universe'], results['best_fused_ring_frame'],
                    output_pdb, args.system, 'fused', results['fused_systems'],
                    results['ligand_ring_centers']
                )
    else:
        logger.info(f"Finding systems in {args.trj_dir}...")
        systems = find_systems(args.trj_dir)
        logger.info(f"Found {len(systems)} systems")
        
        if len(systems) == 0:
            logger.error("No valid systems found. Exiting.")
            return
        
        for system in systems:
            logger.info(f"\n{'='*60}")
            logger.info(f"Processing {system}...")
            logger.info(f"{'='*60}")
            
            fp_path = os.path.join(args.trj_dir, system, 'fp')
            
            try:
                results = process_trajectory_dual(fp_path, system)
                
                if results:
                    all_results[system] = results
                    
                    if results['best_any_ring_frame']:
                        output_pdb = os.path.join(args.output_dir, f"{system}_dual_pi_any_ring.pdb")
                        extract_frames_to_pdb(
                            results['universe'], results['best_any_ring_frame'],
                            output_pdb, system, 'any', results['fused_systems'],
                            results['ligand_ring_centers']
                        )
                    
                    if results['best_fused_ring_frame'] and results['fused_systems']:
                        output_pdb = os.path.join(args.output_dir, f"{system}_dual_pi_4fused.pdb")
                        extract_frames_to_pdb(
                            results['universe'], results['best_fused_ring_frame'],
                            output_pdb, system, 'fused', results['fused_systems'],
                            results['ligand_ring_centers']
                        )
            
            except Exception as e:
                logger.error(f"Error processing {system}: {e}")
                import traceback
                traceback.print_exc()
                continue
    
    if all_results:
        output_csv = os.path.join(args.output_dir, 'dual_arg410_lys414_pi_frames.csv')
        generate_csv_output(all_results, output_csv)
        
        logger.info("\n" + "="*60)
        logger.info("Summary")
        logger.info("="*60)
        for system, results in all_results.items():
            logger.info(f"\n{system}:")
            logger.info(f"  Total frames: {results['total_frames']}")
            logger.info(f"  4-fused ring systems: {len(results['fused_systems'])}")
            logger.info(f"  Dual any-ring frames: {results['any_ring_count']} "
                       f"({100*results['any_ring_count']/results['total_frames']:.2f}%)")
            logger.info(f"  Dual 4-fused frames: {results['fused_ring_count']} "
                       f"({100*results['fused_ring_count']/results['total_frames']:.2f}%)")
            if results['best_any_ring_frame']:
                logger.info(f"  Best any-ring frame: {results['best_any_ring_frame']['frame']} "
                           f"(score: {calculate_frame_score(results['best_any_ring_frame']):.4f})")
            if results['best_fused_ring_frame']:
                logger.info(f"  Best fused-ring frame: {results['best_fused_ring_frame']['frame']} "
                           f"(score: {calculate_frame_score(results['best_fused_ring_frame']):.4f})")
    
    logger.info("\nProcessing complete!")


if __name__ == '__main__':
    main()
