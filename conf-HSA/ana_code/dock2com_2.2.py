#!/usr/bin/env python3
"""
dock2com_2.2.py
===============
One-step workflow: Select best docking pose and create MD-ready system files.
Extended with SASA-based selection metric, halogen support, and CHARMM CGenFF support.

Combines functionality from:
- dock2com_1.py (SDF parsing, model selection, GRO conversion)
- quick_sasa.py (SASA calculation via VDW+Gaussian approximation)

New Features in dock2com_2.2.py:
- CHARMM CGenFF ffbonded.itp support for ligand-specific parameters
- Auto-detection of *_ffbonded.itp from ligand ITP directory
- Proper inclusion of ffbonded.itp in system topology
- Automatic generation of ligand position restraints (posre_lig.itp)
- Ligand heavy atom position restraints for equilibration (under POSRES flag)

New Features in dock2com_2.1.py:
- Fixed halogen handling (Br, Cl, I, F) in atom type recognition
- Proper element inference for CHARMM CGenFF halogen types (BRGR1, CLGR1, etc.)

New Features in dock2com_2.py:
- SASA as a selection metric (--metric sasa)
- Lower SASA indicates better ligand burial within binding site
- Pre-selection filter: models with positive minimizedAffinity are excluded

Pre-selection Filtering
-----------------------
Models with minimizedAffinity > 0 are automatically excluded from selection.
Positive affinity values indicate unfavorable binding and are filtered out
before any metric-based ranking is performed.

Position Restraints
-------------------
Ligand position restraints are generated for heavy atoms only (hydrogens excluded).
The ligand ITP is modified to include posre_lig.itp under the POSRES flag.
Use define = -DPOSRES in MDP files to restrain both receptor and ligand during equilibration.

Example Usage
-------------
    # Basic usage
    python dock2com_2.2.py -i lig.itp -s hsa*-phe_sssD.sdf -r topol.top -t lig.mol2

    # Select by SASA (lower = more buried = better)
    python dock2com_2.2.py -i lig.itp -s hsa*-phe_sssD.sdf -r topol.top -t lig.mol2 --metric sasa

    # List all models with scores
    python dock2com_2.2.py -i lig.itp -s hsa*-phe_sssD.sdf --list-models
    
    # With explicit ffbonded.itp (CHARMM CGenFF)
    python dock2com_2.2.py -i lig.itp -s rec-lig.sdf -r topol.top --lig-ffbonded lig_ffbonded.itp
"""

import argparse
import glob as glob_module
import os
import re
import shutil
import sys
import textwrap
from collections import defaultdict, Counter

import numpy as np
import networkx as nx
from networkx.algorithms import isomorphism

LOWER_IS_BETTER = {"minimizedAffinity", "sasa"}
DEFAULT_METRIC = "minimizedAffinity"

def is_model_valid(model):
    """
    Check if model passes pre-selection filters.
    
    Filters:
    - minimizedAffinity must be <= 0 (positive values indicate unfavorable binding)
    """
    if "minimizedAffinity" in model["scores"]:
        affinity = model["scores"]["minimizedAffinity"]
        if np.isnan(affinity):
            return True
        if affinity > 0:
            return False
    return True
DEFAULT_REC_GRO_PATTERN = "{prefix}.pdb_ali.gro"

DEFAULT_PROBE_RADIUS = 1.4
SASA_GAUSSIAN_LAMBDA = 0.5
SASA_GAUSSIAN_SIGMA = 1.5
SASA_NEIGHBOR_CUTOFF = 10.0

VDW_RADII_ANG = {
    'H': 1.20, 'He': 1.40, 'C': 1.70, 'N': 1.55, 'O': 1.52, 'F': 1.47,
    'Ne': 1.54, 'Na': 2.27, 'Mg': 1.73, 'Si': 2.10, 'P': 1.80, 'S': 1.80,
    'Cl': 1.75, 'Ar': 1.88, 'K': 2.75, 'Ga': 1.87, 'As': 1.85, 'Se': 1.90,
    'Br': 1.85, 'Kr': 2.02, 'In': 1.93, 'Sn': 2.17, 'Te': 2.06, 'I': 1.98,
    'Xe': 2.16, 'Tl': 1.96, 'Pb': 2.02, 'Be': 1.53, 'B': 1.92, 'Al': 1.84,
    'Ca': 2.31, 'Ge': 2.11, 'Rb': 3.03, 'Sr': 2.50, 'Sb': 2.06, 'Cs': 3.43,
    'Ba': 2.68, 'Bi': 2.07, 'Po': 1.97, 'At': 2.02, 'Rn': 2.20, 'Fr': 3.48,
    'Ra': 2.83, 'Li': 1.81, 'Mn': 2.00, 'Fe': 2.00, 'Co': 1.95, 'Ni': 1.63,
    'Cu': 1.40, 'Zn': 1.39, 'Ru': 2.05, 'Rh': 2.00, 'Pd': 1.63, 'Ag': 1.72,
    'Pt': 1.72, 'Au': 1.66, 'Eu': 2.05,
}
VDW_DEFAULT_RADIUS = 1.70

HALOGEN_ATOM_TYPES = {
    'br': 'Br', 'BR': 'Br',
    'cl': 'Cl', 'CL': 'Cl',
    'fl': 'F', 'FL': 'F',
    'f': 'F', 'F': 'F',
    'i': 'I', 'I': 'I',
}


def infer_element_from_atom_name(atom_name):
    special_map = {
        "CA": "C", "CB": "C", "CD": "C", "CD1": "C", "CD2": "C",
        "CG": "C", "CG1": "C", "CG2": "C", "CE": "C", "CE1": "C",
        "CE2": "C", "CE3": "C", "CZ": "C", "CZ2": "C", "CZ3": "C",
        "CH2": "C",
        "NH": "N", "NH1": "N", "NH2": "N", "ND1": "N", "ND2": "N",
        "NE": "N", "NE1": "N", "NE2": "N", "NZ": "N",
        "OH": "O", "OD": "O", "OD1": "O", "OD2": "O", "OE": "O",
        "OE1": "O", "OE2": "O", "OG": "O", "OG1": "O", "OG2": "O",
        "OXT": "O",
        "SD": "S", "SG": "S",
        "EU": "Eu", "FE": "Fe", "ZN": "Zn", "MG": "Mg",
        "MN": "Mn", "CU": "Cu", "NI": "Ni", "CO": "Co",
        "PT": "Pt", "AU": "Au", "AG": "Ag",
        "BR": "Br", "BR": "Br", "CL": "Cl", "CL": "Cl",
    }
    
    upper = atom_name.upper().strip()
    
    if upper in special_map:
        return special_map[upper]
    
    if upper.startswith('H'):
        return "H"
    
    if len(upper) >= 2:
        two_char = upper[:2]
        if two_char.capitalize() in VDW_RADII_ANG:
            return two_char.capitalize()
        if two_char.upper() in VDW_RADII_ANG:
            return two_char.upper()
        if two_char[1].isdigit():
            one_char = two_char[0]
            if one_char in VDW_RADII_ANG:
                return one_char
    
    if upper and upper[0] in VDW_RADII_ANG:
        return upper[0]
    
    return upper[0] if upper else "C"


def calculate_pose_sasa(ligand_atoms, receptor_atoms, probe_radius=1.4):
    if not ligand_atoms:
        return 0.0
    
    all_coords = np.array(
        [[a["x"], a["y"], a["z"]] for a in ligand_atoms + receptor_atoms],
        dtype=np.float32
    )
    
    try:
        from scipy.spatial import cKDTree
        tree = cKDTree(all_coords)
    except ImportError:
        tree = None
    
    n_ligand = len(ligand_atoms)
    sasa_per_atom = np.zeros(n_ligand, dtype=np.float32)
    
    for i, atom in enumerate(ligand_atoms):
        coord = np.array([atom["x"], atom["y"], atom["z"]], dtype=np.float32)
        elem = atom.get("el", "C")
        
        r_vdw = VDW_RADII_ANG.get(elem, VDW_RADII_ANG.get(elem.capitalize(), VDW_DEFAULT_RADIUS))
        r_total = r_vdw + probe_radius
        
        neighbor_sum = 0.0
        
        if tree is not None:
            neighbors = tree.query_ball_point(coord, SASA_NEIGHBOR_CUTOFF)
            for j in neighbors:
                if j == i:
                    continue
                d = np.linalg.norm(all_coords[j] - coord)
                neighbor_sum += np.exp(-(d ** 2) / (2 * SASA_GAUSSIAN_SIGMA ** 2))
        else:
            for j, other_coord in enumerate(all_coords):
                if j == i:
                    continue
                d = np.linalg.norm(other_coord - coord)
                if d < SASA_NEIGHBOR_CUTOFF:
                    neighbor_sum += np.exp(-(d ** 2) / (2 * SASA_GAUSSIAN_SIGMA ** 2))
        
        exposure = np.exp(-SASA_GAUSSIAN_LAMBDA * neighbor_sum)
        sasa_per_atom[i] = 4 * np.pi * (r_total ** 2) * exposure
    
    total_sasa_nm2 = np.sum(sasa_per_atom) / 100.0
    
    return float(total_sasa_nm2)


def parse_gro_for_sasa(gro_path):
    if not os.path.exists(gro_path):
        raise FileNotFoundError(f"GRO file not found: {gro_path}")
    
    atoms = []
    
    with open(gro_path) as f:
        lines = f.readlines()
    
    for line in lines[2:-1]:
        if len(line) < 44:
            continue
        
        atom_name = line[10:15].strip()
        
        try:
            x_nm = float(line[20:28])
            y_nm = float(line[28:36])
            z_nm = float(line[36:44])
        except ValueError:
            continue
        
        x_ang = x_nm * 10.0
        y_ang = y_nm * 10.0
        z_ang = z_nm * 10.0
        
        element = infer_element_from_atom_name(atom_name)
        
        atoms.append({
            "el": element,
            "x": x_ang,
            "y": y_ang,
            "z": z_ang,
        })
    
    return atoms


def parse_sdf(path):
    with open(path) as fh:
        content = fh.read()
    
    models = []
    records = content.split("$$$$")
    for rec_idx, rec in enumerate(records):
        rec = rec.strip()
        if not rec:
            continue
        lines = rec.split("\n")
        if len(lines) < 4:
            continue
        
        title = lines[0].strip()
        counts_line = lines[3]
        try:
            ac = int(counts_line[0:3])
            bc = int(counts_line[3:6])
        except ValueError:
            continue
        
        atoms = []
        for i in range(4, 4 + ac):
            if i >= len(lines):
                break
            ln = lines[i]
            el = ln[31:34].strip()
            if el == "*":
                el = "EU"
            atoms.append({
                "idx": i - 3,
                "el": el,
                "x": float(ln[0:10]),
                "y": float(ln[10:20]),
                "z": float(ln[20:30]),
            })
        
        bonds = []
        for i in range(4 + ac, 4 + ac + bc):
            if i >= len(lines):
                break
            ln = lines[i]
            bonds.append({
                "a1": int(ln[0:3]),
                "a2": int(ln[3:6]),
                "type": int(ln[6:9]),
            })
        
        scores = {}
        for j, ln in enumerate(lines):
            if ln.startswith("> <"):
                key = ln.strip()[3:].rstrip(">").strip()
                if j + 1 < len(lines):
                    try:
                        scores[key] = float(lines[j + 1].strip())
                    except ValueError:
                        pass
        
        models.append({
            "model_num": rec_idx + 1,
            "title": title,
            "atoms": atoms,
            "bonds": bonds,
            "scores": scores,
            "raw_lines": lines,
        })
    
    return models


def collect_sdf_files(sdf_patterns):
    missing = []
    valid_paths = []
    
    for path in sdf_patterns:
        if os.path.isfile(path):
            valid_paths.append(path)
        else:
            missing.append(path)
    
    if missing:
        raise FileNotFoundError(f"SDF file(s) not found: {', '.join(missing)}")
    
    return sorted(valid_paths)


def select_best_across_sdfiles(
    sdf_paths,
    metric=DEFAULT_METRIC,
    rec_gro_pattern=None,
    rec_dir=None,
    probe_radius=DEFAULT_PROBE_RADIUS,
):
    best_model = None
    best_score = None
    best_file = None
    reverse = metric not in LOWER_IS_BETTER
    
    fallback_to_affinity = False
    original_metric = metric
    
    for sdf_path in sdf_paths:
        print(f"Parsing SDF:        {sdf_path}")
        models = parse_sdf(sdf_path)
        print(f"  {len(models)} models found")
        
        receptor_atoms = None
        rec_gro = None
        
        if metric == "sasa":
            if rec_gro_pattern is None:
                rec_gro_pattern = DEFAULT_REC_GRO_PATTERN
            
            rec_gro, prefix = derive_receptor_gro_from_sdf(sdf_path, rec_gro_pattern, rec_dir)
            
            if os.path.exists(rec_gro):
                try:
                    print(f"  Loading receptor for SASA: {rec_gro}")
                    receptor_atoms = parse_gro_for_sasa(rec_gro)
                    print(f"  Receptor atoms: {len(receptor_atoms)}")
                except Exception as e:
                    print(f"  WARNING: Failed to load receptor: {e}")
                    receptor_atoms = None
            else:
                print(f"  WARNING: Receptor GRO not found: {rec_gro}")
        
        if metric == "sasa" and receptor_atoms is None:
            print(f"  WARNING: Cannot calculate SASA without receptor.")
            print(f"  Falling back to 'minimizedAffinity' metric.")
            metric = "minimizedAffinity"
            reverse = metric not in LOWER_IS_BETTER
            fallback_to_affinity = True
        
        for model in models:
            if not is_model_valid(model):
                continue
            
            if original_metric == "sasa" and receptor_atoms is not None:
                if "sasa" not in model["scores"]:
                    try:
                        sasa = calculate_pose_sasa(model["atoms"], receptor_atoms, probe_radius)
                        model["scores"]["sasa"] = sasa
                    except Exception as e:
                        print(f"  WARNING: SASA calculation failed for model {model['model_num']}: {e}")
                        model["scores"]["sasa"] = float('nan')
            
            if metric not in model["scores"]:
                continue
            
            score = model["scores"][metric]
            
            if np.isnan(score):
                continue
            
            if best_score is None:
                best_score = score
                best_model = model
                best_file = sdf_path
            elif reverse:
                if score > best_score:
                    best_score = score
                    best_model = model
                    best_file = sdf_path
            else:
                if score < best_score:
                    best_score = score
                    best_model = model
                    best_file = sdf_path
    
    if best_model is None:
        raise ValueError(f"Metric '{metric}' not found in any model across all files.")
    
    if fallback_to_affinity:
        print(f"\nNOTE: Used 'minimizedAffinity' as fallback metric.")
    
    return best_model, best_file, best_model["model_num"]


def list_models_across_sdfiles(sdf_paths, metrics=None, rec_gro_pattern=None, rec_dir=None, probe_radius=DEFAULT_PROBE_RADIUS):
    if metrics is None:
        all_keys = []
        seen = set()
        for sdf_path in sdf_paths:
            models = parse_sdf(sdf_path)
            for m in models:
                for k in m["scores"]:
                    if k not in seen:
                        all_keys.append(k)
                        seen.add(k)
        metrics = all_keys
    
    max_fn_len = max(len(os.path.basename(p)) for p in sdf_paths)
    fn_width = max(max_fn_len, 8)
    
    header = f"{('File'):>{fn_width}}  {'Model':>6}" + "".join(f"  {k:>16}" for k in metrics)
    print(header)
    print("-" * len(header))
    
    for sdf_path in sdf_paths:
        fn = os.path.basename(sdf_path)
        models = parse_sdf(sdf_path)
        
        receptor_atoms = None
        if "sasa" in metrics:
            if rec_gro_pattern is None:
                rec_gro_pattern = DEFAULT_REC_GRO_PATTERN
            rec_gro, _ = derive_receptor_gro_from_sdf(sdf_path, rec_gro_pattern, rec_dir)
            if os.path.exists(rec_gro):
                try:
                    receptor_atoms = parse_gro_for_sasa(rec_gro)
                except Exception:
                    receptor_atoms = None
        
        for m in models:
            if not is_model_valid(m):
                continue
            
            if "sasa" in metrics and receptor_atoms is not None and "sasa" not in m["scores"]:
                try:
                    sasa = calculate_pose_sasa(m["atoms"], receptor_atoms, probe_radius)
                    m["scores"]["sasa"] = sasa
                except Exception:
                    m["scores"]["sasa"] = float('nan')
            
            row = f"{fn:>{fn_width}}  {m['model_num']:>6}"
            for k in metrics:
                val = m["scores"].get(k, float("nan"))
                row += f"  {val:>16.4f}"
            print(row)


def parse_mol2(path):
    atoms = []
    bonds = []
    header_lines = []
    section = None
    mol_line_count = 0
    
    with open(path) as fh:
        for raw in fh:
            line = raw.rstrip("\n")
            if line.startswith("@<TRIPOS>"):
                section = line.strip()
                mol_line_count = 0
                continue
            
            if section == "@<TRIPOS>MOLECULE":
                mol_line_count += 1
                if mol_line_count <= 5:
                    header_lines.append(line)
            
            elif section == "@<TRIPOS>ATOM":
                parts = line.split()
                if not parts:
                    continue
                atype = parts[5]
                el = _element(atype)
                atoms.append({
                    "idx": int(parts[0]),
                    "name": parts[1],
                    "x": float(parts[2]),
                    "y": float(parts[3]),
                    "z": float(parts[4]),
                    "type": atype,
                    "element": el,
                    "res_num": parts[6] if len(parts) > 6 else "1",
                    "res_name": parts[7] if len(parts) > 7 else "UNL1",
                    "charge": parts[8] if len(parts) > 8 else "0.0000",
                })
            
            elif section == "@<TRIPOS>BOND":
                parts = line.split()
                if len(parts) >= 4:
                    bonds.append({
                        "idx": int(parts[0]),
                        "a1": int(parts[1]),
                        "a2": int(parts[2]),
                        "type": parts[3],
                    })
    
    return atoms, bonds, header_lines


def _element(atom_type):
    el = atom_type.split(".")[0].upper()
    if el.startswith("EU"):
        el = "EU"
    if el in ("BR", "BR"):
        return "Br"
    if el in ("CL", "CL"):
        return "Cl"
    if el == "F":
        return "F"
    if el == "I":
        return "I"
    return el


def parse_itp(path):
    atoms = []
    bonds = []
    in_atoms = False
    in_bonds = False
    
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if line.startswith("["):
                section = line.strip().lower()
                in_atoms = "atoms" in section
                in_bonds = "bonds" in section
                continue
            
            if in_atoms and line and not line.startswith(";") and not line.startswith("#"):
                parts = line.split()
                if len(parts) >= 8:
                    try:
                        atoms.append({
                            "idx": int(parts[0]),
                            "type": parts[1],
                            "res_num": int(parts[2]),
                            "res_name": parts[3],
                            "atom_name": parts[4],
                            "cgnr": int(parts[5]),
                            "charge": float(parts[6]),
                            "mass": float(parts[7]),
                        })
                    except ValueError:
                        continue
            
            if in_bonds and line and not line.startswith(";") and not line.startswith("#"):
                parts = line.split()
                if len(parts) >= 2:
                    try:
                        bonds.append({
                            "a1": int(parts[0]),
                            "a2": int(parts[1]),
                        })
                    except ValueError:
                        continue
    
    return atoms, bonds


def _itp_element(atom_type):
    lower_type = atom_type.lower()
    
    if lower_type.startswith("br"):
        return "Br"
    if lower_type.startswith("cl"):
        return "Cl"
    if lower_type.startswith("fl") or lower_type == "f" or lower_type.startswith("fgr"):
        return "F"
    if lower_type.startswith("ir") or (lower_type == "i" and len(atom_type) <= 2):
        return "I"
    
    if lower_type.startswith("eu"):
        return "EU"
    
    type_map = {
        "ca": "C", "ce": "C", "c": "C", "c3": "C",
        "o": "O", "n": "N", "n2": "N",
        "ha": "H", "h1": "H", "h": "H",
    }
    
    return type_map.get(lower_type, atom_type[0].upper())


def write_gro(path, atoms, title="Generated", box=None, overwrite=False):
    if not overwrite:
        _no_overwrite(path)
    with open(path, "w") as fh:
        fh.write(f"{title}\n")
        fh.write(f"{len(atoms)}\n")
        for a in atoms:
            x_nm = a["x"] * 0.1
            y_nm = a["y"] * 0.1
            z_nm = a["z"] * 0.1
            fh.write(
                f"{a['res_num']:>5d}{a['res_name']:<5s}"
                f"{a['atom_name']:>5s}{a['idx']:>5d}"
                f"{x_nm:>8.3f}{y_nm:>8.3f}{z_nm:>8.3f}\n"
            )
        if box is not None:
            fh.write(f"{box[0]:>10.5f}{box[1]:>10.5f}{box[2]:>10.5f}\n")
        else:
            fh.write("   0.00000   0.00000   0.00000\n")


def parse_pdb(path):
    atoms = []
    with open(path) as fh:
        for line in fh:
            if line.startswith("ATOM") or line.startswith("HETATM"):
                try:
                    atom_name = line[12:16].strip()
                    res_name = line[17:20].strip()
                    res_num = int(line[22:26].strip())
                    x = float(line[30:38])
                    y = float(line[38:46])
                    z = float(line[46:54])
                    element = line[76:78].strip() if len(line) > 76 else atom_name[0]
                    atoms.append({
                        "idx": len(atoms) + 1,
                        "atom_name": atom_name,
                        "res_name": res_name,
                        "res_num": res_num,
                        "x": x,
                        "y": y,
                        "z": z,
                        "element": element,
                    })
                except (ValueError, IndexError):
                    continue
    return atoms


def pdb_to_gro(pdb_path, gro_path, box=None):
    atoms = parse_pdb(pdb_path)
    write_gro(gro_path, atoms, title=f"Converted from {os.path.basename(pdb_path)}", box=box, overwrite=True)
    return len(atoms)


def _build_graph(atoms, bonds, heavy_only=False):
    G = nx.Graph()
    incl = {a["idx"] for a in atoms if not heavy_only or a["element"] != "H"}
    for a in atoms:
        if a["idx"] in incl:
            G.add_node(a["idx"], element=a["element"])
    for b in bonds:
        if b["a1"] in incl and b["a2"] in incl:
            G.add_edge(b["a1"], b["a2"])
    return G


METAL_ELEMENTS = {"M", "EU", "FE", "ZN", "MG", "CA", "MN", "CO", "CU", "NI", "PT", "AU", "AG", "HG", "PB", "CD", "CR", "MO", "W", "V", "TI", "PD", "RH", "RU", "IR", "OS"}
HALOGEN_ELEMENTS = {"F", "CL", "BR", "I", "AT"}


def _is_metal(element):
    return element.upper() in METAL_ELEMENTS


def _is_halogen(element):
    return element.upper() in HALOGEN_ELEMENTS or element.capitalize() in {"Cl", "Br", "I", "F"}


def _node_match(n1, n2):
    e1 = n1["element"].upper()
    e2 = n2["element"].upper()
    if e1 == e2:
        return True
    if _is_metal(e1) and _is_metal(e2):
        return True
    return False


def find_isomorphism(G_ref, G_src):
    gm = isomorphism.GraphMatcher(G_ref, G_src, node_match=_node_match)
    if not gm.is_isomorphic():
        raise ValueError(
            "No graph isomorphism found between the two molecules.\n"
            "Check that both files represent the same molecular graph."
        )
    return next(gm.isomorphisms_iter())


def _vec(a, b):
    v = b - a
    norm = np.linalg.norm(v)
    if norm < 1e-10:
        raise ValueError("Zero-length vector encountered.")
    return v / norm


def _get_coord(atom):
    return np.array([atom["x"], atom["y"], atom["z"]], dtype=float)


def place_h_from_template(ref_heavy_atom, ref_h_atom, ref_heavy_nbrs,
                           new_heavy_coord, new_nbr_coords):
    P_heavy = _get_coord(ref_heavy_atom)
    P_h = _get_coord(ref_h_atom)
    dh = P_h - P_heavy
    
    if not ref_heavy_nbrs:
        new_pos = new_heavy_coord + dh
        return new_pos
    
    nbr_vecs_ref = [_vec(P_heavy, _get_coord(n)) for n in ref_heavy_nbrs]
    nbr_vecs_new = [_vec(new_heavy_coord, c) for c in new_nbr_coords]
    
    e1_ref = np.mean(nbr_vecs_ref, axis=0)
    norm_e1 = np.linalg.norm(e1_ref)
    if norm_e1 < 1e-8:
        e1_ref = nbr_vecs_ref[0]
    else:
        e1_ref /= norm_e1
    
    v2_ref = nbr_vecs_ref[0] - np.dot(nbr_vecs_ref[0], e1_ref) * e1_ref
    if np.linalg.norm(v2_ref) < 1e-8:
        if len(nbr_vecs_ref) > 1:
            v2_ref = nbr_vecs_ref[1] - np.dot(nbr_vecs_ref[1], e1_ref) * e1_ref
        if np.linalg.norm(v2_ref) < 1e-8:
            arb = np.array([1.0, 0.0, 0.0])
            if abs(np.dot(e1_ref, arb)) > 0.9:
                arb = np.array([0.0, 1.0, 0.0])
            v2_ref = arb - np.dot(arb, e1_ref) * e1_ref
    e2_ref = v2_ref / np.linalg.norm(v2_ref)
    e3_ref = np.cross(e1_ref, e2_ref)
    
    h_e1 = np.dot(dh, e1_ref)
    h_e2 = np.dot(dh, e2_ref)
    h_e3 = np.dot(dh, e3_ref)
    
    e1_new = np.mean(nbr_vecs_new, axis=0)
    norm_e1n = np.linalg.norm(e1_new)
    if norm_e1n < 1e-8:
        e1_new = nbr_vecs_new[0]
    else:
        e1_new /= norm_e1n
    
    v2_new = nbr_vecs_new[0] - np.dot(nbr_vecs_new[0], e1_new) * e1_new
    if np.linalg.norm(v2_new) < 1e-8:
        if len(nbr_vecs_new) > 1:
            v2_new = nbr_vecs_new[1] - np.dot(nbr_vecs_new[1], e1_new) * e1_new
        if np.linalg.norm(v2_new) < 1e-8:
            arb = np.array([1.0, 0.0, 0.0])
            if abs(np.dot(e1_new, arb)) > 0.9:
                arb = np.array([0.0, 1.0, 0.0])
            v2_new = arb - np.dot(arb, e1_new) * e1_new
    e2_new = v2_new / np.linalg.norm(v2_new)
    e3_new = np.cross(e1_new, e2_new)
    
    new_pos = new_heavy_coord + h_e1 * e1_new + h_e2 * e2_new + h_e3 * e3_new
    return new_pos


def sdf_pose_to_gro(itp_path, sdf_path, out_path,
                    mol2_template_path=None,
                    metric=DEFAULT_METRIC, model_num=None):
    _no_overwrite(out_path)
    
    print(f"Parsing ITP:         {itp_path}")
    itp_atoms, itp_bonds = parse_itp(itp_path)
    print(f"  {len(itp_atoms)} atoms, {len(itp_bonds)} bonds")
    
    mol2_atoms = None
    mol2_bonds = None
    if mol2_template_path:
        print(f"Parsing mol2 template: {mol2_template_path}")
        mol2_atoms, mol2_bonds, _ = parse_mol2(mol2_template_path)
        print(f"  {len(mol2_atoms)} atoms, {len(mol2_bonds)} bonds")
    
    print(f"Parsing SDF:        {sdf_path}")
    models = parse_sdf(sdf_path)
    print(f"  {len(models)} models found")
    
    model = models[model_num - 1] if model_num else _select_model(models, metric)
    if model_num is None:
        sc = model["scores"]
        print(
            f"  Selected model {model['model_num']}  "
            f"minimizedAffinity={sc.get('minimizedAffinity', float('nan')):.4f}  "
            f"CNNscore={sc.get('CNNscore', float('nan')):.4f}  "
            f"CNNaffinity={sc.get('CNNaffinity', float('nan')):.4f}"
        )
    else:
        print(f"  Using model {model_num}")
    
    sdf_atoms = model["atoms"]
    sdf_bonds = model["bonds"]
    
    sdf_by_idx = {a["idx"]: a for a in sdf_atoms}
    sdf_adj = defaultdict(list)
    for b in sdf_bonds:
        sdf_adj[b["a1"]].append(b["a2"])
        sdf_adj[b["a2"]].append(b["a1"])
    
    new_coords = {}
    
    if mol2_atoms:
        ref_atoms = mol2_atoms
        ref_bonds = mol2_bonds
        
        ref_heavy = [a for a in ref_atoms if a["element"] != "H"]
        ref_heavy_idx = {a["idx"] for a in ref_heavy}
        ref_heavy_bonds = [b for b in ref_bonds
                           if b["a1"] in ref_heavy_idx and b["a2"] in ref_heavy_idx]
        
        sdf_heavy = [{"idx": a["idx"], "element": a["el"]}
                    for a in sdf_atoms if a["el"] != "H"]
        sdf_heavy_bonds = [{"a1": b["a1"], "a2": b["a2"]}
                          for b in sdf_bonds
                          if sdf_by_idx[b["a1"]]["el"] != "H"
                          and sdf_by_idx[b["a2"]]["el"] != "H"]
        
        G_ref_h = _build_graph(ref_heavy, ref_heavy_bonds)
        G_sdf_h = _build_graph(sdf_heavy, sdf_heavy_bonds)
        mol2_to_sdf = find_isomorphism(G_ref_h, G_sdf_h)
        
        ref_by_idx = {a["idx"]: a for a in ref_atoms}
        itp_heavy = [a for a in itp_atoms if _itp_element(a["type"]) != "H"]
        itp_heavy_idx = {a["idx"] for a in itp_heavy}
        itp_heavy_bonds = [b for b in itp_bonds
                           if b["a1"] in itp_heavy_idx and b["a2"] in itp_heavy_idx]
        
        itp_heavy_for_graph = [{"idx": a["idx"], "element": _itp_element(a["type"]).upper()}
                               for a in itp_heavy]
        G_itp_h = _build_graph(itp_heavy_for_graph, itp_heavy_bonds)
        itp_to_mol2 = find_isomorphism(G_itp_h, G_ref_h)
        
        itp_to_sdf = {}
        for itp_idx, mol2_idx in itp_to_mol2.items():
            if mol2_idx in mol2_to_sdf:
                itp_to_sdf[itp_idx] = mol2_to_sdf[mol2_idx]
        
        ref_by_idx = {a["idx"]: a for a in ref_atoms}
        ref_adj = defaultdict(list)
        for b in ref_bonds:
            ref_adj[b["a1"]].append(b["a2"])
            ref_adj[b["a2"]].append(b["a1"])
        
        ref_to_itp = {v: k for k, v in itp_to_mol2.items()}
        
        itp_by_idx = {a["idx"]: a for a in itp_atoms}
        itp_adj = defaultdict(list)
        for b in itp_bonds:
            itp_adj[b["a1"]].append(b["a2"])
            itp_adj[b["a2"]].append(b["a1"])
        
        ref_h_by_parent = defaultdict(list)
        for a in ref_atoms:
            if a["element"] == "H":
                for nb_idx in ref_adj[a["idx"]]:
                    if ref_by_idx[nb_idx]["element"] != "H":
                        ref_h_by_parent[nb_idx].append(a)
        
        sdf_n_h_map = defaultdict(list)
        for a in sdf_atoms:
            if a["el"] == "N":
                for nb_idx in sdf_adj[a["idx"]]:
                    if sdf_by_idx[nb_idx]["el"] == "H":
                        sdf_n_h_map[a["idx"]].append(nb_idx)
        
        used_sdf_h = set()
        
        for itp_idx, sdf_idx in itp_to_sdf.items():
            sdf_a = sdf_by_idx[sdf_idx]
            new_coords[itp_idx] = (sdf_a["x"], sdf_a["y"], sdf_a["z"])
        
        for itp_a in itp_atoms:
            itp_idx = itp_a["idx"]
            if _itp_element(itp_a["type"]) == "H":
                parent_itp_idx = None
                for nb_idx in itp_adj[itp_idx]:
                    if _itp_element(itp_by_idx[nb_idx]["type"]) != "H":
                        parent_itp_idx = nb_idx
                        break
                
                if parent_itp_idx is None:
                    new_coords[itp_idx] = (0.0, 0.0, 0.0)
                    continue
                
                parent_ref_idx = itp_to_mol2.get(parent_itp_idx)
                if parent_ref_idx is None:
                    new_coords[itp_idx] = (0.0, 0.0, 0.0)
                    continue
                
                ref_h_list = ref_h_by_parent.get(parent_ref_idx, [])
                if not ref_h_list:
                    new_coords[itp_idx] = (0.0, 0.0, 0.0)
                    continue
                
                parent_ref = ref_by_idx[parent_ref_idx]
                
                if parent_ref["element"] == "N":
                    sdf_parent_idx = mol2_to_sdf.get(parent_ref_idx)
                    if sdf_parent_idx:
                        avail_sdf_h = [h for h in sdf_n_h_map[sdf_parent_idx]
                                       if h not in used_sdf_h]
                        if avail_sdf_h:
                            sdf_h = sdf_by_idx[avail_sdf_h[0]]
                            used_sdf_h.add(avail_sdf_h[0])
                            new_coords[itp_idx] = (sdf_h["x"], sdf_h["y"], sdf_h["z"])
                            continue
                
                ref_h = ref_h_list[0]
                ref_h_by_parent[parent_ref_idx] = ref_h_list[1:] if len(ref_h_list) > 1 else []
                
                if parent_itp_idx not in new_coords:
                    new_coords[itp_idx] = (ref_h["x"], ref_h["y"], ref_h["z"])
                    continue
                
                new_parent_coord = np.array(new_coords[parent_itp_idx])
                
                other_heavy_ref = [ref_by_idx[n] for n in ref_adj[parent_ref_idx]
                                   if n != ref_h["idx"] and ref_by_idx[n]["element"] != "H"]
                
                other_heavy_new_coords = []
                for oh in other_heavy_ref:
                    oh_itp_idx = ref_to_itp.get(oh["idx"])
                    if oh_itp_idx is not None and oh_itp_idx in new_coords:
                        other_heavy_new_coords.append(np.array(new_coords[oh_itp_idx]))
                    else:
                        other_heavy_new_coords.append(_get_coord(oh))
                
                new_h_xyz = place_h_from_template(
                    parent_ref, ref_h, other_heavy_ref,
                    new_parent_coord, other_heavy_new_coords,
                )
                new_coords[itp_idx] = tuple(new_h_xyz)
    else:
        itp_heavy = [a for a in itp_atoms if _itp_element(a["type"]) != "H"]
        itp_heavy_idx = {a["idx"] for a in itp_heavy}
        
        sdf_heavy = [{"idx": a["idx"], "element": a["el"]}
                    for a in sdf_atoms if a["el"] != "H"]
        sdf_heavy_bonds = [{"a1": b["a1"], "a2": b["a2"]}
                          for b in sdf_bonds
                          if sdf_by_idx[b["a1"]]["el"] != "H"
                          and sdf_by_idx[b["a2"]]["el"] != "H"]
        
        itp_elems = Counter(_itp_element(a["type"]) for a in itp_heavy)
        sdf_elems = Counter(a["element"] for a in sdf_heavy)
        
        if itp_elems != sdf_elems:
            raise ValueError(
                f"Element counts mismatch: ITP={dict(itp_elems)}, SDF={dict(sdf_elems)}"
            )
        
        itp_heavy_for_graph = [{"idx": a["idx"], "element": _itp_element(a["type"]).upper()}
                               for a in itp_heavy]
        itp_heavy_bonds = [{"a1": b["a1"], "a2": b["a2"]}
                          for b in itp_bonds
                          if b["a1"] in itp_heavy_idx and b["a2"] in itp_heavy_idx]
        
        G_itp = _build_graph(itp_heavy_for_graph, itp_heavy_bonds)
        G_sdf = _build_graph(sdf_heavy, sdf_heavy_bonds)
        
        mapping = find_isomorphism(G_itp, G_sdf)
        sdf_to_itp = {v: k for k, v in mapping.items()}
        
        for itp_idx, sdf_idx in mapping.items():
            sdf_a = sdf_by_idx[sdf_idx]
            new_coords[itp_idx] = (sdf_a["x"], sdf_a["y"], sdf_a["z"])
        
        itp_by_idx = {a["idx"]: a for a in itp_atoms}
        itp_adj = defaultdict(list)
        for b in itp_bonds:
            itp_adj[b["a1"]].append(b["a2"])
            itp_adj[b["a2"]].append(b["a1"])
        
        sdf_h_by_parent = defaultdict(list)
        for a in sdf_atoms:
            if a["el"] == "H":
                for nb_idx in sdf_adj[a["idx"]]:
                    if sdf_by_idx[nb_idx]["el"] != "H":
                        sdf_h_by_parent[nb_idx].append(a)
        
        for itp_a in itp_atoms:
            itp_idx = itp_a["idx"]
            if _itp_element(itp_a["type"]) == "H":
                parent_itp_idx = None
                for nb_idx in itp_adj[itp_idx]:
                    if _itp_element(itp_by_idx[nb_idx]["type"]) != "H":
                        parent_itp_idx = nb_idx
                        break
                
                if parent_itp_idx is None:
                    new_coords[itp_idx] = (0.0, 0.0, 0.0)
                    continue
                
                if parent_itp_idx not in new_coords:
                    new_coords[itp_idx] = (0.0, 0.0, 0.0)
                    continue
                
                parent_coord = np.array(new_coords[parent_itp_idx])
                parent_element = _itp_element(itp_by_idx[parent_itp_idx]["type"]).upper()
                
                sdf_parent_idx = mapping.get(parent_itp_idx)
                if sdf_parent_idx and sdf_parent_idx in sdf_h_by_parent:
                    sdf_h_list = sdf_h_by_parent[sdf_parent_idx]
                    if sdf_h_list:
                        sdf_h = sdf_h_list.pop(0)
                        sdf_h_atom = sdf_by_idx[sdf_h["idx"]]
                        new_coords[itp_idx] = (sdf_h_atom["x"], sdf_h_atom["y"], sdf_h_atom["z"])
                        continue
                
                neighbor_coords = []
                for nb_idx in itp_adj[parent_itp_idx]:
                    if nb_idx != itp_idx and _itp_element(itp_by_idx[nb_idx]["type"]) != "H":
                        if nb_idx in new_coords:
                            neighbor_coords.append(np.array(new_coords[nb_idx]))
                
                bond_len = 1.09
                if parent_element == "N":
                    bond_len = 1.01
                elif parent_element == "O":
                    bond_len = 0.96
                elif parent_element == "S":
                    bond_len = 1.34
                
                if neighbor_coords:
                    v_sum = np.zeros(3)
                    for nc in neighbor_coords:
                        v = _vec(parent_coord, nc)
                        v_sum += v
                    h_dir = -v_sum
                    norm = np.linalg.norm(h_dir)
                    if norm > 1e-8:
                        h_dir = h_dir / norm
                    else:
                        h_dir = np.array([1.0, 0.0, 0.0])
                    new_h_pos = parent_coord + bond_len * h_dir
                else:
                    new_h_pos = parent_coord + bond_len * np.array([1.0, 0.0, 0.0])
                
                new_coords[itp_idx] = tuple(new_h_pos)
    
    out_atoms = []
    for i, itp_a in enumerate(itp_atoms, 1):
        if itp_a["idx"] in new_coords:
            x, y, z = new_coords[itp_a["idx"]]
        else:
            x, y, z = 0.0, 0.0, 0.0
        
        out_atoms.append({
            "idx": i,
            "atom_name": itp_a["atom_name"],
            "res_num": itp_a["res_num"],
            "res_name": itp_a["res_name"],
            "x": x,
            "y": y,
            "z": z,
        })
    
    write_gro(out_path, out_atoms)
    print(f"Written: {out_path}")


def _select_model(models, metric=DEFAULT_METRIC):
    valid_models = [m for m in models if is_model_valid(m)]
    candidates = [m for m in valid_models if metric in m["scores"]]
    if not candidates:
        raise ValueError(f"Metric '{metric}' not found in any valid model.")
    reverse = metric not in LOWER_IS_BETTER
    return sorted(candidates, key=lambda m: m["scores"][metric], reverse=reverse)[0]


def extract_water_models(ff_dir):
    watermodels_file = os.path.join(ff_dir, 'watermodels.dat')
    water_models = set()
    
    if os.path.exists(watermodels_file):
        with open(watermodels_file, 'r') as f:
            for line in f:
                line = line.strip()
                if line and not line.startswith('#') and len(line.split()) >= 1:
                    name = line.split()[0]
                    water_models.add(name)
    
    return water_models


def extract_ff_paths_from_top(top_file):
    ff_path = None
    water_itp = None
    ions_itp = None
    ff_dir = None
    
    with open(top_file, 'r') as f:
        for line in f:
            if line.strip().startswith('#include') and '.ff/' in line:
                path = line.strip().split('"')[1] if '"' in line else line.strip()
                if 'forcefield.itp' in path:
                    ff_path = path
                    ff_dir = os.path.dirname(path)
                elif 'ions.itp' in path:
                    ions_itp = path
    
    water_models = extract_water_models(ff_dir) if ff_dir else set()
    
    with open(top_file, 'r') as f:
        for line in f:
            if line.strip().startswith('#include') and '.ff/' in line:
                path = line.strip().split('"')[1] if '"' in line else line.strip()
                if path.endswith('.itp') and 'forcefield' not in path and 'ions' not in path:
                    for wm in water_models:
                        if f'{wm}.itp' in path:
                            water_itp = path
                            break
                    if water_itp:
                        break
    
    return ff_path, water_itp, ions_itp


def get_moleculetype_name(top_file):
    with open(top_file, 'r') as f:
        lines = f.readlines()
    
    in_moleculetype = False
    for line in lines:
        stripped = line.strip()
        if stripped.startswith('[ moleculetype ]') or stripped == '[moleculetype]':
            in_moleculetype = True
            continue
        if in_moleculetype:
            if stripped.startswith('['):
                break
            if stripped and not stripped.startswith(';') and not stripped.startswith('#'):
                return stripped.split()[0]
    
    return None


def get_prefix(filename):
    basename = os.path.basename(filename)
    if '.' in basename:
        return basename.rsplit('.', 1)[0]
    return basename


def combine_coordinates(rec_gro, lig_gro, output_gro):
    with open(rec_gro, 'r') as f:
        rec_lines = f.readlines()
    
    with open(lig_gro, 'r') as f:
        lig_lines = f.readlines()
    
    rec_atom_count = int(rec_lines[1].strip())
    lig_atom_count = int(lig_lines[1].strip())
    total_atoms = rec_atom_count + lig_atom_count
    
    box = rec_lines[-1].strip()
    
    combined = []
    combined.append(rec_lines[0].strip() + '\n')
    combined.append(f' {total_atoms}\n')
    
    for line in rec_lines[2:-1]:
        combined.append(line)
    
    for line in lig_lines[2:-1]:
        combined.append(line)
    
    combined.append(box + '\n')
    
    with open(output_gro, 'w') as f:
        f.writelines(combined)
    
    print(f"Created {output_gro} with {total_atoms} atoms")


def extract_receptor_topology(rec_top, rec_itp):
    with open(rec_top, 'r') as f:
        lines = f.readlines()
    
    start_line = None
    for i, line in enumerate(lines):
        if '[ moleculetype ]' in line:
            start_line = i
            break
    
    if start_line is None:
        raise ValueError(f"Could not find [ moleculetype ] in {rec_top}")
    
    end_idx = None
    for i, line in enumerate(lines):
        if i <= start_line:
            continue
        if line.strip().startswith('#include') and '.ff/' in line:
            end_idx = i
            break
    
    if end_idx is None:
        end_idx = len(lines)
    
    receptor_itp = ''.join(lines[start_line:end_idx])
    
    with open(rec_itp, 'w') as f:
        f.write(receptor_itp)
    
    print(f"Created {rec_itp}")


def clean_itp_for_system(itp_file):
    with open(itp_file, 'r') as f:
        content = f.read()
    
    if '[ atomtypes ]' in content:
        atomtypes_pos = content.find('[ atomtypes ]')
        line_start = content.rfind('\n', 0, atomtypes_pos) + 1
        moleculetype_start = content.find('[ moleculetype ]', atomtypes_pos)
        
        if moleculetype_start != -1:
            content = content[:line_start] + content[moleculetype_start:]
            
            with open(itp_file, 'w') as f:
                f.write(content)
            print(f"Removed [ atomtypes ] section from {itp_file}")


def find_ffbonded_file(lig_itp_path, explicit_path=None):
    """
    Find ligand-specific ffbonded.itp file.
    
    Parameters
    ----------
    lig_itp_path : str
        Path to ligand ITP file
    explicit_path : str, optional
        Explicit path provided by user
    
    Returns
    -------
    str or None
        Path to ffbonded.itp file, or None if not found
    """
    if explicit_path:
        if os.path.exists(explicit_path):
            return os.path.abspath(explicit_path)
        else:
            print(f"WARNING: Specified ffbonded.itp not found: {explicit_path}")
            return None
    
    lig_dir = os.path.dirname(os.path.abspath(lig_itp_path))
    
    patterns = [
        os.path.join(lig_dir, "*_ffbonded.itp"),
        os.path.join(lig_dir, "ffbonded.itp"),
    ]
    
    for pattern in patterns:
        matches = glob_module.glob(pattern)
        if matches:
            return matches[0]
    
    return None


def copy_ffbonded_to_workdir(ffbonded_src, output_name="lig_ffbonded.itp"):
    """
    Copy ffbonded.itp to working directory.
    
    Parameters
    ----------
    ffbonded_src : str
        Source path to ffbonded.itp
    output_name : str
        Output filename in working directory
    
    Returns
    -------
    str
        Path to copied ffbonded.itp
    """
    if not os.path.exists(ffbonded_src):
        raise FileNotFoundError(f"ffbonded.itp not found: {ffbonded_src}")
    
    shutil.copy2(ffbonded_src, output_name)
    print(f"Copied {ffbonded_src} -> {output_name}")
    
    return output_name


def generate_ligand_posre(lig_itp_path, output_path="posre_lig.itp", fc=1000):
    """
    Generate position restraints for ligand heavy atoms.
    
    Parses lig.itp, identifies non-hydrogen atoms, writes posre_lig.itp
    with position restraints for each heavy atom.
    
    Parameters
    ----------
    lig_itp_path : str
        Path to ligand ITP file
    output_path : str
        Output path for posre_lig.itp
    fc : float
        Force constant for position restraints (kJ/mol/nm^2)
    
    Returns
    -------
    str
        Path to generated posre_lig.itp
    """
    atoms, _ = parse_itp(lig_itp_path)
    
    heavy_atoms = [a for a in atoms if _itp_element(a["type"]) != "H"]
    
    if not heavy_atoms:
        print("WARNING: No heavy atoms found in ligand ITP")
        return None
    
    with open(output_path, 'w') as f:
        f.write("[ position_restraints ]\n")
        f.write(";  i funct       fcx        fcy        fcz\n")
        for a in heavy_atoms:
            f.write(f"{a['idx']:>4d}    1  {fc:>10.1f}  {fc:>10.1f}  {fc:>10.1f}\n")
    
    print(f"Generated {output_path} with {len(heavy_atoms)} heavy atom restraints")
    
    return output_path


def modify_lig_itp_posre(lig_itp_path, output_path=None):
    """
    Modify lig.itp to reference posre_lig.itp instead of posre.itp.
    
    Changes the #include statement in the POSRES section from
    "posre.itp" to "posre_lig.itp".
    
    Parameters
    ----------
    lig_itp_path : str
        Path to ligand ITP file (will be modified in place if output_path is None)
    output_path : str, optional
        Output path for modified ITP. If None, modifies in place.
    
    Returns
    -------
    str
        Path to modified lig.itp
    """
    if output_path is None:
        output_path = lig_itp_path
    
    with open(lig_itp_path, 'r') as f:
        content = f.read()
    
    content = content.replace('#include "posre.itp"', '#include "posre_lig.itp"')
    
    if lig_itp_path != output_path:
        shutil.copy2(lig_itp_path, output_path)
    
    with open(output_path, 'w') as f:
        f.write(content)
    
    print(f"Modified {output_path}: posre.itp -> posre_lig.itp")
    
    return output_path


def create_system_topology(args):
    for itp in [args.rec_itp, args.lig_itp]:
        if itp and os.path.exists(itp):
            clean_itp_for_system(itp)
    
    sys_top = f"""; Include forcefield parameters
#include "{args.ff_path}"
"""
    
    if args.lig_ffbonded:
        sys_top += f"""
; Include ligand-specific forcefield parameters (CHARMM CGenFF)
#include "{args.lig_ffbonded}"
"""
    
    sys_top += f"""
; Include receptor topology
#include "{args.rec_itp}"

; Include ligand topology
#include "{args.lig_itp}"
"""
    
    if args.water_itp:
        sys_top += f"""

; Include water topology
#include "{args.water_itp}"

#ifdef POSRES_WATER
; Position restraint for each water oxygen
[ position_restraints ]
;  i funct       fcx        fcy        fcz
   1    1       1000       1000       1000
#endif
"""
    
    if args.ions_itp:
        sys_top += f"""

; Include topology for ions
#include "{args.ions_itp}"
"""
    
    sys_top += f"""

[ system ]
; Name
{args.sys_name}

[ molecules ]
; Compound        #mols
{args.rec_name}             1
{args.lig_name}     1
"""
    
    with open(args.sys_top, 'w') as f:
        f.write(sys_top)
    
    print(f"Created {args.sys_top}")


def derive_receptor_gro_from_sdf(sdf_path, pattern=DEFAULT_REC_GRO_PATTERN, search_dir=None):
    basename = os.path.basename(sdf_path)
    name_without_ext = os.path.splitext(basename)[0]
    
    if '-' in name_without_ext:
        prefix = name_without_ext.split('-')[0]
    else:
        prefix = name_without_ext
    
    rec_gro_name = pattern.format(prefix=prefix)
    
    if search_dir is None:
        search_dir = os.path.dirname(sdf_path)
        if not search_dir:
            search_dir = '.'
    
    rec_gro_path = os.path.join(search_dir, rec_gro_name)
    
    return rec_gro_path, prefix


def _no_overwrite(path):
    if os.path.exists(path):
        print(f"ERROR: '{path}' already exists. Refusing to overwrite.", file=sys.stderr)
        sys.exit(1)


def build_parser():
    parser = argparse.ArgumentParser(
        prog="dock2com_2.2.py",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description=textwrap.dedent(__doc__),
    )
    
    lig_group = parser.add_argument_group("Ligand Inputs")
    lig_group.add_argument("-i", "--lig-itp", required=True, metavar="ITP",
                           help="Ligand ITP file (defines atom ordering)")
    lig_group.add_argument("-s", "--sdf", required=True, nargs="+", metavar="SDF",
                           help="One or more SDF docking files. Supports glob patterns.")
    lig_group.add_argument("-t", "--template", default=None, metavar="MOL2",
                           help="Optional mol2 template for hydrogen reconstruction")
    lig_group.add_argument("--lig-ffbonded", default=None, metavar="ITP",
                           help="Ligand-specific ffbonded.itp (CHARMM CGenFF). Auto-detected if not specified.")
    
    rec_group = parser.add_argument_group("Receptor Inputs")
    rec_group.add_argument("-r", "--rec-top", default=None, metavar="TOP",
                           help="Receptor topology file (.top). Required unless --list-models.")
    rec_group.add_argument("--rec-gro", default=None, metavar="GRO",
                           help="Receptor GRO file. If not provided, auto-detected from SDF filename pattern.")
    rec_group.add_argument("--rec-pdb", default=None, metavar="PDB",
                           help="Receptor PDB file. Converts to GRO if --rec-gro not provided.")
    rec_group.add_argument("--rec-gro-pattern", default=DEFAULT_REC_GRO_PATTERN, metavar="PATTERN",
                           help=f"Pattern for receptor GRO filename (used if --rec-gro not provided). Default: {DEFAULT_REC_GRO_PATTERN}")
    rec_group.add_argument("--rec-dir", default=None, metavar="DIR",
                           help="Directory containing receptor files. Default: same as SDF files.")
    
    out_group = parser.add_argument_group("Outputs")
    out_group.add_argument("--lig-gro", default="best.gro", metavar="GRO",
                           help="Output ligand GRO file. Default: best.gro")
    out_group.add_argument("--com-gro", default="com.gro", metavar="GRO",
                           help="Output combined GRO file. Default: com.gro")
    out_group.add_argument("--rec-itp", default="rec.itp", metavar="ITP",
                           help="Output receptor ITP file. Default: rec.itp")
    out_group.add_argument("--sys-top", default="sys.top", metavar="TOP",
                           help="Output system topology file. Default: sys.top")
    
    sel_group = parser.add_argument_group("Model Selection")
    sel_group.add_argument("--metric", default=DEFAULT_METRIC, metavar="NAME",
                           help=f"Score field for model selection. Options: minimizedAffinity, CNNscore, CNNaffinity, sasa (lower SASA = better burial). Default: {DEFAULT_METRIC}")
    sel_group.add_argument("--model", type=int, default=None, metavar="N",
                           help="Select specific model number (1-based).")
    sel_group.add_argument("--list-models", action="store_true",
                           help="List all models with scores and exit.")
    sel_group.add_argument("--probe-radius", type=float, default=DEFAULT_PROBE_RADIUS, metavar="R",
                           help=f"Solvent probe radius in Angstroms for SASA calculation. Default: {DEFAULT_PROBE_RADIUS}")
    
    top_group = parser.add_argument_group("Topology Options")
    top_group.add_argument("--ff-path", default=None, metavar="PATH",
                           help="Forcefield.itp path. Default: auto-detect from receptor topology.")
    top_group.add_argument("--water-itp", default=None, metavar="PATH",
                           help="Water ITP path. Default: auto-detect from receptor topology.")
    top_group.add_argument("--ions-itp", default=None, metavar="PATH",
                           help="Ions ITP path. Default: auto-detect from receptor topology.")
    top_group.add_argument("--sys-name", default=None, metavar="NAME",
                           help="System name. Default: auto-generate from receptor and ligand names.")
    top_group.add_argument("--rec-name", default=None, metavar="NAME",
                           help="Receptor molecule name. Default: auto-detect from receptor topology.")
    top_group.add_argument("--lig-name", default=None, metavar="NAME",
                           help="Ligand molecule name. Default: auto-detect from ligand ITP.")
    
    return parser


class Config:
    def __init__(self, args):
        self.lig_itp = args.lig_itp
        self.sdf_paths = collect_sdf_files(args.sdf)
        self.mol2_template = args.template
        self.lig_ffbonded_arg = args.lig_ffbonded
        
        self.rec_top = args.rec_top
        self.rec_gro = args.rec_gro
        self.rec_pdb = args.rec_pdb
        self.rec_gro_pattern = args.rec_gro_pattern
        self.rec_dir = args.rec_dir
        
        self.lig_gro = args.lig_gro
        self.com_gro = args.com_gro
        self.rec_itp = args.rec_itp
        self.sys_top = args.sys_top
        
        self.metric = args.metric
        self.model_num = args.model
        self.list_models = args.list_models
        self.probe_radius = args.probe_radius
        
        self.ff_path = args.ff_path
        self.water_itp = args.water_itp
        self.ions_itp = args.ions_itp
        self.sys_name = args.sys_name
        self.rec_name = args.rec_name
        self.lig_name = args.lig_name
        
        self.lig_ffbonded = None


def main():
    parser = build_parser()
    args = parser.parse_args()
    config = Config(args)
    
    if config.list_models:
        list_metrics = None
        if config.metric == "sasa":
            all_keys = []
            seen = set()
            for sdf_path in config.sdf_paths:
                models = parse_sdf(sdf_path)
                for m in models:
                    for k in m["scores"]:
                        if k not in seen:
                            all_keys.append(k)
                            seen.add(k)
            all_keys.append("sasa")
            list_metrics = all_keys
        list_models_across_sdfiles(
            config.sdf_paths,
            metrics=list_metrics,
            rec_gro_pattern=config.rec_gro_pattern,
            rec_dir=config.rec_dir,
            probe_radius=config.probe_radius
        )
        return 0
    
    if config.rec_top is None:
        parser.error("-r/--rec-top is required unless --list-models is specified.")
    
    print("=" * 60)
    print("STEP 1: Selecting best docking model...")
    print("=" * 60)
    
    best_model, best_file, best_model_num = select_best_across_sdfiles(
        config.sdf_paths,
        metric=config.metric,
        rec_gro_pattern=config.rec_gro_pattern,
        rec_dir=config.rec_dir,
        probe_radius=config.probe_radius,
    )
    
    sc = best_model["scores"]
    print(f"\nBest model selected:")
    print(f"  File:   {os.path.basename(best_file)}")
    print(f"  Model:  {best_model_num}")
    print(f"  Scores: minimizedAffinity={sc.get('minimizedAffinity', float('nan')):.4f}  "
          f"CNNscore={sc.get('CNNscore', float('nan')):.4f}  "
          f"CNNaffinity={sc.get('CNNaffinity', float('nan')):.4f}")
    
    if 'sasa' in sc:
        print(f"          SASA={sc['sasa']:.4f} nm²")
    
    print("\n" + "=" * 60)
    print("STEP 2: Locating receptor GRO file...")
    print("=" * 60)
    
    rec_gro = None
    prefix = None
    
    if config.rec_gro:
        rec_gro = config.rec_gro
        prefix = os.path.splitext(os.path.basename(rec_gro))[0]
        if not os.path.exists(rec_gro):
            print(f"ERROR: Receptor GRO file not found: {rec_gro}", file=sys.stderr)
            sys.exit(1)
        print(f"  Using specified GRO: {rec_gro}")
    elif config.rec_pdb:
        if not os.path.exists(config.rec_pdb):
            print(f"ERROR: Receptor PDB file not found: {config.rec_pdb}", file=sys.stderr)
            sys.exit(1)
        prefix = os.path.splitext(os.path.basename(config.rec_pdb))[0]
        rec_gro = f"{prefix}.gro"
        print(f"  Converting PDB to GRO: {config.rec_pdb} -> {rec_gro}")
        n_atoms = pdb_to_gro(config.rec_pdb, rec_gro)
        print(f"  Converted {n_atoms} atoms")
    else:
        rec_gro, prefix = derive_receptor_gro_from_sdf(
            best_file,
            pattern=config.rec_gro_pattern,
            search_dir=config.rec_dir
        )
        
        if not os.path.exists(rec_gro):
            print(f"ERROR: Receptor GRO file not found: {rec_gro}", file=sys.stderr)
            print(f"  Derived from SDF: {best_file}", file=sys.stderr)
            print(f"  Pattern used: {config.rec_gro_pattern}", file=sys.stderr)
            sys.exit(1)
        
        print(f"  Auto-detected: {rec_gro}")
    
    print(f"  Prefix: {prefix}")
    
    print("\n" + "=" * 60)
    print("STEP 3: Converting SDF to ligand GRO...")
    print("=" * 60)
    
    sdf_pose_to_gro(
        itp_path=config.lig_itp,
        sdf_path=best_file,
        out_path=config.lig_gro,
        mol2_template_path=config.mol2_template,
        metric=config.metric,
        model_num=best_model_num
    )
    
    print("\n" + "=" * 60)
    print("STEP 4: Combining receptor and ligand coordinates...")
    print("=" * 60)
    
    combine_coordinates(rec_gro, config.lig_gro, config.com_gro)
    
    print("\n" + "=" * 60)
    print("STEP 5: Extracting receptor topology...")
    print("=" * 60)
    
    extract_receptor_topology(config.rec_top, config.rec_itp)
    
    print("\n" + "=" * 60)
    print("STEP 6: Handling ligand ffbonded.itp (CHARMM CGenFF)...")
    print("=" * 60)
    
    ffbonded_src = find_ffbonded_file(config.lig_itp, config.lig_ffbonded_arg)
    
    if ffbonded_src:
        print(f"  Found ffbonded.itp: {ffbonded_src}")
        config.lig_ffbonded = copy_ffbonded_to_workdir(ffbonded_src, "lig_ffbonded.itp")
    else:
        print("  No ffbonded.itp found (may be needed for CHARMM CGenFF)")
        config.lig_ffbonded = None
    
    print("\n" + "=" * 60)
    print("STEP 7: Generating ligand position restraints...")
    print("=" * 60)
    
    generate_ligand_posre(config.lig_itp, "posre_lig.itp", fc=1000)
    
    modify_lig_itp_posre(config.lig_itp)
    
    print("\n" + "=" * 60)
    print("STEP 8: Creating system topology...")
    print("=" * 60)
    
    ff_path, water_itp, ions_itp = extract_ff_paths_from_top(config.rec_top)
    
    if config.ff_path is None:
        config.ff_path = ff_path
        print(f"  Auto-detected forcefield: {config.ff_path}")
    
    if config.water_itp is None:
        config.water_itp = water_itp if water_itp else ''
        if config.water_itp:
            print(f"  Auto-detected water model: {config.water_itp}")
    
    if config.ions_itp is None:
        config.ions_itp = ions_itp if ions_itp else ''
        if config.ions_itp:
            print(f"  Auto-detected ions ITP: {config.ions_itp}")
    
    if config.rec_name is None:
        config.rec_name = get_moleculetype_name(config.rec_top)
        print(f"  Auto-detected receptor name: {config.rec_name}")
    
    if config.lig_name is None:
        config.lig_name = get_moleculetype_name(config.lig_itp)
        print(f"  Auto-detected ligand name: {config.lig_name}")
    
    if config.sys_name is None:
        config.sys_name = f"{prefix}.pdb_ali_best"
        print(f"  Auto-generated system name: {config.sys_name}")
    
    create_system_topology(config)
    
    print("\n" + "=" * 60)
    print("SUMMARY: Files created")
    print("=" * 60)
    print(f"  Ligand GRO:    {config.lig_gro}")
    print(f"  Combined GRO:  {config.com_gro}")
    print(f"  Receptor ITP:  {config.rec_itp}")
    print(f"  System TOP:    {config.sys_top}")
    print(f"  Ligand posre:  posre_lig.itp")
    if config.lig_ffbonded:
        print(f"  Ligand FF:     {config.lig_ffbonded}")
    print("\nDock2com_2.2 completed successfully!")


if __name__ == "__main__":
    main()
