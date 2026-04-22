#!/usr/bin/env python3
"""
find_similar_poses.py
=====================
Find docking poses similar to a target binding mode.

Workflow:
1. Extract receptor (chain A) and ligand (chain B / MOL) from target.pdb
2. Build a molecular graph for the target ligand using distance-based bonds
3. For each hsa receptor, align it to the target receptor using CA atoms (Kabsch)
4. Apply the same rotation/translation to every docking pose from the SDF files
5. For same-molecule ligands: use graph isomorphism for proper atom mapping, then RMSD
   For different-molecule ligands: use Hungarian element matching for RMSD
6. Output a ranked CSV of all poses sorted by RMSD

Usage:
    python find_similar_poses.py
"""

import csv
import os
import sys
import glob
from collections import defaultdict

import numpy as np
import networkx as nx
from networkx.algorithms import isomorphism

# Typical bond length cutoffs (in Angstroms) for heavy atoms
BOND_CUTOFFS = {
    ("C", "C"): 1.85, ("C", "N"): 1.75, ("C", "O"): 1.70,
    ("N", "N"): 1.70, ("N", "O"): 1.70, ("O", "O"): 1.70,
    ("C", "S"): 2.00, ("N", "S"): 2.00, ("S", "S"): 2.20,
}
DEFAULT_BOND_CUTOFF = 1.85
METAL_BOND_CUTOFF = 2.8  # For EU-X bonds (coordination bonds)

METAL_ELEMENTS = {"EU", "FE", "ZN", "MG", "CA", "MN", "CO", "CU", "NI", "PT", "AU", "AG"}


def _norm_element(el):
    """Normalize element string."""
    e = el.upper().strip()
    if e.startswith("EU"):
        return "EU"
    return e


def _bond_cutoff(el1, el2):
    """Get bond length cutoff for a pair of elements."""
    e1, e2 = _norm_element(el1), _norm_element(el2)
    if e1 in METAL_ELEMENTS or e2 in METAL_ELEMENTS:
        return METAL_BOND_CUTOFF
    pair = tuple(sorted([e1, e2]))
    return BOND_CUTOFFS.get(pair, DEFAULT_BOND_CUTOFF)


# ---------------------------------------------------------------------------
# PDB / SDF parsing
# ---------------------------------------------------------------------------

def parse_pdb_atoms(path):
    """Parse ATOM/HETATM records from a PDB file."""
    atoms = []
    with open(path) as f:
        for line in f:
            if line.startswith("ATOM") or line.startswith("HETATM"):
                try:
                    atom_name = line[12:16].strip()
                    res_name = line[17:20].strip()
                    chain = line[21].strip()
                    res_num = int(line[22:26].strip())
                    x = float(line[30:38])
                    y = float(line[38:46])
                    z = float(line[46:54])
                    element = line[76:78].strip() if len(line) > 76 else ""
                    # Infer element from atom name if missing
                    if not element:
                        name_upper = atom_name.upper()
                        if name_upper.startswith("EU"):
                            element = "EU"
                        elif name_upper.startswith("CL"):
                            element = "CL"
                        elif name_upper.startswith("BR"):
                            element = "BR"
                        elif name_upper.startswith("FE"):
                            element = "FE"
                        elif name_upper.startswith("ZN"):
                            element = "ZN"
                        elif name_upper.startswith("MG"):
                            element = "MG"
                        else:
                            element = name_upper[0] if name_upper else "C"
                    atoms.append({
                        "atom_name": atom_name,
                        "res_name": res_name,
                        "chain": chain,
                        "res_num": res_num,
                        "x": x, "y": y, "z": z,
                        "element": element,
                    })
                except (ValueError, IndexError):
                    continue
    return atoms


def extract_ca_coords(atoms, chain=None):
    """Extract CA atom coordinates, optionally filtered by chain."""
    ca_atoms = []
    for a in atoms:
        if a["atom_name"] == "CA":
            if chain is None or a["chain"] == chain:
                ca_atoms.append(a)
    ca_atoms.sort(key=lambda a: a["res_num"])
    coords = np.array([[a["x"], a["y"], a["z"]] for a in ca_atoms])
    res_nums = [a["res_num"] for a in ca_atoms]
    return coords, res_nums


def extract_ligand_heavy_atoms(atoms, chain="B", res_name="MOL"):
    """Extract heavy (non-H) ligand atoms from target PDB."""
    lig = []
    for a in atoms:
        if a["chain"] == chain and a["res_name"] == res_name:
            el = _norm_element(a["element"])
            if el not in ("H", ""):
                lig.append(a)
    coords = np.array([[a["x"], a["y"], a["z"]] for a in lig])
    elements = [_norm_element(a["element"]) for a in lig]
    return coords, elements


def parse_sdf_models(path):
    """Parse all models from an SDF file. Returns list of dicts with atoms, bonds, scores."""
    with open(path) as f:
        content = f.read()

    models = []
    records = content.split("$$$$")
    for rec_idx, rec in enumerate(records):
        rec = rec.strip()
        if not rec:
            continue
        lines = rec.split("\n")
        if len(lines) < 4:
            continue

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
            try:
                el = ln[31:34].strip()
                if el == "*":
                    el = "EU"
                atoms.append({
                    "idx": i - 3,  # 1-based index
                    "el": _norm_element(el),
                    "x": float(ln[0:10]),
                    "y": float(ln[10:20]),
                    "z": float(ln[20:30]),
                })
            except (ValueError, IndexError):
                continue

        bonds = []
        for i in range(4 + ac, 4 + ac + bc):
            if i >= len(lines):
                break
            ln = lines[i]
            try:
                bonds.append({
                    "a1": int(ln[0:3]),
                    "a2": int(ln[3:6]),
                })
            except (ValueError, IndexError):
                continue

        # Parse scores
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
            "atoms": atoms,
            "bonds": bonds,
            "scores": scores,
        })

    return models


# ---------------------------------------------------------------------------
# Graph construction and isomorphism
# ---------------------------------------------------------------------------

def build_graph_from_coords(coords, elements):
    """Build a molecular graph from coordinates using distance-based bond detection.
    Uses scipy cdist for vectorized pairwise distance computation."""
    from scipy.spatial.distance import cdist

    G = nx.Graph()
    for i, el in enumerate(elements):
        G.add_node(i, element=el)

    n = len(elements)
    if n == 0:
        return G

    # Compute full pairwise distance matrix at once
    dist_matrix = cdist(coords, coords)

    # Build cutoff matrix
    normed = [_norm_element(e) for e in elements]
    cutoff_matrix = np.full((n, n), DEFAULT_BOND_CUTOFF)
    for i in range(n):
        for j in range(i + 1, n):
            c = _bond_cutoff(normed[i], normed[j])
            cutoff_matrix[i, j] = c
            cutoff_matrix[j, i] = c

    # Find all bonded pairs (upper triangle only)
    mask = (dist_matrix < cutoff_matrix) & (dist_matrix > 0)
    i_idx, j_idx = np.where(np.triu(mask, k=1))
    G.add_edges_from(zip(i_idx.tolist(), j_idx.tolist()))

    return G


def build_graph_from_sdf(atoms, bonds, heavy_only=True):
    """Build a molecular graph from SDF atoms and bonds."""
    G = nx.Graph()

    if heavy_only:
        heavy_idx = {a["idx"] for a in atoms if a["el"] != "H"}
    else:
        heavy_idx = {a["idx"] for a in atoms}

    idx_to_node = {}
    for a in atoms:
        if a["idx"] in heavy_idx:
            node_id = len(idx_to_node)
            idx_to_node[a["idx"]] = node_id
            G.add_node(node_id, element=a["el"])

    for b in bonds:
        if b["a1"] in heavy_idx and b["a2"] in heavy_idx:
            n1 = idx_to_node[b["a1"]]
            n2 = idx_to_node[b["a2"]]
            G.add_edge(n1, n2)

    return G, idx_to_node


def _node_match(n1, n2):
    """Node match function for graph isomorphism."""
    e1 = n1["element"]
    e2 = n2["element"]
    if e1 == e2:
        return True
    # Treat all metals as equivalent
    if e1 in METAL_ELEMENTS and e2 in METAL_ELEMENTS:
        return True
    return False


def find_graph_mapping(G_target, G_sdf, target_coords=None, sdf_coords=None):
    """
    Find graph isomorphism mapping from target graph nodes to SDF graph nodes.
    If coordinates are provided, selects the mapping with minimum RMSD.
    Returns dict: target_node -> sdf_node, or None if no isomorphism found.
    """
    gm = isomorphism.GraphMatcher(G_target, G_sdf, node_match=_node_match)
    if not gm.is_isomorphic():
        return None

    if target_coords is None or sdf_coords is None:
        return next(gm.isomorphisms_iter())

    best_mapping = None
    best_rmsd = float("inf")
    for mapping in gm.isomorphisms_iter():
        sq = [np.sum((target_coords[ti] - sdf_coords[di]) ** 2)
              for ti, di in mapping.items()]
        rmsd = np.sqrt(np.mean(sq))
        if rmsd < best_rmsd:
            best_rmsd = rmsd
            best_mapping = mapping

    return best_mapping


# ---------------------------------------------------------------------------
# Alignment (Kabsch algorithm)
# ---------------------------------------------------------------------------

def kabsch_align(P, Q):
    """
    Compute rotation R and translation t to align P onto Q.
    Returns R, t such that Q ~ R @ P + t.
    """
    centroid_P = P.mean(axis=0)
    centroid_Q = Q.mean(axis=0)

    P_c = P - centroid_P
    Q_c = Q - centroid_Q

    H = P_c.T @ Q_c
    U, S, Vt = np.linalg.svd(H)

    d = np.linalg.det(Vt.T @ U.T)
    sign_matrix = np.diag([1, 1, np.sign(d)])

    R = Vt.T @ sign_matrix @ U.T
    t = centroid_Q - R @ centroid_P

    return R, t


def apply_transform(coords, R, t):
    """Apply rotation R and translation t to coordinates."""
    return (R @ coords.T).T + t


# ---------------------------------------------------------------------------
# RMSD computation
# ---------------------------------------------------------------------------

def compute_rmsd_with_mapping(target_coords, docked_coords, mapping):
    """Compute RMSD using a known atom mapping (from graph isomorphism).
    Vectorized version."""
    t_idx = np.array(list(mapping.keys()))
    d_idx = np.array(list(mapping.values()))
    diff = target_coords[t_idx] - docked_coords[d_idx]
    return np.sqrt(np.mean(np.sum(diff ** 2, axis=1))), len(mapping)


def compute_rmsd_with_mapping_multi(target_coords, docked_coords, mappings):
    """Compute RMSD for multiple mappings at once, return minimum.
    mappings: list of dicts (target_node -> sdf_node).
    Returns (best_rmsd, n_matched)."""
    if not mappings:
        return float("inf"), 0
    # Pre-build index arrays for all mappings (they all have same keys for isomorphisms)
    t_idx = np.array(list(mappings[0].keys()))
    n_atoms = len(t_idx)
    n_maps = len(mappings)

    # Build d_idx array: shape (n_maps, n_atoms)
    d_idx_all = np.array([[m[k] for k in t_idx] for m in mappings])

    # target_coords[t_idx]: shape (n_atoms, 3)
    t_coords = target_coords[t_idx]  # (n_atoms, 3)

    # docked_coords[d_idx_all]: shape (n_maps, n_atoms, 3)
    d_coords_all = docked_coords[d_idx_all]  # (n_maps, n_atoms, 3)

    # Compute RMSD for each mapping
    diff = t_coords[None, :, :] - d_coords_all  # (n_maps, n_atoms, 3)
    sq_dists = np.sum(diff ** 2, axis=2)  # (n_maps, n_atoms)
    rmsds = np.sqrt(np.mean(sq_dists, axis=1))  # (n_maps,)

    best_idx = np.argmin(rmsds)
    return float(rmsds[best_idx]), n_atoms


UNMATCHED_PENALTY = 10.0  # Angstroms - penalty distance for each unmatched target atom


def compute_hungarian_rmsd(target_coords, target_elements, docked_coords, docked_elements):
    """
    Compute RMSD using Hungarian (optimal) element-based matching.
    Includes a penalty for unmatched target atoms to avoid bias when
    the docked ligand has fewer atoms than the target.
    
    Each unmatched target atom contributes UNMATCHED_PENALTY^2 to the
    sum of squared distances, ensuring smaller molecules don't get
    artificially lower RMSD values.
    """
    from scipy.optimize import linear_sum_assignment

    if len(target_coords) == 0 or len(docked_coords) == 0:
        return float("inf"), 0

    t_elems = [_norm_element(e) for e in target_elements]
    d_elems = [_norm_element(e) for e in docked_elements]

    target_by_elem = defaultdict(list)
    for i, e in enumerate(t_elems):
        target_by_elem[e].append(i)

    docked_by_elem = defaultdict(list)
    for i, e in enumerate(d_elems):
        docked_by_elem[e].append(i)

    matched_pairs = []
    n_unmatched = 0
    for elem in target_by_elem:
        if elem not in docked_by_elem:
            n_unmatched += len(target_by_elem[elem])
            continue
        t_indices = list(target_by_elem[elem])
        d_indices = list(docked_by_elem[elem])

        t_c = target_coords[t_indices]
        d_c = docked_coords[d_indices]

        diff = t_c[:, None, :] - d_c[None, :, :]
        dist_matrix = np.sqrt(np.sum(diff ** 2, axis=2))

        row_ind, col_ind = linear_sum_assignment(dist_matrix)
        for r, c in zip(row_ind, col_ind):
            matched_pairs.append((t_indices[r], d_indices[c]))

        # Count unmatched target atoms of this element
        n_unmatched += max(0, len(t_indices) - len(d_indices))

    if not matched_pairs:
        return float("inf"), 0

    matched_t = np.array([ti for ti, di in matched_pairs])
    matched_d = np.array([di for ti, di in matched_pairs])
    sq_dists = list(np.sum((target_coords[matched_t] - docked_coords[matched_d]) ** 2, axis=1))
    # Add penalty for unmatched target atoms
    sq_dists.extend([UNMATCHED_PENALTY ** 2] * n_unmatched)

    n_total = len(matched_pairs) + n_unmatched
    return np.sqrt(np.sum(sq_dists) / n_total), len(matched_pairs)


# ---------------------------------------------------------------------------
# Contact distance profile
# ---------------------------------------------------------------------------

CONTACT_CUTOFF = 8.0  # Angstroms - max distance to consider a residue as part of binding pocket


def extract_receptor_heavy_atoms(atoms, chain=None):
    """Extract all heavy (non-H) receptor atoms. If chain is None, use all atoms."""
    heavy = []
    for a in atoms:
        if chain is not None and a["chain"] != chain:
            continue
        el = _norm_element(a["element"])
        if el not in ("H", ""):
            heavy.append(a)
    if not heavy:
        return np.zeros((0, 3)), np.array([], dtype=int)
    coords = np.array([[a["x"], a["y"], a["z"]] for a in heavy])
    res_nums = np.array([a["res_num"] for a in heavy])
    return coords, res_nums


def compute_contact_profile(rec_coords, rec_res_nums, lig_coords, cutoff=CONTACT_CUTOFF):
    """
    Compute contact distance profile: for each residue, the minimum distance
    from any of its heavy atoms to any ligand heavy atom.
    
    Uses KDTree for efficient nearest-neighbor computation.
    Returns dict: {res_num: min_distance} for residues within cutoff.
    """
    from scipy.spatial import cKDTree

    if len(rec_coords) == 0 or len(lig_coords) == 0:
        return {}

    lig_tree = cKDTree(lig_coords)

    # For each receptor atom, find distance to nearest ligand atom
    dists, _ = lig_tree.query(rec_coords, k=1)

    # Group by residue, take minimum distance per residue
    unique_res = np.unique(rec_res_nums)
    profile = {}
    for res in unique_res:
        mask = rec_res_nums == res
        min_dist = np.min(dists[mask])
        if min_dist < cutoff:
            profile[res] = float(min_dist)

    return profile


def contact_similarity(profile1, profile2):
    """
    Compute similarity between two contact distance profiles.
    Uses Pearson correlation on the union of residues (missing residues get cutoff value).
    Returns correlation coefficient (1.0 = identical, 0.0 = uncorrelated, -1.0 = anticorrelated).
    """
    all_res = sorted(set(profile1.keys()) | set(profile2.keys()))
    if len(all_res) < 2:
        return 0.0

    v1 = np.array([profile1.get(r, CONTACT_CUTOFF) for r in all_res])
    v2 = np.array([profile2.get(r, CONTACT_CUTOFF) for r in all_res])

    # Pearson correlation
    v1c = v1 - v1.mean()
    v2c = v2 - v2.mean()
    denom = np.sqrt(np.sum(v1c ** 2) * np.sum(v2c ** 2))
    if denom < 1e-12:
        return 1.0 if np.allclose(v1, v2) else 0.0
    return float(np.sum(v1c * v2c) / denom)


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    script_dir = os.path.dirname(os.path.abspath(__file__))

    # 1. Parse target.pdb
    target_pdb = os.path.join(script_dir, "target.pdb")
    print(f"Parsing target: {target_pdb}")
    target_atoms = parse_pdb_atoms(target_pdb)

    target_ca_coords, target_ca_res = extract_ca_coords(target_atoms, chain="A")
    print(f"  Target receptor CA atoms: {len(target_ca_coords)}")

    target_lig_coords, target_lig_elems = extract_ligand_heavy_atoms(target_atoms, chain="B", res_name="MOL")
    print(f"  Target ligand heavy atoms: {len(target_lig_coords)}")
    elem_counts = defaultdict(int)
    for e in target_lig_elems:
        elem_counts[e] += 1
    print(f"  Target ligand elements: {dict(elem_counts)}")

    # Extract target receptor heavy atoms for contact profiles
    target_rec_coords, target_rec_res = extract_receptor_heavy_atoms(target_atoms, chain="A")
    print(f"  Target receptor heavy atoms: {len(target_rec_coords)}")

    # Compute target contact profile
    target_contact_profile = compute_contact_profile(target_rec_coords, target_rec_res, target_lig_coords)
    print(f"  Target contact profile: {len(target_contact_profile)} residues within {CONTACT_CUTOFF} A")

    # 2. Build molecular graph for target ligand (distance-based bonds)
    target_graph = build_graph_from_coords(target_lig_coords, target_lig_elems)
    print(f"  Target ligand graph: {target_graph.number_of_nodes()} nodes, {target_graph.number_of_edges()} edges")

    # 3. Load hsa receptor PDBs and compute alignments + receptor heavy atoms
    hsa_pdbs = sorted(glob.glob(os.path.join(script_dir, "hsa*.pdb_ali.pdb")))
    print(f"\nFound {len(hsa_pdbs)} receptor PDB files")

    transforms = {}
    hsa_rec_heavy = {}  # {prefix: (aligned_coords, res_nums)} for contact profiles

    for hsa_pdb in hsa_pdbs:
        basename = os.path.basename(hsa_pdb)
        prefix = basename.split(".")[0]

        hsa_atoms = parse_pdb_atoms(hsa_pdb)
        hsa_ca_coords, hsa_ca_res = extract_ca_coords(hsa_atoms)

        target_res_set = set(target_ca_res)
        hsa_res_set = set(hsa_ca_res)
        common_res = sorted(target_res_set & hsa_res_set)

        if len(common_res) < 10:
            print(f"  WARNING: {prefix} has only {len(common_res)} common CA atoms, skipping")
            continue

        t_res_to_idx = {r: i for i, r in enumerate(target_ca_res)}
        h_res_to_idx = {r: i for i, r in enumerate(hsa_ca_res)}

        P = np.array([hsa_ca_coords[h_res_to_idx[r]] for r in common_res])
        Q = np.array([target_ca_coords[t_res_to_idx[r]] for r in common_res])

        R, t = kabsch_align(P, Q)
        transforms[prefix] = (R, t)

        aligned = apply_transform(P, R, t)
        rmsd_align = np.sqrt(np.mean(np.sum((aligned - Q) ** 2, axis=1)))
        print(f"  {prefix}: {len(common_res)} common CA, alignment RMSD = {rmsd_align:.3f} A")

        # Transform receptor heavy atoms to target frame for contact profiles
        rec_coords, rec_res = extract_receptor_heavy_atoms(hsa_atoms)  # no chain filter for hsa PDBs
        rec_aligned = apply_transform(rec_coords, R, t)
        hsa_rec_heavy[prefix] = (rec_aligned, rec_res)

    # 4. Extract starting pose reference: phe_sssD hsa5 model 12
    starting_pose_coords = None
    starting_pose_profile = None
    starting_pose_sdf = os.path.join(script_dir, "phe_sssD", "hsa5-phe_sssD.sdf")
    if os.path.exists(starting_pose_sdf) and "hsa5" in transforms:
        print(f"\nExtracting starting pose reference: phe_sssD hsa5 model 12")
        sp_models = parse_sdf_models(starting_pose_sdf)
        sp_model = None
        for m in sp_models:
            if m["model_num"] == 12:
                sp_model = m
                break
        if sp_model:
            sp_graph, sp_idx_to_node = build_graph_from_sdf(
                sp_model["atoms"], sp_model["bonds"], heavy_only=True
            )
            sp_node_to_orig = {v: k for k, v in sp_idx_to_node.items()}
            atoms_by_idx = {a["idx"]: a for a in sp_model["atoms"]}
            n_sp = sp_graph.number_of_nodes()
            sp_coords = np.zeros((n_sp, 3))
            for node_id in range(n_sp):
                a = atoms_by_idx[sp_node_to_orig[node_id]]
                sp_coords[node_id] = [a["x"], a["y"], a["z"]]

            R_sp, t_sp = transforms["hsa5"]
            starting_pose_coords = apply_transform(sp_coords, R_sp, t_sp)

            # Compute starting pose contact profile using aligned hsa5 receptor
            if "hsa5" in hsa_rec_heavy:
                sp_rec_coords, sp_rec_res = hsa_rec_heavy["hsa5"]
                starting_pose_profile = compute_contact_profile(sp_rec_coords, sp_rec_res, starting_pose_coords)
                print(f"  Starting pose contact profile: {len(starting_pose_profile)} residues")

                # Verify: contact similarity between starting pose and target
                sim = contact_similarity(target_contact_profile, starting_pose_profile)
                print(f"  Starting pose vs target contact similarity: {sim:.4f}")
    else:
        print(f"\nWARNING: Could not find starting pose SDF or hsa5 transform")

    # 5. Find all docking directories
    all_subdirs = sorted([
        d for d in os.listdir(script_dir)
        if os.path.isdir(os.path.join(script_dir, d))
        and not d.startswith("__")
    ])
    print(f"\nDocking directories: {all_subdirs}")

    # 6. Process each directory and SDF file
    results = []

    for subdir in all_subdirs:
        subdir_path = os.path.join(script_dir, subdir)
        sdf_files = sorted(glob.glob(os.path.join(subdir_path, "*.sdf")))

        if not sdf_files:
            continue

        print(f"\nProcessing {subdir}: {len(sdf_files)} SDF files")

        # Parse first model of first SDF to determine molecule type and cache graph/isomorphisms
        first_models = parse_sdf_models(sdf_files[0])
        if not first_models:
            continue

        sdf_graph_cached, sdf_idx_to_node_cached = build_graph_from_sdf(
            first_models[0]["atoms"], first_models[0]["bonds"], heavy_only=True
        )

        # Check graph isomorphism once for this directory
        mapping = find_graph_mapping(target_graph, sdf_graph_cached)
        use_graph_iso = mapping is not None

        # Pre-enumerate ALL isomorphisms once (for same-molecule case)
        all_isomorphisms = None
        if use_graph_iso:
            gm = isomorphism.GraphMatcher(target_graph, sdf_graph_cached, node_match=_node_match)
            all_isomorphisms = list(gm.isomorphisms_iter())
            print(f"  Graph isomorphism: YES (same molecule, {len(mapping)} atoms, {len(all_isomorphisms)} symmetry mappings)")
        else:
            print(f"  Graph isomorphism: NO (different molecule, using Hungarian matching)")

        # Check if starting pose graph iso is available for this directory
        sp_iso = None
        if starting_pose_coords is not None and use_graph_iso:
            # Same molecule as target = same molecule as starting pose
            sp_iso = True

        for sdf_file in sdf_files:
            sdf_basename = os.path.basename(sdf_file)

            if "-" in sdf_basename:
                hsa_prefix = sdf_basename.split("-")[0]
            else:
                hsa_prefix = sdf_basename.split(".")[0]

            if hsa_prefix not in transforms:
                print(f"  WARNING: No alignment for {hsa_prefix}, skipping {sdf_basename}")
                continue

            R, t = transforms[hsa_prefix]
            models = parse_sdf_models(sdf_file)

            # Build SDF graph once per SDF file (all models have same topology)
            if models:
                sdf_graph_file, sdf_idx_to_node_file = build_graph_from_sdf(
                    models[0]["atoms"], models[0]["bonds"], heavy_only=True
                )
                sdf_node_to_orig_file = {v: k for k, v in sdf_idx_to_node_file.items()}
                n_sdf_nodes = sdf_graph_file.number_of_nodes()
            else:
                continue

            # Get aligned receptor heavy atoms for this hsa
            rec_coords_aligned = None
            rec_res_aligned = None
            if hsa_prefix in hsa_rec_heavy:
                rec_coords_aligned, rec_res_aligned = hsa_rec_heavy[hsa_prefix]

            for model in models:
                heavy_atoms = [a for a in model["atoms"] if a["el"] != "H"]
                if not heavy_atoms:
                    continue

                # Reuse cached graph structure — just extract coords for this model
                atoms_by_idx = {a["idx"]: a for a in model["atoms"]}

                sdf_coords = np.zeros((n_sdf_nodes, 3))
                sdf_elems = []
                for node_id in range(n_sdf_nodes):
                    orig_idx = sdf_node_to_orig_file[node_id]
                    a = atoms_by_idx[orig_idx]
                    sdf_coords[node_id] = [a["x"], a["y"], a["z"]]
                    sdf_elems.append(a["el"])

                # Transform to target frame
                sdf_transformed = apply_transform(sdf_coords, R, t)

                # --- RMSD vs MD target ---
                if use_graph_iso and all_isomorphisms:
                    rmsd, n_matched = compute_rmsd_with_mapping_multi(
                        target_lig_coords, sdf_transformed, all_isomorphisms
                    )
                else:
                    rmsd, n_matched = compute_hungarian_rmsd(
                        target_lig_coords, target_lig_elems,
                        sdf_transformed, sdf_elems
                    )

                # --- RMSD vs starting pose ---
                rmsd_sp = float("nan")
                if starting_pose_coords is not None:
                    if sp_iso and all_isomorphisms:
                        rmsd_sp, _ = compute_rmsd_with_mapping_multi(
                            starting_pose_coords, sdf_transformed, all_isomorphisms
                        )
                    else:
                        # Starting pose is phe (80 atoms); for me_* use Hungarian
                        sp_elems = target_lig_elems  # same molecule as target
                        rmsd_sp, _ = compute_hungarian_rmsd(
                            starting_pose_coords, sp_elems,
                            sdf_transformed, sdf_elems
                        )

                # --- Contact distance profile ---
                contact_sim = float("nan")
                contact_sim_sp = float("nan")
                if rec_coords_aligned is not None:
                    pose_profile = compute_contact_profile(
                        rec_coords_aligned, rec_res_aligned, sdf_transformed
                    )
                    contact_sim = contact_similarity(target_contact_profile, pose_profile)
                    if starting_pose_profile is not None:
                        contact_sim_sp = contact_similarity(starting_pose_profile, pose_profile)

                affinity = model["scores"].get("minimizedAffinity", float("nan"))
                cnn_score = model["scores"].get("CNNscore", float("nan"))
                cnn_affinity = model["scores"].get("CNNaffinity", float("nan"))

                results.append({
                    "directory": subdir,
                    "sdf_file": sdf_basename,
                    "hsa_prefix": hsa_prefix,
                    "model_num": model["model_num"],
                    "rmsd": rmsd,
                    "rmsd_sp": rmsd_sp,
                    "contact_sim": contact_sim,
                    "contact_sim_sp": contact_sim_sp,
                    "n_matched": n_matched,
                    "minimizedAffinity": affinity,
                    "CNNscore": cnn_score,
                    "CNNaffinity": cnn_affinity,
                })

    # 7. Sort by contact similarity (descending) and write CSV
    results.sort(key=lambda r: (-r["contact_sim"] if not np.isnan(r["contact_sim"]) else 1))

    out_csv = os.path.join(script_dir, "pose_similarity.csv")
    fieldnames = ["directory", "sdf_file", "hsa_prefix", "model_num",
                  "rmsd", "rmsd_sp", "contact_sim", "contact_sim_sp",
                  "n_matched", "minimizedAffinity", "CNNscore", "CNNaffinity"]

    with open(out_csv, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        for r in results:
            row = dict(r)
            for k in ("rmsd", "rmsd_sp", "contact_sim", "contact_sim_sp",
                       "minimizedAffinity", "CNNscore", "CNNaffinity"):
                row[k] = f"{r[k]:.4f}"
            writer.writerow(row)

    print(f"\n{'='*60}")
    print(f"Results written to: {out_csv}")
    print(f"Total poses evaluated: {len(results)}")

    # --- Rankings by contact similarity vs MD target ---
    print(f"\n{'='*60}")
    print(f"RANKING BY CONTACT SIMILARITY TO MD TARGET (higher = more similar)")
    print(f"{'='*60}")
    results_by_contact = sorted(results, key=lambda r: -r["contact_sim"] if not np.isnan(r["contact_sim"]) else 1)
    print(f"\nTop 20 most similar poses (by contact profile):")
    _print_table_full(results_by_contact[:20])

    # --- Rankings by RMSD vs MD target ---
    print(f"\n{'='*60}")
    print(f"RANKING BY RMSD TO MD TARGET (lower = more similar)")
    print(f"{'='*60}")
    results_by_rmsd = sorted(results, key=lambda r: r["rmsd"])
    print(f"\nTop 20 most similar poses (by RMSD):")
    _print_table_full(results_by_rmsd[:20])

    # --- Rankings by contact similarity vs starting pose ---
    if starting_pose_profile is not None:
        print(f"\n{'='*60}")
        print(f"RANKING BY CONTACT SIMILARITY TO STARTING POSE (phe_sssD hsa5 model 12)")
        print(f"{'='*60}")
        results_by_sp_contact = sorted(results, key=lambda r: -r["contact_sim_sp"] if not np.isnan(r["contact_sim_sp"]) else 1)
        print(f"\nTop 20:")
        _print_table_full(results_by_sp_contact[:20])

    # --- Rankings by RMSD vs starting pose ---
    if starting_pose_coords is not None:
        print(f"\n{'='*60}")
        print(f"RANKING BY RMSD TO STARTING POSE (phe_sssD hsa5 model 12)")
        print(f"{'='*60}")
        results_by_sp_rmsd = sorted(results, key=lambda r: r["rmsd_sp"] if not np.isnan(r["rmsd_sp"]) else float("inf"))
        print(f"\nTop 20:")
        _print_table_full(results_by_sp_rmsd[:20])

    # Per-molecule-type rankings (by contact similarity to target)
    mol_groups = defaultdict(list)
    for r in results:
        dirname = r["directory"]
        mol_type = dirname.split("_")[0]
        mol_groups[mol_type].append(r)

    for mol_type in sorted(mol_groups):
        group = mol_groups[mol_type]
        group_by_contact = sorted(group, key=lambda r: -r["contact_sim"] if not np.isnan(r["contact_sim"]) else 1)

        dir_best = {}
        for r in group_by_contact:
            d = r["directory"]
            if d not in dir_best:
                dir_best[d] = r

        print(f"\n--- Best pose per directory for '{mol_type}' ligands (by contact sim to target) ---")
        _print_table_full(sorted(dir_best.values(), key=lambda r: -r["contact_sim"] if not np.isnan(r["contact_sim"]) else 1))

    threshold = 5.0
    n_similar = sum(1 for r in results if r["rmsd"] < threshold)
    print(f"\nPoses with RMSD < {threshold} A (similar binding mode): {n_similar}")
    print(f"  (Note: me_* RMSD includes {UNMATCHED_PENALTY}A penalty per unmatched target atom)")
    print(f"  contact_sim = Pearson correlation of residue-ligand min-distance profile vs MD target")
    print(f"  contact_sim_sp = same vs starting pose (phe_sssD hsa5 model 12)")
    print(f"  rmsd_sp = heavy atom RMSD vs starting pose")


def _print_table_full(rows):
    """Print a formatted results table with all metrics."""
    print(f"{'Rank':>4}  {'Directory':<20} {'SDF File':<30} {'Model':>5} {'RMSD':>8} {'RMSD_SP':>8} {'ContSim':>7} {'CS_SP':>7} {'Affinity':>10}")
    print("-" * 110)
    for i, r in enumerate(rows):
        fmt = lambda v: f"{v:.4f}" if isinstance(v, float) and not np.isnan(v) else "nan"
        print(f"{i+1:>4}  {r['directory']:<20} {r['sdf_file']:<30} {r['model_num']:>5} {fmt(r['rmsd']):>8} {fmt(r['rmsd_sp']):>8} {fmt(r['contact_sim']):>7} {fmt(r['contact_sim_sp']):>7} {fmt(r['minimizedAffinity']):>10}")


if __name__ == "__main__":
    main()
