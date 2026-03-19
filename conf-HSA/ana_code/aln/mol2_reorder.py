#!/usr/bin/env python3
"""
mol2_reorder.py
===============
Reorder and/or complete atom ordering of molecular structure files to match
a mol2 template, using graph isomorphism (VF2 via networkx).

Two modes
---------
mol2   Reorder a mol2 file so its atom order matches a reference mol2.
       Atom names, types, charges and bond table come from the reference;
       coordinates come from the input mol2.

sdf    Take a docking result SDF (Gnina / Vina style), select one pose,
       place back missing C-H hydrogens from the template mol2 using a
       local-frame reconstruction, and write a properly ordered mol2.

       SDF scoring fields understood:
         minimizedAffinity  Vina score (kcal/mol, lower = better)  [DEFAULT]
         CNN_VS             CNNaffinity * CNNscore (higher = better)
         CNNaffinity        binding affinity prediction in pK units  (higher = better)
         CNNscore           probability pose is within 2 Å of truth  (higher = better)

Usage
-----
  # Reorder mol2 to match template:
  python mol2_reorder.py mol2 -t template.mol2 -i input.mol2 -o output.mol2

  # Convert best SDF pose (by minimizedAffinity) to mol2:
  python mol2_reorder.py sdf -t template.mol2 -i dock_out.sdf -o best_pose.mol2

  # Pick specific model number:
  python mol2_reorder.py sdf -t template.mol2 -i dock_out.sdf -o pose3.mol2 --model 3

  # Pick best by CNN_VS:
  python mol2_reorder.py sdf -t template.mol2 -i dock_out.sdf -o best_cnnvs.mol2 --metric CNN_VS

  # List all models with scores:
  python mol2_reorder.py sdf -t template.mol2 -i dock_out.sdf --list-models
"""

import argparse
import os
import sys
import textwrap
from collections import defaultdict

import numpy as np
import networkx as nx
from networkx.algorithms import isomorphism

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

# Metrics where LOWER value is better (all others: higher is better)
LOWER_IS_BETTER = {"minimizedAffinity"}

# Default metric
DEFAULT_METRIC = "minimizedAffinity"

# ---------------------------------------------------------------------------
# Mol2 I/O
# ---------------------------------------------------------------------------

def parse_mol2(path):
    """
    Parse a mol2 file.

    Returns
    -------
    atoms : list of dict
        Keys: idx, name, x, y, z, type, element, res_num, res_name, charge
    bonds : list of dict
        Keys: idx, a1, a2, type
    header_lines : list of str
        Lines in @<TRIPOS>MOLECULE block (the first 5 data lines).
    """
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
                    "idx":      int(parts[0]),
                    "name":     parts[1],
                    "x":        float(parts[2]),
                    "y":        float(parts[3]),
                    "z":        float(parts[4]),
                    "type":     atype,
                    "element":  el,
                    "res_num":  parts[6] if len(parts) > 6 else "1",
                    "res_name": parts[7] if len(parts) > 7 else "UNL1",
                    "charge":   parts[8] if len(parts) > 8 else "0.0000",
                })

            elif section == "@<TRIPOS>BOND":
                parts = line.split()
                if len(parts) >= 4:
                    bonds.append({
                        "idx":  int(parts[0]),
                        "a1":   int(parts[1]),
                        "a2":   int(parts[2]),
                        "type": parts[3],
                    })

    return atoms, bonds, header_lines


def write_mol2(path, atoms, bonds, header_lines=None, title="*****"):
    """
    Write a mol2 file.

    atoms  : list of dicts with keys idx, name, x, y, z, type,
             res_num, res_name, charge (indices are 1-based, sequential)
    bonds  : list of dicts with keys idx, a1, a2, type
    """
    _no_overwrite(path)
    if header_lines is None:
        header_lines = [title, f" {len(atoms)} {len(bonds)} 0 0 0", "SMALL", "GASTEIGER", ""]

    with open(path, "w") as fh:
        fh.write("@<TRIPOS>MOLECULE\n")
        # Rebuild the counts line correctly
        for i, hl in enumerate(header_lines):
            if i == 1:
                fh.write(f" {len(atoms)} {len(bonds)} 0 0 0\n")
            else:
                fh.write(hl + "\n")
        fh.write("\n")
        fh.write("@<TRIPOS>ATOM\n")
        for a in atoms:
            fh.write(
                f"{a['idx']:>7d} "
                f"{a['name']:<12s}"
                f"{a['x']:>9.4f} "
                f"{a['y']:>9.4f} "
                f"{a['z']:>9.4f} "
                f"{a['type']:<8s}"
                f"{a['res_num']:>2s}  "
                f"{a['res_name']:<12s}"
                f"{a['charge']}\n"
            )
        fh.write("@<TRIPOS>BOND\n")
        for b in bonds:
            fh.write(f"{b['idx']:>6d} {b['a1']:>5d} {b['a2']:>5d}    {b['type']}\n")


# ---------------------------------------------------------------------------
# SDF I/O
# ---------------------------------------------------------------------------

def parse_sdf(path):
    """
    Parse an SDF file into a list of model dicts.

    Each model dict contains:
        model_num   : 1-based index
        title       : first line of the record
        atoms       : list of {idx, el, x, y, z}
                      (idx 1-based; '*' element mapped to 'EU')
        bonds       : list of {a1, a2, type}  (type is int)
        scores      : dict of property_name -> float
        raw_lines   : all raw lines for this record
    """
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
                "el":  el,
                "x":   float(ln[0:10]),
                "y":   float(ln[10:20]),
                "z":   float(ln[20:30]),
            })

        bonds = []
        for i in range(4 + ac, 4 + ac + bc):
            if i >= len(lines):
                break
            ln = lines[i]
            bonds.append({
                "a1":   int(ln[0:3]),
                "a2":   int(ln[3:6]),
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
            "title":     title,
            "atoms":     atoms,
            "bonds":     bonds,
            "scores":    scores,
            "raw_lines": lines,
        })

    return models


def select_model(models, metric=DEFAULT_METRIC, model_num=None):
    """
    Select a model from a list.

    Parameters
    ----------
    model_num : int or None
        If given, return that 1-based model number.
    metric : str
        Score field to optimise when model_num is None.
    """
    if model_num is not None:
        for m in models:
            if m["model_num"] == model_num:
                return m
        raise ValueError(f"Model {model_num} not found (file has {len(models)} models).")

    candidates = [m for m in models if metric in m["scores"]]
    if not candidates:
        raise ValueError(
            f"Metric '{metric}' not found in any model.\n"
            f"Available metrics: {sorted(models[0]['scores'].keys())}"
        )

    reverse = metric not in LOWER_IS_BETTER   # True = higher is better
    best = sorted(candidates, key=lambda m: m["scores"][metric], reverse=reverse)[0]
    return best


def list_models(models, metrics=None):
    """Print a table of all models with their scores."""
    if metrics is None:
        # Collect all score keys present
        all_keys = []
        seen = set()
        for m in models:
            for k in m["scores"]:
                if k not in seen:
                    all_keys.append(k)
                    seen.add(k)
        metrics = all_keys

    header = f"{'Model':>6}" + "".join(f"  {k:>16}" for k in metrics)
    print(header)
    print("-" * len(header))
    for m in models:
        row = f"{m['model_num']:>6}"
        for k in metrics:
            val = m["scores"].get(k, float("nan"))
            row += f"  {val:>16.4f}"
        print(row)


# ---------------------------------------------------------------------------
# Graph utilities
# ---------------------------------------------------------------------------

def _element(atom_type):
    """Extract element symbol from mol2 atom type, normalising Eu variants."""
    el = atom_type.split(".")[0].upper()
    if el.startswith("EU"):
        el = "EU"
    return el


def _build_graph(atoms, bonds, heavy_only=False):
    """
    Build a networkx Graph from atom/bond lists.

    atoms  : list of dicts with 'idx' (int) and 'element' (str)
    bonds  : list of dicts with 'a1' and 'a2'
    heavy_only : if True, skip H atoms and bonds to/from H
    """
    G = nx.Graph()
    incl = {a["idx"] for a in atoms if not heavy_only or a["element"] != "H"}
    for a in atoms:
        if a["idx"] in incl:
            G.add_node(a["idx"], element=a["element"])
    for b in bonds:
        if b["a1"] in incl and b["a2"] in incl:
            G.add_edge(b["a1"], b["a2"])
    return G


def _node_match(n1, n2):
    return n1["element"] == n2["element"]


def find_isomorphism(G_ref, G_src):
    """
    Return the first isomorphism mapping ref_node -> src_node, or raise.
    """
    gm = isomorphism.GraphMatcher(G_ref, G_src, node_match=_node_match)
    if not gm.is_isomorphic():
        raise ValueError(
            "No graph isomorphism found between the two molecules.\n"
            "Check that both files represent the same molecular graph."
        )
    return next(gm.isomorphisms_iter())   # ref_idx -> src_idx


# ---------------------------------------------------------------------------
# Hydrogen placement
# ---------------------------------------------------------------------------

def _vec(a, b):
    """Unit vector from a to b (numpy arrays)."""
    v = b - a
    norm = np.linalg.norm(v)
    if norm < 1e-10:
        raise ValueError("Zero-length vector encountered.")
    return v / norm


def _get_coord(atom):
    return np.array([atom["x"], atom["y"], atom["z"]], dtype=float)


def place_h_from_template(
    ref_heavy_atom,   # mol2 dict for the heavy atom bearing H
    ref_h_atom,       # mol2 dict for the H
    ref_heavy_nbrs,   # list of mol2 dicts: other heavy atoms bonded to ref_heavy_atom
    new_heavy_coord,  # np.array: new (docked) coords of the heavy atom
    new_nbr_coords,   # list of np.array: new coords of the same heavy neighbours
):
    """
    Place a hydrogen in the new coordinate frame using the local geometry
    of the template.

    Strategy
    --------
    Express the H position in a local orthonormal frame built from the
    heavy atom and its neighbours.  Re-apply the same frame to the new
    (docked) heavy atom coordinates.

    The frame is:
        origin = heavy atom
        e1 = mean of (unit vectors to all heavy neighbours)  [if >0 nbrs]
        e2 = perpendicular component of the first neighbour direction
        e3 = e1 x e2

    For isolated atoms (no heavy neighbours) the template H offset vector
    is returned unchanged (rare edge case).
    """
    P_heavy = _get_coord(ref_heavy_atom)
    P_h     = _get_coord(ref_h_atom)
    dh      = P_h - P_heavy          # H offset in template frame

    if not ref_heavy_nbrs:
        # No heavy neighbours: just translate the H offset to new coords
        new_pos = new_heavy_coord + dh
        return new_pos

    # Build orthonormal frame from template
    nbr_vecs_ref = [_vec(P_heavy, _get_coord(n)) for n in ref_heavy_nbrs]
    nbr_vecs_new = [_vec(new_heavy_coord, c) for c in new_nbr_coords]

    # e1: average neighbour direction (template)
    e1_ref = np.mean(nbr_vecs_ref, axis=0)
    norm_e1 = np.linalg.norm(e1_ref)
    if norm_e1 < 1e-8:
        # Degenerate: neighbours cancel; fall back to first neighbour
        e1_ref = nbr_vecs_ref[0]
    else:
        e1_ref /= norm_e1

    # e2: component of first neighbour perp to e1 (template)
    v2_ref = nbr_vecs_ref[0] - np.dot(nbr_vecs_ref[0], e1_ref) * e1_ref
    if np.linalg.norm(v2_ref) < 1e-8:
        # e1 is parallel to first neighbour; try second or construct arbitrary perp
        if len(nbr_vecs_ref) > 1:
            v2_ref = nbr_vecs_ref[1] - np.dot(nbr_vecs_ref[1], e1_ref) * e1_ref
        if np.linalg.norm(v2_ref) < 1e-8:
            # construct arbitrary perpendicular
            arb = np.array([1.0, 0.0, 0.0])
            if abs(np.dot(e1_ref, arb)) > 0.9:
                arb = np.array([0.0, 1.0, 0.0])
            v2_ref = arb - np.dot(arb, e1_ref) * e1_ref
    e2_ref = v2_ref / np.linalg.norm(v2_ref)
    e3_ref = np.cross(e1_ref, e2_ref)

    # Express H offset in template frame
    h_e1 = np.dot(dh, e1_ref)
    h_e2 = np.dot(dh, e2_ref)
    h_e3 = np.dot(dh, e3_ref)

    # Build matching frame in new (docked) coordinates
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


# ---------------------------------------------------------------------------
# Core: mol2 -> mol2 reorder
# ---------------------------------------------------------------------------

def reorder_mol2_to_mol2(ref_path, src_path, out_path):
    """
    Reorder src_path so its atom order matches ref_path.

    - Atom names, types, residue labels, bond table: from ref
    - Coordinates: from src
    - Output: out_path (never overwrites)
    """
    _no_overwrite(out_path)
    print(f"Parsing reference:  {ref_path}")
    ref_atoms, ref_bonds, ref_header = parse_mol2(ref_path)
    print(f"  {len(ref_atoms)} atoms, {len(ref_bonds)} bonds")

    print(f"Parsing source:     {src_path}")
    src_atoms, src_bonds, _ = parse_mol2(src_path)
    print(f"  {len(src_atoms)} atoms, {len(src_bonds)} bonds")

    print("Building molecular graphs and finding isomorphism...")
    G_ref = _build_graph(ref_atoms, ref_bonds)
    G_src = _build_graph(src_atoms, src_bonds)
    mapping = find_isomorphism(G_ref, G_src)   # ref_idx -> src_idx
    print(f"  Isomorphism found ({len(mapping)} atoms).")

    src_by_idx = {a["idx"]: a for a in src_atoms}

    # Assemble output atoms: ref order, src coordinates
    out_atoms = []
    for new_pos, ref_a in enumerate(ref_atoms, 1):
        src_a = src_by_idx[mapping[ref_a["idx"]]]
        out_atoms.append({
            "idx":      new_pos,
            "name":     ref_a["name"],
            "x":        src_a["x"],
            "y":        src_a["y"],
            "z":        src_a["z"],
            "type":     ref_a["type"],
            "element":  ref_a["element"],
            "res_num":  ref_a["res_num"],
            "res_name": ref_a["res_name"],
            "charge":   ref_a["charge"],
        })

    # Bond table unchanged from ref
    write_mol2(out_path, out_atoms, ref_bonds, ref_header)
    print(f"Written: {out_path}")


# ---------------------------------------------------------------------------
# Core: SDF pose -> mol2
# ---------------------------------------------------------------------------

def sdf_pose_to_mol2(
    ref_path, sdf_path, out_path,
    metric=DEFAULT_METRIC, model_num=None
):
    """
    Convert one SDF pose to a mol2, placing missing C-H hydrogens from
    the template mol2.

    Workflow
    --------
    1. Parse reference mol2 (template with full H).
    2. Parse SDF, select best pose (or model_num).
    3. Graph isomorphism on heavy atoms only: find SDF_atom -> ref_heavy_atom.
    4. For each H in template:
       a. If the parent heavy atom has an explicit H in SDF (N-H case):
          use the SDF H coordinates directly.
       b. Otherwise (C-H): reconstruct using local frame from template geometry
          applied to the docked heavy-atom coordinates.
    5. Assemble output in reference mol2 atom order; write mol2.
    """
    _no_overwrite(out_path)

    print(f"Parsing template:   {ref_path}")
    ref_atoms, ref_bonds, ref_header = parse_mol2(ref_path)
    print(f"  {len(ref_atoms)} atoms, {len(ref_bonds)} bonds")

    print(f"Parsing SDF:        {sdf_path}")
    models = parse_sdf(sdf_path)
    print(f"  {len(models)} models found")

    model = select_model(models, metric=metric, model_num=model_num)
    sc = model["scores"]
    print(
        f"  Selected model {model['model_num']}  "
        f"minimizedAffinity={sc.get('minimizedAffinity', float('nan')):.4f}  "
        f"CNNscore={sc.get('CNNscore', float('nan')):.4f}  "
        f"CNNaffinity={sc.get('CNNaffinity', float('nan')):.4f}  "
        f"CNN_VS={sc.get('CNN_VS', float('nan')):.4f}"
    )

    sdf_atoms = model["atoms"]
    sdf_bonds = model["bonds"]

    # --- Build adjacency helpers ---
    ref_by_idx = {a["idx"]: a for a in ref_atoms}
    ref_adj = defaultdict(list)
    for b in ref_bonds:
        ref_adj[b["a1"]].append(b["a2"])
        ref_adj[b["a2"]].append(b["a1"])

    sdf_by_idx = {a["idx"]: a for a in sdf_atoms}
    sdf_adj = defaultdict(list)
    for b in sdf_bonds:
        sdf_adj[b["a1"]].append(b["a2"])
        sdf_adj[b["a2"]].append(b["a1"])

    # --- Graph isomorphism: heavy atoms only ---
    print("Finding heavy-atom isomorphism (SDF -> template)...")
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
    mapping = find_isomorphism(G_ref_h, G_sdf_h)   # ref_heavy_idx -> sdf_idx
    print(f"  Isomorphism found ({len(mapping)} heavy atoms).")

    # --- Collect which SDF N atoms have explicit H ---
    # sdf_n_h_map: sdf_N_idx -> list of sdf_H_idx
    sdf_n_h_map = defaultdict(list)
    for a in sdf_atoms:
        if a["el"] == "N":
            for nb_idx in sdf_adj[a["idx"]]:
                if sdf_by_idx[nb_idx]["el"] == "H":
                    sdf_n_h_map[a["idx"]].append(nb_idx)

    # Track which SDF H atoms have been consumed (N-H)
    used_sdf_h = set()

    # --- Assign coordinates for all ref atoms ---
    # new_coords[ref_idx] = (x, y, z)
    new_coords = {}

    # Heavy atoms: use SDF coordinates
    for ref_idx, sdf_idx in mapping.items():
        sdf_a = sdf_by_idx[sdf_idx]
        new_coords[ref_idx] = (sdf_a["x"], sdf_a["y"], sdf_a["z"])

    # H atoms
    for ref_h in ref_atoms:
        if ref_h["element"] != "H":
            continue

        # Find the unique heavy atom bonded to this H in the template
        heavy_nbrs_ref = [ref_by_idx[n] for n in ref_adj[ref_h["idx"]]
                          if ref_by_idx[n]["element"] != "H"]
        if not heavy_nbrs_ref:
            # Isolated H (shouldn't happen) - skip
            print(f"  WARNING: H atom {ref_h['idx']} has no heavy-atom neighbour in template.")
            new_coords[ref_h["idx"]] = (ref_h["x"], ref_h["y"], ref_h["z"])
            continue

        parent_ref = heavy_nbrs_ref[0]  # the heavy atom bearing this H

        # Is the parent atom a N that has explicit H in SDF?
        if parent_ref["element"] == "N":
            sdf_parent_idx = mapping[parent_ref["idx"]]
            avail_sdf_h = [h for h in sdf_n_h_map[sdf_parent_idx]
                           if h not in used_sdf_h]
            if avail_sdf_h:
                # Use the SDF H position directly
                sdf_h = sdf_by_idx[avail_sdf_h[0]]
                used_sdf_h.add(avail_sdf_h[0])
                new_coords[ref_h["idx"]] = (sdf_h["x"], sdf_h["y"], sdf_h["z"])
                continue

        # General case (C-H or N-H not in SDF): reconstruct from local frame
        sdf_parent_idx = mapping[parent_ref["idx"]]
        new_parent_coord = np.array(new_coords[parent_ref["idx"]])

        # Other heavy neighbours of parent in template
        other_heavy_ref = [ref_by_idx[n] for n in ref_adj[parent_ref["idx"]]
                           if n != ref_h["idx"] and ref_by_idx[n]["element"] != "H"]

        # Their new (docked) coordinates
        other_heavy_new_coords = []
        for oh in other_heavy_ref:
            if oh["idx"] in new_coords:
                other_heavy_new_coords.append(np.array(new_coords[oh["idx"]]))
            else:
                # Heavy neighbour not yet assigned (shouldn't happen)
                other_heavy_new_coords.append(_get_coord(oh))

        new_h_xyz = place_h_from_template(
            parent_ref, ref_h, other_heavy_ref,
            new_parent_coord, other_heavy_new_coords,
        )
        new_coords[ref_h["idx"]] = tuple(new_h_xyz)

    # --- Assemble output in ref atom order ---
    out_atoms = []
    for new_pos, ref_a in enumerate(ref_atoms, 1):
        x, y, z = new_coords[ref_a["idx"]]
        out_atoms.append({
            "idx":      new_pos,
            "name":     ref_a["name"],
            "x":        x,
            "y":        y,
            "z":        z,
            "type":     ref_a["type"],
            "element":  ref_a["element"],
            "res_num":  ref_a["res_num"],
            "res_name": ref_a["res_name"],
            "charge":   ref_a["charge"],
        })

    write_mol2(out_path, out_atoms, ref_bonds, ref_header)
    print(f"Written: {out_path}")


# ---------------------------------------------------------------------------
# Safety helper
# ---------------------------------------------------------------------------

def _no_overwrite(path):
    if os.path.exists(path):
        print(f"ERROR: '{path}' already exists. Refusing to overwrite.", file=sys.stderr)
        sys.exit(1)


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def build_parser():
    parser = argparse.ArgumentParser(
        prog="mol2_reorder.py",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description=textwrap.dedent(__doc__),
    )
    sub = parser.add_subparsers(dest="mode", required=True)

    # --- mol2 subcommand ---
    p_mol2 = sub.add_parser(
        "mol2",
        help="Reorder a mol2 file to match a reference mol2.",
    )
    p_mol2.add_argument("-t", "--template", required=True, metavar="MOL2",
                        help="Reference mol2 (defines atom ordering).")
    p_mol2.add_argument("-i", "--input", required=True, metavar="MOL2",
                        help="Input mol2 to reorder.")
    p_mol2.add_argument("-o", "--output", required=True, metavar="MOL2",
                        help="Output mol2 filename (must not exist).")

    # --- sdf subcommand ---
    p_sdf = sub.add_parser(
        "sdf",
        help="Convert an SDF docking pose to mol2, rebuilding missing H.",
    )
    p_sdf.add_argument("-t", "--template", required=True, metavar="MOL2",
                       help="Reference mol2 with full hydrogens.")
    p_sdf.add_argument("-i", "--input", required=True, metavar="SDF",
                       help="Input SDF docking result.")
    p_sdf.add_argument("-o", "--output", metavar="MOL2",
                       help="Output mol2 filename (must not exist). "
                            "Omit when using --list-models.")
    p_sdf.add_argument("--model", type=int, default=None, metavar="N",
                       help="Select a specific model number (1-based). "
                            "Default: auto-select by --metric.")
    p_sdf.add_argument(
        "--metric",
        default=DEFAULT_METRIC,
        metavar="NAME",
        help=(
            f"Score field to optimise when --model is not given. "
            f"Default: '{DEFAULT_METRIC}' (lower is better). "
            f"Other useful choices: CNN_VS, CNNaffinity, CNNscore "
            f"(all higher is better)."
        ),
    )
    p_sdf.add_argument("--list-models", action="store_true",
                       help="Print all models with scores and exit.")

    return parser


def main():
    parser = build_parser()
    args = parser.parse_args()

    if args.mode == "mol2":
        reorder_mol2_to_mol2(args.template, args.input, args.output)

    elif args.mode == "sdf":
        models = parse_sdf(args.input)
        if args.list_models:
            list_models(models)
            return
        if not args.output:
            parser.error("--output is required unless --list-models is specified.")
        sdf_pose_to_mol2(
            ref_path=args.template,
            sdf_path=args.input,
            out_path=args.output,
            metric=args.metric,
            model_num=args.model,
        )


if __name__ == "__main__":
    main()
