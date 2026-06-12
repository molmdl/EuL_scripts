#!/usr/bin/env python3
"""
T1_atom_mapping.py — Build SAP↔TSAP atom-index mappings.

Reads optimized XYZ reference structures for SAP and TSAP conformers,
uses NetworkX graph isomorphism on heavy-atom bond graphs to enumerate
candidate atom mappings, scores them by Pearson correlation of pairwise
distance matrices, extends to H atoms, validates, and saves the result.

Usage:
    python scripts/T1_atom_mapping.py --data-dir data --out-dir analysis
"""

import argparse
import json
import os
import shutil
import sys
from collections import Counter

import numpy as np
from scipy.stats import pearsonr

# ---------------------------------------------------------------------------
# Bond cutoffs (Å)
# ---------------------------------------------------------------------------
BOND_CUTOFFS = {
    # (sorted element pair): max distance
    ("C", "C"): 1.65,
    ("C", "N"): 1.55,
    ("C", "O"): 1.50,
    ("Eu", "N"): 2.80,
    ("Eu", "O"): 2.70,
}
DEFAULT_HEAVY_CUTOFF = 1.80

H_CUTOFFS = {"C": 1.20, "N": 1.15, "O": 1.15}


# ---------------------------------------------------------------------------
# Step 1: XYZ parser
# ---------------------------------------------------------------------------
def read_xyz(path):
    """Return (elements: list[str], coords: np.ndarray shape (N,3))."""
    with open(path) as f:
        n = int(f.readline().strip())
        f.readline()  # comment line
        elems, coords = [], []
        for _ in range(n):
            parts = f.readline().split()
            elems.append(parts[0])
            coords.append([float(parts[1]), float(parts[2]), float(parts[3])])
    return elems, np.array(coords)


# ---------------------------------------------------------------------------
# Step 3: Build heavy-atom bond graph
# ---------------------------------------------------------------------------
def build_heavy_bond_graph(elems, coords, heavy_idx):
    """
    Build NetworkX graph from heavy atoms.
    Node attrs: element (str)
    Edge attrs: bond_type (sorted tuple of elements)
    """
    import networkx as nx

    G = nx.Graph()
    for i in heavy_idx:
        G.add_node(i, element=elems[i])
    for a in range(len(heavy_idx)):
        for b in range(a + 1, len(heavy_idx)):
            i, j = heavy_idx[a], heavy_idx[b]
            d = np.linalg.norm(coords[i] - coords[j])
            key = tuple(sorted([elems[i], elems[j]]))
            cutoff = BOND_CUTOFFS.get(key, DEFAULT_HEAVY_CUTOFF)
            if d <= cutoff:
                G.add_edge(i, j, bond_type=key)
    return G


# ---------------------------------------------------------------------------
# Step 4: Find all graph isomorphisms
# ---------------------------------------------------------------------------
def find_all_isomorphisms(G_sap, G_tsap):
    """Return list of dicts {tsap_idx: sap_idx} for heavy atoms."""
    from networkx.algorithms.isomorphism import GraphMatcher

    nm = lambda d1, d2: d1["element"] == d2["element"]
    em = lambda d1, d2: d1["bond_type"] == d2["bond_type"]

    gm = GraphMatcher(G_tsap, G_sap, node_match=nm, edge_match=em)
    isos = list(gm.isomorphisms_iter())
    return isos  # each is {tsap_node: sap_node}


# ---------------------------------------------------------------------------
# Step 5: Score isomorphism by distance correlation
# ---------------------------------------------------------------------------
def score_isomorphism(sap_coords, tsap_coords, iso_mapping):
    """
    iso_mapping: {tsap_idx: sap_idx} for heavy atoms.
    Compute Pearson r of all pairwise Euclidean distances between
    the mapped atoms in SAP vs the same atoms in TSAP.
    """
    tsap_indices = list(iso_mapping.keys())
    sap_indices = [iso_mapping[t] for t in tsap_indices]

    sap_pts = sap_coords[sap_indices]
    tsap_pts = tsap_coords[tsap_indices]
    n = len(sap_pts)
    sap_dists = []
    tsap_dists = []
    for i in range(n):
        for j in range(i + 1, n):
            sap_dists.append(np.linalg.norm(sap_pts[i] - sap_pts[j]))
            tsap_dists.append(np.linalg.norm(tsap_pts[i] - tsap_pts[j]))

    r, _ = pearsonr(sap_dists, tsap_dists)
    return r


# ---------------------------------------------------------------------------
# Step 6: Extend heavy-atom mapping to H atoms
# ---------------------------------------------------------------------------
def find_h_parents(elems, coords, heavy_idx, h_idx):
    """Return dict {h_idx: parent_heavy_idx}."""
    parents = {}
    for h in h_idx:
        best_d, best_p = 999.0, None
        for p in heavy_idx:
            d = np.linalg.norm(coords[h] - coords[p])
            cutoff = H_CUTOFFS.get(elems[p], 1.20)
            if d < cutoff and d < best_d:
                best_d, best_p = d, p
        if best_p is not None:
            parents[h] = best_p
    return parents


def extend_mapping_with_hydrogens(
    sap_elems, sap_coords, tsap_elems, tsap_coords, heavy_mapping
):
    """
    heavy_mapping: {tsap_heavy_idx: sap_heavy_idx} (from Step 5)
    Returns full mapping: {tsap_idx: sap_idx} for ALL atoms.
    """
    sap_heavy_set = list(set(heavy_mapping.values()))
    tsap_heavy_set = list(heavy_mapping.keys())

    sap_h_idx = [i for i, e in enumerate(sap_elems) if e == "H"]
    tsap_h_idx = [i for i, e in enumerate(tsap_elems) if e == "H"]

    sap_parents = find_h_parents(sap_elems, sap_coords, sap_heavy_set, sap_h_idx)
    tsap_parents = find_h_parents(tsap_elems, tsap_coords, tsap_heavy_set, tsap_h_idx)

    # Group SAP H atoms by their heavy-atom parent
    sap_H_by_parent = {}
    for h_idx, p_idx in sap_parents.items():
        sap_H_by_parent.setdefault(p_idx, []).append(h_idx)

    full_mapping = dict(heavy_mapping)  # start with heavy atoms

    # For each TSAP H atom, find which TSAP heavy atom it's bonded to,
    # map that heavy atom to SAP, then find the best-matching SAP H bonded to same parent
    for t_h in tsap_h_idx:
        t_parent = tsap_parents.get(t_h)
        if t_parent is None:
            continue  # orphan H — should not happen
        s_parent = heavy_mapping.get(t_parent)
        if s_parent is None:
            continue
        candidate_sap_H = sap_H_by_parent.get(s_parent, [])
        if len(candidate_sap_H) == 0:
            continue
        if len(candidate_sap_H) == 1:
            full_mapping[t_h] = candidate_sap_H[0]
        else:
            # Multiple H on same parent: match by 3D geometry
            t_vec = tsap_coords[t_h] - tsap_coords[t_parent]
            best_match, best_score = None, -2.0
            for s_h in candidate_sap_H:
                if s_h in full_mapping.values():
                    continue  # already assigned
                s_vec = sap_coords[s_h] - sap_coords[s_parent]
                # Cosine similarity of position relative to parent
                norm_t = np.linalg.norm(t_vec)
                norm_s = np.linalg.norm(s_vec)
                if norm_t < 1e-12 or norm_s < 1e-12:
                    continue
                score = np.dot(t_vec, s_vec) / (norm_t * norm_s)
                if score > best_score:
                    best_score = score
                    best_match = s_h
            if best_match is not None:
                full_mapping[t_h] = best_match
                # Remove from candidates so it's not assigned again
                sap_H_by_parent[s_parent].remove(best_match)

    return full_mapping


# ---------------------------------------------------------------------------
# Step 7: Convert mapping dict to array
# ---------------------------------------------------------------------------
def mapping_dict_to_array(full_mapping, n_atoms):
    """Convert {tsap_idx: sap_idx} to list where arr[tsap_idx] = sap_idx."""
    arr = [-1] * n_atoms
    for t, s in full_mapping.items():
        arr[t] = s
    # Verify no unmapped atoms
    unmapped = [i for i, v in enumerate(arr) if v == -1]
    assert -1 not in arr, f"Unmapped atoms: {unmapped}"
    return arr


# ---------------------------------------------------------------------------
# Step 8: Validation
# ---------------------------------------------------------------------------
def validate_mapping(
    sap_elems, sap_coords, tsap_elems, tsap_coords, mapping_array, label
):
    """mapping_array: tsap_to_sap[tsap_idx] = sap_idx"""
    n = len(mapping_array)
    print(f"\n=== Validation: {label} ===")

    # 1. Element counts match
    sap_counts = Counter(sap_elems)
    tsap_counts = Counter(tsap_elems)
    assert sap_counts == tsap_counts, (
        f"Element count mismatch: {dict(sap_counts)} vs {dict(tsap_counts)}"
    )
    print(f"  Element counts: MATCH ({dict(sap_counts)})")

    # 2. Every mapped pair has same element
    mismatches = []
    for t in range(n):
        s = mapping_array[t]
        if tsap_elems[t] != sap_elems[s]:
            mismatches.append((t, tsap_elems[t], s, sap_elems[s]))
    assert len(mismatches) == 0, f"Element mismatches: {mismatches[:10]}"
    print(f"  Element consistency: PASS (all {n} atoms match)")

    # 3. Bond count match
    heavy_sap = [i for i, e in enumerate(sap_elems) if e != "H"]
    heavy_tsap = [i for i, e in enumerate(tsap_elems) if e != "H"]
    G_sap = build_heavy_bond_graph(sap_elems, sap_coords, heavy_sap)
    G_tsap = build_heavy_bond_graph(tsap_elems, tsap_coords, heavy_tsap)
    print(f"  SAP heavy bonds: {G_sap.number_of_edges()}")
    print(f"  TSAP heavy bonds: {G_tsap.number_of_edges()}")
    assert G_sap.number_of_edges() == G_tsap.number_of_edges(), "Bond count mismatch"

    # 4. Bond preservation: every bond in TSAP maps to a bond in SAP
    broken = 0
    for i, j in G_tsap.edges():
        si, sj = mapping_array[i], mapping_array[j]
        if not G_sap.has_edge(si, sj):
            broken += 1
    print(
        f"  Bond preservation: {G_tsap.number_of_edges() - broken}/{G_tsap.number_of_edges()}"
    )
    assert broken == 0, f"{broken} bonds broken by mapping"

    # 5. Euclidean distance preservation (pairwise for heavy atoms)
    tsap_heavy_coords = tsap_coords[heavy_tsap]
    mapped_sap_indices = [mapping_array[t] for t in heavy_tsap]
    sap_heavy_mapped_coords = sap_coords[mapped_sap_indices]

    n_h = len(heavy_tsap)
    sap_d, tsap_d = [], []
    for i in range(n_h):
        for j in range(i + 1, n_h):
            sap_d.append(
                np.linalg.norm(sap_heavy_mapped_coords[i] - sap_heavy_mapped_coords[j])
            )
            tsap_d.append(
                np.linalg.norm(tsap_heavy_coords[i] - tsap_heavy_coords[j])
            )
    r, _ = pearsonr(sap_d, tsap_d)
    rmse = np.sqrt(np.mean([(s - t) ** 2 for s, t in zip(sap_d, tsap_d)]))
    print(f"  Distance correlation (heavy): r = {r:.4f}, RMSE = {rmse:.3f} Å")
    assert r > 0.95, f"Distance correlation too low: {r:.4f}"

    # 6. Bijectivity: mapping is 1-to-1
    mapped_targets = set(mapping_array)
    assert len(mapped_targets) == n, (
        f"Non-bijective: {n} sources -> {len(mapped_targets)} targets"
    )
    print(f"  Bijectivity: PASS (1-to-1)")

    return r


# ---------------------------------------------------------------------------
# Step 9-10: Main
# ---------------------------------------------------------------------------
def map_one_system(system, sap_elems, sap_coords, tsap_elems, tsap_coords):
    """Run the full mapping pipeline for one system (me or phe)."""
    n_atoms = len(sap_elems)
    n_heavy_sap = sum(1 for e in sap_elems if e != "H")
    n_heavy_tsap = sum(1 for e in tsap_elems if e != "H")
    print(f"\n{'='*60}")
    print(f"System: {system}")
    print(f"  SAP atoms: {n_atoms} ({n_heavy_sap} heavy)")
    print(f"  TSAP atoms: {len(tsap_elems)} ({n_heavy_tsap} heavy)")
    assert n_atoms == len(tsap_elems), "Atom count mismatch between SAP and TSAP"

    # Step 2: Separate heavy and H
    heavy_sap = [i for i, e in enumerate(sap_elems) if e != "H"]
    heavy_tsap = [i for i, e in enumerate(tsap_elems) if e != "H"]
    h_sap = [i for i, e in enumerate(sap_elems) if e == "H"]
    h_tsap = [i for i, e in enumerate(tsap_elems) if e == "H"]
    print(f"  Heavy atoms: {len(heavy_sap)} (SAP), {len(heavy_tsap)} (TSAP)")
    print(f"  H atoms: {len(h_sap)} (SAP), {len(h_tsap)} (TSAP)")

    # Step 3: Build bond graphs
    print("  Building heavy-atom bond graphs...")
    G_sap = build_heavy_bond_graph(sap_elems, sap_coords, heavy_sap)
    G_tsap = build_heavy_bond_graph(tsap_elems, tsap_coords, heavy_tsap)
    print(f"  SAP bonds: {G_sap.number_of_edges()}, TSAP bonds: {G_tsap.number_of_edges()}")

    # Step 4: Find isomorphisms
    print("  Enumerating graph isomorphisms...")
    isos = find_all_isomorphisms(G_sap, G_tsap)
    print(f"  Found {len(isos)} isomorphisms")
    if len(isos) == 0:
        print("  ERROR: No isomorphisms found! Check bond cutoffs or element labels.")
        sys.exit(1)

    # Step 5: Score and select best
    print("  Scoring isomorphisms by distance correlation...")
    best_iso = None
    best_r = -1.0
    for iso in isos:
        r = score_isomorphism(sap_coords, tsap_coords, iso)
        if r > best_r:
            best_r = r
            best_iso = iso
    print(f"  Best isomorphism: r = {best_r:.4f}")

    # Step 6: Extend to H atoms
    print("  Extending mapping to H atoms...")
    heavy_mapping = dict(best_iso)
    full_mapping = extend_mapping_with_hydrogens(
        sap_elems, sap_coords, tsap_elems, tsap_coords, heavy_mapping
    )
    n_mapped = len(full_mapping)
    print(f"  Mapped {n_mapped}/{n_atoms} atoms")

    # Step 7: Convert to array
    mapping_array = mapping_dict_to_array(full_mapping, n_atoms)

    # Step 8: Validate
    val_r = validate_mapping(
        sap_elems, sap_coords, tsap_elems, tsap_coords, mapping_array, system
    )

    return {
        "mapping_array": mapping_array,
        "heavy_mapping": heavy_mapping,
        "n_atoms": n_atoms,
        "n_heavy": len(heavy_sap),
        "n_isomorphisms": len(isos),
        "best_r": float(best_r),
        "n_bonds": G_tsap.number_of_edges(),
    }


def main():
    parser = argparse.ArgumentParser(
        description="Build SAP↔TSAP atom index mappings"
    )
    parser.add_argument("--data-dir", default="data", help="Root data directory")
    parser.add_argument("--out-dir", default="analysis", help="Output directory")
    args = parser.parse_args()

    data_dir = args.data_dir
    out_dir = args.out_dir

    # ------------------------------------------------------------------
    # Load structures
    # ------------------------------------------------------------------
    print("Loading XYZ structures...")
    me_sap_elems, me_sap_coords = read_xyz(
        os.path.join(data_dir, "me_rrrD_sap", "xtbopt.xyz")
    )
    me_tsap_elems, me_tsap_coords = read_xyz(
        os.path.join(data_dir, "me_rrrD_tsap", "xtbopt.xyz")
    )
    phe_sap_elems, phe_sap_coords = read_xyz(
        os.path.join(data_dir, "phe_rrrD_sap", "xtbopt.xyz")
    )
    phe_tsap_elems, phe_tsap_coords = read_xyz(
        os.path.join(data_dir, "phe_rrrD_tsap", "xtbopt.xyz")
    )

    print(f"  me:  {len(me_sap_elems)} atoms (SAP), {len(me_tsap_elems)} atoms (TSAP)")
    print(f"  phe: {len(phe_sap_elems)} atoms (SAP), {len(phe_tsap_elems)} atoms (TSAP)")

    # ------------------------------------------------------------------
    # Map both systems
    # ------------------------------------------------------------------
    me_result = map_one_system(
        "me", me_sap_elems, me_sap_coords, me_tsap_elems, me_tsap_coords
    )
    phe_result = map_one_system(
        "phe", phe_sap_elems, phe_sap_coords, phe_tsap_elems, phe_tsap_coords
    )

    # ------------------------------------------------------------------
    # Step 9: Save output
    # ------------------------------------------------------------------
    os.makedirs(out_dir, exist_ok=True)
    out_path = os.path.join(out_dir, "atom_mappings.json")

    # Backup if exists
    if os.path.exists(out_path):
        backup = out_path + ".bak"
        shutil.copy2(out_path, backup)
        print(f"\nBacked up existing mapping to {backup}")

    output = {
        "me_tsap_to_sap": me_result["mapping_array"],
        "phe_tsap_to_sap": phe_result["mapping_array"],
        "metadata": {
            "me": {
                "n_atoms": me_result["n_atoms"],
                "n_heavy": me_result["n_heavy"],
                "n_isomorphisms": me_result["n_isomorphisms"],
                "best_iso_distance_corr": me_result["best_r"],
                "bond_preservation": f"{me_result['n_bonds']}/{me_result['n_bonds']}",
                "heavy_mapping": {
                    str(k): v for k, v in me_result["heavy_mapping"].items()
                },
            },
            "phe": {
                "n_atoms": phe_result["n_atoms"],
                "n_heavy": phe_result["n_heavy"],
                "n_isomorphisms": phe_result["n_isomorphisms"],
                "best_iso_distance_corr": phe_result["best_r"],
                "bond_preservation": f"{phe_result['n_bonds']}/{phe_result['n_bonds']}",
                "heavy_mapping": {
                    str(k): v for k, v in phe_result["heavy_mapping"].items()
                },
            },
            "method": "NetworkX graph isomorphism + geometric distance correlation scoring",
            "bond_cutoffs": {f"{k[0]}-{k[1]}": v for k, v in BOND_CUTOFFS.items()},
        },
    }

    with open(out_path, "w") as f:
        json.dump(output, f, indent=2)
    print(f"\nSaved mapping to {out_path}")
    print("Done!")


if __name__ == "__main__":
    main()
