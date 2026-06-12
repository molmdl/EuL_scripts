#!/usr/bin/env python3
"""
T1_find_indices.py
==================
Determine atom indices for the TSAP conformer that correspond to:
  1. Eu centre
  2. The 9 coordinating donor atoms (O/N)
  3. The N-C-C-N torsion atoms used by metal_geo_analysis.py

Outputs:
  analysis/indices.json          – machine-readable mapping
  analysis/indices_summary.txt   – human-readable report

Usage:
  python scripts/T1_find_indices.py --data-dir data --out-dir analysis
"""

import argparse
import json
import os
import sys

import MDAnalysis as mda
import numpy as np

# ═══════════════════════════════════════════════════════
# Known SAP indices (from metal_geo_analysis.py)
# NOTE: SAP_COORD_IDX is Eu+9 (9 donors including cap N63 at index 63).
# The Phase 02 revised Eu+8 benchmark excludes N63, giving 8 scaffold
# donors + Eu. See SAP_EU8_COORD_IDX below.
# ═══════════════════════════════════════════════════════
SAP_METAL_IDX = 54
SAP_COORD_IDX = [0, 1, 2, 3, 4, 5, 6, 7, 63]
# Eu+8 alignment atoms (Phase 02 benchmark): 8 scaffold donors + Eu,
# excluding cap N63 (index 63) due to asymmetry between SAP/TSAP.
SAP_EU8_COORD_IDX = [0, 1, 2, 3, 4, 5, 6, 7]
SAP_EU8_ALIGN_IDX = sorted(SAP_EU8_COORD_IDX + [SAP_METAL_IDX])  # [0,1,2,3,4,5,6,7,54]
SAP_BOTTOM = [3, 4, 5, 6]
SAP_TOP = [0, 1, 2, 7]
SAP_CAP = [63]
SAP_RING_TORS = [
    ('T1 N3-C8-C9-N4',    3,  8,  9,  4),
    ('T2 N4-C10-C11-N5',  4, 10, 11,  5),
    ('T3 N5-C12-C13-N6',  5, 12, 13,  6),
    ('T4 N6-C14-C15-N3',  6, 14, 15,  3),
]
SAP_CHROM_TOR = ('Tc N3-C16-C23-N7', 3, 16, 23, 7)


def warn(msg):
    print(f'  [WARN] {msg}', file=sys.stderr)


def load_structure(path):
    """Load a structure with MDAnalysis (no topology)."""
    if not os.path.isfile(path):
        raise FileNotFoundError(f'File not found: {path}')
    return mda.Universe(path)


def find_metal(u):
    """Find first atom whose element name is 'Eu' (case-insensitive)."""
    for atom in u.atoms:
        if atom.name.strip().lower() == 'eu':
            return int(atom.index)
    raise ValueError('No Eu atom found in the structure')


def closest_donors(u, metal_idx, n_donors=9, max_distance=3.5):
    """
    Return list of (index, element, distance) for the closest O/N atoms
    to the metal centre.
    """
    metal_pos = u.atoms[metal_idx].position
    donors = []
    for atom in u.atoms:
        name = atom.name.strip()
        if name in ('O', 'N'):
            d = float(np.linalg.norm(atom.position - metal_pos))
            if d <= max_distance:
                donors.append((int(atom.index), name, d))
    donors.sort(key=lambda x: x[2])
    return donors[:n_donors]


def get_neighbors(atoms, idx, cutoff=1.7):
    """Return list of (index, element, distance) for atoms within cutoff."""
    a = atoms[idx]
    out = []
    for b in atoms:
        if a == b:
            continue
        d = np.linalg.norm(a.position - b.position)
        if d <= cutoff:
            out.append((int(b.index), b.name.strip(), float(d)))
    out.sort(key=lambda x: x[2])
    return out


def dihedral(pos, atoms):
    """Compute a single dihedral angle in degrees [-180, 180]."""
    i, j, k, l = atoms
    b1 = pos[j] - pos[i]
    b2 = pos[k] - pos[j]
    b3 = pos[l] - pos[k]
    n1 = np.cross(b1, b2)
    nn1 = np.linalg.norm(n1)
    n2 = np.cross(b2, b3)
    nn2 = np.linalg.norm(n2)
    if nn1 < 1e-9 or nn2 < 1e-9:
        return float(np.nan)
    n1 /= nn1
    n2 /= nn2
    b2n = b2 / np.linalg.norm(b2)
    return float(np.degrees(np.arctan2(np.dot(np.cross(n1, n2), b2n), np.dot(n1, n2))))


def bfs_shortest_path(atoms, start, targets, max_depth=8, bond_cutoff=1.7):
    """
    Breadth-first search for shortest path from start to any target.
    Only considers C,N,O atoms (no H).
    Returns list of indices [start, ..., target] or None.
    """
    from collections import deque

    targets = set(targets)
    queue = deque([(start, [start])])
    visited = {start}
    while queue:
        current, path = queue.popleft()
        if current in targets:
            return path
        if len(path) >= max_depth:
            continue
        for n_idx, name, _ in get_neighbors(atoms, current, bond_cutoff):
            if n_idx in visited:
                continue
            if name not in ('C', 'N', 'O'):
                continue
            visited.add(n_idx)
            queue.append((n_idx, path + [n_idx]))
    return None


def identify_ring_nitrogens(u, n_indices):
    """
    Identify which N donor atoms belong to the macrocyclic ring.

    Strategy:
    - Count C neighbors (≤ 1.8 Å) for each N.
    - Aliphatic ring N's have 3 C neighbors in the DOTA backbone.
    - Aromatic chromophore N's have 2 C neighbors (part of an aromatic ring).
    """
    ring_ns = []
    chrom_ns = []
    n_to_c = {}
    for n_idx in n_indices:
        c_nbrs = [x for x in get_neighbors(u.atoms, n_idx, 1.8) if x[1] == 'C']
        n_to_c[n_idx] = [x[0] for x in c_nbrs]
        if len(c_nbrs) == 2:
            chrom_ns.append(n_idx)
        elif len(c_nbrs) >= 3:
            ring_ns.append(n_idx)
        else:
            warn(f'N {n_idx} has {len(c_nbrs)} C neighbors (unexpected)')
            chrom_ns.append(n_idx)

    return ring_ns, chrom_ns, n_to_c


def build_cc_bridges(u, ring_ns, n_to_c):
    """
    Build C-C bridge tuples that connect adjacent ring N's.
    Returns list of (c1, c2, (n1, n2)).
    """
    c_to_n = {}
    for n_idx, c_list in n_to_c.items():
        if n_idx in ring_ns:
            for c in c_list:
                c_to_n.setdefault(c, []).append(n_idx)

    cc_bridges = []
    all_cs = list(c_to_n.keys())
    for i in range(len(all_cs)):
        for j in range(i + 1, len(all_cs)):
            c1, c2 = all_cs[i], all_cs[j]
            d = float(np.linalg.norm(u.atoms[c1].position - u.atoms[c2].position))
            if d <= 1.65:
                ns = list(set(c_to_n[c1] + c_to_n[c2]))
                if len(ns) == 2:
                    cc_bridges.append((c1, c2, tuple(ns)))
    return cc_bridges


def order_ring_ns_cycle(ring_ns, cc_bridges):
    """
    Given ring N indices and C-C bridge tuples, return a cyclic order.
    """
    adj = {n: [] for n in ring_ns}
    for c1, c2, (a, b) in cc_bridges:
        if a in adj and b in adj:
            adj[a].append((b, c1, c2))
            adj[b].append((a, c1, c2))

    start = ring_ns[0]
    visited = [start]
    current = start
    prev = None
    while len(visited) < len(ring_ns):
        nbrs = [t for t in adj[current] if t[0] != prev]
        if not nbrs:
            break
        next_n, _, _ = nbrs[0]
        visited.append(next_n)
        prev = current
        current = next_n

    return visited


def trace_ring_torsions(u, ordered_ring_ns, n_to_c):
    """
    Build the 4 ring N-C-C-N torsion definitions.
    Returns list of (name, i, j, k, l, angle_deg).
    """
    cc_pairs = {}
    for i in range(len(ordered_ring_ns)):
        n1 = ordered_ring_ns[i]
        n2 = ordered_ring_ns[(i + 1) % len(ordered_ring_ns)]
        c1s = n_to_c[n1]
        c2s = n_to_c[n2]
        pair = None
        for ca in c1s:
            for cb in c2s:
                d = float(np.linalg.norm(u.atoms[ca].position - u.atoms[cb].position))
                if d <= 1.65:
                    pair = (ca, cb)
                    break
            if pair:
                break
        cc_pairs[(n1, n2)] = pair

    torsions = []
    for i in range(len(ordered_ring_ns)):
        n1 = ordered_ring_ns[i]
        n2 = ordered_ring_ns[(i + 1) % len(ordered_ring_ns)]
        pair = cc_pairs.get((n1, n2))
        if not pair:
            continue
        c1, c2 = pair
        c1_nbrs = get_neighbors(u.atoms, c1, 1.8)
        c2_nbrs = get_neighbors(u.atoms, c2, 1.8)
        c1_has_n1 = any(x[0] == n1 and x[1] == 'N' for x in c1_nbrs)
        c2_has_n2 = any(x[0] == n2 and x[1] == 'N' for x in c2_nbrs)
        if c1_has_n1 and c2_has_n2:
            atoms = (n1, c1, c2, n2)
        else:
            atoms = (n1, c2, c1, n2)
        val = dihedral(u.atoms.positions, atoms)
        torsions.append((f'T{i+1}', atoms[0], atoms[1], atoms[2], atoms[3], val))

    return torsions


def trace_chrom_torsion(u, ring_ns, chrom_ns, n_to_c):
    """
    Trace the chromophore N-C-C-N torsion from a ring N to the proximal
    chromophore N, via a C-C bridge.
    Returns ((name, i, j, k, l, angle_deg), prox_n, distal_n) or None.
    """
    proxies = []
    for cn in chrom_ns:
        for rn in ring_ns:
            path = bfs_shortest_path(u.atoms, rn, [cn], max_depth=6)
            if path:
                proxies.append((len(path) - 1, rn, cn, path))

    if not proxies:
        warn('No path found from ring N to chrom N')
        return None

    proxies.sort(key=lambda x: x[0])
    _, ring_prox, prox_n, path = proxies[0]

    distal_n = [n for n in chrom_ns if n != prox_n]
    if not distal_n:
        warn('Could not identify distal chrom N')
        return None
    distal_n = distal_n[0]

    i, j, k, l = path[0], path[1], path[2], path[3]
    val = dihedral(u.atoms.positions, (i, j, k, l))

    best = (i, j, k, l, val)
    for _, alt_rn, _, alt_path in proxies[:4]:
        if len(alt_path) == 4:
            av = dihedral(u.atoms.positions, tuple(alt_path))
            if av is not np.nan:
                if best[4] is np.nan or abs(abs(av) - 38.9) < abs(abs(best[4]) - 38.9):
                    best = (alt_path[0], alt_path[1], alt_path[2], alt_path[3], av)

    return (('Tc', best[0], best[1], best[2], best[3], best[4]), prox_n, distal_n)


def run(data_dir, out_dir):
    """Main analysis routine."""
    os.makedirs(out_dir, exist_ok=True)

    sap_path = os.path.join(data_dir, 'me_rrrD_sap', 'xtbopt.xyz')
    tsap_path = os.path.join(data_dir, 'me_rrrD_tsap', 'xtbopt.xyz')

    if not os.path.isfile(sap_path):
        raise FileNotFoundError(f'SAP structure not found: {sap_path}')
    if not os.path.isfile(tsap_path):
        raise FileNotFoundError(f'TSAP structure not found: {tsap_path}')

    print(f'Loading SAP  : {sap_path}')
    print(f'Loading TSAP : {tsap_path}')

    sap = load_structure(sap_path)
    tsap = load_structure(tsap_path)
    print(f'  SAP atoms  : {sap.atoms.n_atoms}')
    print(f'  TSAP atoms : {tsap.atoms.n_atoms}')

    # ================================================================
    # STEP 1: Find Eu
    # ================================================================
    print()
    print('=' * 60)
    print('STEP 1: Locate Eu centre')
    print('=' * 60)
    sap_metal = find_metal(sap)
    tsap_metal = find_metal(tsap)
    print(f'  SAP  Eu index : {sap_metal}')
    print(f'  TSAP Eu index : {tsap_metal}')

    # ================================================================
    # STEP 2: Find 9 closest O/N donors (TSAP)
    # ================================================================
    print()
    print('=' * 60)
    print('STEP 2: Find 9 closest O/N donors to Eu (TSAP)')
    print('=' * 60)
    tsap_donors = closest_donors(tsap, tsap_metal, n_donors=12, max_distance=6.0)
    print('  Top 12 donors:')
    for idx, elem, dist in tsap_donors:
        print(f'    idx={idx:3d}  {elem:1s}  dist={dist:.3f} Å')

    tsap_top9 = tsap_donors[:9]
    tsap_top9_idx = [x[0] for x in tsap_top9]
    elem_counts = {'O': 0, 'N': 0}
    for _, elem, _ in tsap_top9:
        elem_counts[elem] += 1
    print(f'\n  Top 9 composition: {elem_counts}')
    assert elem_counts['O'] == 3 and elem_counts['N'] == 6, \
        f"Expected 3 O + 6 N, got {elem_counts}"

    if len(tsap_donors) > 9:
        gap = tsap_donors[9][2] - tsap_donors[8][2]
        print(f'  Gap #9 -> #10 : {gap:.3f} Å')
        assert gap > 0.3, f"Expected gap > 0.3 Å, got {gap:.3f} Å"

    n_donors = [(idx, elem, dist) for idx, elem, dist in tsap_top9 if elem == 'N']
    o_donors = [(idx, elem, dist) for idx, elem, dist in tsap_top9 if elem == 'O']

    # ================================================================
    # STEP 3: Identify ring N's vs chromophore N's
    # ================================================================
    print()
    print('=' * 60)
    print('STEP 3: Distinguish ring N vs chromophore N')
    print('=' * 60)
    ring_ns, chrom_ns, n_to_c = identify_ring_nitrogens(
        tsap, [idx for idx, _, _ in n_donors])
    print(f'  Ring N indices      : {ring_ns}')
    print(f'  Chromophore N indices : {chrom_ns}')
    assert len(ring_ns) == 4, f"Expected 4 ring N, found {len(ring_ns)}"
    assert len(chrom_ns) == 2, f"Expected 2 chrom N, found {len(chrom_ns)}"

    # ================================================================
    # STEP 4: Order ring N's into a cycle
    # ================================================================
    print()
    print('=' * 60)
    print('STEP 4: Order ring N cycle + C-C bridges')
    print('=' * 60)
    cc_bridges = build_cc_bridges(tsap, ring_ns, n_to_c)
    ordered_ring = order_ring_ns_cycle(ring_ns, cc_bridges)
    print(f'  Ring N order : {ordered_ring}')
    assert len(ordered_ring) == 4, f"Could not trace full ring cycle, got {ordered_ring}"

    # ================================================================
    # STEP 5: Build ring torsion definitions
    # ================================================================
    print()
    print('=' * 60)
    print('STEP 5: Ring N-C-C-N torsion angles (TSAP)')
    print('=' * 60)
    tsap_ring_torsions = trace_ring_torsions(tsap, ordered_ring, n_to_c)
    for name, i, j, k, l, val in tsap_ring_torsions:
        print(f'  {name} N{i}-C{j}-C{k}-N{l} : {val:7.2f}°')

    print()
    print('SAP ring torsion reference:')
    for name, i, j, k, l in SAP_RING_TORS:
        val = dihedral(sap.atoms.positions, (i, j, k, l))
        print(f'  {name} : {val:7.2f}°')

    # ================================================================
    # STEP 6: Chromophore torsion
    # ================================================================
    print()
    print('=' * 60)
    print('STEP 6: Chromophore torsion')
    print('=' * 60)

    chrom_out = trace_chrom_torsion(tsap, ring_ns, chrom_ns, n_to_c)
    if chrom_out is None:
        raise RuntimeError("Failed to trace chromophore torsion")
    (tc_name, tc_i, tc_j, tc_k, tc_l, tc_val), prox_n, distal_n = chrom_out

    print(f'  Proximal chrom N : {prox_n}')
    print(f'  Distal chrom N   : {distal_n}')
    print(f'  Chrom torsion {tc_name} N{tc_i}-C{tc_j}-C{tc_k}-N{tc_l} : {tc_val:7.2f}°')

    _, si, sj, sk, sl = SAP_CHROM_TOR
    sap_tc_val = dihedral(sap.atoms.positions, (si, sj, sk, sl))
    print(f'  SAP chrom torsion {SAP_CHROM_TOR[0]} : {sap_tc_val:7.2f}°')

    # ================================================================
    # STEP 7: Classification sanity check
    # ================================================================
    print()
    print('=' * 60)
    print('STEP 7: SAP / TSAP classification sanity check')
    print('=' * 60)
    tsap_mean_ring = float(np.mean([val for _, _, _, _, _, val in tsap_ring_torsions]))
    print(f'  TSAP mean ring torsion : {tsap_mean_ring:.2f}°')
    print(f'  TSAP chrom torsion     : {tc_val:.2f}°')
    print(f'  SAP  chrom torsion     : {sap_tc_val:.2f}°')

    if abs(tsap_mean_ring) < 10 or abs(tc_val) < 10:
        print('  -> Classification: UNK (too close to 0°)')
    elif np.sign(tsap_mean_ring) == np.sign(tc_val):
        print('  -> Classification: TSAP (same sign)')
    else:
        print('  -> Classification: SAP (opposite signs)')

    # ================================================================
    # STEP 8: Build coordinate arrays following SAP convention
    # ================================================================
    print()
    print('=' * 60)
    print('STEP 8: Build index tables')
    print('=' * 60)

    o_indices = [idx for idx, _, _ in o_donors]
    assert sorted(o_indices) == [0, 1, 2], f"Unexpected O indices: {o_indices}"

    tsap_coord = o_indices + ordered_ring + [prox_n, distal_n]
    tsap_bottom_idx = list(range(3, 7))
    tsap_top_idx = [0, 1, 2, 7]
    tsap_cap_idx = [8]

    print(f'  TSAP COORD_IDX : {tsap_coord}')
    print(f'  TSAP BOTTOM    : {tsap_bottom_idx} -> {[tsap_coord[i] for i in tsap_bottom_idx]}')
    print(f'  TSAP TOP       : {tsap_top_idx}    -> {[tsap_coord[i] for i in tsap_top_idx]}')
    print(f'  TSAP CAP       : {tsap_cap_idx}     -> {[tsap_coord[i] for i in tsap_cap_idx]}')

    # ================================================================
    # STEP 9: Package and save
    # ================================================================
    print()
    print('=' * 60)
    print('STEP 9: Save outputs')
    print('=' * 60)

    tsap_ring_dicts = [
        {'name': name, 'atoms': [int(i), int(j), int(k), int(l)], 'angle_deg': float(round(val, 2))}
        for name, i, j, k, l, val in tsap_ring_torsions
    ]
    tsap_chrom_dict = {
        'name': tc_name,
        'atoms': [int(tc_i), int(tc_j), int(tc_k), int(tc_l)],
        'angle_deg': float(round(tc_val, 2)),
    }

    data = {
        'sap': {
            'metal': SAP_METAL_IDX,
            'coord': SAP_COORD_IDX,
            'eu8_coord': SAP_EU8_COORD_IDX,
            'eu8_align': SAP_EU8_ALIGN_IDX,
            'bottom': SAP_BOTTOM,
            'top': SAP_TOP,
            'cap': SAP_CAP,
            'ring_torsions': [
                {'name': 'T1', 'atoms': list(SAP_RING_TORS[0][1:])},
                {'name': 'T2', 'atoms': list(SAP_RING_TORS[1][1:])},
                {'name': 'T3', 'atoms': list(SAP_RING_TORS[2][1:])},
                {'name': 'T4', 'atoms': list(SAP_RING_TORS[3][1:])},
            ],
            'chrom_torsion': {'name': 'Tc', 'atoms': list(SAP_CHROM_TOR[1:])},
        },
        'tsap': {
            'metal': tsap_metal,
            'coord': tsap_coord,
            'bottom': tsap_bottom_idx,
            'top': tsap_top_idx,
            'cap': tsap_cap_idx,
            'ring_torsions': tsap_ring_dicts,
            'chrom_torsion': tsap_chrom_dict,
        },
    }

    def _native(x):
        """Recursively convert numpy scalar types to native Python types."""
        import numpy as np
        if isinstance(x, (np.integer, np.floating)):
            return x.item()
        if isinstance(x, (list, tuple)):
            return [_native(i) for i in x]
        if isinstance(x, dict):
            return {k: _native(v) for k, v in x.items()}
        return x

    data = _native(data)
    json_path = os.path.join(out_dir, 'indices.json')
    with open(json_path, 'w') as f:
        json.dump(data, f, indent=2)
    print(f'  Saved JSON: {json_path}')

    summary_path = os.path.join(out_dir, 'indices_summary.txt')
    with open(summary_path, 'w') as f:
        f.write('=' * 60 + '\n')
        f.write('Atom Index Mapping for SAP and TSAP Conformers\n')
        f.write('Generated by T1_find_indices.py\n')
        f.write('=' * 60 + '\n\n')

        # SAP table
        f.write('--- SAP (me_rrrD_sap) ---\n')
        f.write(f'Metal (Eu)      : {SAP_METAL_IDX}\n')
        f.write(f'Coord atoms     : {SAP_COORD_IDX}\n')
        f.write(f'  O atoms       : {SAP_TOP[:3]}\n')
        f.write(f'  Ring N atoms  : {SAP_BOTTOM}\n')
        f.write(f'  Proximal N    : {SAP_TOP[3]}\n')
        f.write(f'  Distal N      : {SAP_CAP[0]}\n')
        f.write('\nRing torsions:\n')
        for name, i, j, k, l in SAP_RING_TORS:
            f.write(f'  {name} : N{i}-C{j}-C{k}-N{l} = {dihedral(sap.atoms.positions, (i,j,k,l)):.2f}°\n')
        f.write(f'  {SAP_CHROM_TOR[0]} : N{si}-C{sj}-C{sk}-N{sl} = {sap_tc_val:.2f}°\n')

        f.write('\n\n--- TSAP (me_rrrD_tsap) ---\n')
        f.write(f'Metal (Eu)      : {tsap_metal}\n')
        f.write(f'Coord atoms     : {tsap_coord}\n')
        f.write(f'  O atoms       : {[tsap_coord[i] for i in tsap_top_idx[:3]]}\n')
        f.write(f'  Ring N atoms  : {[tsap_coord[i] for i in tsap_bottom_idx]}\n')
        f.write(f'  Proximal N    : {prox_n}\n')
        f.write(f'  Distal N      : {distal_n}\n')
        f.write('\nRing torsions:\n')
        for name, i, j, k, l, val in tsap_ring_torsions:
            f.write(f'  {name} : N{i}-C{j}-C{k}-N{l} = {val:.2f}°\n')
        f.write(f'  {tc_name} : N{tc_i}-C{tc_j}-C{tc_k}-N{tc_l} = {tc_val:.2f}°\n')

        f.write('\nDistances from Eu to coordination atoms:\n')
        for pos, abs_idx in enumerate(tsap_coord):
            a = tsap.atoms[abs_idx]
            d = float(np.linalg.norm(a.position - tsap.atoms[tsap_metal].position))
            role = ''
            if pos in tsap_top_idx:
                role = ' (top)'
            elif pos in tsap_bottom_idx:
                role = ' (bottom)'
            elif pos in tsap_cap_idx:
                role = ' (cap)'
            f.write(f'  {abs_idx:3d} {a.name.strip():1s}  {d:.3f} Å{role}\n')

    print(f'  Saved summary: {summary_path}')

    # Validate distances
    print()
    print('Validation: coordination distances')
    for idx in tsap_coord:
        d = float(np.linalg.norm(tsap.atoms[idx].position - tsap.atoms[tsap_metal].position))
        assert 1.5 < d < 3.0, f'Distance {d:.2f} Å for atom {idx} out of range 1.5-3.0'
    print('  All 9 coordination distances within 1.5-3.0 Å  PASS')

    return data


def main():
    parser = argparse.ArgumentParser(
        description='Find TSAP atom indices for coordination geometry analysis.')
    parser.add_argument('--data-dir', default='data',
                        help='Directory containing system subfolders (e.g. data/)')
    parser.add_argument('--out-dir', default='analysis',
                        help='Directory for output files (e.g. analysis/)')
    args = parser.parse_args()

    run(args.data_dir, args.out_dir)


if __name__ == '__main__':
    main()
