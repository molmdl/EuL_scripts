"""
Reorder atoms in sss-delta-EuL2-D.mol2 to match the atom ordering of optfreq.mol2.
Uses graph isomorphism (VF2) on element-labelled graphs.
Atom names, types, residue labels, and bond table come from optfreq.mol2 (reference).
Coordinates come from sss-delta-EuL2-D.mol2 (source).
Output: sss-delta-EuL2-D_reordered.mol2  (never overwrites)
"""

import networkx as nx
from networkx.algorithms import isomorphism
import os, sys

def parse_mol2(filename):
    atoms = []
    bonds = []
    section = None
    with open(filename) as f:
        for line in f:
            line = line.rstrip('\n')
            if line.startswith('@<TRIPOS>'):
                section = line.strip()
                continue
            if section == '@<TRIPOS>ATOM':
                parts = line.split()
                if len(parts) >= 6:
                    atype = parts[5]
                    el = atype.split('.')[0].upper()
                    if el.startswith('EU'):
                        el = 'EU'
                    atoms.append({
                        'idx':      int(parts[0]),
                        'name':     parts[1],
                        'x':        float(parts[2]),
                        'y':        float(parts[3]),
                        'z':        float(parts[4]),
                        'type':     atype,
                        'element':  el,
                        'res_num':  parts[6] if len(parts) > 6 else '1',
                        'res_name': parts[7] if len(parts) > 7 else 'UNL1',
                        'charge':   parts[8] if len(parts) > 8 else '0.0000',
                    })
            elif section == '@<TRIPOS>BOND':
                parts = line.split()
                if len(parts) >= 4:
                    bonds.append({
                        'idx': int(parts[0]),
                        'a1':  int(parts[1]),
                        'a2':  int(parts[2]),
                        'type': parts[3]
                    })
    return atoms, bonds

def build_graph(atoms, bonds):
    G = nx.Graph()
    for a in atoms:
        G.add_node(a['idx'], element=a['element'])
    for b in bonds:
        G.add_edge(b['a1'], b['a2'])
    return G

def node_match(n1, n2):
    return n1['element'] == n2['element']

# ---- Paths ----
ref_file = 'optfreq.mol2'
src_file = 'sss-delta-EuL2-D.mol2'
out_file = 'sss-delta-EuL2-D_reordered.mol2'

if os.path.exists(out_file):
    print(f"ERROR: {out_file} already exists! Refusing to overwrite.")
    sys.exit(1)

print(f"Parsing {ref_file}...")
ref_atoms, ref_bonds = parse_mol2(ref_file)
print(f"  {len(ref_atoms)} atoms, {len(ref_bonds)} bonds")

print(f"Parsing {src_file}...")
src_atoms, src_bonds = parse_mol2(src_file)
print(f"  {len(src_atoms)} atoms, {len(src_bonds)} bonds")

print("Building molecular graphs...")
G_ref = build_graph(ref_atoms, ref_bonds)
G_src = build_graph(src_atoms, src_bonds)

print("Running VF2 graph isomorphism (element-matched)...")
gm = isomorphism.GraphMatcher(G_ref, G_src, node_match=node_match)
if not gm.is_isomorphic():
    print("ERROR: No graph isomorphism found!")
    sys.exit(1)

# mapping: ref_atom_idx -> src_atom_idx
mapping = next(gm.isomorphisms_iter())
print(f"  Isomorphism found ({len(mapping)} atom correspondences).")

# Build src lookup
src_by_idx = {a['idx']: a for a in src_atoms}

# Validation: all ref atoms have a unique corresponding src atom
src_targets = list(mapping.values())
assert len(set(src_targets)) == len(src_targets), "Mapping is not 1-to-1!"

# Write output: atom names/types/bonds from ref, coordinates from src
print(f"Writing {out_file}...")
with open(out_file, 'w') as f:
    f.write('@<TRIPOS>MOLECULE\n')
    f.write('*****\n')
    f.write(f' {len(ref_atoms)} {len(ref_bonds)} 0 0 0\n')
    f.write('SMALL\n')
    f.write('GASTEIGER\n')
    f.write('\n')
    f.write('@<TRIPOS>ATOM\n')
    for new_pos, ref_a in enumerate(ref_atoms, 1):
        src_a = src_by_idx[mapping[ref_a['idx']]]
        f.write(
            f'{new_pos:>7d} '
            f'{ref_a["name"]:<12s}'
            f'{src_a["x"]:>9.4f} '
            f'{src_a["y"]:>9.4f} '
            f'{src_a["z"]:>9.4f} '
            f'{ref_a["type"]:<8s}'
            f'{ref_a["res_num"]:>2s}  '
            f'{ref_a["res_name"]:<12s}'
            f'{ref_a["charge"]}\n'
        )
    f.write('@<TRIPOS>BOND\n')
    for b in ref_bonds:
        f.write(f'{b["idx"]:>6d} {b["a1"]:>5d} {b["a2"]:>5d}    {b["type"]}\n')

print(f"Done. Output saved to: {out_file}")
