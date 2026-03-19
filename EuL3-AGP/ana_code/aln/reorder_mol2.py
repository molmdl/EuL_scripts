"""
Reorder atoms in sss-delta-EuL2-D.mol2 to match the atom ordering of optfreq.mol2.
Uses graph isomorphism (VF2 algorithm via networkx) to find the correspondence.
Never overwrites existing files.
"""

import networkx as nx
from networkx.algorithms import isomorphism
import sys

def parse_mol2(filename):
    """Parse a mol2 file, return atoms list and bonds list."""
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
                    # idx, name, x, y, z, type, ...
                    idx = int(parts[0])
                    name = parts[1]
                    x, y, z = float(parts[2]), float(parts[3]), float(parts[4])
                    atom_type = parts[5]
                    residue_num = parts[6] if len(parts) > 6 else '1'
                    residue_name = parts[7] if len(parts) > 7 else 'UNL1'
                    charge = parts[8] if len(parts) > 8 else '0.0000'
                    atoms.append({
                        'idx': idx, 'name': name,
                        'x': x, 'y': y, 'z': z,
                        'type': atom_type,
                        'element': atom_type.split('.')[0].upper(),
                        'res_num': residue_num,
                        'res_name': residue_name,
                        'charge': charge,
                        'line': line
                    })
            elif section == '@<TRIPOS>BOND':
                parts = line.split()
                if len(parts) >= 4:
                    bond_idx = int(parts[0])
                    a1 = int(parts[1])
                    a2 = int(parts[2])
                    btype = parts[3]
                    bonds.append({'idx': bond_idx, 'a1': a1, 'a2': a2, 'type': btype})
    return atoms, bonds

def build_graph(atoms, bonds):
    """Build a networkx graph with element as node attribute."""
    G = nx.Graph()
    for a in atoms:
        G.add_node(a['idx'], element=a['element'])
    for b in bonds:
        G.add_edge(b['a1'], b['a2'])
    return G

def node_match(n1, n2):
    return n1['element'] == n2['element']

def find_isomorphism(G_ref, G_src):
    """Find a graph isomorphism mapping from G_src node ids to G_ref node ids."""
    gm = isomorphism.GraphMatcher(G_ref, G_src, node_match=node_match)
    if gm.is_isomorphic():
        # iso_map: ref_node -> src_node
        iso = gm.isomorphisms_iter()
        mapping = next(iso)  # ref -> src
        return mapping
    return None

def reorder_and_write(ref_atoms, ref_bonds, src_atoms, src_bonds, mapping, outfile):
    """
    mapping: dict { ref_atom_idx -> src_atom_idx }
    Write a new mol2 where for each position i in ref, we output the src atom
    that corresponds to ref atom i, but using src coordinates/name etc.
    Bond table is rebuilt from ref bond table, with new atom indices.
    """
    # Build src lookup by idx
    src_by_idx = {a['idx']: a for a in src_atoms}

    # For each ref atom (in ref order), get the matching src atom
    reordered = []
    for ref_a in ref_atoms:
        src_idx = mapping[ref_a['idx']]
        src_a = src_by_idx[src_idx]
        reordered.append(src_a)

    # Build a map: new_position (1-based) for each src atom index
    # new_pos[src_idx] = new 1-based index (same as ref order position)
    src_idx_to_new_pos = {}
    for new_pos_0, ref_a in enumerate(ref_atoms):
        src_idx = mapping[ref_a['idx']]
        src_idx_to_new_pos[src_idx] = new_pos_0 + 1

    # Write output
    with open(outfile, 'w') as f:
        f.write('@<TRIPOS>MOLECULE\n')
        f.write('*****\n')
        n_atoms = len(src_atoms)
        n_bonds = len(src_bonds)
        f.write(f' {n_atoms} {n_bonds} 0 0 0\n')
        f.write('SMALL\n')
        f.write('GASTEIGER\n')
        f.write('\n')
        f.write('@<TRIPOS>ATOM\n')
        for new_i, src_a in enumerate(reordered, 1):
            # Use the ref atom's name (to preserve naming convention of optfreq)
            ref_a = ref_atoms[new_i - 1]
            f.write(f'    {new_i:>4d} {ref_a["name"]:<12s}'
                    f'{src_a["x"]:>10.4f}{src_a["y"]:>10.4f}{src_a["z"]:>10.4f} '
                    f'{ref_a["type"]:<8s} {ref_a["res_num"]}  {ref_a["res_name"]:<12s}'
                    f'{ref_a["charge"]}\n')
        f.write('@<TRIPOS>BOND\n')
        # Use ref bond table but atom indices stay the same (they are already ref-order)
        for b in ref_bonds:
            f.write(f'   {b["idx"]:>4d} {b["a1"]:>4d} {b["a2"]:>4d}    {b["type"]}\n')

    print(f"Written {n_atoms} atoms and {n_bonds} bonds to {outfile}")

# ---- Main ----
ref_file = 'optfreq.mol2'
src_file = 'sss-delta-EuL2-D.mol2'
out_file = 'sss-delta-EuL2-D_reordered.mol2'

import os
if os.path.exists(out_file):
    print(f"ERROR: {out_file} already exists! Refusing to overwrite.")
    sys.exit(1)

print(f"Parsing {ref_file}...")
ref_atoms, ref_bonds = parse_mol2(ref_file)
print(f"  {len(ref_atoms)} atoms, {len(ref_bonds)} bonds")

print(f"Parsing {src_file}...")
src_atoms, src_bonds = parse_mol2(src_file)
print(f"  {len(src_atoms)} atoms, {len(src_bonds)} bonds")

print("Building graphs...")
G_ref = build_graph(ref_atoms, ref_bonds)
G_src = build_graph(src_atoms, src_bonds)

print("Running graph isomorphism (VF2)...")
mapping = find_isomorphism(G_ref, G_src)

if mapping is None:
    print("ERROR: No graph isomorphism found! The molecules may not be the same connectivity.")
    sys.exit(1)

print(f"Isomorphism found. Writing reordered file...")
reorder_and_write(ref_atoms, ref_bonds, src_atoms, src_bonds, mapping, out_file)
print("Done.")
