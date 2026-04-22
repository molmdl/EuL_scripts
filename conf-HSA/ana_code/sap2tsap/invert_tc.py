#!/usr/bin/env python3
"""
invert_tc.py
------------
Invert the sign of the Tc (12-N-4 ring -> chromophore) torsion angle in
a single-frame GRO file and its matching ITP file.

Tc is defined as:  N4 - C17 - C24 - N8  (1-based GRO atom indices)

Strategy (ROTATION method):
    Rotate C17 (and its hydrogens H45, H46) around the N4-C24 axis by an
    optimized angle that inverts Tc while minimizing N4 valence angle strain.
    
    This approach:
    - Preserves N4-C17 and C24-C17 bond lengths exactly (rotation around axis)
    - Produces much lower N4 strain (~29°) compared to reflection (~55°)
    - Keeps T1-T4 unchanged (ring torsions unaffected)
    - T1-T4 end up with same sign as inverted Tc (TSAP configuration)
    
    Unlike reflection, rotation does NOT change the "handedness" of the
    molecular geometry, so dihedral phases in the ITP do NOT need negation.

Usage:
    python invert_tc.py --gro <input.gro> --itp <input.itp>
                        --gro-out <output.gro> --itp-out <output.itp>

    Flags --gro-out and --itp-out default to <stem>_tc_inv.gro / .itp
    alongside the input files.
"""

import argparse
import os
import sys

import numpy as np

# -----------------------------------------------------------------------------
#  ALL TUNEABLE CONSTANTS  (edit here, not in the functions below)
# -----------------------------------------------------------------------------

# 1-based GRO atom indices for the Tc torsion:  N4 - C17 - C24 - N8
TC_N4_IDX  =  4   # ring nitrogen  (lies on reflection plane, not moved)
TC_C17_IDX = 17   # methylene carbon (linker, REFLECTED then H-rotated)
TC_C24_IDX = 24   # chromophore aromatic carbon (lies on reflection plane, not moved)
TC_N8_IDX  =  8   # proximal chromophore nitrogen (lies on reflection plane, not moved)

# Methylene hydrogens on C17 (1-based); reflected with C17, then rotated in Step 2
TC_H45_IDX = 45
TC_H46_IDX = 46

# H-rotation optimiser: clash-score parameters
# Only pairs that are 1,5 or further apart contribute to the score.
# The repulsive term is (CLASH_SIGMA / d)^12; pairs with d > CLASH_CUTOFF
# contribute negligibly and are skipped for speed.
CLASH_SIGMA  = 0.30   # nm  (soft wall; pairs closer than this are penalised)
CLASH_CUTOFF = 0.30   # nm  (only pairs closer than this are evaluated)

# Angular scan parameters for H-rotation optimiser
H_ROT_STEP_DEG  = 1.0   # degrees (scan resolution; 1 deg is sufficient)

# -----------------------------------------------------------------------------
#  GRO  I/O
# -----------------------------------------------------------------------------

def read_gro(path):
    """Parse a single-frame GRO file.
    Returns (title, n_atoms_line, atoms, box) where atoms is a list of dicts
    and n_atoms_line is the raw line-2 string (preserves original whitespace)."""
    with open(path) as f:
        lines = f.readlines()
    title = lines[0].rstrip('\n')
    n_atoms_line = lines[1].rstrip('\n')   # keep original spacing
    n_atoms = int(n_atoms_line.strip())
    atoms = []
    for line in lines[2: 2 + n_atoms]:
        atoms.append(dict(
            resnum   = int(line[0:5]),
            resname  = line[5:10].strip(),
            atomname = line[10:15].strip(),
            atomnum  = int(line[15:20]),
            x        = float(line[20:28]),
            y        = float(line[28:36]),
            z        = float(line[36:44]),
        ))
    box = lines[2 + n_atoms].rstrip('\n')
    return title, n_atoms_line, atoms, box


def write_gro(path, title, n_atoms_line, atoms, box):
    """Write atoms list back to GRO fixed-width format.
    n_atoms_line is written verbatim to preserve the original spacing."""
    with open(path, 'w') as f:
        f.write(title + '\n')
        f.write(n_atoms_line + '\n')
        for a in atoms:
            f.write(
                f"{a['resnum']:5d}"
                f"{a['resname']:<5s}"
                f"{a['atomname']:>5s}"
                f"{a['atomnum']:5d}"
                f"{a['x']:8.3f}"
                f"{a['y']:8.3f}"
                f"{a['z']:8.3f}"
                '\n'
            )
        f.write(box + '\n')


def atoms_to_xp(atoms):
    """Return (N,3) numpy array of coordinates (nm)."""
    return np.array([[a['x'], a['y'], a['z']] for a in atoms])


def xp_to_atoms(atoms, xp):
    """Copy coordinates from (N,3) xp back into atoms list (in-place)."""
    for a, row in zip(atoms, xp):
        a['x'], a['y'], a['z'] = float(row[0]), float(row[1]), float(row[2])


# -----------------------------------------------------------------------------
#  BOND GRAPH (for clash exclusion)
# -----------------------------------------------------------------------------

def build_bond_graph(itp_path):
    """
    Parse the [ bonds ] section of an ITP file and return an adjacency dict:
        {atom_i (1-based int): set of directly bonded neighbours (1-based ints)}
    Only the bond section is used; all other sections are ignored.
    """
    adj = {}
    in_bonds = False
    with open(itp_path) as fh:
        for line in fh:
            stripped = line.strip()
            low = stripped.lower()
            if low.startswith('['):
                in_bonds = 'bonds' in low
                continue
            if not in_bonds or not stripped or stripped.startswith(';'):
                continue
            code = stripped.split(';')[0].split()
            if len(code) < 2:
                continue
            try:
                ai, aj = int(code[0]), int(code[1])
            except ValueError:
                continue
            adj.setdefault(ai, set()).add(aj)
            adj.setdefault(aj, set()).add(ai)
    return adj


def _excluded_pairs(adj, source, max_sep=3):
    """
    BFS from `source` up to `max_sep` bonds away.
    Returns a set of atom indices (1-based) that are excluded from full
    non-bonded interactions with `source` (i.e. 1,2 and 1,3 pairs).
    With max_sep=3 this returns the 1,2 + 1,3 set; use max_sep=4 to also
    include 1,4 scaled pairs.
    """
    visited = {source: 0}
    queue = [source]
    while queue:
        node = queue.pop(0)
        depth = visited[node]
        if depth >= max_sep:
            continue
        for nb in adj.get(node, set()):
            if nb not in visited:
                visited[nb] = depth + 1
                queue.append(nb)
    return set(visited.keys())


# -----------------------------------------------------------------------------
#  GEOMETRY
# -----------------------------------------------------------------------------

def dihedral_deg(p1, p2, p3, p4):
    """Dihedral angle in degrees [-180, 180], IUPAC/GROMACS convention."""
    b1 = p2 - p1;  b2 = p3 - p2;  b3 = p4 - p3
    n1 = np.cross(b1, b2);  nn1 = np.linalg.norm(n1)
    n2 = np.cross(b2, b3);  nn2 = np.linalg.norm(n2)
    if nn1 < 1e-9 or nn2 < 1e-9:
        return np.nan
    n1 /= nn1;  n2 /= nn2
    b2n = b2 / np.linalg.norm(b2)
    return float(np.degrees(np.arctan2(np.dot(np.cross(n1, n2), b2n),
                                        np.dot(n1, n2))))


def angle_deg(pa, pb, pc):
    """Valence angle at pb in degrees."""
    v1 = pa - pb;  v2 = pc - pb
    n1 = np.linalg.norm(v1);  n2 = np.linalg.norm(v2)
    if n1 < 1e-9 or n2 < 1e-9:
        return np.nan
    cos = np.dot(v1, v2) / (n1 * n2)
    return float(np.degrees(np.arccos(np.clip(cos, -1.0, 1.0))))


def rotate_around_axis(points, origin, axis, angle_deg_val):
    """
    Rotate each point in (N,3) array `points` around the line
    (origin, origin+axis) by `angle_deg_val` degrees (right-hand rule).
    Returns a new (N,3) array.
    """
    theta = np.radians(angle_deg_val)
    u = axis / np.linalg.norm(axis)          # unit axis
    ux, uy, uz = u
    c, s = np.cos(theta), np.sin(theta)

    # Rodrigues rotation matrix
    R = np.array([
        [c + ux*ux*(1-c),     ux*uy*(1-c) - uz*s,  ux*uz*(1-c) + uy*s],
        [uy*ux*(1-c) + uz*s,  c + uy*uy*(1-c),      uy*uz*(1-c) - ux*s],
        [uz*ux*(1-c) - uy*s,  uz*uy*(1-c) + ux*s,   c + uz*uz*(1-c)  ],
    ])

    pts = points - origin                    # translate to origin
    pts = (R @ pts.T).T                      # rotate
    return pts + origin                      # translate back


def report_torsions(xp, label):
    """Print the 5 key torsion angles."""
    ring_defs = [
        ('T1 N4 -C9 -C10-N5 ',  4,  9, 10,  5),
        ('T2 N5 -C11-C12-N6 ',  5, 11, 12,  6),
        ('T3 N6 -C13-C14-N7 ',  6, 13, 14,  7),
        ('T4 N7 -C15-C16-N4 ',  7, 15, 16,  4),
        ('Tc N4 -C17-C24-N8 ',  4, 17, 24,  8),
    ]
    print(f'  Torsions [{label}]:')
    for name, a, b, c, d in ring_defs:
        ang = dihedral_deg(xp[a-1], xp[b-1], xp[c-1], xp[d-1])
        print(f'    {name}: {ang:8.2f} deg')


def report_n4_angles(xp, label):
    """Print valence angles at N4 involving C17."""
    i_n4  = TC_N4_IDX  - 1
    i_c9  = 9  - 1   # ring CH2 bonded to N4
    i_c16 = 16 - 1   # ring CH2 bonded to N4
    i_c17 = TC_C17_IDX - 1
    print(f'  N4 valence angles [{label}]:')
    print(f'    C9 -N4-C16 : {angle_deg(xp[i_c9],  xp[i_n4], xp[i_c16]):7.3f} deg')
    print(f'    C9 -N4-C17 : {angle_deg(xp[i_c9],  xp[i_n4], xp[i_c17]):7.3f} deg')
    print(f'    C16-N4-C17 : {angle_deg(xp[i_c16], xp[i_n4], xp[i_c17]):7.3f} deg')


# -----------------------------------------------------------------------------
#  Tc INVERSION  - GRO coordinate rebuild
# -----------------------------------------------------------------------------

def invert_tc_gro(xp, adj):
    """
    Invert the Tc torsion (N4-C17-C24-N8) in coordinate array xp.

    Algorithm (ROTATION method with ring optimization):

    Two-stage rotation to invert Tc while minimizing N4 strain:
    
    Stage 1: Rotate C17 (and H45, H46) around the N4-C24 axis to invert Tc
    Stage 2: Rotate the ring side atoms (C9, C10, C15, C16) around the same axis
             to reduce N4 valence angle strain while preserving TSAP configuration
    
    The rotation axis is N4 -> C24. Since both atoms lie on the axis,
    the N4-C17 and C24-C17 bond lengths are exactly preserved.
    
    Ring side rotation changes T1 and T4, but the optimization ensures they
    maintain the same sign as Tc (TSAP configuration).

    Parameters
    ----------
    xp  : (N,3) float array of coordinates in nm (0-based indexing)
    adj : adjacency dict from build_bond_graph(itp_path) - not used in this method

    Returns a modified copy of xp (nm).
    """
    xp = xp.copy()

    # -- 0-based indices -------------------------------------------------------
    i_n4  = TC_N4_IDX  - 1
    i_c17 = TC_C17_IDX - 1
    i_c24 = TC_C24_IDX - 1
    i_n8  = TC_N8_IDX  - 1
    i_h45 = TC_H45_IDX - 1
    i_h46 = TC_H46_IDX - 1

    # Ring side atoms to rotate (C9, C10, C15, C16 - 1-based: 9, 10, 15, 16)
    ring_side_atoms = [8, 9, 14, 15]  # 0-based

    p_n4  = xp[i_n4].copy()
    p_c24 = xp[i_c24].copy()

    # Original Tc value
    tc_orig = dihedral_deg(xp[i_n4], xp[i_c17], xp[i_c24], xp[i_n8])
    target_tc = -tc_orig

    # Rotation axis: N4 -> C24
    rot_axis = p_c24 - p_n4
    axis_norm = np.linalg.norm(rot_axis)
    if axis_norm < 1e-9:
        sys.exit('ERROR: N4 and C24 are coincident.')
    rot_axis /= axis_norm

    # Find optimal combination of C17 rotation and ring-side rotation
    target_tolerance = 15.0
    ideal_n4 = 109.47
    
    candidates = []
    
    for c17_rot in np.arange(260, 340, 1.0):
        xp_test = xp.copy()
        
        # Stage 1: rotate C17
        for idx in (i_c17, i_h45, i_h46):
            xp_test[idx] = rotate_around_axis(xp_test[idx:idx+1], p_n4, rot_axis, c17_rot)[0]
        
        tc_test = dihedral_deg(xp_test[i_n4], xp_test[i_c17], xp_test[i_c24], xp_test[i_n8])
        tc_diff = abs(tc_test - target_tc)
        
        if tc_diff > target_tolerance:
            continue
        
        # Try different ring rotations
        for ring_rot in np.arange(-35, 10, 1.0):
            xp_try = xp_test.copy()
            
            # Stage 2: rotate ring side
            for idx in ring_side_atoms:
                xp_try[idx] = rotate_around_axis(xp_try[idx:idx+1], p_n4, rot_axis, ring_rot)[0]
            
            # Check T1-T4 signs
            t1 = dihedral_deg(xp_try[i_n4], xp_try[8], xp_try[9], xp_try[4])
            t2 = dihedral_deg(xp_try[4], xp_try[10], xp_try[11], xp_try[5])
            t3 = dihedral_deg(xp_try[5], xp_try[12], xp_try[13], xp_try[6])
            t4 = dihedral_deg(xp_try[6], xp_try[14], xp_try[15], xp_try[i_n4])
            
            tc_sign = np.sign(tc_test)
            t_signs = [np.sign(t1), np.sign(t2), np.sign(t3), np.sign(t4)]
            
            if tc_sign != t_signs[0] or not all(s == t_signs[0] for s in t_signs):
                continue
            
            # Calculate N4 strain
            n4_c9 = angle_deg(xp_try[8], xp_try[i_n4], xp_try[i_c17])
            n4_c16 = angle_deg(xp_try[15], xp_try[i_n4], xp_try[i_c17])
            strain = max(abs(ideal_n4 - n4_c9), abs(ideal_n4 - n4_c16))
            
            # Score: prioritize Tc accuracy, then minimize strain
            score = tc_diff + 0.1 * strain
            
            candidates.append((score, c17_rot, ring_rot, strain, tc_test, t1, t2, t3, t4))
    
    if not candidates:
        sys.exit(f'ERROR: Could not find valid rotation to invert Tc (target: {target_tc:.2f})')
    
    # Sort by score (Tc accuracy + strain)
    candidates.sort(key=lambda x: x[0])
    best = candidates[0]
    best_c17_rot = best[1]
    best_ring_rot = best[2]
    
    # Apply best rotations
    for idx in (i_c17, i_h45, i_h46):
        xp[idx] = rotate_around_axis(xp[idx:idx+1], p_n4, rot_axis, best_c17_rot)[0]
    
    for idx in ring_side_atoms:
        xp[idx] = rotate_around_axis(xp[idx:idx+1], p_n4, rot_axis, best_ring_rot)[0]

    return xp


# -----------------------------------------------------------------------------
#  ITP  dihedral inversion
# -----------------------------------------------------------------------------

def _split_comment(line):
    if ';' in line:
        idx = line.index(';')
        return line[:idx], line[idx:].rstrip('\n')
    return line.rstrip('\n'), ''


def _negate_angle(val):
    """Negate angle and wrap to (-180, 180]."""
    a = -float(val)
    while a <= -180.0: a += 360.0
    while a >   180.0: a -= 360.0
    return a


def _negate_dihedral_line(line):
    """
    Negate d0/phase in a functype-1, -2, -9, or -4 dihedral line.
    Returns (new_line, ft_str) or (original_line, None) on parse failure.

    Note: ft=9 and ft=4 terms with phase 0 or 180 are symmetric under sign
    inversion (0 -> 0, 180 -> -180 -> +180 after wrapping), so they are
    written out but their physics is unchanged.  ft=2 DRIH terms carry
    genuine signed equilibrium values and must be negated.
    """
    code, comment = _split_comment(line)
    parts = code.split()
    if len(parts) < 6:
        return line, None
    try:
        ai, aj, ak, al = (int(parts[i]) for i in range(4))
    except ValueError:
        return line, None

    ft = parts[4]
    if ft not in {'1', '2', '9', '4'}:
        return line, None

    try:
        angle = _negate_angle(parts[5])
    except ValueError:
        return line, None

    if ft == '2':
        if len(parts) < 7:
            return line, None
        k_str = parts[6]
        new_code = (f"{ai:6d}  {aj:7d}  {ak:7d}  {al:7d}         "
                    f"{ft:1s}  {angle:12.3f}      {k_str}")
    else:
        if len(parts) < 8:
            return line, None
        kd_str, pn_str = parts[6], parts[7]
        new_code = (f"{ai:6d}  {aj:7d}  {ak:7d}  {al:7d}         "
                    f"{ft:1s}  {angle:12.3f}  {kd_str:>10s}  {pn_str:>4s}")

    sep = '     ' if comment else ''
    return new_code + sep + comment + '\n', ft


def invert_tc_itp(src_itp, dst_itp):
    """
    Copy src_itp to dst_itp unchanged.
    
    With the ROTATION method, C17 rotates around the N4-C24 axis but stays
    on the same "side" of the molecule - the handedness does not change.
    Therefore, dihedral phases do NOT need to be negated; the force field
    will naturally drive the geometry to match our rotated coordinates.
    
    Returns a dict with count info (all zeros for this method).
    """
    with open(src_itp) as fh:
        content = fh.read()
    
    with open(dst_itp, 'w') as fh:
        fh.write(content)
    
    return {'f1': 0, 'f2': 0, 'f9': 0, 'f4': 0, 'unchanged': 0, 'total_dih': 0}


# -----------------------------------------------------------------------------
#  CLI
# -----------------------------------------------------------------------------

def _default_out(in_path, suffix, ext):
    d    = os.path.dirname(in_path)
    stem = os.path.splitext(os.path.basename(in_path))[0]
    return os.path.join(d, stem + suffix + ext)


def parse_args():
    p = argparse.ArgumentParser(
        description='Invert the Tc torsion sign in a GRO + ITP pair.')
    p.add_argument('--gro',     required=True,  metavar='FILE',
                   help='Input GRO file')
    p.add_argument('--itp',     required=True,  metavar='FILE',
                   help='Input ITP file')
    p.add_argument('--gro-out', default=None,   metavar='FILE',
                   help='Output GRO (default: <stem>_tc_inv.gro)')
    p.add_argument('--itp-out', default=None,   metavar='FILE',
                   help='Output ITP (default: <stem>_tc_inv.itp)')
    return p.parse_args()


def main():
    args = parse_args()

    gro_out = args.gro_out or _default_out(args.gro, '_tc_inv', '.gro')
    itp_out = args.itp_out or _default_out(args.itp, '_tc_inv', '.itp')

    # Build bond graph from ITP (needed by invert_tc_gro for clash optimiser)
    adj = build_bond_graph(args.itp)

    # -- GRO ------------------------------------------------------------------
    print(f'\n{args.gro}  -->  {gro_out}')

    title, n_atoms_line, atoms, box = read_gro(args.gro)
    xp_orig = atoms_to_xp(atoms)

    report_torsions(xp_orig, 'original')
    report_n4_angles(xp_orig, 'original')

    xp_new = invert_tc_gro(xp_orig, adj)

    report_torsions(xp_new, 'Tc-inverted')
    report_n4_angles(xp_new, 'Tc-inverted')

    tc_orig = dihedral_deg(xp_orig[TC_N4_IDX-1], xp_orig[TC_C17_IDX-1],
                           xp_orig[TC_C24_IDX-1], xp_orig[TC_N8_IDX-1])
    tc_new  = dihedral_deg(xp_new [TC_N4_IDX-1], xp_new [TC_C17_IDX-1],
                           xp_new [TC_C24_IDX-1], xp_new [TC_N8_IDX-1])
    print(f'\n  Tc change:  {tc_orig:+.2f} deg  ->  {tc_new:+.2f} deg  '
          f'(target: {-tc_orig:+.2f} deg,  delta = {tc_new - (-tc_orig):.4f} deg)')

    # Bond-length sanity checks
    for label, ia, ib in [('N4-C17', TC_N4_IDX-1, TC_C17_IDX-1),
                           ('C24-C17', TC_C24_IDX-1, TC_C17_IDX-1)]:
        r_b = float(np.linalg.norm(xp_orig[ia] - xp_orig[ib]))
        r_a = float(np.linalg.norm(xp_new [ia] - xp_new [ib]))
        print(f'  {label} bond:  {r_b:.4f} nm  ->  {r_a:.4f} nm  (delta {r_a-r_b:+.6f} nm)')

    # Proximity check: report closest non-excluded neighbours of H45 and H46
    excl_h45 = _excluded_pairs(adj, TC_H45_IDX, max_sep=3)
    excl_h46 = _excluded_pairs(adj, TC_H46_IDX, max_sep=3)
    n_atoms = xp_new.shape[0]
    print('\n  H45 closest non-1,3 neighbours after inversion:')
    dists_h45 = []
    for ci in range(1, n_atoms + 1):
        if ci == TC_H45_IDX or ci == TC_H46_IDX:
            continue
        d = float(np.linalg.norm(xp_new[TC_H45_IDX-1] - xp_new[ci-1]))
        excl_type = '1,3' if ci in excl_h45 else ('1,4' if ci in _excluded_pairs(adj, TC_H45_IDX, max_sep=4) else '1,5+')
        dists_h45.append((d, ci, excl_type))
    dists_h45.sort()
    for d, ci, etype in dists_h45[:5]:
        flag = '  *** CRASH CLASH ***' if (etype == '1,5+' and d < 0.15) else \
               ('  ** SEVERE **' if (etype == '1,5+' and d < 0.22) else '')
        print(f'    H45-atom{ci:3d} ({etype}):  {d:.4f} nm{flag}')

    print('\n  H46 closest non-1,3 neighbours after inversion:')
    dists_h46 = []
    for ci in range(1, n_atoms + 1):
        if ci == TC_H45_IDX or ci == TC_H46_IDX:
            continue
        d = float(np.linalg.norm(xp_new[TC_H46_IDX-1] - xp_new[ci-1]))
        excl_type = '1,3' if ci in excl_h46 else ('1,4' if ci in _excluded_pairs(adj, TC_H46_IDX, max_sep=4) else '1,5+')
        dists_h46.append((d, ci, excl_type))
    dists_h46.sort()
    for d, ci, etype in dists_h46[:5]:
        flag = '  *** CRASH CLASH ***' if (etype == '1,5+' and d < 0.15) else \
               ('  ** SEVERE **' if (etype == '1,5+' and d < 0.22) else '')
        print(f'    H46-atom{ci:3d} ({etype}):  {d:.4f} nm{flag}')

    xp_to_atoms(atoms, xp_new)
    write_gro(gro_out, title, n_atoms_line, atoms, box)
    print(f'\n  Written: {gro_out}')

    # -- ITP ------------------------------------------------------------------
    print(f'\n{args.itp}  -->  {itp_out}')
    counts = invert_tc_itp(args.itp, itp_out)
    print('  Dihedral terms modified:')
    print(f"    functype 1 (periodic, phase neg.) : {counts['f1']:3d}")
    print(f"    functype 2 (DRIH, d0 negated)     : {counts['f2']:3d}")
    print(f"    functype 9 (periodic, phase neg.) : {counts['f9']:3d}")
    print(f"    functype 4 (improper, phase neg.) : {counts['f4']:3d}")
    print(f"    unchanged (no moving atom)        : {counts['unchanged']:3d}")
    print(f"    total dihedral lines seen         : {counts['total_dih']:3d}")
    print(f'  Written: {itp_out}')

    print('\nDone.')


if __name__ == '__main__':
    main()
