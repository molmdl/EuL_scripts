#!/usr/bin/env python3
"""
invert_tc.py
------------
Invert the sign of the Tc (12-N-4 ring -> chromophore) torsion angle in
a single-frame GRO file and its matching ITP file.

Tc is defined as:  N4 - C17 - C24 - N8  (1-based GRO atom indices)

Strategy:
  GRO (two-step):
    Step 1 — Reflect C17, H45, H46 across the plane defined by N4, C24, N8.
      Because N4 and C24 both lie on the reflection plane, the N4-C17 and
      C24-C17 bond lengths are exactly preserved.  This operation flips Tc
      from its original value to -Tc exactly.

    Step 2 — Rotate H45 and H46 around the N4->C17(reflected) axis to
      minimise 1,4+ non-bonded clashes introduced by the reflection.
      H45 and H46 are not part of the Tc definition (N4-C17-C24-N8), so
      this rotation does not change Tc.  Only H positions change; all
      heavy-atom geometry (in particular N4-C17-C24-N8) is unchanged.

    The N4 valence angles C9-N4-C17 and C16-N4-C17 do change after the
    reflection (this is physically real for TSAP geometry).

  ITP:
    Every dihedral term (functype 1, 2, 9, or 4) whose four atom indices
    include any of the MOVING atoms (C17, H45, H46) has its d0 / phase negated.
    C24 is NOT in the pivot set because it does not move during the reflection.
    All other dihedral lines are copied verbatim.

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
H_ROT_STEP_DEG  = 0.1   # degrees (scan resolution)

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

    Algorithm (two steps):

    Step 1 — Reflect C17, H45, H46 across the plane defined by N4, C24, N8.
      N4 and C24 lie on the plane, so:
        * N4-C17 bond length is exactly preserved
        * C24-C17 bond length is exactly preserved
        * Tc flips from tc0 to -tc0 exactly

    Step 2 — Rotate H45, H46 around the N4->C17(reflected) axis to minimise
      1,5+ non-bonded clashes.  The H atoms are not part of the Tc definition
      (N4-C17-C24-N8), so rotating them cannot change Tc.

    Parameters
    ----------
    xp  : (N,3) float array of coordinates in nm (0-based indexing)
    adj : adjacency dict from build_bond_graph(itp_path)

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

    p_n4  = xp[i_n4].copy()
    p_c24 = xp[i_c24].copy()
    p_n8  = xp[i_n8].copy()

    # ------------------------------------------------------------------
    # Step 1: reflect C17, H45, H46 across plane(N4, C24, N8)
    # ------------------------------------------------------------------
    plane_n = np.cross(p_c24 - p_n4, p_n8 - p_n4)
    pn_norm = np.linalg.norm(plane_n)
    if pn_norm < 1e-9:
        sys.exit('ERROR: N4, C24, N8 are collinear - cannot define reflection plane.')
    plane_n /= pn_norm

    for idx in (i_c17, i_h45, i_h46):
        rel = xp[idx] - p_n4
        xp[idx] = xp[idx] - 2.0 * np.dot(rel, plane_n) * plane_n

    # ------------------------------------------------------------------
    # Step 2: rotate H45 and H46 around N4->C17(reflected) to minimise
    #         1,5+ non-bonded clashes
    # ------------------------------------------------------------------
    # Rotation axis: N4 -> C17 (after reflection)
    c17_refl = xp[i_c17].copy()
    rot_axis = c17_refl - p_n4
    axis_norm = np.linalg.norm(rot_axis)
    if axis_norm < 1e-9:
        sys.exit('ERROR: N4 and C17(reflected) are coincident.')
    rot_axis /= axis_norm

    # Build excluded-pair sets (1,2 + 1,3) for H45 and H46 (1-based)
    excl_h45 = _excluded_pairs(adj, TC_H45_IDX, max_sep=3)
    excl_h46 = _excluded_pairs(adj, TC_H46_IDX, max_sep=3)

    # Gather all atoms that can clash with H45 or H46 (1,4 and further = not excluded)
    n_atoms = xp.shape[0]
    all_indices_1based = set(range(1, n_atoms + 1))
    candidates_h45 = sorted(all_indices_1based - excl_h45)
    candidates_h46 = sorted(all_indices_1based - excl_h46)

    # Score function: repulsive (sigma/d)^12 for pairs closer than CLASH_CUTOFF
    def clash_score(h45_pos, h46_pos):
        score = 0.0
        for ci in candidates_h45:
            idx0 = ci - 1
            if idx0 == i_h45 or idx0 == i_h46:
                continue
            d = float(np.linalg.norm(h45_pos - xp[idx0]))
            if d < CLASH_CUTOFF:
                score += (CLASH_SIGMA / d) ** 12
        for ci in candidates_h46:
            idx0 = ci - 1
            if idx0 == i_h45 or idx0 == i_h46:
                continue
            d = float(np.linalg.norm(h46_pos - xp[idx0]))
            if d < CLASH_CUTOFF:
                score += (CLASH_SIGMA / d) ** 12
        return score

    # Angular scan over [-180, 180) at H_ROT_STEP_DEG resolution
    angles = np.arange(-180.0, 180.0, H_ROT_STEP_DEG)
    best_angle = 0.0
    best_score = None

    h_pts_orig = np.array([xp[i_h45], xp[i_h46]])

    for ang in angles:
        h_pts_trial = rotate_around_axis(h_pts_orig, p_n4, rot_axis, ang)
        s = clash_score(h_pts_trial[0], h_pts_trial[1])
        if best_score is None or s < best_score:
            best_score = s
            best_angle = ang

    # Apply best rotation to H45 and H46
    h_pts_best = rotate_around_axis(h_pts_orig, p_n4, rot_axis, best_angle)
    xp[i_h45] = h_pts_best[0]
    xp[i_h46] = h_pts_best[1]

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
    Copy src_itp to dst_itp, negating d0/phase for every dihedral term
    whose four atom indices include any MOVING atom: C17, H45, or H46.

    C24 is intentionally excluded from the pivot set because it does not
    move during the reflection.  Including C24 would incorrectly modify
    aromatic-ring dihedral terms that are unaffected by the Tc change
    (e.g. ring-planarity impropers, ca-ca-ca-X terms, EU-N-C-C DRIH terms
    involving N8-C24 but not C17).

    Returns a dict of modification counts.
    """
    # Only the atoms that physically move are in the pivot set
    pivot_atoms = {TC_C17_IDX, TC_H45_IDX, TC_H46_IDX}

    with open(src_itp) as fh:
        lines = fh.readlines()

    out_lines  = []
    in_dih     = False
    counts     = {'f1': 0, 'f2': 0, 'f9': 0, 'f4': 0, 'unchanged': 0, 'total_dih': 0}

    for line in lines:
        stripped = line.strip()
        low      = stripped.lower()

        # Track section: any '[' line updates the flag
        if low.startswith('['):
            in_dih = 'dihedrals' in low
            out_lines.append(line)
            continue

        if not in_dih or not stripped or stripped.startswith(';'):
            out_lines.append(line)
            continue

        # Parse atom indices
        code, _ = _split_comment(line)
        parts   = code.split()
        if len(parts) < 5:
            out_lines.append(line)
            continue
        try:
            ai, aj, ak, al = (int(parts[i]) for i in range(4))
        except ValueError:
            out_lines.append(line)
            continue

        counts['total_dih'] += 1

        # Only touch terms that involve at least one moving atom
        if not ({ai, aj, ak, al} & pivot_atoms):
            counts['unchanged'] += 1
            out_lines.append(line)
            continue

        new_line, ft = _negate_dihedral_line(line)
        out_lines.append(new_line)
        if   ft == '1': counts['f1'] += 1
        elif ft == '2': counts['f2'] += 1
        elif ft == '9': counts['f9'] += 1
        elif ft == '4': counts['f4'] += 1
        else:
            counts['unchanged'] += 1   # parse failed, kept original

    with open(dst_itp, 'w') as fh:
        fh.writelines(out_lines)

    return counts


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
