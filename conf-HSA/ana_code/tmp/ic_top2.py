#!/usr/bin/env python3
"""
invert_chirality.py
-------------------
For each *sssL*.gro in ./invert/, produce the sssD conformer by:

  1. Mirror all coordinates:  z -> -z
  2. For each of the 3 stereocentre carbons (C74, C75, C76):
       - swap the C->H and C->methyl bond directions
         (place methyl-C along old C-H direction at C-C bond length,
          place stereocentre-H along old C-methyl direction at C-H bond length)
       - rebuild the 3 methyl hydrogens at tetrahedral geometry
         using calxyz (angles 180, -60, 60 degrees)

This mirrors the approach of legacy_make_mirror.py, extended to all three
stereocentres (sss) and adapted to GRO format / current atom numbering.

Torsions inverted (metal_geo_analysis.py definitions, 0-based):
  T1 N3-C8-C9-N4   (GRO: N4-C9-C10-N5)
  T2 N4-C10-C11-N5 (GRO: N5-C11-C12-N6)
  T3 N5-C12-C13-N6 (GRO: N6-C13-C14-N7)
  T4 N6-C14-C15-N3 (GRO: N7-C15-C16-N4)
  Tc N3-C16-C23-N7 (GRO: N4-C17-C24-N8)

Stereocentre carbons (S config preserved):
  C74 : stereocentre-H=H127, methyl=C85 (H86 H87 H88), aryl=C92,  arm-N=N73
  C75 : stereocentre-H=H126, methyl=C81 (H82 H83 H84), aryl=C103, arm-N=N72
  C76 : stereocentre-H=H125, methyl=C77 (H78 H79 H80), aryl=C114, arm-N=N71

Output:
  *sssL*.gro  ->  *sssD*.gro
  phe_sssL.itp -> phe_sssD.itp (all dihedrals inverted except stereocentre-local)
"""

import glob
import os

import numpy as np


# ── GRO I/O ──────────────────────────────────────────────────────────────────

def read_gro(path):
    """Parse a single-frame GRO file.
    Returns (title, atoms, box) where atoms is a list of dicts."""
    with open(path) as f:
        lines = f.readlines()
    title = lines[0].rstrip('\n')
    n_atoms = int(lines[1].strip())
    atoms = []
    for line in lines[2 : 2 + n_atoms]:
        atoms.append(
            dict(
                resnum=int(line[0:5]),
                resname=line[5:10].strip(),
                atomname=line[10:15].strip(),
                atomnum=int(line[15:20]),
                x=float(line[20:28]),
                y=float(line[28:36]),
                z=float(line[36:44]),
            )
        )
    box = lines[2 + n_atoms].rstrip('\n')
    return title, atoms, box


def write_gro(path, title, atoms, box):
    """Write atoms list back to GRO fixed-width format."""
    with open(path, 'w') as f:
        f.write(title + '\n')
        f.write(f'{len(atoms):5d}\n')
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


# ── Internal coordinate builder (calxyz equivalent) ──────────────────────────

def calxyz(a, b, c, r, theta_deg, phi_deg):
    """
    Place atom D at bond length r from C, angle B-C-D = theta, dihedral A-B-C-D = phi.
    Equivalent to the calxyz function used in the legacy script.

    a, b, c : (3,) positions of reference atoms A, B, C
    r       : C-D bond length
    theta   : B-C-D angle in degrees
    phi     : A-B-C-D dihedral in degrees
    Returns (3,) position of D.
    """
    theta = np.radians(theta_deg)
    phi = np.radians(phi_deg)

    bc = c - b
    bc /= np.linalg.norm(bc)
    ab = b - a
    ab /= np.linalg.norm(ab)

    n = np.cross(ab, bc)
    nm = np.linalg.norm(n)
    if nm < 1e-9:
        tmp = np.array([1.0, 0.0, 0.0]) if abs(bc[0]) < 0.9 else np.array([0.0, 1.0, 0.0])
        n = np.cross(bc, tmp)
    n /= np.linalg.norm(n)

    m = np.cross(n, bc)

    d = r * (
        -np.cos(theta) * bc
        + np.sin(theta) * np.cos(phi) * m
        + np.sin(theta) * np.sin(phi) * n
    )
    return c + d


# ── Dihedral angle ────────────────────────────────────────────────────────────

def dihedral_deg(p1, p2, p3, p4):
    b1 = p2 - p1
    b2 = p3 - p2
    b3 = p4 - p3
    n1 = np.cross(b1, b2)
    n2 = np.cross(b2, b3)
    nn1 = np.linalg.norm(n1)
    nn2 = np.linalg.norm(n2)
    if nn1 < 1e-9 or nn2 < 1e-9:
        return np.nan
    n1 /= nn1
    n2 /= nn2
    b2n = b2 / np.linalg.norm(b2)
    return float(np.degrees(np.arctan2(np.dot(np.cross(n1, n2), b2n), np.dot(n1, n2))))


def dihedral_from_xp(xp, i, j, k, l):
    """Compute dihedral for 1-based indices i,j,k,l from coordinate array xp."""
    return dihedral_deg(xp[i - 1], xp[j - 1], xp[k - 1], xp[l - 1])


def report_torsions(xp, label):
    """Print the 5 chirality-defining torsions (GRO 1-based atom numbers)."""

    def p(i):
        return xp[i - 1]

    ring_defs = [
        ('T1 N4 -C9 -C10-N5 ', 4, 9, 10, 5),
        ('T2 N5 -C11-C12-N6 ', 5, 11, 12, 6),
        ('T3 N6 -C13-C14-N7 ', 6, 13, 14, 7),
        ('T4 N7 -C15-C16-N4 ', 7, 15, 16, 4),
        ('Tc N4 -C17-C24-N8 ', 4, 17, 24, 8),
    ]
    print(f'  Torsions [{label}]:')
    for name, a, b, c, d in ring_defs:
        ang = dihedral_deg(p(a), p(b), p(c), p(d))
        print(f'    {name}: {ang:8.2f}°')


# ── Stereocentre definitions ──────────────────────────────────────────────────

STEREOCENTRES = [
    # (centre, H_idx, methyl_C, methyl_Hs, aryl_C, arm_N)
    (74, 127, 85, [86, 87, 88], 92, 73),
    (75, 126, 81, [82, 83, 84], 103, 72),
    (76, 125, 77, [78, 79, 80], 114, 71),
]


# ── Main inversion routine ────────────────────────────────────────────────────

def invert_molecule(xp):
    """
    Apply the legacy mirror + stereocentre fix to coordinate array xp (nm).

    Steps (following legacy_make_mirror.py):
      1. Mirror z:  xp[:,2] *= -1
      2. For each stereocentre:
           a. Capture bond lengths C-H and C-methyl from the post-mirror coords.
           b. Capture unit vectors:  dH = (H - C)/|H - C|
                                     dC = (methyl_C - C)/|methyl_C - C|
           c. Swap: place methyl_C at C + dH * bond_CC
                    place H         at C + dC * bond_CH
           d. Rebuild methyl Hs via calxyz at tetrahedral angles
              (dihedral 180°, -60°, +60° relative to the arm-N reference).

    Returns a modified copy of xp.
    """
    xp = xp.copy()

    xp[:, 2] *= -1.0

    for centre, h_idx, methyl_c, methyl_hs, aryl_c, arm_n in STEREOCENTRES:
        i_c = centre - 1
        i_h = h_idx - 1
        i_me = methyl_c - 1
        i_n = arm_n - 1

        p_c = xp[i_c].copy()
        p_h = xp[i_h].copy()
        p_me = xp[i_me].copy()

        d_h = p_h - p_c
        bond_ch = np.linalg.norm(d_h)
        d_h /= bond_ch

        d_me = p_me - p_c
        bond_cc = np.linalg.norm(d_me)
        d_me /= bond_cc

        xp[i_me] = p_c + d_h * bond_cc
        xp[i_h] = p_c + d_me * bond_ch

        p_n = xp[i_n]
        p_c_node = xp[i_c]
        p_me_new = xp[i_me]

        # Use pre-swap methyl position for methyl C-H length estimate.
        ch_len = np.mean([np.linalg.norm(xp[hj - 1] - p_me) for hj in methyl_hs])

        for hj, phi in zip(methyl_hs, [180.0, -60.0, 60.0]):
            xp[hj - 1] = calxyz(p_n, p_c_node, p_me_new, ch_len, 109.0, phi)

    return xp


# ── ITP utilities ────────────────────────────────────────────────────────────

# Atom sets that are purely local to each stereocentre (1-based GRO indices).
# A dihedral (i,j,k,l) is excluded from inversion iff all 4 indices fall
# within one of these sets, because the GRO inversion explicitly restores the
# S-configuration at C74, C75, C76 (the three stereocentres).
#
#   C74: centre=74, H=127, methyl_C=85, methyl_Hs=86/87/88, aryl=92,  arm_N=73
#   C75: centre=75, H=126, methyl_C=81, methyl_Hs=82/83/84, aryl=103, arm_N=72
#   C76: centre=76, H=125, methyl_C=77, methyl_Hs=78/79/80, aryl=114, arm_N=71

STEREOCENTRE_LOCAL_SETS = [
    {74, 73, 85, 86, 87, 88, 92, 127},   # C74
    {75, 72, 81, 82, 83, 84, 103, 126},  # C75
    {76, 71, 77, 78, 79, 80, 114, 125},  # C76
]

# Atoms physically relocated by invert_molecule(): the stereocentre-H and
# methyl-C at each stereocentre are swapped to new positions.  Any ft=2 DRIH
# term that includes one of these atoms has its geometry changed by the swap in
# a way that simple negation of d0_L cannot reproduce.  For these terms we read
# the actual dihedral angle from the D-GRO coordinates and use that as d0_D.
#
#   C74: methyl_C=85, stereocentre-H=127
#   C75: methyl_C=81, stereocentre-H=126
#   C76: methyl_C=77, stereocentre-H=125
SWAPPED_ATOMS = {85, 127, 81, 126, 77, 125}


def is_stereocentre_local(a, b, c, d):
    """Return True if all 4 atom indices belong to one stereocentre's local set."""
    quad = {a, b, c, d}
    return any(quad <= s for s in STEREOCENTRE_LOCAL_SETS)


def has_swapped_atom(a, b, c, d):
    """Return True if any of the 4 indices is a physically-swapped stereocentre atom."""
    return bool({a, b, c, d} & SWAPPED_ATOMS)


def split_comment(line):
    if ';' in line:
        idx = line.index(';')
        return line[:idx], line[idx:].rstrip('\n')
    return line.rstrip('\n'), ''


def negate_dihedral_angle(line):
    """Negate dihedral d0/phase for functype 2/9/4.

    Returns (new_line, ft_str) or (original_line, None) on parse failure.
    """
    code, comment = split_comment(line)
    stripped = code.strip()
    if not stripped:
        return line, None

    parts = stripped.split()
    if len(parts) < 6:
        return line, None

    try:
        ai, aj, ak, al = (int(parts[i]) for i in range(4))
    except ValueError:
        return line, None

    ft = parts[4]
    if ft not in {'2', '9', '4'}:
        return line, None

    try:
        angle = -float(parts[5])
    except ValueError:
        return line, None

    while angle <= -180.0:
        angle += 360.0
    while angle > 180.0:
        angle -= 360.0

    if ft == '2':
        if len(parts) < 7:
            return line, None
        k_str = parts[6]
        new_code = (
            f"{ai:6d}  {aj:7d}  {ak:7d}  {al:7d}         "
            f"{ft:1s}  {angle:12.3f}      {k_str}"
        )
    else:
        if len(parts) < 8:
            return line, None
        kd_str = parts[6]
        pn_str = parts[7]
        new_code = (
            f"{ai:6d}  {aj:7d}  {ak:7d}  {al:7d}         "
            f"{ft:1s}  {angle:12.3f}  {kd_str:>10s}  {pn_str:>4s}"
        )

    if comment:
        return new_code + '     ' + comment + '\n', ft
    return new_code + '\n', ft


def replace_ft2_d0(line, new_angle):
    """Replace d0 in a functype-2 DRIH line with new_angle (degrees).

    Returns (new_line, '2') or (original_line, None) on parse failure.
    """
    code, comment = split_comment(line)
    stripped = code.strip()
    if not stripped:
        return line, None
    parts = stripped.split()
    if len(parts) < 7:
        return line, None
    try:
        ai, aj, ak, al = (int(parts[i]) for i in range(4))
    except ValueError:
        return line, None
    if parts[4] != '2':
        return line, None

    angle = float(new_angle)
    while angle <= -180.0:
        angle += 360.0
    while angle > 180.0:
        angle -= 360.0

    k_str = parts[6]
    new_code = (
        f"{ai:6d}  {aj:7d}  {ak:7d}  {al:7d}         "
        f"2  {angle:12.3f}      {k_str}"
    )
    if comment:
        return new_code + '     ' + comment + '\n', '2'
    return new_code + '\n', '2'



def invert_all_itp_dihedrals(src_itp, dst_itp, xp_D=None):
    """Invert ALL dihedral terms in src_itp, writing result to dst_itp.

    Strategy (with stereocentre awareness):

    For functype-9 and functype-4 terms:
      - Negate phase for every entry.
      - EXCEPTION: skip if all 4 indices belong to one stereocentre's local set.

    For functype-2 DRIH terms:
      - Default: negate d0 (d0_D = -d0_L).
        Rationale: ITP encodes d0 = act (equilibrium dihedral). Under mirror,
        act_D = -act_L, so d0_D = -d0_L.
      - EXCEPTION 1: skip if purely stereocentre-local (S-config preserved by GRO
        inversion -- these terms describe geometry that is explicitly unchanged).
      - EXCEPTION 2 (swapped-atom correction): if any of the 4 indices is a
        physically-relocated atom (methyl_C or stereocentre-H at C74/C75/C76,
        which are swapped in invert_molecule()), then act_D != -act_L because the
        atom was moved to a new position (~118 deg offset). For these 12 terms,
        read the actual dihedral angle from the D-GRO coordinate array (xp_D)
        and use that as d0_D directly (d0_D = act_D from D-GRO).

    Parameters
    ----------
    src_itp : str
        Path to the L-form ITP.
    dst_itp : str
        Path to write the D-form ITP.
    xp_D : np.ndarray, shape (N,3), optional
        D-form coordinate array (0-based, nm).  Required for swapped-atom
        correction; if None the correction is skipped.
    """
    with open(src_itp) as fh:
        lines = fh.readlines()

    out_lines = []
    in_dihedrals = False
    counts = {
        'f2_neg': 0, 'f2_geo': 0, 'f9': 0, 'f4': 0,
        'skipped': 0, 'total': 0,
    }

    for line in lines:
        stripped = line.strip()
        low = stripped.lower()

        if low.startswith('['):
            in_dihedrals = ('dihedrals' in low)
            out_lines.append(line)
            continue

        if not in_dihedrals or not stripped or stripped.startswith(';'):
            out_lines.append(line)
            continue

        code, _comment = split_comment(line)
        parts = code.split()
        if len(parts) < 5:
            out_lines.append(line)
            continue

        try:
            ai, aj, ak, al = (int(parts[i]) for i in range(4))
        except ValueError:
            out_lines.append(line)
            continue

        # Skip dihedrals that are purely local to a stereocentre.
        if is_stereocentre_local(ai, aj, ak, al):
            counts['skipped'] += 1
            out_lines.append(line)
            continue

        ft = parts[4] if len(parts) > 4 else ''

        # Swapped-atom correction for ft=2: use actual D-GRO dihedral as d0.
        if ft == '2' and xp_D is not None and has_swapped_atom(ai, aj, ak, al):
            act_D = dihedral_from_xp(xp_D, ai, aj, ak, al)
            if not np.isnan(act_D):
                new_line, ftr = replace_ft2_d0(line, act_D)
                if ftr is not None:
                    out_lines.append(new_line)
                    counts['f2_geo'] += 1
                    counts['total'] += 1
                    continue
            # Fall through to standard negation if geometry read failed.

        new_line, ftr = negate_dihedral_angle(line)
        if ftr is None:
            out_lines.append(line)
            continue
        out_lines.append(new_line)

        if ftr == '2':
            counts['f2_neg'] += 1
            counts['total'] += 1
        elif ftr == '9':
            counts['f9'] += 1
            counts['total'] += 1
        elif ftr == '4':
            counts['f4'] += 1
            counts['total'] += 1

    with open(dst_itp, 'w') as fh:
        fh.writelines(out_lines)

    return counts


# ── Main ─────────────────────────────────────────────────────────────────────

def main():
    script_dir = os.path.dirname(os.path.abspath(__file__))
    pattern = os.path.join(script_dir, '*sssL*.gro')
    gro_files = sorted(glob.glob(pattern))

    if not gro_files:
        print(f'No *sssL*.gro files found in {script_dir}')
        return

    # Keep D-GRO coords for phe so we can pass them to the ITP inverter.
    phe_xp_D = None

    for src in gro_files:
        basename = os.path.basename(src)
        dst_name = basename.replace('L', 'D')
        dst = os.path.join(script_dir, dst_name)

        print(f'\n{basename}  -->  {dst_name}')

        title, atoms, box = read_gro(src)
        xp_orig = atoms_to_xp(atoms)

        report_torsions(xp_orig, 'original (L)')

        xp_new = invert_molecule(xp_orig)

        report_torsions(xp_new, 'inverted (D)')

        xp_to_atoms(atoms, xp_new)
        write_gro(dst, title, atoms, box)
        print(f'  Written: {dst}')

        if 'phe' in basename:
            phe_xp_D = xp_new

    src_itp = os.path.join(script_dir, 'phe_sssL.itp')
    dst_itp = os.path.join(script_dir, 'phe_sssD.itp')
    if os.path.exists(src_itp):
        if phe_xp_D is None:
            # Fallback: read from already-written D-GRO file.
            _, atoms_d, _ = read_gro(os.path.join(script_dir, 'phe_sssD.gro'))
            phe_xp_D = atoms_to_xp(atoms_d)
        counts = invert_all_itp_dihedrals(src_itp, dst_itp, xp_D=phe_xp_D)
        print(f'\n{os.path.basename(src_itp)}  -->  {os.path.basename(dst_itp)}')
        print('  Topology modifications:')
        print(f"    functype 2 (d0 negated)      : {counts['f2_neg']:3d}")
        print(f"    functype 2 (d0 from D-GRO)   : {counts['f2_geo']:3d}  [swapped-atom correction]")
        print(f"    functype 9 (phase negated)    : {counts['f9']:3d}")
        print(f"    functype 4 (phase negated)    : {counts['f4']:3d}")
        print(f"    skipped (stereocentre-local)  : {counts['skipped']:3d}")
        print(f"    total modified                : {counts['total']:3d}")
    else:
        print(f'\nTopology source not found, skipping: {src_itp}')

    print('\nDone.')


if __name__ == '__main__':
    main()
