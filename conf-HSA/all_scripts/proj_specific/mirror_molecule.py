#!/usr/bin/env python3
"""
mirror_molecule.py
==================
Mirror a GROMACS .gro coordinate file and adjust the .itp topology accordingly.

PHYSICS RATIONALE
-----------------
A mirror reflection (e.g. x -> -x) is an improper rotation.  Under such a
transformation the following terms change:

  .gro coordinates
    The chosen axis component of every atom coordinate is negated.

  Bonds (harmonic, functype 1)
    Bond lengths are scalars -> UNCHANGED.

  Angles (harmonic, functype 1 and bond-bond cross terms functype 3)
    Bond angles are defined by the cosine of the angle between two bond
    vectors.  Cosines are even functions and the angle is between 0 and 180
    deg, so they are invariant under reflection -> UNCHANGED.

  Proper dihedrals - DRIH style (functype 2, harmonic in d0)
    The dihedral angle is a *pseudoscalar*: it changes sign under any
    improper rotation (including mirror reflection).  Therefore d0 -> -d0.

  Proper dihedrals - periodic cosine (functype 9)
    E = kd * (1 + cos(pn*phi - phase))
    A reflection sends phi -> -phi, so the potential becomes
    E' = kd * (1 + cos(-pn*phi - phase)) = kd * (1 + cos(pn*phi + phase))
    This is equivalent to phase -> -phase.
    In this file the only phase values are 0.0 and 180.0, for which
    -0 = 0 and -180 mod 360 = 180, so these are numerically invariant.
    The code negates them anyway (correctly, as a general rule).

  Improper dihedrals (functype 4, periodic cosine)
    Same transformation as functype 9: phase -> -phase.
    All improper phases here are 180.0, so -180 mod 360 = 180 -> UNCHANGED.
    The code negates them anyway for correctness.

USAGE
-----
  python mirror_molecule.py -g lig.gro -i lig_g.itp
  python mirror_molecule.py -g lig.gro -i lig_g.itp --axis y
  python mirror_molecule.py -g lig.gro -i lig_g.itp -og mirror.gro -oi mirror.itp
"""

import argparse
import re
import sys
import os


# ── GRO parsing / writing ─────────────────────────────────────────────────────

def parse_gro(path):
    """Return (title, n_atoms, atom_lines, box_line)."""
    with open(path) as fh:
        lines = fh.readlines()
    title = lines[0]
    n_atoms = int(lines[1].strip())
    atom_lines = lines[2 : 2 + n_atoms]
    box_line = lines[2 + n_atoms]
    return title, n_atoms, atom_lines, box_line


def mirror_gro_line(line, axis):
    """Negate one coordinate component in a fixed-format .gro atom line.

    .gro fixed-format:
      [0:20]  residue num + name, atom name, atom num
      [20:28] x (8.3f)
      [28:36] y (8.3f)
      [36:44] z (8.3f)
      [44:]   optional velocities + newline
    """
    prefix = line[:20]
    x = float(line[20:28])
    y = float(line[28:36])
    z = float(line[36:44])
    suffix = line[44:]

    if   axis == 'x': x = -x
    elif axis == 'y': y = -y
    elif axis == 'z': z = -z
    else: raise ValueError(f"Unknown axis '{axis}'")

    return f"{prefix}{x:8.3f}{y:8.3f}{z:8.3f}{suffix}"


def write_gro(path, title, n_atoms, atom_lines, box_line):
    with open(path, 'w') as fh:
        fh.write(title)
        fh.write(f"{n_atoms:5d}\n")
        fh.writelines(atom_lines)
        fh.write(box_line)


# ── ITP line modification ─────────────────────────────────────────────────────

def split_comment(line):
    """Return (code_part, comment_str) where comment_str includes the leading
    semicolon and preserves the original comment text."""
    if ';' in line:
        idx = line.index(';')
        return line[:idx], line[idx:].rstrip('\n')
    return line.rstrip('\n'), ''


def negate_dihedral_angle(line, functype_expected):
    """Return the line with the relevant angle field negated.

    functype 2  (DRIH harmonic):  field index 5 = d0 (degrees)
    functype 9  (periodic proper): field index 5 = phase (degrees)
    functype 4  (periodic improper): field index 5 = phase (degrees)

    The function returns the line unchanged if the functype does not match.
    Column widths are preserved as closely as possible by using the original
    whitespace between tokens.
    """
    code, comment = split_comment(line)
    stripped = code.strip()
    if not stripped:
        return line

    parts = stripped.split()
    if len(parts) < 6:
        return line

    functype = parts[4]
    if functype not in functype_expected:
        return line

    angle = float(parts[5])
    angle = -angle

    # Normalise to (-180, 180] for cleanliness
    while angle <= -180.0:
        angle += 360.0
    while angle > 180.0:
        angle -= 360.0

    # Rebuild preserving original field layout as much as possible.
    # Original GROMACS/Sobtop format for dihedral lines:
    #   functype 2:  "%6d %8d %8d %8d %9d %12.3f %16s" + comment
    #   functype 9:  "%6d %8d %8d %8d %9d %12.3f %12s %4s" + comment
    # We reconstruct with the same field widths.
    ai, aj, ak, al = (int(parts[i]) for i in range(4))

    if functype == '2':
        # d0  k
        k_str = parts[6]
        new_code = (f"{ai:6d}  {aj:7d}  {ak:7d}  {al:7d}         "
                    f"{functype:1s}  {angle:12.3f}      {k_str}")
    elif functype in ('9', '4'):
        # phase  kd  pn
        kd_str = parts[6]
        pn_str = parts[7] if len(parts) > 7 else ''
        new_code = (f"{ai:6d}  {aj:7d}  {ak:7d}  {al:7d}         "
                    f"{functype:1s}  {angle:12.3f}  {kd_str:>10s}  {pn_str:>4s}")
    else:
        return line

    if comment:
        return new_code + '     ' + comment + '\n'
    return new_code + '\n'


# ── ITP section parsing / writing ─────────────────────────────────────────────

def parse_itp(path):
    """Return list of [tag, lines] pairs.

    tag is the section header text (e.g. 'dihedrals ; propers') or None for
    lines before any header.  lines includes the header line itself.
    """
    with open(path) as fh:
        raw = fh.readlines()

    sections = []
    current_tag = None
    current_lines = []
    header_re = re.compile(r'^\s*\[([^\]]+)\]')

    for line in raw:
        m = header_re.match(line)
        if m:
            sections.append([current_tag, current_lines])
            current_tag = m.group(1).strip()
            current_lines = [line]
        else:
            current_lines.append(line)

    sections.append([current_tag, current_lines])
    return sections


def classify_dihedral_section(tag):
    """Return 'proper', 'improper', or None."""
    tag_l = tag.lower()
    if 'dihedrals' not in tag_l:
        return None
    if 'improper' in tag_l:
        return 'improper'
    return 'proper'


def mirror_itp_sections(sections):
    """Return new sections list with dihedral phases/d0 negated."""
    # Count statistics
    counts = {'f2': 0, 'f9': 0, 'f4': 0}
    new_sections = []

    for tag, lines in sections:
        if tag is None:
            new_sections.append([tag, lines])
            continue

        kind = classify_dihedral_section(tag)
        if kind is None:
            new_sections.append([tag, lines])
            continue

        new_lines = []
        for line in lines:
            stripped = line.strip()
            # Skip blank lines, comments, and section headers
            if not stripped or stripped.startswith(';') or stripped.startswith('['):
                new_lines.append(line)
                continue

            parts = stripped.split()
            if len(parts) < 5:
                new_lines.append(line)
                continue

            ft = parts[4]

            if kind == 'proper' and ft == '2':
                new_lines.append(negate_dihedral_angle(line, {'2'}))
                counts['f2'] += 1
            elif kind == 'proper' and ft == '9':
                new_lines.append(negate_dihedral_angle(line, {'9'}))
                counts['f9'] += 1
            elif kind == 'improper' and ft == '4':
                new_lines.append(negate_dihedral_angle(line, {'4'}))
                counts['f4'] += 1
            else:
                new_lines.append(line)

        new_sections.append([tag, new_lines])

    return new_sections, counts


def write_itp(path, sections):
    with open(path, 'w') as fh:
        for tag, lines in sections:
            fh.writelines(lines)


# ── Main ──────────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(
        description='Mirror a GROMACS molecule (.gro + .itp) through a coordinate plane.',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__)
    parser.add_argument('-g', '--gro', required=True, help='Input .gro file')
    parser.add_argument('-i', '--itp', required=True, help='Input .itp file')
    parser.add_argument('--axis', default='x', choices=['x', 'y', 'z'],
                        help='Mirror through the plane perpendicular to this axis '
                             '(default: x, i.e. negate x-coordinates)')
    parser.add_argument('-og', '--out_gro', default=None,
                        help='Output .gro filename (default: <stem>_mirror.gro)')
    parser.add_argument('-oi', '--out_itp', default=None,
                        help='Output .itp filename (default: <stem>_mirror.itp)')
    args = parser.parse_args()

    gro_stem = os.path.splitext(args.gro)[0]
    itp_stem = os.path.splitext(args.itp)[0]
    out_gro  = args.out_gro or f"{gro_stem}_mirror.gro"
    out_itp  = args.out_itp or f"{itp_stem}_mirror.itp"

    # ── Mirror coordinates ────────────────────────────────────────────────
    print(f"Reading coordinates : {args.gro}")
    title, n_atoms, atom_lines, box_line = parse_gro(args.gro)
    mirrored_atoms = [mirror_gro_line(l, args.axis) for l in atom_lines]
    orig_title = title.rstrip('\n')
    new_title = orig_title + f"  [mirror-{args.axis}]\n"
    write_gro(out_gro, new_title, n_atoms, mirrored_atoms, box_line)
    print(f"Wrote               : {out_gro}  ({n_atoms} atoms, {args.axis} negated)")

    # ── Mirror topology ───────────────────────────────────────────────────
    print(f"Reading topology    : {args.itp}")
    sections = parse_itp(args.itp)
    new_sections, counts = mirror_itp_sections(sections)
    write_itp(out_itp, new_sections)
    print(f"Wrote               : {out_itp}")

    # ── Report ────────────────────────────────────────────────────────────
    print()
    print("Topology modifications (dihedral angles are pseudoscalars -> negate):")
    print(f"  [ dihedrals ; propers  ]  functype 2  (DRIH harmonic d0) : {counts['f2']:4d} entries  d0 -> -d0")
    print(f"  [ dihedrals ; propers  ]  functype 9  (periodic phase)   : {counts['f9']:4d} entries  phase -> -phase")
    f4_note = "(180 -> 180, numerically unchanged)" if counts['f4'] > 0 else ""
    print(f"  [ dihedrals ; impropers]  functype 4  (periodic phase)   : {counts['f4']:4d} entries  phase -> -phase  {f4_note}")
    print()
    print("Topology terms NOT changed (scalars, invariant under reflection):")
    print("  Bonds (functype 1)  |  Angles (functype 1 & 3)  |  charges  |  masses  |  LJ params")
    print()
    print("Done.")


if __name__ == '__main__':
    main()
