#!/usr/bin/env python3
"""
Print 1-based atom indices for selections used by cal_o1_abs.py.

Usage:
  python get_index.py [path/to/structure1 [path/to/structure2 ...]]

If no paths are provided the script will try `E64/ion.gro` and `E92/ion.gro`.
"""
import sys

try:
    import MDAnalysis as mda
except Exception as e:
    print(f"ERROR: failed to import MDAnalysis: {e}")
    sys.exit(2)


# Selections copied from cal_o1_abs.py (avoid importing that module because
# it executes analysis on import).
SELECTIONS = {
    "CA34": "protein and resid 34 and name CA",
    "CA42": "protein and resid 42 and name CA",
    "COM (E3P N1..N4)": "resname E3P and name N1 N2 N3 N4",
    "EU (EU3)": "name EU3",
}


def print_selection_info(u, sel_name, sel_str):
    atoms = u.select_atoms(sel_str)
    print(f"  - {sel_name}: '{sel_str}'  -> {atoms.n_atoms} atom(s)")
    if atoms.n_atoms == 0:
        print("      (no atoms selected)")
        return
    one_based = (atoms.indices + 1).tolist()
    print("      1-based indices:", ", ".join(map(str, one_based)))
    for a in atoms:
        print(f"      idx={a.index + 1:6d}  resid={a.resid:6d}  resname={a.resname:6s}  name={a.name:>4s}")


def main(paths):
    if not paths:
        paths = ["E64/ion.gro", "E92/ion.gro"]
    for path in paths:
        print(f"File: {path}")
        try:
            u = mda.Universe(path)
        except Exception as e:
            print(f"  ERROR: Failed to load '{path}': {e}")
            continue
        for key, sel in SELECTIONS.items():
            print_selection_info(u, key, sel)
        print()


if __name__ == "__main__":
    main(sys.argv[1:])
