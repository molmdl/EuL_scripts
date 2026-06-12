"""B0: Inspect prod_0.tpr atoms and compare with xTB xtbopt.xyz.

Read a GROMACS TPR (with solvent), extract ligand (resname != water/ions),
write atom names to JSON, and compare with xTB reference.
"""

import argparse
import json
import sys
from pathlib import Path

try:
    import MDAnalysis as mda
except ImportError:
    print("ERROR: MDAnalysis not available in current conda env.")
    sys.exit(1)


def read_xtb_atoms(xyz_path: Path):
    """Read atom symbols from xTB XYZ (n atoms, then comment, then n lines of symbol x y z)."""
    lines = xyz_path.read_text().strip().splitlines()
    n = int(lines[0].strip())
    atoms = []
    for line in lines[2:2+n]:
        tok = line.strip().split()
        atoms.append(tok[0])
    assert len(atoms) == n, f"Expected {n} atoms, got {len(atoms)}"
    return atoms


def read_tpr_atoms(tpr_path: Path):
    """Read TPR, identify ligand residue(s), return atom names and residue info."""
    u = mda.Universe(str(tpr_path))
    resnames = set(u.residues.resnames)
    # Identify non-standard residues (likely ligand)
    solv_res = {"SOL", "WAT", "HOH", "TIP3P", "OPC3", "NA", "CL", "K", "CA", "MG"}
    ligand_resnames = sorted(r for r in resnames if r not in solv_res)
    print(f"  TPR resnames: {sorted(resnames)}")
    print(f"  Candidate ligand resname(s): {ligand_resnames}")

    # Select ligand atoms
    if len(ligand_resnames) == 0:
        print("  WARNING: no non-solvent residue found")
        return None
    sel_str = " or ".join(f"resname {r}" for r in ligand_resnames)
    lig = u.select_atoms(sel_str)
    print(f"  Ligand selection: '{sel_str}' -> {len(lig)} atoms")
    atoms = list(lig.names)
    residues = list(lig.resnames)
    return {
        "ligand_resnames": ligand_resnames,
        "n_atoms": len(atoms),
        "atoms": atoms,
        "residues": residues,
        "n_residues_total": len(u.residues),
        "n_atoms_total": len(u.atoms),
    }


def main():
    parser = argparse.ArgumentParser(description="Inspect prod_0.tpr atoms vs xTB xyz")
    parser.add_argument("--tpr", required=True, type=Path)
    parser.add_argument("--xyz", required=True, type=Path)
    parser.add_argument("--out", type=Path, default=None)
    args = parser.parse_args()

    print(f"\n=== {args.tpr.parent.name} ===")
    print(f"TPR: {args.tpr}")
    tpr_info = read_tpr_atoms(args.tpr)

    print(f"xTB: {args.xyz}")
    xtb_atoms = read_xtb_atoms(args.xyz)
    print(f"  xTB atoms: {len(xtb_atoms)}")

    if tpr_info is None:
        print("  No comparison possible.")
        return

    tpr_atoms = tpr_info["atoms"]
    print(f"\n  Comparison: TPR ligand {len(tpr_atoms)} atoms vs xTB {len(xtb_atoms)} atoms")
    if len(tpr_atoms) == len(xtb_atoms):
        matches = sum(1 for a, b in zip(tpr_atoms, xtb_atoms) if a == b)
        print(f"  Direct match: {matches}/{len(tpr_atoms)} ({100*matches/len(tpr_atoms):.1f}%)")
    else:
        print(f"  LENGTH MISMATCH: cannot compare 1:1")

    # Try element mapping (strip numbers from GROMACS atom names)
    tpr_elem = [ ''.join(c for c in a if c.isalpha()) for a in tpr_atoms ]
    xtb_elem = [ ''.join(c for c in a if c.isalpha()) for a in xtb_atoms ]
    print(f"  Element sequences: TPR={tpr_elem[:10]}..., xTB={xtb_elem[:10]}...")
    if len(tpr_elem) == len(xtb_elem):
        elem_matches = sum(1 for a, b in zip(tpr_elem, xtb_elem) if a == b)
        print(f"  Element match: {elem_matches}/{len(tpr_elem)} ({100*elem_matches/len(tpr_elem):.1f}%)")
        # Also try reverse or offset? Report first 20 detail
        print(f"  First 20 TPR: {tpr_elem[:20]}")
        print(f"  First 20 xTB: {xtb_elem[:20]}")
    else:
        print(f"  Cannot compare elements (length mismatch)")

    # Write JSON
    if args.out:
        out = {
            "system": args.tpr.parent.name,
            "tpr_path": str(args.tpr),
            "xyz_path": str(args.xyz),
            "tpr": tpr_info,
            "xtb_atoms": xtb_atoms,
        }
        args.out.parent.mkdir(parents=True, exist_ok=True)
        args.out.write_text(json.dumps(out, indent=2))
        print(f"  Saved: {args.out}")


if __name__ == "__main__":
    main()
