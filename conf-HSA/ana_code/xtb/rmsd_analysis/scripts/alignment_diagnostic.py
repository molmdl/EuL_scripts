#!/usr/bin/env python3
"""
alignment_diagnostic.py — Verify alignment indices for SAP and TSAP systems.

For one SAP system and one TSAP system per method (solv_md, xtb),
print the atom names at computed alignment indices and confirm they
match expected donor atoms.

Usage:
  python rmsd_analysis/scripts/alignment_diagnostic.py
"""

import json
import sys
from pathlib import Path

import MDAnalysis as mda
import numpy as np

# ───────────────────────────────────────────────────────────────
# Constants (mirrored from rmsd_compute.py)
# ───────────────────────────────────────────────────────────────

ALIGN_EU8_SAP = sorted([0, 1, 2, 3, 4, 5, 6, 7, 54])
ALIGN_EU9_SAP = sorted([0, 1, 2, 3, 4, 5, 6, 7, 54, 63])


def load_tsap_mappings(mappings_path: Path):
    """Load TSAP→SAP atom index mappings and build inverse (SAP→TSAP)."""
    with open(mappings_path) as f:
        data = json.load(f)

    sap_to_tsap = {}
    for species in ("me", "phe"):
        key = f"{species}_tsap_to_sap"
        if key not in data:
            continue
        tsap_to_sap = data[key]
        n = len(tsap_to_sap)
        inv = [0] * n
        for tsap_idx, sap_idx in enumerate(tsap_to_sap):
            inv[sap_idx] = tsap_idx
        sap_to_tsap[species] = inv

    return sap_to_tsap


def get_eu8_align_idx(conformer, species, sap_to_tsap):
    """Get Eu+8 alignment indices for conformer and species."""
    if conformer == "sap":
        return ALIGN_EU8_SAP
    else:
        inv = sap_to_tsap[species]
        return sorted([inv[i] for i in ALIGN_EU8_SAP])


def get_eu9_align_idx(conformer, species, sap_to_tsap):
    """Get Eu+9 alignment indices for conformer and species."""
    if conformer == "sap":
        return ALIGN_EU9_SAP
    else:
        inv = sap_to_tsap[species]
        return sorted([inv[i] for i in ALIGN_EU9_SAP])


def check_solv_md(system_name, solv_dir, sap_to_tsap):
    """Check solv_md alignment indices for one system."""
    species = system_name.split("_")[0]
    conformer = system_name.split("_")[2]

    tpr_path = solv_dir / system_name / "prod_0.tpr"
    xtc_path = solv_dir / system_name / "solv_all.xtc"

    if not tpr_path.exists():
        print(f"    SKIP: {tpr_path} not found")
        return False

    u = mda.Universe(str(tpr_path), str(xtc_path))
    mol = u.select_atoms("resname MOL")

    eu8_idx = get_eu8_align_idx(conformer, species, sap_to_tsap)

    print(f"  solv_md: {system_name} ({conformer})")
    print(f"    Eu+8 align indices: {eu8_idx}")
    print(f"    Atom names at indices:")
    for idx in eu8_idx:
        atom = mol[idx]
        print(f"      [{idx}] = {atom.name} (element: {atom.element})")

    # Expected: 3 O donors, 5 N donors, 1 Eu (EU3 in GROMACS)
    eu_count = sum(1 for idx in eu8_idx if mol[idx].element in ("Eu", "EU"))
    o_count = sum(1 for idx in eu8_idx if mol[idx].element == "O")
    n_count = sum(1 for idx in eu8_idx if mol[idx].element == "N")

    expected_eu = 1
    expected_o = 3
    expected_n = 5

    ok = (eu_count == expected_eu) and (o_count == expected_o) and (n_count == expected_n)
    print(f"    Element counts: Eu={eu_count} (exp {expected_eu}), "
          f"O={o_count} (exp {expected_o}), N={n_count} (exp {expected_n})")
    print(f"    Result: {'PASS' if ok else 'FAIL'}")

    del u
    return ok


def check_xtb(system_name, data_dir, sap_to_tsap):
    """Check xtb alignment indices for one system."""
    species = system_name.split("_")[0]
    conformer = system_name.split("_")[2]

    traj_path = data_dir / system_name / "traj.xyz"

    if not traj_path.exists():
        print(f"    SKIP: {traj_path} not found")
        return False

    u = mda.Universe(str(traj_path), in_memory=True)

    eu8_idx = get_eu8_align_idx(conformer, species, sap_to_tsap)
    eu9_idx = get_eu9_align_idx(conformer, species, sap_to_tsap)

    print(f"  xtb: {system_name} ({conformer})")
    print(f"    Eu+8 align indices: {eu8_idx}")
    print(f"    Atom names at Eu+8 indices:")
    for idx in eu8_idx:
        atom = u.atoms[idx]
        print(f"      [{idx}] = {atom.name} (element: {atom.element})")
    print(f"    Eu+9 align indices: {eu9_idx}")
    print(f"    Atom names at Eu+9 indices:")
    for idx in eu9_idx:
        atom = u.atoms[idx]
        print(f"      [{idx}] = {atom.name} (element: {atom.element})")

    # Expected Eu+8: 3 O + 5 N + 1 Eu = 9 atoms
    # Expected Eu+9: 3 O + 6 N + 1 Eu = 10 atoms
    eu8_eu = sum(1 for idx in eu8_idx if u.atoms[idx].element == "Eu")
    eu8_o = sum(1 for idx in eu8_idx if u.atoms[idx].element == "O")
    eu8_n = sum(1 for idx in eu8_idx if u.atoms[idx].element == "N")
    eu9_eu = sum(1 for idx in eu9_idx if u.atoms[idx].element == "Eu")
    eu9_o = sum(1 for idx in eu9_idx if u.atoms[idx].element == "O")
    eu9_n = sum(1 for idx in eu9_idx if u.atoms[idx].element == "N")

    ok_eu8 = (eu8_eu == 1) and (eu8_o == 3) and (eu8_n == 5)
    ok_eu9 = (eu9_eu == 1) and (eu9_o == 3) and (eu9_n == 6)

    print(f"    Eu+8 element counts: Eu={eu8_eu} (exp 1), O={eu8_o} (exp 3), N={eu8_n} (exp 5)")
    print(f"    Eu+9 element counts: Eu={eu9_eu} (exp 1), O={eu9_o} (exp 3), N={eu9_n} (exp 6)")
    print(f"    Eu+8: {'PASS' if ok_eu8 else 'FAIL'}, Eu+9: {'PASS' if ok_eu9 else 'FAIL'}")

    del u
    return ok_eu8 and ok_eu9


def main():
    mappings_path = Path("analysis/atom_mappings.json")
    data_dir = Path("data")
    solv_dir = Path("solv_md")

    sap_to_tsap = load_tsap_mappings(mappings_path)
    print(f"Loaded mappings for: {list(sap_to_tsap.keys())}")
    print()

    # Print mapping verification
    for species, inv in sap_to_tsap.items():
        print(f"  {species}: SAP[54](Eu) -> TSAP[{inv[54]}], "
              f"SAP[7](5th N) -> TSAP[{inv[7]}], "
              f"SAP[63](cap N) -> TSAP[{inv[63]}]")
    print()

    # Test systems: one SAP, one TSAP, per species
    all_pass = True

    # solv_md checks (SAP and TSAP)
    print("=== solv_md Alignment Check ===")
    for sys_name in ["phe_sssD_sap", "phe_sssD_tsap"]:
        ok = check_solv_md(sys_name, solv_dir, sap_to_tsap)
        all_pass = all_pass and ok
    print()

    # xtb checks (SAP and TSAP)
    print("=== xtb Alignment Check ===")
    for sys_name in ["phe_sssD_sap", "phe_sssD_tsap"]:
        ok = check_xtb(sys_name, data_dir, sap_to_tsap)
        all_pass = all_pass and ok
    print()

    # com_md check (uses same Eu+8 indices as solv_md)
    print("=== com_md Alignment Check ===")
    com_dir = Path("com_md")
    for sys_name in ["phe_sssD_sap", "phe_sssD_tsap"]:
        pdb_path = com_dir / f"{sys_name}_fp" / "v1.pdb"
        xtc_path = com_dir / f"{sys_name}_fp" / "v1.xtc"
        if not pdb_path.exists():
            print(f"  SKIP {sys_name}: com_md not found")
            continue

        species = sys_name.split("_")[0]
        conformer = sys_name.split("_")[2]
        u = mda.Universe(str(pdb_path), str(xtc_path))
        mol = u.select_atoms("resname MOL")
        eu8_idx = get_eu8_align_idx(conformer, species, sap_to_tsap)

        print(f"  com_md: {sys_name} ({conformer})")
        print(f"    Eu+8 align indices: {eu8_idx}")
        print(f"    Atom names at indices:")
        for idx in eu8_idx:
            atom = mol[idx]
            print(f"      [{idx}] = {atom.name} (element: {atom.element})")

        eu_count = sum(1 for idx in eu8_idx if mol[idx].element in ("Eu", "EU"))
        o_count = sum(1 for idx in eu8_idx if mol[idx].element == "O")
        n_count = sum(1 for idx in eu8_idx if mol[idx].element == "N")

        ok = (eu_count == 1) and (o_count == 3) and (n_count == 5)
        print(f"    Element counts: Eu={eu_count} (exp 1), "
              f"O={o_count} (exp 3), N={n_count} (exp 5)")
        print(f"    Result: {'PASS' if ok else 'FAIL'}")
        all_pass = all_pass and ok
        del u
    print()

    # Final verdict
    if all_pass:
        print("ALIGNMENT_CHECK: PASS")
    else:
        print("ALIGNMENT_CHECK: FAIL (fix required)")


if __name__ == "__main__":
    main()
