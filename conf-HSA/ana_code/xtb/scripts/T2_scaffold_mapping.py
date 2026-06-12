#!/usr/bin/env python3
"""
T2_scaffold_mapping.py — Build me ↔ phe scaffold atom mapping.

Identifies which atom indices form the **common scaffold** (shared between
me and phe, excluding substituents) and provides cross-species index
mappings for joint PCA on common atoms.

Key insight (verified): The two molecules differ only at 2 attachment
points. In me_SAP, C[127] and C[128] are methyl carbons (CH₃ groups).
In phe_SAP, the same positions are the first two carbons of the phenyl
ring. Everything else (68 heavy + 59 H = 127 atoms) is the common
scaffold.

Usage:
    python scripts/T2_scaffold_mapping.py --data-dir data --out-dir analysis
"""

import argparse
import json
import os
import shutil
import sys
from collections import Counter

import numpy as np


# ---------------------------------------------------------------------------
# XYZ parser (single-frame)
# ---------------------------------------------------------------------------
def read_xyz(path):
    """Return (elements: list[str], coords: np.ndarray shape (N,3))."""
    with open(path) as f:
        n = int(f.readline().strip())
        f.readline()  # comment line
        elems, coords = [], []
        for _ in range(n):
            parts = f.readline().split()
            elems.append(parts[0])
            coords.append([float(parts[1]), float(parts[2]), float(parts[3])])
    return elems, np.array(coords)


# ---------------------------------------------------------------------------
# Helper: find H atoms bonded to specified heavy atoms
# ---------------------------------------------------------------------------
def find_h_bonded_to(elems, coords, heavy_indices, cutoff=1.2):
    """Find all H atoms bonded to any of the specified heavy atoms.

    Parameters
    ----------
    elems : list[str]
        Element symbols for each atom.
    coords : np.ndarray shape (N, 3)
        Cartesian coordinates.
    heavy_indices : list[int]
        Indices of heavy atoms whose bonded H's to find.
    cutoff : float
        Maximum H-heavy distance to consider bonded (Å).

    Returns
    -------
    list[int]
        Sorted list of H-atom indices bonded to the specified heavy atoms.
    """
    h_bonded = []
    for h_idx in range(len(elems)):
        if elems[h_idx] != "H":
            continue
        for p_idx in heavy_indices:
            d = np.linalg.norm(coords[h_idx] - coords[p_idx])
            if d < cutoff:
                h_bonded.append(h_idx)
                break
    return sorted(h_bonded)


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def main():
    parser = argparse.ArgumentParser(
        description="Build me ↔ phe scaffold atom mapping"
    )
    parser.add_argument(
        "--data-dir", default="data",
        help="Directory containing species subdirectories (default: data)"
    )
    parser.add_argument(
        "--out-dir", default="analysis",
        help="Directory for output files (default: analysis)"
    )
    args = parser.parse_args()

    data_dir = args.data_dir
    out_dir = args.out_dir

    # ===================================================================
    # Step 1: Load reference XYZ structures and atom mappings
    # ===================================================================
    print("=" * 60)
    print("Step 1: Loading reference structures and atom mappings")
    print("=" * 60)

    me_sap_e, me_sap_c = read_xyz(os.path.join(data_dir, "me_rrrD_sap", "xtbopt.xyz"))
    me_tsap_e, me_tsap_c = read_xyz(os.path.join(data_dir, "me_rrrD_tsap", "xtbopt.xyz"))
    phe_sap_e, phe_sap_c = read_xyz(os.path.join(data_dir, "phe_rrrD_sap", "xtbopt.xyz"))
    phe_tsap_e, phe_tsap_c = read_xyz(os.path.join(data_dir, "phe_rrrD_tsap", "xtbopt.xyz"))

    print(f"  me_SAP:  {len(me_sap_e)} atoms  ({Counter(me_sap_e)['C']}C + "
          f"{Counter(me_sap_e)['N']}N + {Counter(me_sap_e)['O']}O + "
          f"{Counter(me_sap_e)['Eu']}Eu + {Counter(me_sap_e)['H']}H)")
    print(f"  me_TSAP: {len(me_tsap_e)} atoms")
    print(f"  phe_SAP: {len(phe_sap_e)} atoms  ({Counter(phe_sap_e)['C']}C + "
          f"{Counter(phe_sap_e)['N']}N + {Counter(phe_sap_e)['O']}O + "
          f"{Counter(phe_sap_e)['Eu']}Eu + {Counter(phe_sap_e)['H']}H)")
    print(f"  phe_TSAP:{len(phe_tsap_e)} atoms")

    # Load atom mappings from T1
    mappings_path = os.path.join(out_dir, "atom_mappings.json")
    with open(mappings_path) as f:
        mappings = json.load(f)

    me_tsap_to_sap = mappings["me_tsap_to_sap"]
    phe_tsap_to_sap = mappings["phe_tsap_to_sap"]

    # Build inverse mappings: SAP → TSAP
    me_sap_to_tsap = [0] * len(me_sap_e)
    for t, s in enumerate(me_tsap_to_sap):
        me_sap_to_tsap[s] = t

    phe_sap_to_tsap = [0] * len(phe_sap_e)
    for t, s in enumerate(phe_tsap_to_sap):
        phe_sap_to_tsap[s] = t

    print(f"  Loaded atom_mappings.json ({len(me_tsap_to_sap)} me, "
          f"{len(phe_tsap_to_sap)} phe mappings)")

    # ===================================================================
    # Step 2: Identify scaffold vs substituent in SAP ordering
    # ===================================================================
    print("\n" + "=" * 60)
    print("Step 2: Identifying scaffold vs substituent (SAP ordering)")
    print("=" * 60)

    # Find element differences in the first min(len(me), len(phe)) positions
    n_common = min(len(me_sap_e), len(phe_sap_e))
    sap_diffs = []
    for i in range(n_common):
        if me_sap_e[i] != phe_sap_e[i]:
            sap_diffs.append(i)

    print(f"  Element differences in first {n_common} positions: {sap_diffs}")
    for idx in sap_diffs:
        print(f"    idx {idx}: me={me_sap_e[idx]}, phe={phe_sap_e[idx]}")

    # Identify me substituent atoms
    me_subst_heavy_sap = [127, 128]  # the 2 methyl C atoms
    me_subst_h_sap = find_h_bonded_to(me_sap_e, me_sap_c, me_subst_heavy_sap, cutoff=1.2)
    print(f"  me methyl H atoms (bonded to C[127,128]): {me_subst_h_sap}")
    me_subst_sap = sorted(me_subst_heavy_sap + me_subst_h_sap)
    print(f"  me substituent (SAP): {me_subst_sap} ({len(me_subst_sap)} atoms)")

    # phe substituent atoms: everything from index 127 onward
    phe_subst_sap = list(range(127, len(phe_sap_e)))
    print(f"  phe substituent (SAP): indices {127}–{len(phe_sap_e)-1} "
          f"({len(phe_subst_sap)} atoms)")

    # Common scaffold SAP indices
    me_scaffold_sap = sorted(set(range(len(me_sap_e))) - set(me_subst_sap))
    phe_scaffold_sap = sorted(set(range(len(phe_sap_e))) - set(phe_subst_sap))

    print(f"  me scaffold (SAP): {len(me_scaffold_sap)} atoms "
          f"(indices 0–{max(me_scaffold_sap)})")
    print(f"  phe scaffold (SAP): {len(phe_scaffold_sap)} atoms "
          f"(indices 0–{max(phe_scaffold_sap)})")

    # ===================================================================
    # Step 3: Validate scaffold identification
    # ===================================================================
    print("\n" + "=" * 60)
    print("Step 3: Validating scaffold identification")
    print("=" * 60)

    # 3a. Same number of scaffold atoms
    assert len(me_scaffold_sap) == len(phe_scaffold_sap), (
        f"Scaffold size mismatch: me={len(me_scaffold_sap)} "
        f"vs phe={len(phe_scaffold_sap)}"
    )
    print(f"  ✓ Same scaffold size: {len(me_scaffold_sap)}")

    # 3b. Same elements at scaffold positions
    for sap_idx in me_scaffold_sap:
        assert me_sap_e[sap_idx] == phe_sap_e[sap_idx], (
            f"Element mismatch at scaffold SAP idx {sap_idx}: "
            f"me={me_sap_e[sap_idx]} vs phe={phe_sap_e[sap_idx]}"
        )
    print(f"  ✓ Same elements at all scaffold positions")

    # 3c. Same number of heavy scaffold atoms
    me_scaffold_heavy = [i for i in me_scaffold_sap if me_sap_e[i] != "H"]
    phe_scaffold_heavy = [i for i in phe_scaffold_sap if phe_sap_e[i] != "H"]
    assert len(me_scaffold_heavy) == len(phe_scaffold_heavy), (
        f"Heavy scaffold size mismatch: {len(me_scaffold_heavy)} "
        f"vs {len(phe_scaffold_heavy)}"
    )
    print(f"  ✓ Same heavy scaffold count: {len(me_scaffold_heavy)}")

    # 3d. All scaffold heavy indices are identical
    assert me_scaffold_heavy == phe_scaffold_heavy, "Scaffold heavy indices differ"
    print(f"  ✓ Scaffold heavy indices identical between species")

    n_h_scaffold = len(me_scaffold_sap) - len(me_scaffold_heavy)
    print(f"  Scaffold: {len(me_scaffold_sap)} atoms "
          f"({len(me_scaffold_heavy)} heavy + {n_h_scaffold} H)")
    print(f"  me substituent: {len(me_subst_sap)} atoms")
    print(f"  phe substituent: {len(phe_subst_sap)} atoms")

    # ===================================================================
    # Step 4: Map scaffold atoms to TSAP ordering
    # ===================================================================
    print("\n" + "=" * 60)
    print("Step 4: Mapping scaffold atoms to TSAP ordering")
    print("=" * 60)

    # Map SAP scaffold → TSAP scaffold
    me_scaffold_tsap = [me_sap_to_tsap[i] for i in me_scaffold_sap]
    phe_scaffold_tsap = [phe_sap_to_tsap[i] for i in phe_scaffold_sap]

    # Map SAP substituent → TSAP substituent
    me_subst_tsap = [me_sap_to_tsap[i] for i in me_subst_sap]
    phe_subst_tsap = [phe_sap_to_tsap[i] for i in phe_subst_sap]

    print(f"  me scaffold TSAP indices: {len(me_scaffold_tsap)} atoms, "
          f"range [{min(me_scaffold_tsap)}, {max(me_scaffold_tsap)}]")
    print(f"  phe scaffold TSAP indices: {len(phe_scaffold_tsap)} atoms, "
          f"range [{min(phe_scaffold_tsap)}, {max(phe_scaffold_tsap)}]")
    print(f"  me substituent TSAP: {sorted(me_subst_tsap)}")
    print(f"  phe substituent TSAP: {sorted(phe_subst_tsap)}")

    # Validate: same elements at mapped positions
    for j, sap_idx in enumerate(me_scaffold_sap):
        me_t = me_scaffold_tsap[j]
        phe_t = phe_scaffold_tsap[j]
        assert me_tsap_e[me_t] == phe_tsap_e[phe_t], (
            f"Element mismatch at scaffold pos {j} (SAP {sap_idx}): "
            f"me_TSAP[{me_t}]={me_tsap_e[me_t]} "
            f"vs phe_TSAP[{phe_t}]={phe_tsap_e[phe_t]}"
        )
    print(f"  ✓ Same elements at all mapped TSAP scaffold positions")

    # ===================================================================
    # Step 5: Build phe→me cross-species mappings
    # ===================================================================
    print("\n" + "=" * 60)
    print("Step 5: Building phe→me cross-species mappings")
    print("=" * 60)

    # SAP mapping: identity for scaffold atoms (same indices in both species)
    phe_to_me_sap = {str(phe_sap_idx): me_sap_idx
                     for phe_sap_idx, me_sap_idx
                     in zip(phe_scaffold_sap, me_scaffold_sap)}
    n_identity_sap = sum(1 for k, v in phe_to_me_sap.items() if int(k) == v)
    print(f"  phe_to_me_sap: {len(phe_to_me_sap)} scaffold atoms, "
          f"{n_identity_sap} identity mappings "
          f"({len(phe_to_me_sap) - n_identity_sap} non-identity)")

    # TSAP mapping: link through SAP anchor
    phe_to_me_tsap = {}
    for j, sap_idx in enumerate(me_scaffold_sap):
        phe_t = phe_scaffold_tsap[j]
        me_t = me_scaffold_tsap[j]
        phe_to_me_tsap[str(phe_t)] = me_t

    tsap_diffs = sum(1 for k, v in phe_to_me_tsap.items() if int(k) != v)
    print(f"  phe_to_me_tsap: {len(phe_to_me_tsap)} scaffold atoms, "
          f"{tsap_diffs} differ from identity")

    # Identify the non-identity mappings for reporting
    non_identity = [(int(k), v) for k, v in phe_to_me_tsap.items()
                    if int(k) != v]
    if non_identity:
        print(f"  Non-identity TSAP mappings (phe_TSAP → me_TSAP):")
        for phe_t, me_t in non_identity:
            # Find which SAP index this corresponds to
            sap_idx = None
            for j, s in enumerate(me_scaffold_sap):
                if me_scaffold_tsap[j] == me_t:
                    sap_idx = s
                    break
            print(f"    phe_TSAP[{phe_t}] ({phe_tsap_e[phe_t]}) → "
                  f"me_TSAP[{me_t}] ({me_tsap_e[me_t]})  "
                  f"[SAP idx {sap_idx}]")

    # ===================================================================
    # Step 6: Save output
    # ===================================================================
    print("\n" + "=" * 60)
    print("Step 6: Saving scaffold_mapping.json")
    print("=" * 60)

    output = {
        "scaffold_atoms_sap": {
            "me": me_scaffold_sap,
            "phe": phe_scaffold_sap,
        },
        "scaffold_atoms_tsap": {
            "me": sorted(me_scaffold_tsap),
            "phe": sorted(phe_scaffold_tsap),
        },
        "substituent_atoms": {
            "me": me_subst_sap,
            "phe": phe_subst_sap,
        },
        "phe_to_me_sap": phe_to_me_sap,
        "phe_to_me_tsap": phe_to_me_tsap,
        "n_scaffold": len(me_scaffold_sap),
        "n_substituent": {
            "me": len(me_subst_sap),
            "phe": len(phe_subst_sap),
        },
    }

    os.makedirs(out_dir, exist_ok=True)
    out_path = os.path.join(out_dir, "scaffold_mapping.json")

    # Back up existing file if present
    if os.path.exists(out_path):
        backup_path = out_path + ".bak"
        shutil.copy2(out_path, backup_path)
        print(f"  Backed up existing file to {backup_path}")

    with open(out_path, "w") as f:
        json.dump(output, f, indent=2)
    print(f"  Saved to {out_path}")

    # ===================================================================
    # Step 7: Final validation and summary
    # ===================================================================
    print("\n" + "=" * 60)
    print("Step 7: Final validation and summary")
    print("=" * 60)

    # Verify scaffold atoms partition the full atom list
    assert set(me_scaffold_sap) | set(me_subst_sap) == set(range(len(me_sap_e))), \
        "me scaffold + substituent do not partition complete atom set"
    assert set(me_scaffold_sap) & set(me_subst_sap) == set(), \
        "me scaffold and substituent overlap"
    assert set(phe_scaffold_sap) | set(phe_subst_sap) == set(range(len(phe_sap_e))), \
        "phe scaffold + substituent do not partition complete atom set"
    assert set(phe_scaffold_sap) & set(phe_subst_sap) == set(), \
        "phe scaffold and substituent overlap"
    print(f"  ✓ Scaffold + substituent = complete atom set (both species)")

    # Verify TSAP scaffold indices are valid
    assert all(0 <= i < len(me_tsap_e) for i in me_scaffold_tsap), \
        "me TSAP scaffold indices out of range"
    assert all(0 <= i < len(phe_tsap_e) for i in phe_scaffold_tsap), \
        "phe TSAP scaffold indices out of range"
    print(f"  ✓ All TSAP scaffold indices valid")

    # Verify phe_to_me_tsap values are in me's atom range
    assert all(0 <= v < len(me_tsap_e) for v in phe_to_me_tsap.values()), \
        "phe_to_me_tsap values out of me atom range"
    print(f"  ✓ All phe_to_me_tsap values valid")

    # Verify phe_to_me_tsap mapping count matches scaffold size
    assert len(phe_to_me_tsap) == len(me_scaffold_sap), \
        f"phe_to_me_tsap has {len(phe_to_me_tsap)} entries, expected {len(me_scaffold_sap)}"
    print(f"  ✓ phe_to_me_tsap has correct number of entries ({len(phe_to_me_tsap)})")

    # Summary
    print(f"\n{'='*60}")
    print(f"=== Scaffold Mapping Summary ===")
    print(f"{'='*60}")
    print(f"Common scaffold: {output['n_scaffold']} atoms "
          f"({len(me_scaffold_heavy)} heavy + {n_h_scaffold} H)")
    print(f"me substituent: {output['n_substituent']['me']} atoms "
          f"(2 methyl C + 6 methyl H)")
    print(f"phe substituent: {output['n_substituent']['phe']} atoms "
          f"(12 phenyl C + 10 phenyl H)")
    print(f"TSAP mapping: {len(phe_to_me_tsap)} scaffold atoms, "
          f"{tsap_diffs} differ between phe_TSAP and me_TSAP indices")

    print(f"\nALL VALIDATIONS PASSED ✓")


if __name__ == "__main__":
    main()
