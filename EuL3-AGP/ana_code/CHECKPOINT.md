# Chirality Inversion Checkpoint
*Last updated: 2026-03-18*

---

## Goal

Produce `phe_sssD.itp` and `*sssD*.gro` from their L-form counterparts, inverting
the topology/geometry for the D-enantiomer of a Eu(III) complex.

**Configuration rule:**
- All backbone/ring torsions → inverted (D-form)
- The 3 chiral carbons (C74, C75, C76) → **S-configuration preserved**

---

## Status: PRIMARY GOAL COMPLETE ✅

All three output files are correct and verified:

| File | Lines/Bytes | Status |
|------|-------------|--------|
| `invert/phe_sssD.gro` | 6810 bytes | ✅ verified |
| `invert/me_sssD.gro`  | 6180 bytes | ✅ verified |
| `invert/phe_sssD.itp` | 155255 bytes | ✅ gold-standard check passes |

The main script is `invert/invert_chirality.py` (553 lines). Run it from `invert/`:

```bash
cd invert/
python3 invert_chirality.py
```

Expected output summary:
```
functype 2 (d0 negated)      : 218
functype 2 (d0 from D-GRO)   :  12  [swapped-atom correction]
functype 9 (phase negated)    : 350
functype 4 (phase negated)    :  45
skipped (stereocentre-local)  :  27
total modified                : 625
```

---

## Key Discoveries (hard-won, do not forget)

### 1. DRIH convention: d0 = act (NOT -act)
The L-form ITP encodes `d0_L = act_L` (the actual geometric dihedral from the GRO
file). All 230 ft=2 terms satisfy `|d0_L - act_L| < 2°`. This was confirmed by
explicit cross-check against the GRO coordinates.

### 2. Mirror rule: d0_D = -d0_L
Under z-mirror: `act_D = -act_L`, so `d0_D = act_D = -act_L = -d0_L`. Negation
is correct for all non-exception ft=2 terms (218 terms).

### 3. Sign-convention pitfall (fixed)
An earlier pure-Python dihedral function had the wrong sign convention (returned
`-act` instead of `act`), causing false passes/failures. The correct function is
`dihedral_deg()` in `invert_chirality.py` (numpy arctan2, returns +act).

### 4. Swapped-atom correction (12 terms)
`invert_molecule()` physically swaps the stereocentre-H and methyl-C at each of
the 3 stereocentres. Any ft=2 term containing one of these atoms has
`act_D ≠ -act_L` (~118° offset). Fix: read `act_D` directly from D-GRO coords
and use that as `d0_D`.

Swapped atoms (1-based GRO indices):
```
SWAPPED_ATOMS = {85, 127, 81, 126, 77, 125}
# C74: methyl_C=85, H=127
# C75: methyl_C=81, H=126
# C76: methyl_C=77, H=125
```

### 5. Stereocentre-local terms skipped (27 terms)
A dihedral is skipped if all 4 indices fall within one stereocentre's local set:
```python
STEREOCENTRE_LOCAL_SETS = [
    {74, 73, 85, 86, 87, 88, 92, 127},   # C74
    {75, 72, 81, 82, 83, 84, 103, 126},  # C75
    {76, 71, 77, 78, 79, 80, 114, 125},  # C76
]
```
These are kept at L-form values because `invert_molecule()` explicitly restores
the S-configuration geometry at each stereocentre.

### 6. ft=9 and ft=4 phases
All phases are exactly 0° or 180°, so negation is numerically a no-op, but it is
done for correctness (350 ft=9 terms, 45 ft=4 terms).

### 7. Section-awareness not needed
The ITP has no cross-section functype mismatches; section-tracking is unnecessary.

---

## Stereocentre definitions

```python
STEREOCENTRES = [
    # (centre, H_idx, methyl_C, methyl_Hs, aryl_C, arm_N)
    (74, 127, 85, [86, 87, 88], 92,  73),   # C74
    (75, 126, 81, [82, 83, 84], 103, 72),   # C75
    (76, 125, 77, [78, 79, 80], 114, 71),   # C76
]
```

GRO inversion steps (per stereocentre, after z-mirror):
1. Capture unit vectors `d_H = (H - C)/|H-C|` and `d_Me = (Me - C)/|Me-C|`
2. Swap: place Me at `C + d_H * bond_CC`; place H at `C + d_Me * bond_CH`
3. Rebuild 3 methyl Hs via `calxyz` at tetrahedral angles (φ = 180°, −60°, +60°
   relative to arm-N reference)

---

## ITP Inversion Rule Summary

| Term type | Condition | Action |
|-----------|-----------|--------|
| ft=2 | stereocentre-local | **Skip** (keep d0_L) |
| ft=2 | has swapped atom | `d0_D = act_D` from D-GRO |
| ft=2 | all other | `d0_D = -d0_L` (negate) |
| ft=9 | any | `phase_D = -phase_L` (negate) |
| ft=4 | any | `phase_D = -phase_L` (negate) |

---

## Gold-Standard Verification (passes ✅)

Check: for all 230 ft=2 terms in `phe_sssD.itp`, `|d0_D - act_D| < 2°`
where `act_D` is computed from `phe_sssD.gro` using `dihedral_deg()`.

Also verified: L-form check passes (all 230 ft=2 terms: `|d0_L - act_L| < 2°`).

---

## File Structure

```
ana_code/
├── CHECKPOINT.md                        ← this file
├── invert/
│   ├── invert_chirality.py              ← MAIN SCRIPT (553 lines) ✅
│   ├── phe_sssL.gro                     ← input GRO, L-form, 149 atoms
│   ├── phe_sssD.gro                     ← output GRO, D-form ✅
│   ├── me_sssL.gro                      ← input GRO, me-variant L, 135 atoms
│   ├── me_sssD.gro                      ← output GRO, me-variant D ✅
│   ├── phe_sssL.itp                     ← input topology, L-form (1795 lines)
│   └── phe_sssD.itp                     ← output topology, D-form ✅
├── mirror/
│   ├── mirror_molecule.py               ← reference: full-molecule ITP inversion
│   ├── lig_g.itp                        ← reference L-form ITP
│   └── lig_mirror.gro                   ← reference mirror GRO (x-mirror)
├── phe_sssL_sap/
│   └── pr_0.gro                         ← MD-simulated L-form (14104 atoms, lig=first 149)
├── legacy_make_mirror.py                ← GRO-only mirror (old, 1 stereocentre, PDB numbering)
└── metal_geo_analysis.py                ← defines T1-T4, Tc torsion indices (0-based)
```

---

## Deferred / Open Items

1. **Amide-C / Eu distance issue** (deferred by user): an older GRO output was
   reported to have the amide carbon bonded to O1 too close to Eu. Investigation
   has NOT been done. The current `phe_sssD.gro` was not explicitly checked for
   this; it may or may not be present.

2. **`me_sssD.itp`**: not produced (no `me_sssL.itp` was provided). If the me
   variant needs an ITP, it would require the same procedure applied to a
   `me_sssL.itp` source file.

