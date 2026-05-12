# Amide-Aromatic Interaction Research

## 1. Definition of Amide-Aromatic (Amide-π) Interactions

Amide-aromatic interactions are noncovalent interactions between:
- **Amide groups** (from backbone peptide bonds or side chains) and
- **Aromatic rings** (Phe, Tyr, Trp, His side chains)

Two main types exist:
- **NH-π interactions**: The amide N-H hydrogen points toward the aromatic π-electron cloud (hydrogen bond-like)
- **π-π stacking**: The planar amide group stacks parallel to the aromatic ring
- **Lone pair-π (n-π*)**: Interaction between amide oxygen lone pairs and aromatic π-system

---

## 2. Which Amide Groups Are Involved

**Both backbone and side chain amides participate:**

| Amide Source | Examples | Notes |
|-------------|----------|-------|
| **Backbone amides** | Peptide CONH | Most commonly studied; NH-π interactions with aromatic rings |
| **Side chain amides** | Asn, Gln | Can participate but less frequently documented |
| **Acetylated lysine** | Ac-Lys | Forms amide-π interactions replacing cation-π upon acetylation |

**Key finding** (Hughes & Waters, 2006, JACS): Acetylated lysine forms amide-π interactions with Trp that are energetically comparable to the cation-π interaction formed by unmodified lysine.

---

## 3. Geometric Criteria from Literature

### Distance Cutoffs:

| Study | Distance Criterion | Measurement |
|-------|-------------------|-------------|
| **Goyal et al. 2017 (JACS)** | ~3.3 Å | Distance between amide N/H and aromatic ring center |
| **Baskaran et al. 2021 (Magnetic Resonance)** | < 8 Å (for analysis) | H to ring center; strong effects within 3-4 Å |
| **Levitt & Perutz 1988** | ~3.3 Å | N to ring center (local minimum in VdW energy) |
| **van der Spoel et al. 1996** | 1.0 nm cutoff used | Based on geometric centers |

**Preferred distances:**
- **NH-π interaction**: ~3.3-3.4 Å (H to ring center)
- **π-π stacking**: ~3.3 Å (ring-to-ring centroid distance)
- **Hydrogen bond-like**: ~2.8 Å peak (but NH-π is slightly longer)

### Angle Criteria:

| Study | Angle Definition | Preferred Range |
|-------|-----------------|-----------------|
| **Goyal et al. 2017** | N-H...π angle | 155-160° (most preferred) |
| **Baskaran et al. 2021** | Azimuthal angle θ | ~25° (Tyr/Phe outliers); Trp shows angles above ring plane; His shows in-plane |

### Specific Atom Measurements:

- **For NH-π**: Measure from amide **hydrogen** to aromatic **ring centroid**
- **For n-π***: Measure from amide **oxygen** to aromatic **ring centroid**
- **Ring centroid**: Geometric center of aromatic ring atoms
- **Ring normal**: Vector perpendicular to ring plane (used for angle calculations)

---

## 4. Computational Detection Methods and Tools

### Published Methods:

| Tool/Method | Reference | Notes |
|-------------|-----------|-------|
| **Custom geometric criteria** | Goyal et al. 2017 | Distance + angle thresholds implemented in MD analysis |
| **BMRB/PDB mining** | Baskaran et al. 2021 | Database surveys using chemical shift anomalies |
| **PDB structural surveys** | Burley & Petsko 1986; Brandl et al. 2001 | Crystal structure analysis |
| **MD simulation analysis** | van der Spoel et al. 1996 | Molecular dynamics with geometric cutoffs |

### Detection Approach (from literature):

```python
# Typical criteria for NH-π detection (based on Goyal et al. 2017):
# 1. Distance: H(ring_centroid) < 4.0 Å (peak at ~3.3 Å)
# 2. Angle: N-H...ring_centroid angle > 150° (optimal 155-160°)
# 3. Azimuthal angle: < 45° for above-ring interactions
```

---

## 5. Energy/Strength Comparison

### Interaction Energies:

| Interaction Type | Energy (kcal/mol) | Source |
|-----------------|-------------------|--------|
| **Amide-π (NH-π)** | ~1.5-3.0 kcal/mol | Perutz 1993; experimental estimates |
| **Amide-π** | ~0.5-1.5 kcal/mol | Goyal et al. 2017 (QM calculations) |
| **Cation-π** | ~2-5 kcal/mol | General literature |
| **Conventional H-bond** | ~4-6 kcal/mol | Standard values |

### Key Findings:

**Hughes & Waters 2006 (JACS):**
- Amide-π interaction (acetyl-Lys...Trp) is **energetically comparable** to cation-π (Lys...Trp)
- Despite losing positive charge, the interaction energy is "not significantly perturbed"
- **Enthalpically driven** (polar-π interaction nature)

**Zytkiewicz et al. 2023 (Biochemistry):**
- **sp²O-aromatic sp²C**: Enthalpy-driven, entropically unfavorable (direct chemical interaction/lone pair-π)
- **sp³C/sp²C-aromatic sp²C**: Entropy-driven, enthalpically unfavorable (hydrophobic effects)
- Enthalpic and entropic contributions are much larger in magnitude than the free energy itself

**Cheng et al. 2017 (JACS):**
- Amide-O interaction is **more favorable** than amide-N interaction with aromatics

---

## 6. Key Literature References

| Ref | Year | Key Parameters | Journal |
|-----|------|----------------|---------|
| **Levitt & Perutz** | 1988 | d(N-ring) = 3.3 Å (VdW minimum); ~3 kcal/mol | J Mol Biol |
| **Burley & Petsko** | 1986 | First systematic survey; 60% aromatics in π-π | FEBS Lett |
| **Brandl et al.** | 2001 | 14,087 C-H...π in 1154 PDB structures | J Mol Biol |
| **Hughes & Waters** | 2006 | Amide-π vs cation-π comparison | JACS |
| **Goyal et al.** | 2017 | d = 3.3 Å; angle = 155-160°; QM + MD | JACS |
| **Baskaran et al.** | 2021 | d < 8 Å analysis; θ analysis; NMR evidence | Magnetic Resonance |
| **Zytkiewicz et al.** | 2023 | Thermodynamic decomposition; enthalpy vs entropy | Biochemistry |
| **Cheng et al.** | 2017 | Amide-O vs amide-N comparison | JACS |

---

## Summary for Implementation

For detecting amide-aromatic interactions in protein structures, use:

1. **Distance cutoff**: H or O to ring centroid ≤ 4.0 Å (with peak probability at ~3.3 Å)
2. **Angle criterion**: N-H...ring angle ≥ 150° for NH-π
3. **Atoms to measure**: Amide H (for NH-π) or O (for n-π*) to aromatic ring centroid
4. **Aromatic rings**: Phe, Tyr, Trp, His side chains
5. **Amide sources**: Both backbone CONH and side chain Asn/Gln
