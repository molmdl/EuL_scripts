#!/usr/bin/env python3
"""
Metal coordination geometry analysis for 9-coordinate Eu complexes.

Analyses SAP (square antiprism) and TSAP (twisted square antiprism)
geometries in MD trajectories OR single-frame GRO structures.

Section I   – Structure/trajectory loading, PBC handling, alignment
Section II  – Geometry analysis A: RMSD to ideal SAP/TSAP polyhedra
Section III – Geometry analysis B: N-C-C-N torsion angles

── Input modes ──────────────────────────────────────────────────────────────
Mode 1  TPR + XTC   (original, unchanged)
    python metal_geo_analysis_v2.py \\
        --system me_sssL_sap  --tpr me_sssL_sap/pr_0.tpr  --xtc me_sssL_sap/pr_0.xtc

Mode 2  Single GRO  (new – ligand-only OR receptor–ligand complex)
    python metal_geo_analysis_v2.py \\
        --system me_sssL_sap  --gro solv/me_sssL_sap.gro

Mode 3  Batch from list.txt  (new – runs a whole directory at once)
    python metal_geo_analysis_v2.py \\
        --list solv/list.txt  --dir solv  --outdir results_solv
    python metal_geo_analysis_v2.py \\
        --list com/list.txt   --dir com   --outdir results_com

    Each line of list.txt is a system name; the script looks for
    <dir>/<name>.gro and writes output to <outdir>/<name>/.

Modes can be freely mixed in one invocation, e.g.:
    python metal_geo_analysis_v2.py \\
        --system A --tpr A.tpr --xtc A.xtc \\
        --system B --gro solv/B.gro \\
        --list com/list.txt --dir com --outdir results_com

── Notes ────────────────────────────────────────────────────────────────────
• For GRO inputs the MOL residue is identified by `resname MOL`; the rest of
  the atoms (receptor, solvent) are ignored.  This works for both a
  ligand-only GRO (all atoms are MOL) and a complex GRO.
• Topology (atom indices for EU, coordinating atoms, torsions) is detected
  automatically from the first frame: two known topologies are supported —
  "SAP-type" (EU3 at mol index 54) and "TSAP-type" (Eu at mol index 84).
  See TOPOLOGIES dict below for full index definitions.
• Part A uses per-frame local idealisation: each frame keeps its own
  coordination-shell size/height parameters while only the angular SAP/TSAP
  arrangement is idealised before RMSD evaluation.
"""

import argparse
import os
import textwrap
import warnings
import csv

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import MDAnalysis as mda

warnings.filterwarnings('ignore')

# ═══════════════════════════════════════════════════════
#  TOPOLOGY DEFINITIONS
# ═══════════════════════════════════════════════════════
#
# Two distinct MOL atom orderings exist across the GRO files:
#
#  SAP-type  – used by all me_* and all phe_*sap* files
#              EU3 (type E) is at MOL index 54
#
#  TSAP-type – used by all phe_*tsap* files
#              Eu  (type E) is at MOL index 84
#              The chromophore arm is rotated relative to the ring, shifting
#              many atom indices.
#
# For each topology the dict stores:
#   metal_idx : int          – MOL-local index of the Eu atom
#   coord_idx : list[int]    – 9 MOL-local indices of coordinating atoms
#                              order: [O0,O1,O2, N_bottom×4, N_top, N_cap]
#                              i.e. local positions within the 9-array map to
#                              BOTTOM=[3,4,5,6]  TOP=[0,1,2,7]  CAP=[8]
#   ring_tors : list of (name,i,j,k,l)  MOL-local indices, N-C-C-N pattern
#   chrom_tor : (name,i,j,k,l)          MOL-local indices, chromophore torsion
#
# The BOTTOM/TOP/CAP sub-group indices into the 9-element coord sub-array
# are the SAME for both topologies (determined by position in coord_idx):
#   BOTTOM = [3,4,5,6]   ring N atoms
#   TOP    = [0,1,2,7]   O-arms + proximal chromophore N
#   CAP    = [8]         distal chromophore N

BOTTOM = [3, 4, 5, 6]
TOP    = [0, 1, 2, 7]
CAP    = [8]

TOPOLOGIES = {
    # ── SAP-type ──────────────────────────────────────────────────
    # EU3 at mol index 54  (atom name "EU3", type "E")
    # Coord: O0(0) O1(1) O2(2) N3(3) N4(4) N5(5) N6(6) N7(7) N63(63)
    # Ring torsions bridge adjacent ring-N pairs via two carbons each.
    # Chrom torsion: N3-C16-C23-N7 (N3→chromophore arm→proximal coord N)
    'sap': dict(
        metal_idx = 54,
        coord_idx = [0, 1, 2, 3, 4, 5, 6, 7, 63],
        ring_tors = [
            ('T1 N3-C8-C9-N4',    3,  8,  9,  4),
            ('T2 N4-C10-C11-N5',  4, 10, 11,  5),
            ('T3 N5-C12-C13-N6',  5, 12, 13,  6),
            ('T4 N6-C14-C15-N3',  6, 14, 15,  3),
        ],
        chrom_tor = ('Tc N3-C16-C23-N7', 3, 16, 23, 7),
    ),

    # ── TSAP-type ─────────────────────────────────────────────────
    # Eu at mol index 84  (atom name "Eu", type "E")
    # Coord: O0(0) O1(1) O2(2) N3(3) N4(4) N5(5) N6(6) N106(106) N107(107)
    #   → local positions 0-8 map identically to BOTTOM/TOP/CAP as SAP-type.
    # Ring torsions determined by bond topology (different C bridge atoms).
    # Chrom torsion: N6-C81-C89-N106 (ring N6→chromophore arm→proximal coord N)
    'tsap': dict(
        metal_idx = 84,
        coord_idx = [0, 1, 2, 3, 4, 5, 6, 106, 107],
        ring_tors = [
            ('T1 N3-C13-C16-N4',   3, 13, 16,  4),
            ('T2 N4-C19-C22-N5',   4, 19, 22,  5),
            ('T3 N5-C25-C28-N6',   5, 25, 28,  6),
            ('T4 N6-C31-C86-N3',   6, 31, 86,  3),
        ],
        chrom_tor = ('Tc N6-C81-C89-N106', 6, 81, 89, 106),
    ),
}

# Convenience alias kept for geometry helpers that use module-level constants.
# These are overridden per-system at runtime via the `topo` dict.
_DEFAULT_TOPO = TOPOLOGIES['sap']
METAL_IDX = _DEFAULT_TOPO['metal_idx']
COORD_IDX  = _DEFAULT_TOPO['coord_idx']
RING_TORS  = _DEFAULT_TOPO['ring_tors']
CHROM_TOR  = _DEFAULT_TOPO['chrom_tor']

# Target twist angles for ideal polyhedra
# Twist = mean_azimuth(TOP) − mean_azimuth(BOTTOM), wrapped to [−90, 90]
SAP_TWIST_DEG  = -45.0   # perfect square antiprism  (measured ≈ −40° in frame 0)
TSAP_TWIST_DEG =   0.0   # eclipsed / TSAP

# Plotting colours (SAP blue, TSAP purple)
C_SAP   = 'steelblue'
C_TSAP  = 'mediumorchid'
C_UNK   = 'lightgray'
# Cool teal-to-violet palette: T1 T2 T3 T4 Tc
TOR_COLS = ['#00b4d8', '#0077b6', '#48cae4', '#7b2d8b', '#c77dff']  # T1 T2 T3 T4 Tc

# ═══════════════════════════════════════════════════════
#  GEOMETRY / MATH HELPERS  (fully vectorised)
# ═══════════════════════════════════════════════════════

def calc_dihedral(pos, i, j, k, l):
    """Single dihedral angle in degrees [−180, 180]."""
    b1 = pos[j] - pos[i]
    b2 = pos[k] - pos[j]
    b3 = pos[l] - pos[k]
    n1 = np.cross(b1, b2);  nn1 = np.linalg.norm(n1)
    n2 = np.cross(b2, b3);  nn2 = np.linalg.norm(n2)
    if nn1 < 1e-9 or nn2 < 1e-9:
        return np.nan
    n1 /= nn1;  n2 /= nn2
    b2n = b2 / np.linalg.norm(b2)
    # IUPAC/GROMACS convention: positive = clockwise when viewed along j→k
    return float(np.degrees(np.arctan2(np.dot(np.cross(n1, n2), b2n), np.dot(n1, n2))))


def calc_dihedrals_batch(pos, ijkl):
    """
    Compute T dihedral angles in one vectorised pass.

    pos  : (N, 3) positions
    ijkl : (T, 4) integer index array, each row = [i, j, k, l]

    Returns (T,) array of angles in degrees [−180, 180].
    NaN is returned for degenerate cases.
    """
    i, j, k, l = ijkl[:, 0], ijkl[:, 1], ijkl[:, 2], ijkl[:, 3]
    b1 = pos[j] - pos[i]          # (T,3)
    b2 = pos[k] - pos[j]          # (T,3)
    b3 = pos[l] - pos[k]          # (T,3)
    n1 = np.cross(b1, b2)         # (T,3)
    n2 = np.cross(b2, b3)         # (T,3)
    nn1 = np.linalg.norm(n1, axis=1, keepdims=True)   # (T,1)
    nn2 = np.linalg.norm(n2, axis=1, keepdims=True)
    b2n = b2 / np.linalg.norm(b2, axis=1, keepdims=True)
    # avoid division by zero
    valid = (nn1[:, 0] > 1e-9) & (nn2[:, 0] > 1e-9)
    nn1 = np.where(nn1 > 1e-9, nn1, 1.0)
    nn2 = np.where(nn2 > 1e-9, nn2, 1.0)
    n1 /= nn1;  n2 /= nn2
    # IUPAC/GROMACS convention: sin term = (n1 × n2) · b2_hat
    n1xn2 = np.cross(n1, n2)       # (T,3)
    dot_sin = np.einsum('ti,ti->t', n1xn2, b2n)
    dot_n1_n2 = np.einsum('ti,ti->t', n1, n2)
    angles = np.degrees(np.arctan2(dot_sin, dot_n1_n2))
    angles[~valid] = np.nan
    return angles


def kabsch(P, Q):
    """
    Kabsch algorithm.  Find rotation R minimising ||P @ R.T − Q||.
    P, Q : (N, 3) arrays assumed centred (or weighted-centred).
    Returns (rmsd, R).
    """
    H = P.T @ Q
    U, _, Vt = np.linalg.svd(H)
    d = float(np.sign(np.linalg.det(Vt.T @ U.T)))
    R = Vt.T @ np.diag([1.0, 1.0, d]) @ U.T
    diff = P @ R.T - Q
    rmsd = float(np.sqrt(np.mean(np.einsum('ij,ij->i', diff, diff))))
    return rmsd, R


def _kabsch_batch(Ps, Q):
    """
    Batched Kabsch: compute RMSD for K permutations of P against fixed Q.

    Ps : (K, N, 3)  – K permuted copies of the ideal coordinates
    Q  : (N, 3)     – target (frame) coordinates, centred

    Returns rmsds : (K,) array of RMSD values, and Rs : (K, 3, 3) rotations.
    Uses a single np.linalg.svd call over the batch dimension for speed.
    """
    H = np.einsum('kni,nj->kij', Ps, Q)          # (K,3,3)
    U, S, Vt = np.linalg.svd(H)
    dets  = np.linalg.det(np.einsum('kij,kjl->kil',
                                     Vt.swapaxes(1,2), U.swapaxes(1,2)))
    signs = np.sign(dets)                         # (K,)
    diag  = np.zeros((len(Ps), 3, 3))
    diag[:, 0, 0] = 1.0;  diag[:, 1, 1] = 1.0;  diag[:, 2, 2] = signs
    Rs = np.einsum('kij,kjl,klm->kim',
                   Vt.swapaxes(1,2), diag, U.swapaxes(1,2))  # (K,3,3)
    rotated = np.einsum('kni,kji->knj', Ps, Rs)   # (K,N,3)
    diff    = rotated - Q[np.newaxis]              # (K,N,3)
    rmsds   = np.sqrt(np.mean(np.sum(diff**2, axis=2), axis=1))  # (K,)
    return rmsds, Rs


def align_pos(mob_pos, ref_pos, mask):
    """
    Kabsch-align mob_pos to ref_pos on atoms selected by boolean mask.
    Returns new positions for ALL atoms.
    """
    mob_h = mob_pos[mask];  ref_h = ref_pos[mask]
    cm = mob_h.mean(0);     cr = ref_h.mean(0)
    _, R = kabsch(mob_h - cm, ref_h - cr)
    return (mob_pos - cm) @ R.T + cr


def rodrigues(axis, deg):
    """Rotation matrix about unit axis by deg degrees (Rodrigues)."""
    a = np.radians(deg)
    c, s = np.cos(a), np.sin(a)
    kx, ky, kz = axis
    K = np.array([[0, -kz, ky], [kz, 0, -kx], [-ky, kx, 0]])
    return c * np.eye(3) + s * K + (1 - c) * np.outer(axis, axis)


def circ_mean_deg(angles):
    """Circular mean of a list of angles in degrees."""
    rad = np.radians(angles)
    return float(np.degrees(np.arctan2(np.mean(np.sin(rad)),
                                        np.mean(np.cos(rad)))))


def get_c4_frame(coord9):
    """
    Compute orthonormal frame (x, y, z) where z is the C4-like axis:
    normal to the BOTTOM (ring-N) best-fit plane, pointing toward CAP.

    coord9 : (9, 3) – Eu-centred coordinating atoms.
    Returns (x_hat, y_hat, z_hat), each shape (3,).
    """
    ring_N   = coord9[BOTTOM]                   # (4, 3)
    centroid = ring_N.mean(0)
    _, _, Vt = np.linalg.svd(ring_N - centroid)
    z = Vt[-1]
    if np.dot(z, coord9[CAP[0]]) < 0:
        z = -z

    # Use the CAP in-plane projection to define a physically meaningful
    # azimuth origin whenever possible; otherwise fall back to an arbitrary
    # perpendicular direction.
    cap_vec  = coord9[CAP[0]]
    cap_proj = cap_vec - np.dot(cap_vec, z) * z
    cap_norm = np.linalg.norm(cap_proj)
    if cap_norm > 1e-8:
        x = cap_proj / cap_norm
    else:
        tmp = np.array([1.0, 0.0, 0.0]) if abs(z[0]) < 0.9 else np.array([0.0, 1.0, 0.0])
        x = np.cross(z, tmp)
        x /= np.linalg.norm(x)
    y = np.cross(z, x)
    y /= np.linalg.norm(y)
    return x, y, z


def _azimuths(vecs, x, y, z):
    """
    Return azimuthal angles (degrees) for each row-vector in vecs.

    vecs : (M, 3)
    x, y, z : (3,) orthonormal C4 frame axes
    """
    # project onto the xy plane of the C4 frame
    vp = vecs - (vecs @ z)[:, np.newaxis] * z   # (M,3)
    return np.degrees(np.arctan2(vp @ y, vp @ x))  # (M,)


def get_twist_angle(coord9):
    """
    Twist angle between TOP and BOTTOM layers, wrapped to [−90, 90].

    Uses sorted-azimuth inter-layer offset to avoid the circular-mean
    singularity that occurs when atoms lie exactly at 0/90/180/270 deg.

    Fully vectorised: no Python loops over atoms.
    """
    x, y, z = get_c4_frame(coord9)
    phi_bot = np.sort(_azimuths(coord9[BOTTOM], x, y, z))   # (4,)
    phi_top = np.sort(_azimuths(coord9[TOP],    x, y, z))   # (4,)
    offsets = (phi_top - phi_bot + 90) % 180 - 90           # (4,)
    return float(offsets.mean())


def build_local_ideal_directions(coord9):
    """
    Build frame-specific ideal SAP and TSAP *direction* templates (unit vectors).

    Implements the same approach as UCSF Chimera Metal Geometry 'Distance RMSD':
    only the angular (in-plane) arrangement is idealised.  Per-atom radial
    scaling is applied in min_rmsd_cyclic using the real metal-ligand distances,
    so no radial spread contributes to the RMSD — only angular distortion does.

    The C4 frame, layer heights, and layer-averaged in-plane radii are all taken
    from the instantaneous frame so the 3D polar angles are preserved correctly.
    The vectors are normalised to unit length after construction.

    coord9 : (9, 3) Eu-centred coordinating-atom coordinates for one frame.
    Returns (dirs_sap, dirs_tsap) – both (9, 3) arrays of unit vectors.
    """
    x, y, z = get_c4_frame(coord9)
    zh    = coord9 @ z
    vp    = coord9 - zh[:, np.newaxis] * z
    radii = np.linalg.norm(vp, axis=1)

    r_bot = float(radii[BOTTOM].mean())
    z_bot = float(zh[BOTTOM].mean())
    r_top = float(radii[TOP].mean())
    z_top = float(zh[TOP].mean())
    r_cap = float(radii[CAP[0]])
    z_cap = float(zh[CAP[0]])

    def make_directions(twist_deg):
        out = np.zeros((9, 3))
        phi_b = np.radians(np.arange(4) * 90.0)
        phi_t = phi_b + np.radians(twist_deg)
        # Use layer-averaged in-plane radii to preserve correct 3D polar angles
        out[BOTTOM] = (
            r_bot * np.cos(phi_b)[:, np.newaxis] * x[np.newaxis, :]
            + r_bot * np.sin(phi_b)[:, np.newaxis] * y[np.newaxis, :]
            + z_bot * z[np.newaxis, :]
        )
        out[TOP] = (
            r_top * np.cos(phi_t)[:, np.newaxis] * x[np.newaxis, :]
            + r_top * np.sin(phi_t)[:, np.newaxis] * y[np.newaxis, :]
            + z_top * z[np.newaxis, :]
        )
        out[CAP[0]] = r_cap * x + z_cap * z
        norms = np.linalg.norm(out, axis=1, keepdims=True)
        norms = np.where(norms > 1e-12, norms, 1.0)
        return out / norms   # (9, 3) unit vectors

    return make_directions(SAP_TWIST_DEG), make_directions(TSAP_TWIST_DEG)


# Pre-compute all 16 cyclic permutation index arrays (4 BOTTOM x 4 TOP)
_CYCLIC_PERMS = []
for _db in range(4):
    for _dt in range(4):
        _p = list(range(9))
        for _ii, _orig in enumerate(BOTTOM):
            _p[_orig] = BOTTOM[(_ii + _db) % 4]
        for _ii, _orig in enumerate(TOP):
            _p[_orig] = TOP[(_ii + _dt) % 4]
        _CYCLIC_PERMS.append(_p)
_CYCLIC_PERMS = np.array(_CYCLIC_PERMS)   # (16, 9)


def min_rmsd_cyclic(coord9, dirs9, return_transform=False):
    """
    Chimera-style RMSD to ideal geometry, minimised over all 16 cyclic
    permutations of BOTTOM (4 cyclic) and TOP (4 cyclic) subgroups.

    Algorithm (matches UCSF Chimera Metal Geometry 'Distance RMSD'):
      1. Kabsch-align ideal unit directions (permuted) to real unit directions.
      2. Place each ideal atom along the rotated direction at the real atom's
         distance from the metal (per-atom radial scaling).
      3. RMSD between scaled ideal positions and real positions.

    This removes all radial deviation from the RMSD; only the angular distortion
    from the ideal SAP/TSAP pattern is measured.

    coord9 : (9, 3) real coordinates, centred on Eu.
    dirs9  : (9, 3) unit direction vectors for the ideal geometry
             (from build_local_ideal_directions).
    """
    real_r = np.linalg.norm(coord9, axis=1)                   # (9,)
    Q_unit = coord9 / np.where(real_r > 1e-12, real_r, 1.0)[:, np.newaxis]  # (9,3)

    Ps_unit = dirs9[_CYCLIC_PERMS]                             # (16, 9, 3) permuted unit dirs

    # Step 1: Kabsch rotation on unit-sphere coordinates
    _, Rs = _kabsch_batch(Ps_unit, Q_unit)                     # Rs: (16, 3, 3)

    # Step 2: rotate ideal unit dirs, then scale by real distances
    rotated_unit  = np.einsum('kni,kji->knj', Ps_unit, Rs)    # (16, 9, 3)
    ideal_placed  = rotated_unit * real_r[np.newaxis, :, np.newaxis]  # (16, 9, 3)

    # Step 3: RMSD
    diff   = ideal_placed - coord9[np.newaxis]                 # (16, 9, 3)
    rmsds  = np.sqrt(np.mean(np.sum(diff**2, axis=2), axis=1))  # (16,)

    best_k    = int(np.argmin(rmsds))
    best_rmsd = float(rmsds[best_k])
    if return_transform:
        return best_rmsd, ideal_placed[best_k]
    return best_rmsd


def classify_geom(ring_tors, chrom_tor, thresh=10.0):
    """
    SAP  if sign(mean_ring) ≠ sign(chrom_tor)
    TSAP if sign(mean_ring) = sign(chrom_tor)
    UNK  if either value is within thresh of 0
    """
    ring_arr = np.asarray(ring_tors, dtype=float)
    valid    = ring_arr[~np.isnan(ring_arr)]
    if len(valid) == 0 or np.isnan(chrom_tor):
        return 'UNK'
    mean_ring = valid.mean()
    if abs(mean_ring) < thresh or abs(chrom_tor) < thresh:
        return 'UNK'
    return 'SAP' if np.sign(mean_ring) != np.sign(chrom_tor) else 'TSAP'


# ═══════════════════════════════════════════════════════
#  STRUCTURE / TRAJECTORY LOADING  (Section I)
# ═══════════════════════════════════════════════════════

EU_MASS = 151.964   # standard atomic weight of Eu (IUPAC 2021)

def load_mol_gro(gro):
    """
    Load a single-frame GRO file and return (u, mol_ag).

    Works for both ligand-only GROs (all atoms are MOL) and complex GROs
    that also contain receptor/solvent atoms — only the MOL residue is
    returned and used for analysis.

    Mass correction
    ---------------
    MDAnalysis guesses masses from atom type when reading a bare GRO file.
    The Eu atom has type "E", which maps to mass 0.0 by default.  This
    function patches any MOL atom whose mass is 0.0 **and** whose type
    starts with "e" (case-insensitive) to EU_MASS (151.964 Da).

    Returns
    -------
    u      : MDAnalysis Universe  (1-frame trajectory)
    mol_ag : AtomGroup for resname MOL
    """
    u = mda.Universe(gro)
    mol_ag = u.select_atoms('resname MOL')
    if mol_ag.n_atoms == 0:
        raise ValueError(
            f"No atoms with resname MOL found in {gro}.\n"
            "Check that the ligand residue name is MOL in the GRO file."
        )
    # Patch EU mass: MDAnalysis gives type "E" a mass of 0.0
    for atom in mol_ag.atoms:
        if atom.mass == 0.0 and atom.type.lower().startswith('e'):
            atom.mass = EU_MASS
    return u, mol_ag


def _find_matching_gro(tpr_path, n_atoms):
    """
    Search for a GRO file with exactly n_atoms atoms near the given TPR.

    Looks in:
      1. The directory containing the TPR  (e.g. fp/)
      2. The parent of that directory       (e.g. the system root)

    Returns the path of the first matching GRO, or None.
    """
    import glob
    tpr_dir = os.path.dirname(os.path.abspath(tpr_path))
    search_dirs = [tpr_dir, os.path.dirname(tpr_dir)]
    for sdir in search_dirs:
        for gro in sorted(glob.glob(os.path.join(sdir, '*.gro'))):
            try:
                with open(gro) as fh:
                    fh.readline()           # title line
                    n = int(fh.readline())  # atom count line
                if n == n_atoms:
                    return gro
            except Exception:
                continue
    return None


def load_mol_universe(tpr, xtc):
    """
    Load the MOL residue from tpr+xtc, handling three cases transparently:

    Case A – stripped XTC (atoms == MOL only):
        Build a MOL-only universe via Merge + load_new.  mol_ag = u.atoms.

    Case B – full XTC (atoms == whole simulation box):
        Load Universe(tpr, xtc) directly.  mol_ag is the persistent
        'resname MOL' AtomGroup inside that universe.

    Case C – stripped XTC (atoms == complex minus solvent, i.e. between
              MOL count and full-system count):
        Find a GRO file in the same or parent directory whose atom count
        matches the XTC.  Load Universe(gro, xtc) to get per-frame
        positions, then copy TPR masses onto the MOL atoms so that Eu
        gets its correct mass (~152 Da) rather than the 0.0 that
        MDAnalysis guesses from the bare GRO atom type "E".

    Returns (u, mol_ag) where:
        u       – Universe whose .trajectory is iterated in the main loop
        mol_ag  – AtomGroup of exactly the MOL atoms (positions updated
                  automatically each frame because it is a live selection)
    """
    import MDAnalysis.coordinates.XTC as _XTC
    xtc_natoms = _XTC.XTCReader(xtc).n_atoms

    u_tpr   = mda.Universe(tpr)
    mol_tpr = u_tpr.select_atoms('resname MOL')
    n_mol   = mol_tpr.n_atoms

    if xtc_natoms == n_mol:
        # Case A: stripped XTC – original Merge approach
        u      = mda.Merge(mol_tpr)
        u.load_new(xtc)
        mol_ag = u.atoms
    elif xtc_natoms == u_tpr.atoms.n_atoms:
        # Case B: full XTC – load directly, select MOL as a live AtomGroup
        u      = mda.Universe(tpr, xtc)
        mol_ag = u.select_atoms('resname MOL')
    else:
        # Case C: stripped XTC that is a subset of the full system
        # (e.g. complex without bulk solvent).  The TPR atom count does not
        # match the XTC, so we locate a companion GRO with the right count
        # and use it as the topology for loading the trajectory.
        # We then transplant TPR masses onto the MOL atoms so Eu mass is
        # correct (TPR stores ~152 Da; bare GRO gives 0.0 for type "E").
        gro_path = _find_matching_gro(tpr, xtc_natoms)
        if gro_path is None:
            raise ValueError(
                f"XTC has {xtc_natoms} atoms, but TPR has {u_tpr.atoms.n_atoms} "
                f"total and {n_mol} MOL atoms.  Could not find a GRO file with "
                f"{xtc_natoms} atoms near {tpr} to use as topology."
            )
        print(f'  [load] Case C (stripped complex XTC): using {gro_path} as topology')
        u      = mda.Universe(gro_path, xtc)
        mol_ag = u.select_atoms('resname MOL')
        if mol_ag.n_atoms == 0:
            raise ValueError(
                f"No resname MOL atoms found in {gro_path}."
            )
        # Copy masses from TPR so Eu gets ~152 Da instead of 0.0
        if mol_ag.n_atoms == n_mol:
            mol_ag.masses = mol_tpr.masses
        else:
            # Atom count mismatch between TPR MOL and GRO MOL — patch only Eu
            for atom in mol_ag.atoms:
                if atom.mass == 0.0 and atom.type.lower().startswith('e'):
                    atom.mass = EU_MASS

    return u, mol_ag


def heavy_mask(mol_ag):
    """Boolean mask over MOL atoms: True for non-hydrogen atoms (mass > 1.5 u)."""
    return mol_ag.masses > 1.5


def make_mol_whole(u, mol_ag=None):
    """
    Apply MDAnalysis unwrap to make the covalent complex whole across PBC.

    If mol_ag is supplied (the MOL AtomGroup), only that group is unwrapped —
    this avoids the huge cost of unwrapping the full solvent box when the
    trajectory contains all simulation atoms.

    Silently skips if box info is unavailable.
    """
    try:
        ts = u.trajectory.ts
        if ts.dimensions is not None and ts.dimensions[0] > 0:
            ag = mol_ag if mol_ag is not None else u.atoms
            ag.unwrap(compound='residues')
    except Exception:
        pass


# ═══════════════════════════════════════════════════════
#  PDB WRITING HELPER
# ═══════════════════════════════════════════════════════

def write_pdb_with_dummies(mol_atoms, pos, ideal_placed, eu_pos,
                           outfile, frame_idx, label):
    """
    Write a PDB for one frame with the actual complex (chain A) and
    9 dummy atoms representing the ideal coordination polyhedron (chain B).

    mol_atoms   : MDAnalysis AtomGroup  (for name / resname / type metadata)
    pos         : (N_mol, 3) aligned positions for this frame
    ideal_placed: (9, 3) ideal coord atom positions in the SAME frame (Eu-centred)
    eu_pos      : (3,) absolute Eu position
    """
    ELEM = {'c': 'C', 'n': 'N', 'o': 'O', 'h': 'H',
            'e': 'Eu', 'eu': 'Eu', 's': 'S', 'p': 'P'}
    DUM_NAMES = ['DO0 ', 'DO1 ', 'DO2 ',
                 'DN3 ', 'DN4 ', 'DN5 ', 'DN6 ',
                 'DN7 ', 'DN63']

    def elem_symbol(atype):
        t = atype.lower()
        for k in sorted(ELEM, key=len, reverse=True):
            if t.startswith(k):
                return ELEM[k]
        return atype[0].upper()

    with open(outfile, 'w') as f:
        f.write(f'REMARK  frame={frame_idx}  label={label}\n')
        f.write( 'REMARK  Chain A = complex  |  Chain B = ideal coord polyhedron\n')

        # Complex atoms
        for i, atom in enumerate(mol_atoms):
            p    = pos[i]
            name = atom.name[:4].ljust(4)
            elem = elem_symbol(atom.type)
            f.write(
                f'ATOM  {i+1:5d} {name} {atom.resname:<3s} A{1:4d}    '
                f'{p[0]:8.3f}{p[1]:8.3f}{p[2]:8.3f}  1.00  0.00'
                f'          {elem:>2s}\n'
            )

        # Dummy atoms for ideal geometry
        abs_ideal = eu_pos + ideal_placed   # absolute Å
        for j, (dname, dp) in enumerate(zip(DUM_NAMES, abs_ideal)):
            serial = len(mol_atoms) + j + 1
            f.write(
                f'HETATM{serial:5d} {dname} DUM B{2:4d}    '
                f'{dp[0]:8.3f}{dp[1]:8.3f}{dp[2]:8.3f}  1.00  0.00'
                f'           X\n'
            )

        f.write('END\n')
    print(f'    Wrote {outfile}')


# ═══════════════════════════════════════════════════════
#  PEAK FINDING HELPER
# ═══════════════════════════════════════════════════════

def histogram_peaks(values, n_bins=50):
    """Return list of bin-centre positions for histogram peaks."""
    vals = [v for v in values if not np.isnan(v)]
    if not vals:
        return [0.0]
    counts, edges = np.histogram(vals, bins=n_bins)
    centres = 0.5 * (edges[:-1] + edges[1:])
    try:
        from scipy.signal import find_peaks as sp_peaks
        idx, _ = sp_peaks(counts,
                          height=max(counts) * 0.15,
                          prominence=max(counts) * 0.05)
    except Exception:
        idx = []
    if len(idx) == 0:
        idx = [int(np.argmax(counts))]
    return [float(centres[i]) for i in idx]


def nearest_frame(series, target):
    """Index of the element in series closest to target (ignores NaN)."""
    arr  = np.asarray(series, dtype=float)
    diffs = np.where(np.isnan(arr), np.inf, np.abs(arr - target))
    return int(np.argmin(diffs))


# ═══════════════════════════════════════════════════════
#  SINGLE-PASS DATA COLLECTION
# ═══════════════════════════════════════════════════════

# ═══════════════════════════════════════════════════════
#  TOPOLOGY DETECTION
# ═══════════════════════════════════════════════════════

def detect_topology(mol_ag):
    """
    Identify which topology ('sap' or 'tsap') applies to this MOL AtomGroup
    by inspecting the atom type at the expected EU positions.

    Detection logic (checked in priority order):
      1. If mol index 84 has element type 'E' (Eu)  → 'tsap'
      2. If mol index 54 has element type 'E' (EU3) → 'sap'
      3. Fallback: search for any atom with type 'E' and pick the closest
         topology by index.

    Raises ValueError if no Eu atom is found.
    Returns the topology key ('sap' or 'tsap') and the topology dict.
    """
    atoms = mol_ag.atoms
    n = mol_ag.n_atoms

    def _is_eu(idx):
        return idx < n and atoms[idx].type.lower().startswith('e')

    if _is_eu(84):
        key = 'tsap'
    elif _is_eu(54):
        key = 'sap'
    else:
        # Fallback: find first atom whose type starts with 'e'
        eu_candidates = [a.index - atoms[0].index
                         for a in atoms if a.type.lower().startswith('e')]
        if not eu_candidates:
            raise ValueError(
                "Cannot detect topology: no Eu atom found in MOL residue. "
                "Check that the metal atom type starts with 'E' in the input."
            )
        eu_idx = eu_candidates[0]
        # Pick whichever known metal_idx is closest
        dist_sap  = abs(eu_idx - TOPOLOGIES['sap']['metal_idx'])
        dist_tsap = abs(eu_idx - TOPOLOGIES['tsap']['metal_idx'])
        key = 'sap' if dist_sap <= dist_tsap else 'tsap'
        print(f'  [topology] WARNING: EU not at expected index; '
              f'found at {eu_idx}, using "{key}" topology.')

    topo = TOPOLOGIES[key]
    # Sanity-check that coordinating atom indices are within range
    max_idx = max(topo['coord_idx'])
    if max_idx >= n:
        raise ValueError(
            f'Topology "{key}" requires MOL atom index {max_idx}, '
            f'but MOL only has {n} atoms.'
        )
    return key, topo


def _make_torsion_ijkl(topo):
    """Build (5, 4) int array of torsion indices from a topology dict."""
    ct = topo['chrom_tor']
    return np.array(
        [[i, j, k, l] for _, i, j, k, l in topo['ring_tors']]
        + [[ct[1], ct[2], ct[3], ct[4]]],
        dtype=int
    )   # (5, 4) – rows: T1 T2 T3 T4 Tc


def collect_trajectory_data(u, mol_ag, ref_pos, hmask, topo):
    """
    Iterate over every frame, compute Part-A and Part-B quantities.

    Parameters
    ----------
    u       : MDAnalysis Universe
    mol_ag  : AtomGroup for the MOL residue
    ref_pos : (N_mol, 3) reference positions for Kabsch alignment
    hmask   : boolean mask for heavy atoms (used in alignment)
    topo    : topology dict from TOPOLOGIES (contains metal_idx, coord_idx,
              ring_tors, chrom_tor)

    Returns a list of per-frame dicts.
    """
    metal_idx   = topo['metal_idx']
    coord_idx   = topo['coord_idx']
    torsion_ijkl = _make_torsion_ijkl(topo)

    records = []
    for ts in u.trajectory:
        make_mol_whole(u, mol_ag)
        raw_pos = mol_ag.positions.copy()            # (N_mol, 3) unwrapped
        aln_pos = align_pos(raw_pos, ref_pos, hmask) # (N_mol, 3) aligned

        eu    = aln_pos[metal_idx]
        coord = aln_pos[coord_idx] - eu              # (9, 3) centred on Eu

        # ── Part A : per-frame local idealisation + batched RMSD ────
        dirs_sap, dirs_tsap = build_local_ideal_directions(coord)
        rmsd_sap  = min_rmsd_cyclic(coord, dirs_sap)
        rmsd_tsap = min_rmsd_cyclic(coord, dirs_tsap)

        # ── Part B : all 5 torsions in one vectorised call ──────────
        all_tors = calc_dihedrals_batch(raw_pos, torsion_ijkl)  # (5,)
        ring_t   = list(all_tors[:4])
        tc       = float(all_tors[4])
        geom     = classify_geom(ring_t, tc)

        records.append({
            'frame'     : int(ts.frame),
            'time'      : float(ts.time),
            'pos'       : aln_pos.copy(),
            'coord9'    : coord.copy(),
            'eu'        : eu.copy(),
            'dirs_sap'  : dirs_sap.copy(),
            'dirs_tsap' : dirs_tsap.copy(),
            'rmsd_sap'  : rmsd_sap,
            'rmsd_tsap' : rmsd_tsap,
            'ring_tors' : ring_t,
            'chrom_tor' : tc,
            'geom'      : geom,
        })

    return records


def run_part_A(records, mol_atoms, outdir, sysname):
    os.makedirs(outdir, exist_ok=True)

    times     = np.array([r['time']      for r in records])
    rs_sap    = np.array([r['rmsd_sap']  for r in records])
    rs_tsap   = np.array([r['rmsd_tsap'] for r in records])

    # ── Time-series CSV ────────────────────────────────────────────
    ts_csv = os.path.join(outdir, 'A_rmsd_timeseries.csv')
    with open(ts_csv, 'w', newline='') as f:
        w = csv.writer(f)
        w.writerow(['frame', 'time_ps', 'RMSD_to_ideal_SAP_A',
                    'RMSD_to_ideal_TSAP_A'])
        for r, rs, rt in zip(records, rs_sap, rs_tsap):
            w.writerow([r['frame'], f"{r['time']:.1f}",
                        f'{rs:.4f}', f'{rt:.4f}'])
    print(f'  [A] {ts_csv}')

    # ── Time-series PNG ────────────────────────────────────────────
    fig, ax = plt.subplots(figsize=(10, 3.8))
    ax.plot(times / 1000, rs_sap,  color=C_SAP,  lw=1.2, alpha=0.85,
            label='RMSD → ideal SAP')
    ax.plot(times / 1000, rs_tsap, color=C_TSAP, lw=1.2, alpha=0.85,
            label='RMSD → ideal TSAP')
    ax.set_xlabel('Time (ns)', fontsize=11)
    ax.set_ylabel('RMSD  (Å)', fontsize=11)
    ax.set_title(f'{sysname}  –  Coordination-sphere RMSD to ideal geometry',
                 fontsize=11)
    ax.legend(framealpha=0.8)
    plt.tight_layout()
    fig.savefig(os.path.join(outdir, 'A_rmsd_timeseries.png'), dpi=150)
    plt.close(fig)
    print(f'  [A] {os.path.join(outdir, "A_rmsd_timeseries.png")}')

    # ── Histogram CSV + PNG (fraction) ────────────────────────────
    n_bins = 50
    cnt_sap,  edg_sap  = np.histogram(rs_sap,  bins=n_bins)
    cnt_tsap, edg_tsap = np.histogram(rs_tsap, bins=n_bins)
    ctr_sap  = 0.5 * (edg_sap[:-1]  + edg_sap[1:])
    ctr_tsap = 0.5 * (edg_tsap[:-1] + edg_tsap[1:])
    nf = len(rs_sap)
    frac_sap  = cnt_sap  / nf
    frac_tsap = cnt_tsap / nf

    hist_csv = os.path.join(outdir, 'A_rmsd_histogram.csv')
    with open(hist_csv, 'w', newline='') as f:
        w = csv.writer(f)
        w.writerow(['bin_SAP_A', 'fraction_SAP',
                    'bin_TSAP_A', 'fraction_TSAP'])
        for row in zip(ctr_sap, frac_sap, ctr_tsap, frac_tsap):
            w.writerow([f'{x:.4f}' for x in row])
    print(f'  [A] {hist_csv}')

    fig, ax = plt.subplots(figsize=(7, 4))
    ax.plot(ctr_sap,  frac_sap,  color=C_SAP,  lw=2,
            label='RMSD → ideal SAP')
    ax.plot(ctr_tsap, frac_tsap, color=C_TSAP, lw=2,
            label='RMSD → ideal TSAP')
    ax.fill_between(ctr_sap,  frac_sap,  alpha=0.25, color=C_SAP)
    ax.fill_between(ctr_tsap, frac_tsap, alpha=0.25, color=C_TSAP)
    ax.set_xlabel('RMSD  (Å)', fontsize=11)
    ax.set_ylabel('Fraction of frames', fontsize=11)
    ax.set_title(f'{sysname}  –  RMSD histogram', fontsize=11)
    ax.legend(framealpha=0.8)
    plt.tight_layout()
    fig.savefig(os.path.join(outdir, 'A_rmsd_histogram.png'), dpi=150)
    plt.close(fig)
    print(f'  [A] {os.path.join(outdir, "A_rmsd_histogram.png")}')

    # ── Representative frames + PDB ───────────────────────────────
    peaks_sap  = histogram_peaks(rs_sap.tolist())
    peaks_tsap = histogram_peaks(rs_tsap.tolist())

    rep_jobs = []
    for pk in peaks_sap:
        fi = nearest_frame(rs_sap.tolist(), pk)
        rep_jobs.append((fi, 'SAP', pk))
    for pk in peaks_tsap:
        fi = nearest_frame(rs_tsap.tolist(), pk)
        rep_jobs.append((fi, 'TSAP', pk))

    written = set()
    for fi, geom, pk in rep_jobs:
        key = (fi, geom)
        if key in written:
            continue
        written.add(key)

        r = records[fi]
        dirs_local = r['dirs_sap'] if geom == 'SAP' else r['dirs_tsap']
        _, ideal_placed = min_rmsd_cyclic(r['coord9'], dirs_local,
                                          return_transform=True)
        label   = f'{geom}_peak_{pk:.2f}A_frame{fi}'
        pdbfile = os.path.join(outdir, f'A_rep_{geom}_{fi:04d}.pdb')
        write_pdb_with_dummies(mol_atoms, r['pos'], ideal_placed,
                                r['eu'], pdbfile, fi, label)

    return rs_sap, rs_tsap, times


# ═══════════════════════════════════════════════════════
#  PART B  – Torsion analysis
# ═══════════════════════════════════════════════════════

def run_part_B(records, outdir, sysname):
    os.makedirs(outdir, exist_ok=True)

    times  = np.array([r['time']      for r in records])
    geoms  = [r['geom']               for r in records]
    tnames = [n for n, *_ in RING_TORS] + [CHROM_TOR[0]]
    # shape (n_tors, n_frames)
    tor_mat = np.array(
        [r['ring_tors'] + [r['chrom_tor']] for r in records]
    ).T

    cls_color = {'SAP': C_SAP, 'TSAP': C_TSAP, 'UNK': C_UNK}

    # ── Torsion time-series CSV ────────────────────────────────────
    ts_csv = os.path.join(outdir, 'B_torsions_timeseries.csv')
    with open(ts_csv, 'w', newline='') as f:
        w = csv.writer(f)
        w.writerow(['frame', 'time_ps'] + tnames + ['geom'])
        for i, r in enumerate(records):
            vals = [f'{v:.2f}' for v in tor_mat[:, i]]
            w.writerow([r['frame'], f"{r['time']:.1f}"] + vals + [geoms[i]])
    print(f'  [B] {ts_csv}')

    # ── Torsion time-series PNG  (all torsions, one overlaid plot) ──
    fig, ax = plt.subplots(figsize=(11, 4.5))
    for name, col, vals in zip(tnames, TOR_COLS, tor_mat):
        ax.plot(times / 1000, vals, color=col, lw=1.4, alpha=0.85, label=name)
    ax.axhline(0, color='k', lw=0.6, ls='--', alpha=0.5)
    ax.set_xlabel('Time (ns)', fontsize=11)
    ax.set_ylabel('Torsion angle (°)', fontsize=11)
    ax.set_ylim(-185, 185)
    ax.set_yticks([-180, -90, 0, 90, 180])
    ax.legend(fontsize=9, framealpha=0.85, loc='upper right',
              ncol=len(tnames))
    ax.set_title(f'{sysname}  –  Torsion time series', fontsize=11)
    plt.tight_layout()
    fig.savefig(os.path.join(outdir, 'B_torsions_timeseries.png'), dpi=150)
    plt.close(fig)
    print(f'  [B] {os.path.join(outdir, "B_torsions_timeseries.png")}')

    # ── Torsion histogram CSV + PNG (SAP blue, TSAP purple) ───────
    n_bins = 60
    bin_edges = np.linspace(-180, 180, n_bins + 1)
    bin_ctrs  = 0.5 * (bin_edges[:-1] + bin_edges[1:])

    sap_mask  = np.array([g == 'SAP'  for g in geoms])
    tsap_mask = np.array([g == 'TSAP' for g in geoms])

    nf_total = len(geoms)
    n_sap_f  = int(np.sum(sap_mask))
    n_tsap_f = int(np.sum(tsap_mask))

    hist_csv = os.path.join(outdir, 'B_torsion_histograms.csv')
    header   = ['bin_deg']
    for n in tnames:
        header += [f'{n}_SAP', f'{n}_TSAP', f'{n}_all']

    all_hists = {}   # name → (frac_sap, frac_tsap, frac_all)
    for name, vals_row in zip(tnames, tor_mat):
        def cnt(mask):
            v = vals_row[mask]
            return np.histogram(v[~np.isnan(v)], bins=bin_edges)[0]
        cs  = cnt(sap_mask)
        ct  = cnt(tsap_mask)
        ca  = cnt(np.ones(nf_total, dtype=bool))
        # normalise each curve by its own frame count so fractions sum to 1
        dbin = bin_edges[1] - bin_edges[0]
        fs  = cs  / max(n_sap_f,  1) / dbin if n_sap_f  > 0 else cs * 0.0
        ft  = ct  / max(n_tsap_f, 1) / dbin if n_tsap_f > 0 else ct * 0.0
        fa  = ca  / nf_total          / dbin
        all_hists[name] = (fs, ft, fa)

    with open(hist_csv, 'w', newline='') as f:
        w = csv.writer(f)
        w.writerow(header)
        for j, bc in enumerate(bin_ctrs):
            row = [f'{bc:.2f}']
            for name in tnames:
                fs, ft, fa = all_hists[name]
                row += [f'{fs[j]:.6f}', f'{ft[j]:.6f}', f'{fa[j]:.6f}']
            w.writerow(row)
    print(f'  [B] {hist_csv}')

    # PNG: all torsions overlaid in one plot (solid = all frames)
    fig, ax = plt.subplots(figsize=(8, 4.5))
    for name, col in zip(tnames, TOR_COLS):
        fa = all_hists[name][2]
        ax.plot(bin_ctrs, fa, color=col, lw=1.8, alpha=0.90, label=name)
        ax.fill_between(bin_ctrs, fa, alpha=0.12, color=col)
    ax.axvline(0, color='k', lw=0.6, ls='--', alpha=0.5)
    ax.set_xlabel('Torsion angle (°)', fontsize=11)
    ax.set_ylabel('Probability density (deg⁻¹)', fontsize=11)
    ax.set_xlim(-180, 180)
    ax.set_xticks([-180, -90, 0, 90, 180])
    ax.legend(fontsize=9, framealpha=0.85, ncol=len(tnames))
    ax.set_title(f'{sysname}  –  Torsion angle distributions', fontsize=11)
    plt.tight_layout()
    fig.savefig(os.path.join(outdir, 'B_torsion_histograms.png'), dpi=150)
    plt.close(fig)
    print(f'  [B] {os.path.join(outdir, "B_torsion_histograms.png")}')

    # ── Classification time-series CSV + PNG ──────────────────────
    cls_csv = os.path.join(outdir, 'B_classification.csv')
    with open(cls_csv, 'w', newline='') as f:
        w = csv.writer(f)
        w.writerow(['frame', 'time_ps', 'geometry'])
        for r, g in zip(records, geoms):
            w.writerow([r['frame'], f"{r['time']:.1f}", g])
    print(f'  [B] {cls_csv}')

    # Classification panel: top = torsion scatter, bottom = geom bar
    mean_ring = np.array([np.nanmean(r['ring_tors']) for r in records])
    chrom_t   = np.array([r['chrom_tor']             for r in records])
    pt_colors = [cls_color[g] for g in geoms]

    cls_int   = {'SAP': 1, 'TSAP': -1, 'UNK': 0}
    cls_vals  = np.array([cls_int[g] for g in geoms])

    n_sap  = geoms.count('SAP')
    n_tsap = geoms.count('TSAP')
    n_unk  = geoms.count('UNK')
    nf     = len(geoms)

    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 5.5), sharex=True,
                                    gridspec_kw={'height_ratios': [2, 1]})

    ax1.scatter(times / 1000, mean_ring, c=pt_colors, s=12, alpha=0.8,
                linewidths=0, label='mean ring tors.')
    ax1.scatter(times / 1000, chrom_t,  c=pt_colors, s=12, alpha=0.8,
                linewidths=0, marker='^', label='chrom. tors.')
    ax1.axhline(0, color='k', lw=0.5, ls='--')
    ax1.set_ylabel('Torsion (°)', fontsize=10)
    ax1.set_ylim(-185, 185)
    ax1.set_yticks([-180, -90, 0, 90, 180])
    ax1.legend(fontsize=8, loc='upper right')

    ax2.scatter(times / 1000, cls_vals, c=pt_colors, s=18, alpha=0.85,
                linewidths=0)
    ax2.set_yticks([-1, 0, 1])
    ax2.set_yticklabels(['TSAP', 'UNK', 'SAP'], fontsize=9)
    ax2.set_xlabel('Time (ns)', fontsize=10)
    ax2.set_ylabel('Geometry', fontsize=10)

    fig.suptitle(
        f'{sysname}  –  Torsion-based classification\n'
        f'SAP: {n_sap} ({100*n_sap/nf:.1f}%)   '
        f'TSAP: {n_tsap} ({100*n_tsap/nf:.1f}%)   '
        f'UNK: {n_unk} ({100*n_unk/nf:.1f}%)',
        fontsize=11)
    plt.tight_layout()
    fig.savefig(os.path.join(outdir, 'B_classification.png'), dpi=150)
    plt.close(fig)
    print(f'  [B] {os.path.join(outdir, "B_classification.png")}')

    # ── Classification histogram ───────────────────────────────────
    cls_hist_csv = os.path.join(outdir, 'B_classification_histogram.csv')
    with open(cls_hist_csv, 'w', newline='') as f:
        w = csv.writer(f)
        w.writerow(['geometry', 'count', 'fraction'])
        for g, n in [('SAP', n_sap), ('TSAP', n_tsap), ('UNK', n_unk)]:
            w.writerow([g, n, f'{n/nf:.4f}'])
    print(f'  [B] {cls_hist_csv}')

    fig, ax = plt.subplots(figsize=(4, 3.5))
    bars = ax.bar(['SAP', 'TSAP', 'UNK'],
                  [n_sap / nf, n_tsap / nf, n_unk / nf],
                  color=[C_SAP, C_TSAP, C_UNK],
                  edgecolor='k', linewidth=0.7)
    ax.set_ylabel('Fraction of frames', fontsize=10)
    ax.set_ylim(0, 1.05)
    ax.set_title(f'{sysname}  –  Geometry population', fontsize=10)
    for bar, val in zip(bars, [n_sap, n_tsap, n_unk]):
        ax.text(bar.get_x() + bar.get_width() / 2, bar.get_height() + 0.02,
                str(val), ha='center', va='bottom', fontsize=9)
    plt.tight_layout()
    fig.savefig(os.path.join(outdir, 'B_classification_histogram.png'), dpi=150)
    plt.close(fig)
    print(f'  [B] {os.path.join(outdir, "B_classification_histogram.png")}')

    return geoms, tor_mat, times


# ═══════════════════════════════════════════════════════
#  ARGUMENT PARSING
# ═══════════════════════════════════════════════════════

def parse_args():
    """
    Build the SYSTEMS list from command-line flags.

    Three input modes (freely mixable in one call):

    Mode 1  --system <name> --tpr <file> --xtc <file>
        Original TPR+XTC mode.  Repeat the triplet for each system.
        Output goes to <name>/analysis/ (or <name>/<outdir> with --outdir).

    Mode 2  --system <name> --gro <file>
        Single-frame GRO mode.  The MOL residue is extracted from the GRO
        regardless of whether the file is ligand-only or a complex.
        Output goes to <name>/analysis/ (or <name>/<outdir> with --outdir).
        --outdir applies to the immediately preceding --system just like in
        Mode 1.

    Mode 3  --list <list.txt> --dir <directory> [--batch-outdir <root>]
        Batch mode: each line of list.txt is a system name; the script
        looks in <dir> for input files and writes output to
        <root>/<name>/  (default root = "analysis").

        File resolution per name (first match wins):
          1. <dir>/<name>.tpr + .xtc     → TPR+XTC trajectory mode (preferred:
                                           TPR carries accurate Eu mass ~152 Da)
          2. <dir>/<name>.gro            → single-frame GRO mode (fallback;
                                           Eu mass patched to 151.964 Da)

        --list/--dir may be repeated to process multiple batches in one run.
    """
    parser = argparse.ArgumentParser(
        description='Metal coordination geometry analysis (SAP/TSAP) for 9-coord Eu complexes.',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=textwrap.dedent("""
        Examples
        --------
        # Mode 1 – TPR + XTC (unchanged from v1):
        python metal_geo_analysis_v2.py \\
            --system me_sssL_sap  --tpr me_sssL_sap/pr_0.tpr  --xtc me_sssL_sap/pr_0.xtc \\
            --system phe_sssL_sap --tpr phe_sssL_sap/pr_0.tpr --xtc phe_sssL_sap/pr_0.xtc

        # Mode 2 – single GRO (ligand-only or complex):
        python metal_geo_analysis_v2.py \\
            --system me_sssL_sap  --gro solv/me_sssL_sap.gro \\
            --system me_sssL_sap  --gro com/me_sssL_sap.gro  --outdir com_analysis

        # Mode 3 – batch from list.txt:
        python metal_geo_analysis_v2.py \\
            --list solv/list.txt --dir solv --batch-outdir results_solv \\
            --list com/list.txt  --dir com  --batch-outdir results_com

        # For each name in list.txt the script looks for (in priority order):
        #   <dir>/<name>.tpr + .xtc   (TPR+XTC trajectory — preferred, accurate Eu mass)
        #   <dir>/<name>.gro          (single-frame GRO — fallback, Eu mass patched)

        # Mixed:
        python metal_geo_analysis_v2.py \\
            --system ref --tpr ref.tpr --xtc ref.xtc \\
            --list com/list.txt --dir com --outdir results_com
        """)
    )

    # ── Per-system flags (Mode 1 & 2) ────────────────────────────
    parser.add_argument('--system', dest='names',  action='append',
                        metavar='NAME',
                        help='System name (output directory / plot title).')
    parser.add_argument('--outdir', dest='outdirs', action='append',
                        metavar='DIR',
                        help='Output sub-directory for the preceding --system '
                             '(default: "analysis").')
    parser.add_argument('--tpr',    dest='tprs',   action='append',
                        metavar='FILE',
                        help='GROMACS TPR file (Mode 1).')
    parser.add_argument('--xtc',    dest='xtcs',   action='append',
                        metavar='FILE',
                        help='XTC trajectory file (Mode 1).')
    parser.add_argument('--gro',    dest='gros',   action='append',
                        metavar='FILE',
                        help='Single-frame GRO file (Mode 2). '
                             'Works for ligand-only and complex GROs.')

    # ── Batch flags (Mode 3) ──────────────────────────────────────
    parser.add_argument('--list',   dest='lists',  action='append',
                        metavar='FILE',
                        help='Path to list.txt for batch mode (Mode 3). '
                             'Each line is a system name.')
    parser.add_argument('--dir',    dest='dirs',   action='append',
                        metavar='DIR',
                        help='Directory containing input files for batch mode. '
                             'For each system name the script looks for '
                             '<dir>/<name>.gro first, then <dir>/<name>.tpr + .xtc. '
                             'Must be paired with --list.')
    parser.add_argument('--batch-outdir', dest='batch_outdirs',
                        action='append', metavar='DIR',
                        help='Root output directory for a --list/--dir batch '
                             '(default: "analysis"). '
                             'Output for each system goes to <batch-outdir>/<name>/.')

    args = parser.parse_args()

    # ─────────────────────────────────────────────────────────────
    # Build the ordered list of systems to process.
    # Each entry: {'name': str, 'mode': 'tpr_xtc'|'gro', ...cfg...}
    # ─────────────────────────────────────────────────────────────
    systems = []   # list of dicts, order preserved for output

    # ── Mode 1 & 2: per-system flags ─────────────────────────────
    names  = args.names  or []
    outdirs = args.outdirs or []
    tprs   = args.tprs   or []
    xtcs   = args.xtcs   or []
    gros   = args.gros   or []

    # We allow --outdir to be omitted; pad to match --system count
    while len(outdirs) < len(names):
        outdirs.append('analysis')

    # Validate: each --system needs either (--tpr + --xtc) or --gro
    # We consume tprs/xtcs/gros in order; user must supply them
    # in the same order as --system.
    tpr_idx = 0
    gro_idx = 0
    for idx, (name, odir) in enumerate(zip(names, outdirs)):
        # Decide mode by checking which pool has remaining entries
        # at the current position.  We rely on the user ordering them
        # consistently; if both pools are exhausted, raise an error.
        has_tpr = tpr_idx < len(tprs)
        has_gro = gro_idx < len(gros)
        has_xtc = tpr_idx < len(xtcs)

        if has_tpr and has_xtc:
            # Check if the next available tpr actually belongs to this system
            # by consuming it. (User must pair them in order.)
            systems.append({
                'name'  : name,
                'mode'  : 'tpr_xtc',
                'tpr'   : tprs[tpr_idx],
                'xtc'   : xtcs[tpr_idx],
                'outdir': os.path.join(name, odir),
            })
            tpr_idx += 1
        elif has_gro:
            systems.append({
                'name'  : name,
                'mode'  : 'gro',
                'gro'   : gros[gro_idx],
                'outdir': os.path.join(name, odir),
            })
            gro_idx += 1
        else:
            parser.error(
                f'System "{name}": no --tpr/--xtc pair or --gro file found.\n'
                'Each --system must be followed by either '
                '--tpr <f> --xtc <f>  or  --gro <f>.'
            )

    # Warn about unused tpr/xtc/gro entries
    if tpr_idx < len(tprs):
        print(f'WARNING: {len(tprs)-tpr_idx} unused --tpr entries.')
    if gro_idx < len(gros):
        print(f'WARNING: {len(gros)-gro_idx} unused --gro entries.')

    # ── Mode 3: batch ─────────────────────────────────────────────
    lists        = args.lists        or []
    batch_dirs   = args.dirs         or []
    batch_oroots = args.batch_outdirs or []

    if len(lists) != len(batch_dirs):
        parser.error(
            f'--list and --dir must be paired: got {len(lists)} --list '
            f'and {len(batch_dirs)} --dir.'
        )

    # Pad batch output roots
    while len(batch_oroots) < len(lists):
        batch_oroots.append('analysis')

    for lst_file, gro_dir, oroot in zip(lists, batch_dirs, batch_oroots):
        if not os.path.isfile(lst_file):
            parser.error(f'list file not found: {lst_file}')
        if not os.path.isdir(gro_dir):
            parser.error(f'--dir not found: {gro_dir}')
        with open(lst_file) as fh:
            batch_names = [line.strip() for line in fh if line.strip()]
        for bname in batch_names:
            gro_path = os.path.join(gro_dir, bname + '.gro')
            tpr_path = os.path.join(gro_dir, bname + '.tpr')
            xtc_path = os.path.join(gro_dir, bname + '.xtc')
            # TPR+XTC takes priority over GRO: TPR carries accurate masses
            # (EU ~152 Da), whereas bare GRO relies on MDAnalysis guessing
            # from atom type and gives EU mass = 0.0.
            if os.path.isfile(tpr_path) and os.path.isfile(xtc_path):
                systems.append({
                    'name'  : bname,
                    'mode'  : 'tpr_xtc',
                    'tpr'   : tpr_path,
                    'xtc'   : xtc_path,
                    'outdir': os.path.join(oroot, bname),
                })
            elif os.path.isfile(gro_path):
                systems.append({
                    'name'  : bname,
                    'mode'  : 'gro',
                    'gro'   : gro_path,
                    'outdir': os.path.join(oroot, bname),
                })
            else:
                parser.error(
                    f'No input found for batch entry "{bname}" in {gro_dir}.\n'
                    f'  Looked for: {tpr_path} + {xtc_path}\n'
                    f'         or : {gro_path}'
                )

    if not systems:
        parser.error(
            'No systems specified.\n'
            'Use one of:\n'
            '  --system <name> --tpr <file> --xtc <file>\n'
            '  --system <name> --gro <file>\n'
            '  --list <list.txt> --dir <directory> [--batch-outdir <root>]\n'
            '    (each name resolved to <dir>/<name>.gro  or  <dir>/<name>.tpr + .xtc)'
        )

    # ── File existence checks ──────────────────────────────────────
    for cfg in systems:
        if cfg['mode'] == 'tpr_xtc':
            if not os.path.isfile(cfg['tpr']):
                parser.error(f'TPR not found: {cfg["tpr"]}')
            if not os.path.isfile(cfg['xtc']):
                parser.error(f'XTC not found: {cfg["xtc"]}')
        else:
            if not os.path.isfile(cfg['gro']):
                parser.error(f'GRO not found: {cfg["gro"]}')

    return systems


# ═══════════════════════════════════════════════════════
#  MAIN
# ═══════════════════════════════════════════════════════

def main():
    systems = parse_args()

    for cfg in systems:
        sysname = cfg['name']
        print(f'\n{"═"*60}')
        print(f'  System : {sysname}  [{cfg["mode"]}]')
        if cfg['mode'] == 'tpr_xtc':
            print(f'  TPR    : {cfg["tpr"]}')
            print(f'  XTC    : {cfg["xtc"]}')
        else:
            print(f'  GRO    : {cfg["gro"]}')
        print(f'  Output : {cfg["outdir"]}')
        print(f'{"═"*60}')

        # ── Load structure / trajectory ───────────────────────────
        if cfg['mode'] == 'tpr_xtc':
            u, mol_ag = load_mol_universe(cfg['tpr'], cfg['xtc'])
        else:
            u, mol_ag = load_mol_gro(cfg['gro'])

        print(f'  MOL atoms : {mol_ag.n_atoms}   Frames : {u.trajectory.n_frames}')

        # ── Topology detection ────────────────────────────────────
        topo_key, topo = detect_topology(mol_ag)
        print(f'  Topology  : {topo_key}  '
              f'(EU mol-idx={topo["metal_idx"]}, '
              f'coord={topo["coord_idx"]})')

        # Reference positions for alignment = first frame
        u.trajectory[0]
        make_mol_whole(u, mol_ag)
        ref_pos = mol_ag.positions.copy()
        hmask   = heavy_mask(mol_ag)

        print('  Part A reference: per-frame local ideal SAP/TSAP polyhedra')

        # Single trajectory pass – collect all data
        print('  Processing frames …')
        records = collect_trajectory_data(u, mol_ag, ref_pos, hmask, topo)
        print(f'  Done.  {len(records)} frames processed.')

        # Part A output
        print('\n  ── Part A : RMSD to ideal SAP / TSAP ──')
        rs_sap, rs_tsap, _ = run_part_A(records, mol_ag, cfg['outdir'], sysname)

        # Part B output
        print('\n  ── Part B : Torsion analysis ──')
        geoms, _, _ = run_part_B(records, cfg['outdir'], sysname)

        # Summary
        n_sap  = geoms.count('SAP')
        n_tsap = geoms.count('TSAP')
        n_unk  = geoms.count('UNK')
        nf     = len(geoms)
        print(f'\n  Summary  {sysname}:')
        print(f'    Geometry (torsion): SAP={n_sap} ({100*n_sap/nf:.1f}%)  '
              f'TSAP={n_tsap} ({100*n_tsap/nf:.1f}%)  UNK={n_unk} ({100*n_unk/nf:.1f}%)')
        print(f'    RMSD→SAP  : {rs_sap.mean():.3f} ± {rs_sap.std():.3f} Å')
        print(f'    RMSD→TSAP : {rs_tsap.mean():.3f} ± {rs_tsap.std():.3f} Å')

    print('\nAll done.')


if __name__ == '__main__':
    main()
