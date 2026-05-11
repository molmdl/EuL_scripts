#!/usr/bin/env python3
"""
Metal coordination geometry analysis for 9-coordinate Eu complexes.

Analyzes SAP (square antiprism) and TSAP (twisted square antiprism)
geometries in MD trajectories.

Section I   – Trajectory loading, PBC handling, alignment
Section II  – Geometry analysis A: RMSD to ideal SAP/TSAP polyhedra
Section III – Geometry analysis B: N-C-C-N torsion angles

Usage:
    # Run one or more systems by specifying name / TPR / XTC triplets:
    python metal_geo_analysis.py \
        --system me_sssL_sap  --tpr me_sssL_sap/pr_0.tpr  --xtc me_sssL_sap/pr_0.xtc \
        --system phe_sssL_sap --tpr phe_sssL_sap/pr_0.tpr --xtc phe_sssL_sap/pr_0.xtc

    Output is written to <system_name>/<outdir>/ for each system.
    Part A uses per-frame local idealisation: each frame keeps its own
    coordination-shell size/height parameters, while only the angular
    SAP/TSAP arrangement is idealised before RMSD evaluation.
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
#  CONFIGURATION
# ═══════════════════════════════════════════════════════

# Index of EU3 metal centre in the MOL atom array
METAL_IDX = 54

# Indices of the 9 coordinating atoms in the MOL atom array:
#   O0(0), O1(1), O2(2), N3(3), N4(4), N5(5), N6(6), N7(7), N63(8-th element)
COORD_IDX = [0, 1, 2, 3, 4, 5, 6, 7, 63]

# Sub-group indices *within* the 9-element coord sub-array
BOTTOM = [3, 4, 5, 6]   # ring N   : N3 N4 N5 N6
TOP    = [0, 1, 2, 7]   # O-arms + proximal chrom N  : O0 O1 O2 N7
CAP    = [8]            # distal chrom N              : N63

# Torsion definitions  (name, i, j, k, l)  using MOL atom indices
RING_TORS = [
    ('T1 N3-C8-C9-N4',    3,  8,  9,  4),
    ('T2 N4-C10-C11-N5',  4, 10, 11,  5),
    ('T3 N5-C12-C13-N6',  5, 12, 13,  6),
    ('T4 N6-C14-C15-N3',  6, 14, 15,  3),
]
CHROM_TOR = ('Tc N3-C16-C23-N7', 3, 16, 23, 7)

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
#  TRAJECTORY LOADING  (Section I)
# ═══════════════════════════════════════════════════════

def load_mol_universe(tpr, xtc):
    """
    Load the MOL residue from tpr+xtc, handling two cases transparently:

    Case A – stripped XTC (atoms == MOL only):
        Build a MOL-only universe via Merge + load_new.  mol_ag = u.atoms.

    Case B – full XTC (atoms == whole simulation box):
        Load Universe(tpr, xtc) directly.  mol_ag is the persistent
        'resname MOL' AtomGroup inside that universe.

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
        raise ValueError(
            f"XTC has {xtc_natoms} atoms, but TPR has {u_tpr.atoms.n_atoms} total "
            f"and {n_mol} MOL atoms.  Cannot determine loading strategy."
        )

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

# Pre-build torsion index array once at module level
_TORSION_IJKL = np.array(
    [[i, j, k, l] for _, i, j, k, l in RING_TORS]
    + [[CHROM_TOR[1], CHROM_TOR[2], CHROM_TOR[3], CHROM_TOR[4]]],
    dtype=int
)   # (5, 4)  – rows: T1 T2 T3 T4 Tc


def nearest_frame(series, target):
    """Index of the element in series closest to target (ignores NaN)."""
    arr  = np.asarray(series, dtype=float)
    diffs = np.where(np.isnan(arr), np.inf, np.abs(arr - target))
    return int(np.argmin(diffs))


def collect_trajectory_data(u, mol_ag, ref_pos, hmask):
    """
    Iterate over every frame, compute Part-A and Part-B quantities.

    Alignment (Kabsch on heavy atoms) for Part A; raw unwrapped positions
    for Part B torsions.  All per-frame geometry is vectorised:

      • 5 torsions computed with a single calc_dihedrals_batch call.
      • frame-specific SAP/TSAP direction templates (Chimera-style, unit vectors).
      • RMSD via per-atom radial scaling + batched Kabsch over 16 permutations.

    Returns a list of per-frame dicts.
    """
    records = []
    for ts in u.trajectory:
        make_mol_whole(u, mol_ag)
        raw_pos = mol_ag.positions.copy()            # (N_mol, 3) unwrapped
        aln_pos = align_pos(raw_pos, ref_pos, hmask) # (N_mol, 3) aligned

        eu    = aln_pos[METAL_IDX]
        coord = aln_pos[COORD_IDX] - eu              # (9, 3) centred on Eu

        # ── Part A : per-frame local idealisation + batched RMSD ────
        # Chimera-style: direction templates, per-atom radial scaling in min_rmsd_cyclic
        dirs_sap, dirs_tsap = build_local_ideal_directions(coord)
        rmsd_sap  = min_rmsd_cyclic(coord, dirs_sap)
        rmsd_tsap = min_rmsd_cyclic(coord, dirs_tsap)

        # ── Part B : all 5 torsions in one vectorised call ──────────
        all_tors = calc_dihedrals_batch(raw_pos, _TORSION_IJKL)  # (5,)
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
    Build the SYSTEMS dict from command-line flags.

    Each system requires three consecutive flags:
        --system <name>  --tpr <path>  --xtc <path>

    These may be repeated for as many systems as desired.
    """
    parser = argparse.ArgumentParser(
        description='Metal coordination geometry analysis (SAP/TSAP) for 9-coord Eu complexes.',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=textwrap.dedent("""
        Examples
        --------
        # Two systems:
        python metal_geo_analysis.py \
            --system me_sssL_sap  --tpr me_sssL_sap/pr_0.tpr  --xtc me_sssL_sap/pr_0.xtc \
            --system phe_sssL_sap --tpr phe_sssL_sap/pr_0.tpr --xtc phe_sssL_sap/pr_0.xtc

        Output is written to <system_name>/<outdir>/ for each system.
        Part A uses frame-specific local ideal SAP/TSAP references.
        """)
    )
    parser.add_argument('--system', dest='names',   action='append', metavar='NAME',
                        help='System name (used for output directory and plot titles).')
    parser.add_argument('--outdir', dest='outdir',   action='append', metavar='NAME',
                        help='Output directory (used for output directory).')
    parser.add_argument('--tpr',    dest='tprs',    action='append', metavar='FILE',
                        help='GROMACS TPR file for this system.')
    parser.add_argument('--xtc',    dest='xtcs',    action='append', metavar='FILE',
                        help='XTC trajectory file for this system (may be solvent-stripped).')

    args = parser.parse_args()

    names = args.names or []
    outdir = args.outdir or ['analysis']
    tprs  = args.tprs  or []
    xtcs  = args.xtcs  or []

    if not names:
        parser.error(
            'No systems specified.\n'
            'Use: --system <name> --tpr <file> --xtc <file>  (repeat for each system)'
        )

    if not (len(names) == len(tprs) == len(xtcs)):
        parser.error(
            f'Mismatch: got {len(names)} --system, {len(tprs)} --tpr, {len(xtcs)} --xtc.\n'
            'Each --system must be followed by exactly one --tpr and one --xtc.'
        )

    systems = {}
    for name, outdir, tpr, xtc in zip(names, outdir, tprs, xtcs):
        if name in systems:
            parser.error(f'Duplicate system name: "{name}"')
        if not os.path.isfile(tpr):
            parser.error(f'TPR not found: {tpr}')
        if not os.path.isfile(xtc):
            parser.error(f'XTC not found: {xtc}')
        systems[name] = {
            'tpr':    tpr,
            'xtc':    xtc,
            'outdir': os.path.join(name, outdir),
        }

    return systems


# ═══════════════════════════════════════════════════════
#  MAIN
# ═══════════════════════════════════════════════════════

def main():
    systems = parse_args()
    sys_items = list(systems.items())

    # ── Per-system analysis ───────────────────────────────────────
    for sysname, cfg in sys_items:
        print(f'\n{"═"*60}')
        print(f'  System : {sysname}')
        print(f'{"═"*60}')

        u, mol_ag = load_mol_universe(cfg['tpr'], cfg['xtc'])
        print(f'  MOL atoms : {mol_ag.n_atoms}   Frames : {u.trajectory.n_frames}')

        # Reference positions for alignment = first frame of this system
        u.trajectory[0]
        make_mol_whole(u, mol_ag)
        ref_pos = mol_ag.positions.copy()
        hmask   = heavy_mask(mol_ag)

        print('  Part A reference: per-frame local ideal SAP/TSAP polyhedra')

        # Single trajectory pass – collect all data
        print('  Processing frames …')
        records = collect_trajectory_data(u, mol_ag, ref_pos, hmask)
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
