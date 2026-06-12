"""
Shared utilities for PCA analysis pipeline.

Consolidates duplicated constants, FES computation, DCCM computation,
and binding-site definitions that were previously scattered across
multiple scripts.

Created during P1 fix cycle (Groups 5, 6, 9).
"""

import numpy as np

# ---------------------------------------------------------------------------
# Gas constant
# ---------------------------------------------------------------------------
R_GAS = 0.008314463  # kJ/(mol·K), IUPAC 2019 CODATA value

# ---------------------------------------------------------------------------
# Binding site definitions (three conventions used across analyses)
# ---------------------------------------------------------------------------
# 1. PCA residue contribution analysis (broad pocket)
BINDING_SITE_RESID_RANGE = (377, 490)       # 114 residues
BINDING_SITE_N_RESIDUES = 114

# 2. DCCM visualization rectangle (CA atom indices in the 582×582 matrix)
DCCM_BINDING_SITE_START = 376               # CA matrix index for resid 379
DCCM_BINDING_SITE_END = 487                 # CA matrix index for resid 490 (inclusive)
DCCM_BINDING_SITE_RESID_OFFSET = 3          # resid = CA_index + RESID_OFFSET

# 3. Contact analysis (specific ligand-contacting residues)
BINDING_SITE_CONTACT_RESIDUES = set(range(377, 389)) | set(range(483, 491))  # 20 residues

# ---------------------------------------------------------------------------
# DCCM computation
# ---------------------------------------------------------------------------

def compute_dccm(ca_coords, n_ca_atoms=582, n_coords_per_ca=3):
    """Compute the dynamic cross-correlation matrix (DCCM) from CA coordinates.

    Parameters
    ----------
    ca_coords : ndarray, shape (n_frames, n_ca_atoms * 3)
        Flattened CA coordinates for all frames.
    n_ca_atoms : int
        Number of CA atoms (default 582).
    n_coords_per_ca : int
        Coordinates per CA atom (default 3).

    Returns
    -------
    dccm : ndarray, shape (n_ca_atoms, n_ca_atoms)
        Normalized cross-correlation matrix. Values in [-1, 1].
    """
    nframes = ca_coords.shape[0]
    mean = ca_coords.mean(axis=0).astype(np.float64)
    delta = (ca_coords.astype(np.float64) - mean)
    cov = (delta.T @ delta) / (nframes - 1)
    cov_4d = cov.reshape(n_ca_atoms, n_coords_per_ca, n_ca_atoms, n_coords_per_ca)
    numerator = np.einsum("imjm->ij", cov_4d)
    var_i = numerator.diagonal().copy()
    denom = np.sqrt(np.outer(var_i, var_i))
    denom[denom == 0] = 1.0
    dccm = numerator / denom
    return dccm

# ---------------------------------------------------------------------------
# Free energy surface (FES) computation
# ---------------------------------------------------------------------------

def compute_fes(pc1, pc2, n_bins=20, kT=2.4942, x_range=None, y_range=None):
    """Compute the free energy surface from two PC projection arrays.

    Uses histogram binning with automatic empty-bin masking and
    minimum-to-zero shift (standard convention for FES).

    Parameters
    ----------
    pc1, pc2 : array-like
        PC projection scores (1-D arrays of equal length).
    n_bins : int
        Number of bins per dimension (default 20).
    kT : float
        Thermal energy in kJ/mol (default 2.4942 at 300 K).
    x_range, y_range : tuple of (min, max), optional
        Custom bin range for each axis. If None, uses data range (auto).

    Returns
    -------
    G : ndarray, shape (n_bins, n_bins)
        Free energy grid in kJ/mol. Empty bins are NaN. Minimum = 0.
    xedges : ndarray
        Bin edges along PC1.
    yedges : ndarray
        Bin edges along PC2.
    """
    if x_range is not None and y_range is not None:
        xedges = np.linspace(x_range[0], x_range[1], n_bins + 1)
        yedges = np.linspace(y_range[0], y_range[1], n_bins + 1)
        H, _, _ = np.histogram2d(pc1, pc2, bins=(xedges, yedges))
    else:
        H, xedges, yedges = np.histogram2d(pc1, pc2, bins=n_bins)
    H = H.T
    prob = H / H.sum()
    G = np.full_like(prob, np.nan, dtype=float)
    mask = prob > 0
    G[mask] = -kT * np.log(prob[mask])
    G -= G[mask].min()  # Shift so minimum = 0
    return G, xedges, yedges
