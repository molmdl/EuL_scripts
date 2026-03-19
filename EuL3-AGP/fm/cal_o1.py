#!/usr/bin/env python

import logging
import numpy as np
import MDAnalysis as mda

# Default variables
TOPOLOGY_FILE = 'v1.pdb'
TRAJECTORY_FILE = 'v1.xtc'
OUTPUT_FILE = 'o1.dat'

HELIX_SELECTION = 'protein and resid 34:42 and backbone'
N_TERM_SELECTION = 'protein and resid 34:37 and backbone'
C_TERM_SELECTION = 'protein and resid 39:42 and backbone'
COM_SELECTION = 'resname E3P and name N1 N2 N3 N4'
EU_SELECTION = 'name EU3'  # Assuming unique; adjust if needed

# Setup logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s %(levelname)s: %(message)s')

try:
    logging.info('Loading universe...')
    u = mda.Universe(TOPOLOGY_FILE, TRAJECTORY_FILE)

    helix = u.select_atoms(HELIX_SELECTION)
    n_term = u.select_atoms(N_TERM_SELECTION)
    c_term = u.select_atoms(C_TERM_SELECTION)
    com_atoms = u.select_atoms(COM_SELECTION)
    eu_atom = u.select_atoms(EU_SELECTION)

    if helix.n_atoms == 0:
        raise ValueError('No atoms found in helix selection')
    if n_term.n_atoms == 0 or c_term.n_atoms == 0:
        raise ValueError('No atoms found in N/C-term selections')
    if com_atoms.n_atoms != 4:
        raise ValueError(f'Expected 4 atoms for COM, found {com_atoms.n_atoms}')
    if eu_atom.n_atoms != 1:
        raise ValueError(f'Expected 1 EU atom, found {eu_atom.n_atoms}')

    angles = []

    logging.info(f'Starting analysis over {len(u.trajectory)} frames...')
    for ts in u.trajectory:
        # Principal axes: shape (3,3), rows sorted decreasing eigenvalue
        p_axes = helix.principal_axes()
        v1 = p_axes[2]  # Smallest eigenvalue → helix long axis

        # Orient v1 consistently (N to C direction)
        com_n = n_term.center_of_mass()
        com_c = c_term.center_of_mass()
        dir_vec = com_c - com_n
        if np.dot(v1, dir_vec) < 0:
            v1 = -v1

        # Compute v2
        com = com_atoms.center_of_mass()
        eu_pos = eu_atom.positions[0]
        v2 = eu_pos - com

        # Normalize
        v1_norm = v1 / np.linalg.norm(v1)
        v2_norm = v2 / np.linalg.norm(v2)

        # Compute dot and cross
        dot = np.dot(v1_norm, v2_norm)
        cross = np.cross(v1_norm, v2_norm)

        # Unsigned angle (0-180)
        angle = np.degrees(np.arccos(np.clip(dot, -1.0, 1.0)))

        # Reference normal for sign: from helix COM to ligand COM
        helix_com = helix.center_of_mass()
        n_ref = com - helix_com
        sign = np.sign(np.dot(n_ref, cross))

        # Signed angle (-180 to 180)
        theta_signed = angle * sign

        angles.append(theta_signed)

    # Save output: header + single column
    logging.info(f'Saving results to {OUTPUT_FILE}')
    np.savetxt(OUTPUT_FILE, angles, fmt='%.6f', header='o1', comments='', delimiter=' ')

    logging.info('Analysis completed successfully.')

except Exception as e:
    logging.error(f'Error during analysis: {e}')
    raise
