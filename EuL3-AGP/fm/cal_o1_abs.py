#!/usr/bin/env python

import logging
import numpy as np
import MDAnalysis as mda

# Default variables
TOPOLOGY_FILE = 'v1.pdb'
TRAJECTORY_FILE = 'v1.xtc'
OUTPUT_FILE = 'o1.dat'

SIGNED_ANGLE = True  # False for unsigned (0-360°)

CA34_SELECTION = 'protein and resid 34 and name CA'
CA42_SELECTION = 'protein and resid 42 and name CA'
COM_SELECTION = 'resname E3P and name N1 N2 N3 N4'
EU_SELECTION = 'name EU3'

# Setup logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s %(levelname)s: %(message)s')

try:
    logging.info('Loading universe...')
    u = mda.Universe(TOPOLOGY_FILE, TRAJECTORY_FILE)

    ca34 = u.select_atoms(CA34_SELECTION)
    ca42 = u.select_atoms(CA42_SELECTION)
    com_atoms = u.select_atoms(COM_SELECTION)
    eu_atom = u.select_atoms(EU_SELECTION)

    if ca34.n_atoms != 1 or ca42.n_atoms != 1:
        raise ValueError('Invalid CA selections')
    if com_atoms.n_atoms != 4:
        raise ValueError(f'Expected 4 atoms for COM, found {com_atoms.n_atoms}')
    if eu_atom.n_atoms != 1:
        raise ValueError(f'Expected 1 EU atom, found {eu_atom.n_atoms}')

    angles = []

    logging.info(f'Starting analysis over {len(u.trajectory)} frames...')
    for ts in u.trajectory:
        # Get positions for dihedral points
        pA = ca42.positions[0]
        pB = ca34.positions[0]
        pC = eu_atom.positions[0]
        pD = com_atoms.center_of_mass()

        # Compute dihedral in radians
        angle_rad = mda.lib.distances.calc_dihedrals(pA, pB, pC, pD)

        if SIGNED_ANGLE:
            # Stable signed angle (-180 to 180)
            theta = np.degrees(np.arctan2(np.sin(angle_rad), np.cos(angle_rad)))
        else:
            # Unsigned (0 to 360)
            theta = np.degrees(angle_rad) % 360

        angles.append(theta)

    # Save output: header + single column
    logging.info(f'Saving results to {OUTPUT_FILE}')
    np.savetxt(OUTPUT_FILE, angles, fmt='%.6f', header='o1', comments='', delimiter=' ')

    logging.info('Analysis completed successfully.')

except Exception as e:
    logging.error(f'Error during analysis: {e}')
    raise
