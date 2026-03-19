#!/usr/bin/env python

import logging
import numpy as np
import MDAnalysis as mda

# Default variables
TOPOLOGY_FILE = 'v1.pdb'
TRAJECTORY_FILE = 'v1.xtc'
OUTPUT_FILE = 'o1.dat'

#SIGNED_ANGLE = True          # False for unsigned (0-180°, smoother near 0°/180°)
SIGNED_ANGLE = False          # False for unsigned (0-180°, smoother near 0°/180°)
USE_PRINCIPAL_AXIS = False   # True to use inertia principal axis instead of CA vector

HELIX_COM_SELECTION = 'protein and resid 34:42 and backbone'
CA34_SELECTION = 'protein and resid 34 and name CA'
CA42_SELECTION = 'protein and resid 42 and name CA'
PROTEIN_SELECTION = 'protein'
COM_SELECTION = 'resname E3P and name N1 N2 N3 N4'
EU_SELECTION = 'name EU3'

# Setup logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s %(levelname)s: %(message)s')

try:
    logging.info('Loading universe...')
    u = mda.Universe(TOPOLOGY_FILE, TRAJECTORY_FILE)

    helix_com_ag = u.select_atoms(HELIX_COM_SELECTION)
    ca34 = u.select_atoms(CA34_SELECTION)
    ca42 = u.select_atoms(CA42_SELECTION)
    protein = u.select_atoms(PROTEIN_SELECTION)
    com_atoms = u.select_atoms(COM_SELECTION)
    eu_atom = u.select_atoms(EU_SELECTION)

    if helix_com_ag.n_atoms == 0 or ca34.n_atoms != 1 or ca42.n_atoms != 1:
        raise ValueError('Invalid helix or CA selections')
    if protein.n_atoms == 0:
        raise ValueError('No protein atoms found')
    if com_atoms.n_atoms != 4:
        raise ValueError(f'Expected 4 atoms for ligand COM, found {com_atoms.n_atoms}')
    if eu_atom.n_atoms != 1:
        raise ValueError(f'Expected 1 EU atom, found {eu_atom.n_atoms}')

    angles = []

    logging.info(f'Starting analysis over {len(u.trajectory)} frames...')
    for ts in u.trajectory:
        # Define v1 (helix axis)
        if USE_PRINCIPAL_AXIS:
            p_axes = helix_com_ag.principal_axes()
            v1 = p_axes[2]  # smallest eigenvalue
        else:
            pos34 = ca34.positions[0]
            pos42 = ca42.positions[0]
            v1 = pos42 - pos34

        v1_norm_raw = v1 / np.linalg.norm(v1)

        # Consistent orientation: point towards protein COM from CA34
        protein_com = protein.center_of_mass()
        orient_vec = protein_com - pos34
        if np.dot(v1_norm_raw, orient_vec) < 0:
            v1 = -v1

        v1_norm = v1 / np.linalg.norm(v1)

        # v2: ligand COM to EU3
        ligand_com = com_atoms.center_of_mass()
        eu_pos = eu_atom.positions[0]
        v2 = eu_pos - ligand_com
        v2_norm = v2 / np.linalg.norm(v2)

        # Dot and cross
        dot = np.dot(v1_norm, v2_norm)
        cross = np.cross(v1_norm, v2_norm)

        # Base angle (0 to 180°)
        base_angle = np.degrees(np.arctan2(np.linalg.norm(cross), dot))

        if not SIGNED_ANGLE:
            theta = base_angle
        else:
            # Stable reference normal: protein COM - helix COM
            helix_com = helix_com_ag.center_of_mass()
            n_ref = protein_com - helix_com
            n_ref_norm = n_ref / np.linalg.norm(n_ref) if np.linalg.norm(n_ref) > 0 else np.array([0,0,1])

            sign = np.sign(np.dot(n_ref_norm, cross))
            if sign == 0:
                sign = 1.0
            theta = base_angle * sign

        angles.append(theta)

    # Save output
    logging.info(f'Saving results to {OUTPUT_FILE}')
    np.savetxt(OUTPUT_FILE, angles, fmt='%.6f', header='o1', comments='', delimiter=' ')

    logging.info('Analysis completed successfully.')

except Exception as e:
    logging.error(f'Error during analysis: {e}')
    raise
