"""
Standardized MDAnalysis selection strings for HREX analysis tools.

This module provides consistent default selections across rmsd.py and com_dist.py
for receptor-ligand analysis. These constants ensure uniform atom group definitions
when multiple analysis tools are combined in production workflows.

Usage:
    from selection_defaults import DEFAULT_LIGAND_SEL, DEFAULT_RECEPTOR_SEL
    
    parser.add_argument('--ligand_sel', default=DEFAULT_LIGAND_SEL,
                        help=f"Ligand selection (default: {DEFAULT_LIGAND_SEL})")

Selection Syntax:
    MDAnalysis uses a powerful selection language for defining atom groups:
    
    - 'protein': Built-in keyword matching standard protein residues
    - 'resname LIG': Matches residues with residue name "LIG"
    - 'not name H*': Excludes atoms with names starting with H (hydrogens)
    - 'and', 'or', 'not': Boolean operators (always use explicit operators)
    - 'name CA': Matches atoms named CA (e.g., C-alpha carbons)
    - 'nucleic': Built-in keyword for nucleic acid residues
    - 'name P': Matches phosphorus atoms (backbone of nucleic acids)
    
    Full syntax reference:
    https://docs.mdanalysis.org/stable/documentation_pages/selections.html

Why heavy atoms only?
    Most HREX analysis focuses on heavy atoms (non-hydrogen) because:
    - Hydrogen positions are often unreliable in MD trajectories
    - Heavy atom distances are more meaningful for binding interactions
    - RMSD calculations without hydrogens are standard practice
    - Reduces computational cost without sacrificing accuracy

Note on fp.py compatibility:
    The interaction fingerprint tool (fp.py) uses resid-based selections 
    (e.g., 'resid 1-80') rather than these resname-based defaults. This is
    intentional - fingerprints track receptor residue positions for contact
    mapping, while these constants identify molecular components by type.
    See fp.py documentation for its selection paradigm.

Constants:
    DEFAULT_LIGAND_SEL: Standard ligand heavy atoms (resname LIG)
    DEFAULT_RECEPTOR_SEL: Standard receptor heavy atoms (protein residues)
    DEFAULT_ALIGNMENT_SEL: C-alpha atoms for protein alignment (RMSD)
    DEFAULT_ALIGNMENT_NUCLEIC: Phosphorus atoms for nucleic acid alignment
    DEFAULT_RMSD_SEL: Alias for DEFAULT_LIGAND_SEL (RMSD calculation target)
"""

# Ligand selections (resname-based)
# Select ligand heavy atoms by residue name "LIG", excluding hydrogens
DEFAULT_LIGAND_SEL = "resname LIG and not name H*"

# Receptor selections (keyword-based)
# Select all protein heavy atoms, excluding hydrogens
DEFAULT_RECEPTOR_SEL = "protein and not name H*"

# Alignment selections for RMSD calculations
# Protein systems: Use C-alpha backbone atoms for structural alignment
DEFAULT_ALIGNMENT_SEL = "protein and name CA"

# Nucleic acid systems: Use phosphorus backbone atoms for structural alignment
DEFAULT_ALIGNMENT_NUCLEIC = "nucleic and name P"

# RMSD calculation selection (by convention, same as ligand heavy atoms)
# This measures ligand conformational changes relative to reference structure
DEFAULT_RMSD_SEL = DEFAULT_LIGAND_SEL
