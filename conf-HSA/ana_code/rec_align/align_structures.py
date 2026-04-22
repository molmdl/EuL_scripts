#!/usr/bin/env python3
"""
Structural alignment tool for PDB and GRO files.
Uses MDAnalysis for structural superposition with sequence alignment for residue matching.
Input format is preserved in output (PDB -> PDB, GRO -> GRO).

Usage:
    python align_structures.py -r reference.pdb -t target1.pdb target2.gro
    python align_structures.py -r reference.pdb -t hsa*.pdb hsa*.gro --output-dir ./output
"""

import argparse
import warnings
from pathlib import Path
from typing import Tuple, List, Dict

import numpy as np
import MDAnalysis as mda
from MDAnalysis.analysis import align
from MDAnalysis.analysis.rms import rmsd
from Bio.Align import PairwiseAligner


# =============================================================================
# Configuration Constants
# =============================================================================

DEFAULT_OUTPUT_SUFFIX_PDB = '_ali.pdb'
DEFAULT_OUTPUT_SUFFIX_GRO = '_ali.gro'

warnings.filterwarnings('ignore', message='.*discontinuous.*')


# =============================================================================
# File Loading Functions
# =============================================================================

def load_structure(filepath: Path) -> mda.Universe:
    """Load structure file using MDAnalysis (handles PDB and GRO)."""
    return mda.Universe(str(filepath))


# =============================================================================
# Alignment Functions
# =============================================================================

def get_aligned_residue_pairs(ref_universe: mda.Universe, 
                               target_universe: mda.Universe) -> List[Tuple[int, int]]:
    """
    Use sequence alignment to find corresponding residues between structures.
    
    Returns list of (ref_residue_index, target_residue_index) pairs.
    """
    ref_protein = ref_universe.select_atoms("protein")
    target_protein = target_universe.select_atoms("protein")
    
    # Get sequences as one-letter codes
    ref_seq = ref_protein.residues.sequence(format="Seq")
    target_seq = target_protein.residues.sequence(format="Seq")
    
    # Perform global sequence alignment
    aligner = PairwiseAligner()
    aligner.mode = "global"
    alignment = aligner.align(ref_seq, target_seq)[0]
    
    # Parse the alignment to get corresponding positions
    ref_aln = str(alignment[0])
    target_aln = str(alignment[1])
    
    # Build the correspondence by iterating through alignment
    ref_idx = 0
    target_idx = 0
    pairs = []
    
    for r_char, t_char in zip(ref_aln, target_aln):
        if r_char != '-' and t_char != '-':
            pairs.append((ref_idx, target_idx))
            ref_idx += 1
            target_idx += 1
        elif r_char != '-':
            ref_idx += 1  # gap in target
        elif t_char != '-':
            target_idx += 1  # gap in reference
    
    return pairs


def align_structure(ref_universe: mda.Universe, target_universe: mda.Universe,
                    output_path: Path) -> float:
    """
    Align target structure to reference using MDAnalysis alignto().
    
    Uses sequence alignment to identify corresponding residues for structures
    with different residue numbers or insertions/deletions.
    
    Parameters
    ----------
    ref_universe : mda.Universe
        Reference structure
    target_universe : mda.Universe
        Target structure to align
    output_path : Path
        Output file path (format determined by extension)
    
    Returns
    -------
    float
        RMSD after alignment
    """
    # Get corresponding residues via sequence alignment
    residue_pairs = get_aligned_residue_pairs(ref_universe, target_universe)
    
    if len(residue_pairs) < 3:
        raise ValueError(f"Insufficient aligned residues ({len(residue_pairs)}) for superposition")
    
    # Get residue objects
    ref_protein = ref_universe.select_atoms("protein")
    target_protein = target_universe.select_atoms("protein")
    ref_residues = list(ref_protein.residues)
    target_residues = list(target_protein.residues)
    
    # Build selection by residue ID
    ref_resids = [ref_residues[i].resid for i, j in residue_pairs]
    target_resids = [target_residues[j].resid for i, j in residue_pairs]
    
    # Create selection strings
    ref_sel = "resid " + " ".join(map(str, ref_resids)) + " and name CA"
    target_sel = "resid " + " ".join(map(str, target_resids)) + " and name CA"
    
    select_dict = {
        'reference': ref_sel,
        'mobile': target_sel
    }
    
    # Perform alignment using MDAnalysis alignto()
    old_rmsd, new_rmsd = align.alignto(
        target_universe, 
        ref_universe, 
        select=select_dict,
        match_atoms=True
    )
    
    # Write output (format determined by extension)
    target_universe.atoms.write(str(output_path))
    
    return new_rmsd


# =============================================================================
# Verification Functions
# =============================================================================

def verify_alignment(output_path: Path, expected_path: Path, 
                     tolerance: float = 0.1) -> Tuple[float, bool]:
    """
    Verify alignment by comparing CA coordinates.
    
    Uses superposition to compare internal geometry, not absolute coordinates.
    """
    output_u = mda.Universe(str(output_path))
    expected_u = mda.Universe(str(expected_path))
    
    output_ca = output_u.select_atoms("protein and name CA")
    expected_ca = expected_u.select_atoms("protein and name CA")
    
    n = min(len(output_ca), len(expected_ca))
    if n < 3:
        return float('inf'), False
    
    # Compute RMSD with superposition (compares internal geometry)
    rmsd_val = rmsd(output_ca.positions[:n], expected_ca.positions[:n], 
                    superposition=True)
    
    return rmsd_val, rmsd_val < tolerance


# =============================================================================
# CLI Interface
# =============================================================================

def parse_arguments():
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(
        description='Align PDB and GRO structures to a reference PDB. Output format matches input format.',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  Single PDB file:
    python align_structures.py -r reference.pdb -t target.pdb
  
  Multiple files (PDB outputs PDB, GRO outputs GRO):
    python align_structures.py -r reference.pdb -t target1.pdb target2.gro
  
  Batch processing:
    python align_structures.py -r reference.pdb -t hsa*.pdb hsa*.gro --output-dir ./aligned
"""
    )
    
    parser.add_argument(
        '-r', '--reference',
        type=Path,
        required=True,
        help='Reference PDB file for alignment'
    )
    
    parser.add_argument(
        '-t', '--targets',
        type=Path,
        nargs='+',
        required=True,
        help='Target structure files (PDB or GRO) to align'
    )
    
    parser.add_argument(
        '-o', '--output-dir',
        type=Path,
        default=Path('.'),
        help='Output directory for aligned structures (default: current directory)'
    )
    
    parser.add_argument(
        '--suffix-pdb',
        type=str,
        default=DEFAULT_OUTPUT_SUFFIX_PDB,
        help=f'Output filename suffix for PDB files (default: {DEFAULT_OUTPUT_SUFFIX_PDB})'
    )
    
    parser.add_argument(
        '--suffix-gro',
        type=str,
        default=DEFAULT_OUTPUT_SUFFIX_GRO,
        help=f'Output filename suffix for GRO files (default: {DEFAULT_OUTPUT_SUFFIX_GRO})'
    )
    
    parser.add_argument(
        '--verify',
        action='store_true',
        help='Verify outputs against expected files after alignment'
    )
    
    parser.add_argument(
        '--expected-dir',
        type=Path,
        help='Directory containing expected outputs for verification'
    )
    
    return parser.parse_args()


def main():
    """Main entry point for structural alignment."""
    args = parse_arguments()
    
    if not args.reference.exists():
        raise FileNotFoundError(f"Reference file not found: {args.reference}")
    
    args.output_dir.mkdir(parents=True, exist_ok=True)
    
    print(f"Loading reference: {args.reference}")
    ref_universe = load_structure(args.reference)
    
    results = []
    
    for target_path in args.targets:
        if not target_path.exists():
            print(f"Warning: Target file not found, skipping: {target_path}")
            continue
        
        extension = target_path.suffix.lower()
        is_gro = extension == '.gro'
        suffix = args.suffix_gro if is_gro else args.suffix_pdb
        
        output_name = target_path.stem + suffix
        output_path = args.output_dir / output_name
        
        print(f"\nProcessing: {target_path}")
        
        try:
            target_universe = load_structure(target_path)
            rmsd = align_structure(ref_universe, target_universe, output_path)
            
            results.append((target_path.name, rmsd, output_path.name, extension))
            print(f"  -> {output_path.name} (RMSD: {rmsd:.3f} A)")
            
        except Exception as e:
            print(f"  Error: {e}")
            results.append((target_path.name, None, None, extension))
    
    print("\n" + "=" * 60)
    print("ALIGNMENT SUMMARY")
    print("=" * 60)
    
    successful = [r for r in results if r[1] is not None]
    failed = [r for r in results if r[1] is None]
    
    if successful:
        rmsd_values = np.array([r[1] for r in successful])
        print(f"Successfully aligned: {len(successful)}/{len(results)}")
        print(f"Mean RMSD: {np.mean(rmsd_values):.3f} A")
        print(f"Min RMSD:  {np.min(rmsd_values):.3f} A")
        print(f"Max RMSD:  {np.max(rmsd_values):.3f} A")
    
    if failed:
        print(f"\nFailed: {len(failed)}")
        for name, _, _, _ in failed:
            print(f"  - {name}")
    
    print("\nOutput files:")
    for name, rmsd, output_name, ext in successful:
        fmt = "GRO" if ext == '.gro' else "PDB"
        print(f"  {output_name} ({fmt})")
    
    if args.verify and args.expected_dir:
        print("\n" + "=" * 60)
        print("VERIFICATION AGAINST EXPECTED FILES")
        print("=" * 60)
        
        verify_results = []
        for name, rmsd, output_name, ext in successful:
            if ext == '.gro':
                print(f"{name}: Skipping GRO verification (expected files are PDB)")
                continue
            
            # Expected files are named: hsa0.pdb_ali.pdb
            expected_name = name + '_ali.pdb'
            expected_path = args.expected_dir / expected_name
            
            if expected_path.exists():
                output_path = args.output_dir / output_name
                verify_rmsd, passed = verify_alignment(output_path, expected_path)
                verify_results.append((name, verify_rmsd, passed))
                status = 'PASS' if passed else 'FAIL'
                print(f"{name}: RMSD = {verify_rmsd:.6f} A [{status}]")
            else:
                print(f"{name}: Expected file not found: {expected_path}")
        
        if verify_results:
            passed_count = sum(1 for r in verify_results if r[2])
            print(f"\nVerification: {passed_count}/{len(verify_results)} passed")


if __name__ == '__main__':
    main()
