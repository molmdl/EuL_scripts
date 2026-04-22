#!/usr/bin/env python3

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))

from modules.load_align import load_and_align_trajectory
from modules.rmsd import calculate_rmsd
from modules.contacts import calculate_contacts
from modules.hydrogen_bonds import calculate_hydrogen_bonds_wrapper
from modules.mmpbsa_plots import plot_mmpbsa_results
from utils.cli import parse_arguments
from utils.io import setup_output_directory


def main():
    args = parse_arguments()
    
    output_dir = setup_output_directory(args.output_dir)
    
    if args.aim in ['all', '0']:
        universe, reference_atoms, receptor, ligand = load_and_align_trajectory(
            topology=args.tpr,
            trajectory=args.xtc,
            index_file=args.ndx,
            align_selection=args.receptor_backbone_sel,
            receptor_sel=args.receptor_sel,
            ligand_sel=args.ligand_sel
        )
        print(f"Loaded and aligned trajectory: {universe.trajectory.n_frames} frames")
    
    if args.aim in ['all', '1']:
        rmsd_results = calculate_rmsd(
            universe=universe,
            receptor_sel=args.receptor_backbone_sel,
            ligand_sel=args.ligand_sel + ' and not name H*',
            output_dir=str(output_dir / 'rmsd')
        )
        print(f"RMSD analysis complete: {len(rmsd_results[0])} frames")
    
    if args.aim in ['all', '2']:
        contact_results = calculate_contacts(
            universe=universe,
            receptor_sel=args.receptor_sel,
            ligand_sel=args.ligand_sel,
            cutoff=args.contact_cutoff,
            output_dir=str(output_dir / 'contacts')
        )
        print(f"Contact analysis complete")
    
    if args.aim in ['all', '3']:
        hbond_results = calculate_hydrogen_bonds_wrapper(
            universe=universe,
            receptor_sel=args.receptor_sel,
            ligand_sel=args.ligand_sel,
            distance_cutoff=args.hbond_distance_cutoff,
            angle_cutoff=args.hbond_angle_cutoff,
            output_dir=str(output_dir / 'hydrogen_bonds')
        )
        print(f"Hydrogen bond analysis complete")
    
    if args.aim in ['all', '4']:
        plot_mmpbsa_results(
            mmpbsa_dir=args.mmpbsa_dir,
            output_dir=str(output_dir / 'mmpbsa')
        )
        print(f"MMPBSA plots generated")
    
    print(f"\nAll analyses complete. Output saved to: {output_dir}")


if __name__ == '__main__':
    main()