#!/usr/bin/env python3
"""
Extract PDB frames for best contacts, hbonds, and RMSD from trajectory files.
Uses gmx trjconv with index file for frame extraction.
"""

import os
import pandas as pd
import subprocess


def extract_frame(system_name, trial_id, frame_num, criterion, trj_dir):
    """
    Extract a specific frame from trajectory and save as PDB using gmx trjconv.
    
    Parameters:
    -----------
    system_name : str
        Name of the system
    trial_id : str
        Trial ID (e.g., 'mmpbsa_0')
    frame_num : int
        Frame number to extract
    criterion : str
        Type of criterion (contacts, hbonds, rmsd)
    trj_dir : str
        Path to trj directory
    """
    # Build paths
    system_trj_path = os.path.join(trj_dir, system_name)
    trial_path = os.path.join(system_trj_path, trial_id)
    tpr_file = os.path.join(system_trj_path, 'prod_0.tpr')  # Use production TPR
    xtc_file = os.path.join(trial_path, 'com_traj.xtc')
    ndx_file = os.path.join(system_trj_path, 'index.ndx')
    
    # Check if files exist
    if not os.path.exists(tpr_file):
        print(f"  Warning: TPR file not found: {tpr_file}")
        return False
    if not os.path.exists(xtc_file):
        print(f"  Warning: XTC file not found: {xtc_file}")
        return False
    
    # Output filename
    output_pdb = f"{system_name}_best_{criterion}.pdb"
    
    # Convert frame number to time in ps (assuming 100ps per frame)
    time_ps = frame_num * 100
    
    # Build gmx trjconv command
    cmd = [
        'gmx', 'trjconv',
        '-s', tpr_file,
        '-f', xtc_file,
        '-n', ndx_file,
        '-dump', str(time_ps),
        '-o', output_pdb
    ]
    
    try:
        # Run gmx trjconv with Protein_Other group
        # Protein_Other contains protein + ligand (non-solvent)
        result = subprocess.run(
            cmd,
            input='Protein_Other\n',  # Select Protein_Other group
            capture_output=True,
            text=True,
            timeout=30
        )
        
        if result.returncode != 0:
            print(f"  Warning: gmx trjconv failed")
            if result.stderr:
                print(f"  Error: {result.stderr[:200]}")
            return False
        
        if os.path.exists(output_pdb):
            print(f"  Saved {output_pdb} (trial: {trial_id}, frame: {frame_num})")
            return True
        else:
            print(f"  Warning: Output file not created")
            return False
            
    except subprocess.TimeoutExpired:
        print(f"  Warning: gmx trjconv timed out")
        return False
    except Exception as e:
        print(f"  Warning: Error running gmx trjconv: {e}")
        return False


def main():
    """Main function to extract all best frames."""
    # Read summary CSV
    csv_file = 'best_frames_summary.csv'
    if not os.path.exists(csv_file):
        print(f"Error: {csv_file} not found. Run select_best_frames.py first.")
        return
    
    df = pd.read_csv(csv_file)
    
    trj_dir = 'trj'
    
    print(f"Processing {len(df)} systems...")
    
    success_count = 0
    fail_count = 0
    
    for idx, row in df.iterrows():
        system = row['system']
        print(f"\nProcessing {system}...")
        
        # Check if system exists in trj directory
        system_trj_path = os.path.join(trj_dir, system)
        if not os.path.exists(system_trj_path):
            print(f"  Warning: System {system} not found in {trj_dir} directory")
            fail_count += 3
            continue
        
        # Extract frame with most contacts
        if pd.notna(row['most_contacts_trial_id']) and row['most_contacts_trial_id'] != 'N/A':
            if extract_frame(
                system,
                row['most_contacts_trial_id'],
                int(row['most_contacts_frame_num']),
                'contacts',
                trj_dir
            ):
                success_count += 1
            else:
                fail_count += 1
        
        # Extract frame with most hbonds
        if pd.notna(row['most_hbonds_trial_id']) and row['most_hbonds_trial_id'] != 'N/A':
            if extract_frame(
                system,
                row['most_hbonds_trial_id'],
                int(row['most_hbonds_frame_num']),
                'hbonds',
                trj_dir
            ):
                success_count += 1
            else:
                fail_count += 1
        
        # Extract frame with lowest RMSD
        if pd.notna(row['lowest_rmsd_trial_id']) and row['lowest_rmsd_trial_id'] != 'N/A':
            if extract_frame(
                system,
                row['lowest_rmsd_trial_id'],
                int(row['lowest_rmsd_frame_num']),
                'rmsd',
                trj_dir
            ):
                success_count += 1
            else:
                fail_count += 1
    
    print(f"\n{'='*60}")
    print(f"Extraction complete!")
    print(f"  Success: {success_count} frames")
    print(f"  Failed: {fail_count} frames")
    print(f"{'='*60}")


if __name__ == '__main__':
    main()
