#!/usr/bin/env python3
"""
Script to select frames with:
- Most contacts
- Most hydrogen bonds
- Lowest ligand RMSD

for each system in symbolic links to the current directory.
"""

import os
import glob
import pandas as pd
from pathlib import Path


def find_systems():
    """Find all symbolic links in current directory"""
    systems = []
    for item in os.listdir('.'):
        if os.path.islink(item):
            systems.append(item)
    return sorted(systems)


def process_contacts(contacts_file):
    """
    Process contacts timeseries file.
    Sum all residue contacts to get total contacts per frame.
    Returns dataframe with frame, time_ps, total_contacts
    """
    try:
        df = pd.read_csv(contacts_file)
        
        if df.empty:
            return pd.DataFrame()
        
        # Sum all other columns (residue columns are numbered)
        residue_cols = [col for col in df.columns if col != 'time_ps']
        df['total_contacts'] = df[residue_cols].sum(axis=1)
        
        # Create frame number (assuming 100ps spacing)
        df['frame'] = df.index
        
        return df[['frame', 'time_ps', 'total_contacts']]
    except pd.errors.EmptyDataError:
        return pd.DataFrame()


def process_hbonds(hbond_file):
    """
    Process hydrogen bonds timeseries file.
    Returns dataframe with frame, time_ps, hbond_count
    """
    try:
        df = pd.read_csv(hbond_file)
        
        if df.empty:
            return pd.DataFrame()
        
        return df[['frame', 'time_ps', 'hbond_count']]
    except pd.errors.EmptyDataError:
        return pd.DataFrame()


def process_rmsd(rmsd_file):
    """
    Process ligand RMSD file.
    Returns dataframe with frame, time_ps, rmsd
    """
    try:
        df = pd.read_csv(rmsd_file)
        
        if df.empty:
            return pd.DataFrame()
        
        return df[['frame', 'time_ps', 'rmsd']]
    except pd.errors.EmptyDataError:
        return pd.DataFrame()


def process_system(system_path):
    """
    Process all trials in a system directory.
    Returns combined dataframe with contacts, hbonds, and rmsd.
    """
    all_contacts = []
    all_hbonds = []
    all_rmsd = []
    
    # Find all trial directories
    trial_dirs = sorted(glob.glob(os.path.join(system_path, 'per_trial', 'mmpbsa_*')))
    
    for trial_idx, trial_dir in enumerate(trial_dirs):
        trial_name = os.path.basename(trial_dir)
        
        # Process contacts
        contacts_file = os.path.join(trial_dir, 'contacts', 'contacts_timeseries.csv')
        if os.path.exists(contacts_file):
            contacts_df = process_contacts(contacts_file)
            if not contacts_df.empty:
                contacts_df['trial'] = trial_name
                contacts_df['global_frame'] = f"{trial_name}_" + contacts_df['frame'].astype(str)
                all_contacts.append(contacts_df)
        
        # Process hydrogen bonds
        hbond_file = os.path.join(trial_dir, 'hydrogen_bonds', 'hbond_timeseries.csv')
        if os.path.exists(hbond_file):
            hbond_df = process_hbonds(hbond_file)
            if not hbond_df.empty:
                hbond_df['trial'] = trial_name
                hbond_df['global_frame'] = f"{trial_name}_" + hbond_df['frame'].astype(str)
                all_hbonds.append(hbond_df)
        
        # Process RMSD
        rmsd_file = os.path.join(trial_dir, 'rmsd', 'rmsd_ligand.csv')
        if os.path.exists(rmsd_file):
            rmsd_df = process_rmsd(rmsd_file)
            if not rmsd_df.empty:
                rmsd_df['trial'] = trial_name
                rmsd_df['global_frame'] = f"{trial_name}_" + rmsd_df['frame'].astype(str)
                all_rmsd.append(rmsd_df)
    
    # Concatenate all trials
    combined_contacts = pd.concat(all_contacts, ignore_index=True) if all_contacts else pd.DataFrame()
    combined_hbonds = pd.concat(all_hbonds, ignore_index=True) if all_hbonds else pd.DataFrame()
    combined_rmsd = pd.concat(all_rmsd, ignore_index=True) if all_rmsd else pd.DataFrame()
    
    return combined_contacts, combined_hbonds, combined_rmsd


def find_best_frames(system_name, contacts_df, hbonds_df, rmsd_df):
    """
    Find frames with most contacts, most hbonds, and lowest rmsd.
    Returns dictionary with results.
    """
    results = {'system': system_name}
    
    # Most contacts
    if not contacts_df.empty:
        best_contact_idx = contacts_df['total_contacts'].idxmax()
        best_contact = contacts_df.loc[best_contact_idx]
        results['most_contacts_trial_id'] = best_contact['trial']
        results['most_contacts_frame_num'] = best_contact['frame']
        results['most_contacts_time_ps'] = best_contact['time_ps']
        results['most_contacts_value'] = best_contact['total_contacts']
    else:
        results['most_contacts_trial_id'] = 'N/A'
        results['most_contacts_frame_num'] = 'N/A'
        results['most_contacts_time_ps'] = 'N/A'
        results['most_contacts_value'] = 'N/A'
    
    # Most hydrogen bonds
    if not hbonds_df.empty:
        best_hbond_idx = hbonds_df['hbond_count'].idxmax()
        best_hbond = hbonds_df.loc[best_hbond_idx]
        results['most_hbonds_trial_id'] = best_hbond['trial']
        results['most_hbonds_frame_num'] = best_hbond['frame']
        results['most_hbonds_time_ps'] = best_hbond['time_ps']
        results['most_hbonds_value'] = best_hbond['hbond_count']
    else:
        results['most_hbonds_trial_id'] = 'N/A'
        results['most_hbonds_frame_num'] = 'N/A'
        results['most_hbonds_time_ps'] = 'N/A'
        results['most_hbonds_value'] = 'N/A'
    
    # Lowest RMSD (excluding frame 0 which is typically 0.0)
    if not rmsd_df.empty:
        # Filter out frame 0 which usually has RMSD of 0.0
        rmsd_filtered = rmsd_df[rmsd_df['frame'] > 0]
        if not rmsd_filtered.empty:
            best_rmsd_idx = rmsd_filtered['rmsd'].idxmin()
            best_rmsd = rmsd_filtered.loc[best_rmsd_idx]
            results['lowest_rmsd_trial_id'] = best_rmsd['trial']
            results['lowest_rmsd_frame_num'] = best_rmsd['frame']
            results['lowest_rmsd_time_ps'] = best_rmsd['time_ps']
            results['lowest_rmsd_value'] = best_rmsd['rmsd']
        else:
            results['lowest_rmsd_trial_id'] = 'N/A'
            results['lowest_rmsd_frame_num'] = 'N/A'
            results['lowest_rmsd_time_ps'] = 'N/A'
            results['lowest_rmsd_value'] = 'N/A'
    else:
        results['lowest_rmsd_trial_id'] = 'N/A'
        results['lowest_rmsd_frame_num'] = 'N/A'
        results['lowest_rmsd_time_ps'] = 'N/A'
        results['lowest_rmsd_value'] = 'N/A'
    
    return results


def main():
    """Main function to process all systems and output results."""
    print("Finding systems...")
    systems = find_systems()
    print(f"Found {len(systems)} systems: {', '.join(systems)}")
    
    all_results = []
    
    for system in systems:
        print(f"\nProcessing {system}...")
        
        # Process system
        contacts_df, hbonds_df, rmsd_df = process_system(system)
        
        # Find best frames
        results = find_best_frames(system, contacts_df, hbonds_df, rmsd_df)
        all_results.append(results)
        
        # Print summary
        print(f"  Most contacts: {results['most_contacts_value']} (trial: {results['most_contacts_trial_id']}, frame: {results['most_contacts_frame_num']})")
        print(f"  Most hbonds: {results['most_hbonds_value']} (trial: {results['most_hbonds_trial_id']}, frame: {results['most_hbonds_frame_num']})")
        if results['lowest_rmsd_value'] != 'N/A':
            print(f"  Lowest RMSD: {results['lowest_rmsd_value']:.3f} Å (trial: {results['lowest_rmsd_trial_id']}, frame: {results['lowest_rmsd_frame_num']})")
        else:
            print(f"  Lowest RMSD: N/A")
    
    # Create output dataframe
    output_df = pd.DataFrame(all_results)
    
    # Save to CSV
    output_file = 'best_frames_summary.csv'
    output_df.to_csv(output_file, index=False)
    print(f"\nResults saved to {output_file}")
    
    return output_df


if __name__ == '__main__':
    main()
