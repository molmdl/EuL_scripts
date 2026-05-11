#!/usr/bin/env python3

import argparse
import warnings
from pathlib import Path
from collections import defaultdict
from io import StringIO

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
import MDAnalysis as mda
from MDAnalysis.analysis import distances
from MDAnalysis.analysis.hydrogenbonds.hbond_analysis import HydrogenBondAnalysis as HBA
from MDAnalysis.coordinates.XTC import XTCReader

DEFAULT_TOPOLOGY = 'com.tpr'
DEFAULT_TRAJECTORY = 'com_traj.xtc'
DEFAULT_INDEX_FILE = 'index.ndx'
DEFAULT_MMPBSA_DIR = '.'
DEFAULT_RECEPTOR_BACKBONE_SEL = 'protein and name N CA C O'
DEFAULT_RECEPTOR_SEL = 'protein'
DEFAULT_LIGAND_SEL = 'resname MOL'
DEFAULT_CONTACT_CUTOFF = 4.5
DEFAULT_HBOND_DISTANCE_CUTOFF = 3.5
DEFAULT_HBOND_ANGLE_CUTOFF = 150
DEFAULT_OUTPUT_DIR = 'output'
DEFAULT_PLOT_DPI = 300
DEFAULT_PLOT_FORMAT = 'png'
DEFAULT_WORKSPACE = None
DEFAULT_TRIALS_FILE = 'trials.txt'

sns.set_theme(style='darkgrid')

GMX_MMPBSA_STYLE = {
    'Line Plot': {
        'line_width': 0.7,
        'line_color': '#000000',
        'rolling_avg_color': '#FF0000',
        'rolling_avg_style': 'dashed',
        'rolling_avg_width': 0.8,
        'figure_width': 8,
        'figure_height': 4,
        'fontsize': {'x_label': 12, 'y_label': 12, 'title': 14, 'legend': 9, 'x_ticks': 10, 'y_ticks': 10},
    },
    'Bar Plot': {
        'palette': 'husl',
        'color': '#2C69B0',
        'error_line_width': 0.7,
        'error_line_color': '#000000',
        'error_cap_size': 0,
        'figure_width': 6,
        'figure_height': 5,
        'fontsize': {'x_label': 12, 'y_label': 12, 'title': 14, 'legend': 9, 'x_ticks': 10, 'y_ticks': 10},
        'x_rotation': 45,
        'y_inverted': True,
    },
    'Heatmap Plot': {
        'palette': 'seismic',
        'figure_width': 10,
        'figure_height': 7,
        'fontsize': {'x_label': 12, 'y_label': 12, 'title': 14, 'legend': 9, 'x_ticks': 10, 'y_ticks': 10},
        'x_rotation': 0,
    },
}

PALETTE = {
    'primary': '#2C69B0',
    'secondary': '#E63946',
    'receptor': '#2C69B0',
    'ligand': '#E63946',
    'favorable': '#2C69B0',
    'unfavorable': '#E63946',
}

def parse_arguments():
    parser = argparse.ArgumentParser(
        description='Comprehensive protein-ligand trajectory analysis',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  Single analysis:
    python com_ana_trj.py --tpr com.tpr --xtc com_traj.xtc
    python com_ana_trj.py --ana rmsd --tpr com.tpr --xtc com_traj.xtc
    
  Multi-trial (single ligand):
    python com_ana_trj.py --workspace ./ligand1 --output-dir ./results
    
  Multi-ligand:
    python com_ana_trj.py --workspace ./workspace --output-dir ./results
    
Directory structure for multi-trial/multi-ligand:
  workspace/
    ligand1/
      trials.txt     # Contains: trial1, trial2, trial3
      trial1/
        com.tpr, com_traj.xtc, FINAL_RESULTS_MMPBSA.csv
      trial2/
        ...
    ligand2/
      trials.txt
      ...
        """
    )
    
    parser.add_argument('--tpr', type=str, default=DEFAULT_TOPOLOGY,
                        help=f"Path to TPR topology file (default: {DEFAULT_TOPOLOGY})")
    parser.add_argument('--xtc', type=str, default=DEFAULT_TRAJECTORY,
                        help=f"Path to XTC trajectory file (default: {DEFAULT_TRAJECTORY})")
    parser.add_argument('--ndx', type=str, default=None,
                        help="Path to GROMACS index file (default: None)")
    
    parser.add_argument('--workspace', type=str, default=DEFAULT_WORKSPACE,
                        help="Path to workspace directory containing ligand subdirectories with trials.txt")
    parser.add_argument('--trials-file', type=str, default=DEFAULT_TRIALS_FILE,
                        help=f"Name of trials file in each ligand directory (default: {DEFAULT_TRIALS_FILE})")
    
    parser.add_argument('--receptor-backbone-sel', type=str,
                        default=DEFAULT_RECEPTOR_BACKBONE_SEL,
                        help=f"Selection for receptor backbone (default: {DEFAULT_RECEPTOR_BACKBONE_SEL})")
    parser.add_argument('--receptor-sel', type=str, default=DEFAULT_RECEPTOR_SEL,
                        help=f"Selection for receptor (default: {DEFAULT_RECEPTOR_SEL})")
    parser.add_argument('--ligand-sel', type=str, default=DEFAULT_LIGAND_SEL,
                        help=f"Selection for ligand (default: {DEFAULT_LIGAND_SEL})")
    
    parser.add_argument('--contact-cutoff', type=float, default=DEFAULT_CONTACT_CUTOFF,
                        help=f"Contact distance cutoff in Angstroms (default: {DEFAULT_CONTACT_CUTOFF})")
    parser.add_argument('--hbond-distance-cutoff', type=float, default=DEFAULT_HBOND_DISTANCE_CUTOFF,
                        help=f"Hydrogen bond distance cutoff in Angstroms (default: {DEFAULT_HBOND_DISTANCE_CUTOFF})")
    parser.add_argument('--hbond-angle-cutoff', type=float, default=DEFAULT_HBOND_ANGLE_CUTOFF,
                        help=f"Hydrogen bond angle cutoff in degrees (default: {DEFAULT_HBOND_ANGLE_CUTOFF})")
    
    parser.add_argument('--output-dir', type=str, default=DEFAULT_OUTPUT_DIR,
                        help=f"Output directory (default: {DEFAULT_OUTPUT_DIR})")
    parser.add_argument('--mmpbsa-dir', type=str, default=DEFAULT_MMPBSA_DIR,
                        help=f"Directory containing MMPBSA output files (default: {DEFAULT_MMPBSA_DIR})")
    
    parser.add_argument('--ana', type=str, default='all',
                        choices=['all', 'load', 'rmsd', 'contacts', 'hbonds', 'mmpbsa'],
                        help='Analysis to run: load=load/align, rmsd=RMSD, contacts=per-residue contacts, hbonds=hydrogen bonds, mmpbsa=MMPBSA plots (default: all)')
    
    parser.add_argument('--plot-dpi', type=int, default=DEFAULT_PLOT_DPI,
                        help=f"DPI for output plots (default: {DEFAULT_PLOT_DPI})")
    parser.add_argument('--plot-format', type=str, default=DEFAULT_PLOT_FORMAT,
                        choices=['png', 'pdf', 'svg'],
                        help=f"Output format for plots (default: {DEFAULT_PLOT_FORMAT})")
    
    parser.add_argument('--skip-align', action='store_true',
                        help='Skip alignment (trajectory already aligned)')
    parser.add_argument('--n-jobs', type=int, default=1,
                        help='Number of parallel jobs (default: 1)')
    
    return parser.parse_args()


def setup_output_directory(output_dir):
    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)
    return output_path


def load_universe(topology, trajectory, index_file=None):
    try:
        if index_file and Path(index_file).exists():
            u = mda.Universe(topology, trajectory, index=index_file)
        else:
            u = mda.Universe(topology, trajectory)
        return u
    except ValueError as e:
        if "don't have the same number of atoms" in str(e):
            print(f"Atom count mismatch detected. Creating universe from topology first...")
            
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                
                if index_file and Path(index_file).exists():
                    u_topo = mda.Universe(topology, index=index_file)
                else:
                    u_topo = mda.Universe(topology)
                
                n_atoms_topo = u_topo.atoms.n_atoms
                
                traj_reader = XTCReader(trajectory)
                n_atoms_traj = traj_reader.n_atoms
                traj_reader.close()
                
                print(f"Topology has {n_atoms_topo} atoms, trajectory has {n_atoms_traj} atoms")
                print(f"Selecting first {n_atoms_traj} atoms from topology...")
                
                u_select = mda.Universe(topology, in_memory=False)
                u_select.atoms = u_select.atoms[:n_atoms_traj]
                u_select.load_new(trajectory)
                
                print(f"Successfully loaded universe with {len(u_select.atoms)} atoms")
                return u_select
        raise


def _kabsch_rotate(P, Q):
    H = P.T @ Q
    U, S, Vt = np.linalg.svd(H)
    d = np.linalg.det(Vt.T @ U.T)
    sign_matrix = np.eye(3)
    sign_matrix[2, 2] = d
    R = Vt.T @ sign_matrix @ U.T
    return R


def align_trajectory(universe, reference_sel='protein and name N CA C O'):
    reference_atoms = universe.select_atoms(reference_sel)
    
    if len(reference_atoms) == 0:
        print(f"Warning: No atoms found for selection '{reference_sel}'")
        print("Attempting alternative backbone selection...")
        reference_atoms = universe.select_atoms('name N CA C O')
    
    print(f"Using {len(reference_atoms)} atoms for alignment")
    print(f"Selection: {reference_sel}")
    
    ref_positions = reference_atoms.positions.copy()
    ref_com = ref_positions.mean(axis=0)
    ref_centered = ref_positions - ref_com
    
    def align_to_reference(ts):
        mobile_positions = reference_atoms.positions
        mobile_com = mobile_positions.mean(axis=0)
        mobile_centered = mobile_positions - mobile_com
        R = _kabsch_rotate(mobile_centered, ref_centered)
        aligned = (ts.positions - mobile_com) @ R + ref_com
        ts.positions[:] = aligned
        return ts
    
    universe.trajectory.add_transformations(align_to_reference)
    
    return universe, reference_atoms


def get_selections(universe, receptor_sel='protein', ligand_sel='resname MOL'):
    receptor = universe.select_atoms(receptor_sel)
    ligand = universe.select_atoms(ligand_sel)
    
    print(f"Receptor: {len(receptor)} atoms ({receptor.n_residues} residues)")
    print(f"Ligand: {len(ligand)} atoms")
    
    return receptor, ligand


def load_and_align_trajectory(topology, trajectory, index_file=None,
                              align_selection='protein and name N CA C O',
                              receptor_sel='protein', ligand_sel='resname MOL'):
    universe = load_universe(topology, trajectory, index_file)
    universe, reference_atoms = align_trajectory(universe, align_selection)
    receptor, ligand = get_selections(universe, receptor_sel, ligand_sel)
    
    return universe, reference_atoms, receptor, ligand


def calculate_receptor_rmsd(universe, selection='protein and name N CA C O'):
    backbone = universe.select_atoms(selection)
    reference = backbone.positions.copy()
    
    results = []
    for ts in universe.trajectory:
        rmsd = np.sqrt(np.mean(np.sum((backbone.positions - reference)**2, axis=1)))
        results.append({
            'frame': ts.frame,
            'time_ps': ts.time,
            'rmsd': rmsd
        })
    
    return pd.DataFrame(results)


def calculate_ligand_rmsd(universe, selection='resname MOL and not name H*'):
    ligand = universe.select_atoms(selection)
    
    if len(ligand) == 0:
        print("Warning: No ligand atoms found with selection, trying all ligand atoms")
        ligand = universe.select_atoms('resname MOL')
    
    reference = ligand.positions.copy()
    
    results = []
    for ts in universe.trajectory:
        rmsd = np.sqrt(np.mean(np.sum((ligand.positions - reference)**2, axis=1)))
        results.append({
            'frame': ts.frame,
            'time_ps': ts.time,
            'rmsd': rmsd
        })
    
    return pd.DataFrame(results)


def plot_rmsd_overlay(df_receptor, df_ligand, output_path, dpi=300):
    style = GMX_MMPBSA_STYLE['Line Plot']
    fig, ax = plt.subplots(figsize=(style['figure_width'], style['figure_height']))
    
    sns.lineplot(data=df_receptor, x='time_ps', y='rmsd',
                 color=PALETTE['receptor'], linewidth=style['line_width'],
                 label='Receptor Backbone', ax=ax)
    sns.lineplot(data=df_ligand, x='time_ps', y='rmsd',
                 color=PALETTE['ligand'], linewidth=style['line_width'],
                 label='Ligand Heavy Atoms', ax=ax)
    
    window = min(5, len(df_receptor))
    if window > 1:
        ma_rec = df_receptor['rmsd'].rolling(window=window, center=True).mean()
        ma_lig = df_ligand['rmsd'].rolling(window=window, center=True).mean()
        sns.lineplot(x=df_receptor['time_ps'], y=ma_rec,
                     color=PALETTE['receptor'], linewidth=style['rolling_avg_width'],
                     linestyle=style['rolling_avg_style'], label='Receptor (MA)', ax=ax)
        sns.lineplot(x=df_ligand['time_ps'], y=ma_lig,
                     color=PALETTE['ligand'], linewidth=style['rolling_avg_width'],
                     linestyle=style['rolling_avg_style'], label='Ligand (MA)', ax=ax)
    
    ax.set_xlabel('Time (ps)', fontsize=style['fontsize']['x_label'])
    ax.set_ylabel('RMSD (A)', fontsize=style['fontsize']['y_label'])
    ax.set_title('RMSD: Receptor Backbone vs Ligand Heavy Atoms', fontsize=style['fontsize']['title'])
    ax.legend(loc='upper left', fontsize=style['fontsize']['legend'], bbox_to_anchor=(1.02, 1), borderaxespad=0)
    
    text = (f"Receptor: mean={df_receptor['rmsd'].mean():.2f}A, std={df_receptor['rmsd'].std():.2f}A\n"
            f"Ligand: mean={df_ligand['rmsd'].mean():.2f}A, std={df_ligand['rmsd'].std():.2f}A")
    ax.text(1.02, 0.02, text, transform=ax.transAxes, fontsize=9,
            verticalalignment='bottom', bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))
    
    plt.tight_layout()
    plt.savefig(output_path, dpi=dpi, bbox_inches='tight')
    plt.close()


def calculate_rmsd(universe, receptor_sel='protein and name N CA C O',
                   ligand_sel='resname MOL and not name H*',
                   output_dir='output/rmsd'):
    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)
    
    print("Calculating receptor backbone RMSD...")
    df_receptor = calculate_receptor_rmsd(universe, receptor_sel)
    
    print("Calculating ligand heavy atom RMSD...")
    df_ligand = calculate_ligand_rmsd(universe, ligand_sel)
    
    df_receptor.to_csv(output_path / 'rmsd_receptor.csv', index=False)
    df_ligand.to_csv(output_path / 'rmsd_ligand.csv', index=False)
    
    df_combined = pd.merge(df_receptor, df_ligand, on=['frame', 'time_ps'], 
                           suffixes=('_receptor', '_ligand'))
    df_combined.to_csv(output_path / 'rmsd_combined.csv', index=False)
    
    plot_rmsd_overlay(df_receptor, df_ligand, str(output_path / 'rmsd_plot.png'))
    
    return df_receptor, df_ligand


def calculate_per_residue_contacts(universe, receptor_sel='protein',
                                   ligand_sel='resname MOL', cutoff=4.5):
    receptor = universe.select_atoms(receptor_sel)
    ligand = universe.select_atoms(ligand_sel)
    residues = receptor.residues
    
    res_atom_indices = {}
    for res in residues:
        res_atom_indices[res.resid] = res.atoms.indices
    
    contact_timeseries = {res.resid: [] for res in residues}
    times = []
    
    for ts in universe.trajectory:
        times.append(ts.time)
        
        dist_mat = distances.distance_array(
            receptor.positions,
            ligand.positions,
            box=universe.dimensions
        )
        
        for resid, atom_indices in res_atom_indices.items():
            res_distances = dist_mat[atom_indices, :]
            n_contacts = (res_distances <= cutoff).sum()
            contact_timeseries[resid].append(n_contacts)
    
    return contact_timeseries, times


def save_contacts_csv(contact_timeseries, times, output_dir):
    df = pd.DataFrame(contact_timeseries, index=times)
    df.index.name = 'time_ps'
    df.to_csv(output_dir / 'contacts_timeseries.csv')
    
    summary_data = []
    for resid, contacts in contact_timeseries.items():
        contacts_arr = np.array(contacts)
        summary_data.append({
            'residue_id': resid,
            'mean_contacts': np.mean(contacts_arr),
            'std_contacts': np.std(contacts_arr),
            'max_contacts': np.max(contacts_arr),
            'min_contacts': np.min(contacts_arr),
            'fraction_in_contact': np.mean(contacts_arr > 0)
        })
    
    df_summary = pd.DataFrame(summary_data)
    df_summary = df_summary.sort_values('mean_contacts', ascending=False)
    df_summary.to_csv(output_dir / 'contacts_summary.csv', index=False)
    
    df_total = pd.DataFrame({'time_ps': times, 
                             'total_contacts': df.sum(axis=1).values})
    df_total.to_csv(output_dir / 'contacts_total.csv', index=False)


def plot_contact_heatmap(contact_timeseries, times, output_path, top_n=50, dpi=300):
    style = GMX_MMPBSA_STYLE['Heatmap Plot']
    df = pd.DataFrame(contact_timeseries, index=times)
    
    avg_contacts = df.mean().sort_values(ascending=False).head(top_n)
    df_top = df[avg_contacts.index]
    
    fig, ax = plt.subplots(figsize=(style['figure_width'], style['figure_height']))
    
    vmax = max(df_top.values.max(), abs(df_top.values.min()))
    im = ax.imshow(df_top.values.T, aspect='auto', origin='lower',
                   cmap=style['palette'], interpolation='nearest', vmin=0, vmax=vmax)
    
    ax.set_xlabel('Frame', fontsize=style['fontsize']['x_label'])
    ax.set_ylabel('Residue ID', fontsize=style['fontsize']['y_label'])
    ax.set_title(f'Receptor-Ligand Contact Heatmap (Top {top_n} Residues)', fontsize=style['fontsize']['title'])
    
    cbar = plt.colorbar(im, ax=ax, shrink=0.8)
    cbar.set_label('Number of Contacts', fontsize=style['fontsize']['y_label'])
    
    n_xticks = min(10, len(times))
    xtick_positions = np.linspace(0, len(times)-1, n_xticks).astype(int)
    xtick_labels = [f'{times[i]:.0f}' for i in xtick_positions]
    ax.set_xticks(xtick_positions)
    ax.set_xticklabels(xtick_labels, fontsize=style['fontsize']['x_ticks'], rotation=0)
    
    n_yticks = min(20, len(df_top.columns))
    ytick_positions = np.linspace(0, len(df_top.columns)-1, n_yticks).astype(int)
    ax.set_yticks(ytick_positions)
    ax.set_yticklabels([str(df_top.columns[i]) for i in ytick_positions], fontsize=style['fontsize']['y_ticks'])
    
    plt.tight_layout()
    plt.savefig(output_path, dpi=dpi, bbox_inches='tight')
    plt.close()


def plot_contact_bar(contact_timeseries, output_path, top_n=20, dpi=300):
    style = GMX_MMPBSA_STYLE['Bar Plot']
    df = pd.DataFrame(contact_timeseries)
    avg_contacts = df.mean().sort_values(ascending=False).head(top_n)
    
    fig, ax = plt.subplots(figsize=(style['figure_width'], style['figure_height']))
    
    y_pos = np.arange(len(avg_contacts))
    palette = sns.color_palette(style['palette'], n_colors=len(avg_contacts))
    
    sns.barplot(x=avg_contacts.values, y=y_pos, palette=palette, ax=ax, orient='h')
    
    ax.set_yticks(y_pos)
    ax.set_yticklabels(avg_contacts.index, fontsize=style['fontsize']['y_ticks'])
    ax.set_xlabel('Average Contact Number', fontsize=style['fontsize']['x_label'])
    ax.set_title(f'Top {top_n} Residues by Average Contacts', fontsize=style['fontsize']['title'])
    if style['y_inverted']:
        ax.invert_yaxis()
    
    plt.tight_layout()
    plt.savefig(output_path, dpi=dpi, bbox_inches='tight')
    plt.close()


def calculate_contacts(universe, receptor_sel='protein', ligand_sel='resname MOL',
                       cutoff=4.5, output_dir='output/contacts'):
    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)
    
    print(f"Calculating contacts with {cutoff}A cutoff...")
    contact_timeseries, times = calculate_per_residue_contacts(
        universe, receptor_sel, ligand_sel, cutoff
    )
    
    save_contacts_csv(contact_timeseries, times, output_path)
    
    plot_contact_heatmap(contact_timeseries, times, 
                         str(output_path / 'contacts_heatmap.png'))
    plot_contact_bar(contact_timeseries, 
                     str(output_path / 'contacts_bar.png'))
    
    return contact_timeseries


def calculate_hydrogen_bonds(universe, receptor_sel='protein', ligand_sel='resname MOL',
                             distance_cutoff=3.5, angle_cutoff=150):
    print(f"Running hydrogen bond analysis (dist={distance_cutoff}A, angle={angle_cutoff}deg)...")
    
    ligand_name = ligand_sel.replace('resname ', '').replace(' and not name H*', '')
    
    hbonds = HBA(
        universe=universe,
        between=[receptor_sel, ligand_sel],
        d_a_cutoff=distance_cutoff,
        d_h_a_angle_cutoff=angle_cutoff
    )
    
    hbonds.run(verbose=False)
    
    hbond_data = hbonds.results.hbonds
    
    if len(hbond_data) == 0:
        print("Warning: No hydrogen bonds found!")
        return pd.DataFrame(), pd.DataFrame()
    
    atom_to_residue = {}
    for atom in universe.atoms:
        atom_to_residue[atom.index] = {
            'resid': atom.resid,
            'resname': atom.resname,
            'resnum': atom.resnum
        }
    
    frames = np.unique(hbond_data[:, 0])
    time_series = []
    
    for frame in frames:
        frame_hbonds = hbond_data[hbond_data[:, 0] == frame]
        n_hbonds = len(frame_hbonds)
        
        universe.trajectory[int(frame)]
        time_ps = universe.trajectory.time
        
        time_series.append({
            'frame': int(frame),
            'time_ps': time_ps,
            'hbond_count': n_hbonds
        })
    
    df_time = pd.DataFrame(time_series)
    
    residue_data = defaultdict(lambda: {
        'count': 0,
        'distances': [],
        'angles': []
    })
    
    for hbond in hbond_data:
        frame, donor_idx, hydrogen_idx, acceptor_idx, distance, angle = hbond
        
        for idx in [donor_idx, acceptor_idx]:
            res_info = atom_to_residue.get(idx)
            if res_info and res_info['resname'] != ligand_name:
                key = f"{res_info['resname']}{res_info['resnum']}"
                residue_data[key]['count'] += 1
                residue_data[key]['distances'].append(distance)
                residue_data[key]['angles'].append(angle)
    
    residue_rows = []
    for key, data in residue_data.items():
        residue_rows.append({
            'residue': key,
            'hbond_count': data['count'],
            'avg_distance_A': np.mean(data['distances']) if data['distances'] else 0,
            'avg_angle_deg': np.mean(data['angles']) if data['angles'] else 0
        })
    
    df_residue = pd.DataFrame(residue_rows)
    df_residue = df_residue.sort_values('hbond_count', ascending=False)
    
    return df_time, df_residue


def plot_hydrogen_bonds(df_time, df_residue, output_path, dpi=300):
    style = GMX_MMPBSA_STYLE['Bar Plot']
    line_style = GMX_MMPBSA_STYLE['Line Plot']
    
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(style['figure_width'], style['figure_height'] * 1.5))
    
    sns.lineplot(data=df_time, x='time_ps', y='hbond_count',
                 color=line_style['line_color'], linewidth=line_style['line_width'],
                 ax=ax1, label='H-bond count')
    
    ax1.fill_between(df_time['time_ps'], 0, df_time['hbond_count'], alpha=0.3, color=PALETTE['primary'])
    
    window = min(5, len(df_time))
    if window > 1:
        ma = df_time['hbond_count'].rolling(window=window, center=True).mean()
        sns.lineplot(x=df_time['time_ps'], y=ma,
                     color=line_style['rolling_avg_color'], linewidth=line_style['rolling_avg_width'],
                     linestyle=line_style['rolling_avg_style'], label='Mov. Av.', ax=ax1)
    
    ax1.set_xlabel('Time (ps)', fontsize=line_style['fontsize']['x_label'])
    ax1.set_ylabel('Hydrogen Bond Count', fontsize=line_style['fontsize']['y_label'])
    ax1.set_title('Receptor-Ligand Hydrogen Bonds Over Time', fontsize=line_style['fontsize']['title'])
    ax1.legend(loc='upper right', fontsize=line_style['fontsize']['legend'])
    
    top_n = 15
    top_res = df_residue.head(top_n)
    
    y_pos = np.arange(len(top_res))
    palette = sns.color_palette(style['palette'], n_colors=len(top_res))
    
    sns.barplot(x=top_res['hbond_count'], y=y_pos, palette=palette, ax=ax2, orient='h')
    
    ax2.set_yticks(y_pos)
    ax2.set_yticklabels(top_res['residue'], fontsize=style['fontsize']['y_ticks'])
    ax2.set_xlabel('Total Hydrogen Bond Count', fontsize=style['fontsize']['x_label'])
    ax2.set_title(f'Top {top_n} Residues by Hydrogen Bond Frequency', fontsize=style['fontsize']['title'])
    if style['y_inverted']:
        ax2.invert_yaxis()
    
    plt.tight_layout()
    plt.savefig(output_path, dpi=dpi, bbox_inches='tight')
    plt.close()


def calculate_hydrogen_bonds_wrapper(universe, receptor_sel='protein', ligand_sel='resname MOL',
                                      distance_cutoff=3.5, angle_cutoff=150,
                                      output_dir='output/hydrogen_bonds'):
    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)
    
    df_time, df_residue = calculate_hydrogen_bonds(
        universe, receptor_sel, ligand_sel, distance_cutoff, angle_cutoff
    )
    
    df_time.to_csv(output_path / 'hbond_timeseries.csv', index=False)
    df_residue.to_csv(output_path / 'hbond_per_residue.csv', index=False)
    
    if len(df_time) > 0:
        plot_hydrogen_bonds(df_time, df_residue, str(output_path / 'hbond_plot.png'))
    
    return df_time, df_residue


def parse_results_csv(csv_path):
    with open(csv_path, 'r') as f:
        content = f.read()
    
    sections = content.split('\n\n')
    data = {}
    
    for section in sections:
        lines = [l.strip() for l in section.strip().split('\n') if l.strip()]
        if not lines:
            continue
            
        section_title = lines[0]
        if section_title == 'GENERALIZED BORN:' and len(lines) > 1:
            section_title = lines[1]
        
        if section_title == 'Complex Energy Terms':
            header_idx = 0
            for i, line in enumerate(lines):
                if line.startswith('Frame #'):
                    header_idx = i
                    break
            data['complex'] = pd.read_csv(StringIO('\n'.join(lines[header_idx:])))
        elif section_title == 'Receptor Energy Terms':
            data['receptor'] = pd.read_csv(StringIO('\n'.join(lines[1:])))
        elif section_title == 'Ligand Energy Terms':
            data['ligand'] = pd.read_csv(StringIO('\n'.join(lines[1:])))
        elif section_title == 'Delta Energy Terms':
            data['delta'] = pd.read_csv(StringIO('\n'.join(lines[1:])))
    
    return data


def parse_decomposition_csv(csv_path):
    with open(csv_path, 'r') as f:
        lines = f.readlines()
    
    data = {'complex': [], 'delta': []}
    current_section = None
    header_found = False
    
    for line in lines:
        line = line.strip()
        
        if line == 'Complex:':
            current_section = 'complex'
            header_found = False
            continue
        elif line == 'Receptor:':
            current_section = 'receptor'
            header_found = False
            continue
        elif line == 'Ligand:':
            current_section = 'ligand'
            header_found = False
            continue
        elif 'Delta' in line and ':' in line:
            current_section = 'delta'
            header_found = False
            continue
        
        if line.startswith('Frame #') or line.startswith('"Frame'):
            header_found = True
            continue
        
        if not line or not header_found:
            continue
        
        parts = line.split(',')
        if len(parts) >= 8:
            try:
                frame = int(parts[0])
                residue = parts[1]
                internal = float(parts[2])
                vdw = float(parts[3])
                electrostatic = float(parts[4])
                polar = float(parts[5])
                nonpolar = float(parts[6])
                total = float(parts[7])
                
                if current_section:
                    if current_section not in data:
                        data[current_section] = []
                    data[current_section].append({
                        'frame': frame,
                        'residue': residue,
                        'internal': internal,
                        'vdw': vdw,
                        'electrostatic': electrostatic,
                        'polar': polar,
                        'nonpolar': nonpolar,
                        'total': total
                    })
            except (ValueError, IndexError):
                continue
    
    if 'complex' in data and len(data['complex']) > 0:
        if 'delta' not in data or len(data.get('delta', [])) == 0:
            data['delta'] = data['complex']
    
    return {k: pd.DataFrame(v) for k, v in data.items() if v}


def plot_energy_timecourse(df_delta, output_path, dpi=300):
    style = GMX_MMPBSA_STYLE['Line Plot']
    fig, ax = plt.subplots(figsize=(style['figure_width'], style['figure_height']))
    
    frame_col = None
    for col in df_delta.columns:
        if 'frame' in col.lower():
            frame_col = col
            break
    if frame_col is None:
        frame_col = df_delta.columns[0]
    
    frames = df_delta[frame_col]
    terms = ['VDWAALS', 'EEL', 'EGB', 'ESURF', 'TOTAL']
    
    palette = sns.color_palette('husl', n_colors=len(terms))
    
    for i, term in enumerate(terms):
        if term in df_delta.columns:
            sns.lineplot(x=frames, y=df_delta[term], label=term,
                        color=palette[i], linewidth=style['line_width'], ax=ax)
    
    if 'TOTAL' in df_delta.columns:
        window = min(5, len(df_delta))
        ma = df_delta['TOTAL'].rolling(window=window, center=True).mean()
        sns.lineplot(x=frames, y=ma, color=style['rolling_avg_color'],
                    linewidth=style['rolling_avg_width'], linestyle=style['rolling_avg_style'],
                    label='Mov. Av.', ax=ax)
    
    ax.set_xlabel('Frame', fontsize=style['fontsize']['x_label'])
    ax.set_ylabel('Energy (kcal/mol)', fontsize=style['fontsize']['y_label'])
    ax.set_title('Binding Free Energy Decomposition', fontsize=style['fontsize']['title'])
    ax.legend(loc='upper left', fontsize=style['fontsize']['legend'], bbox_to_anchor=(1.02, 1), borderaxespad=0)
    ax.axhline(y=0, color='gray', linestyle='-', lw=0.5)
    
    plt.tight_layout()
    plt.savefig(output_path, dpi=dpi, bbox_inches='tight')
    plt.close()


def plot_energy_summary(df_delta, output_path, dpi=300):
    style = GMX_MMPBSA_STYLE['Bar Plot']
    terms = ['VDWAALS', 'EEL', 'EGB', 'ESURF', 'GGAS', 'GSOLV', 'TOTAL']
    means = []
    stds = []
    labels = []
    
    for term in terms:
        col = term if term in df_delta.columns else f'D{term}'
        if col in df_delta.columns:
            means.append(df_delta[col].mean())
            stds.append(df_delta[col].std())
            labels.append(term)
    
    fig, ax = plt.subplots(figsize=(style['figure_width'], style['figure_height']))
    
    x = np.arange(len(labels))
    palette = sns.color_palette('husl', n_colors=len(labels))
    
    bars = ax.bar(x, means, yerr=stds, capsize=style['error_cap_size'],
                  color=palette, alpha=0.85,
                  error_kw=dict(elinewidth=style['error_line_width'],
                               ecolor=style['error_line_color'],
                               capsize=style['error_cap_size']))
    
    ax.set_xticks(x)
    ax.set_xticklabels(labels, rotation=style['x_rotation'], ha='right', fontsize=style['fontsize']['x_ticks'])
    ax.set_ylabel('Energy (kcal/mol)', fontsize=style['fontsize']['y_label'])
    ax.set_title('Binding Energy Summary', fontsize=style['fontsize']['title'])
    ax.axhline(y=0, color='black', lw=0.5)
    
    plt.tight_layout()
    plt.savefig(output_path, dpi=dpi, bbox_inches='tight')
    plt.close()


def plot_decomposition_bar(df_delta, output_path, top_n=15, dpi=300):
    if df_delta.empty:
        print("No decomposition data to plot")
        return
    
    style = GMX_MMPBSA_STYLE['Bar Plot']
    residue_avg = df_delta.groupby('residue')['total'].mean().sort_values()
    top = residue_avg.head(top_n)
    
    labels = []
    for res in top.index:
        parts = res.split(':')
        if len(parts) >= 4:
            labels.append(f"{parts[2][0]}{parts[3]}")
        else:
            labels.append(res)
    
    fig, ax = plt.subplots(figsize=(style['figure_width'], style['figure_height']))
    
    colors = [PALETTE['favorable'] if v < 0 else PALETTE['unfavorable'] for v in top.values]
    y = np.arange(len(labels))
    
    sns.barplot(x=top.values, y=y, palette=colors, ax=ax, orient='h')
    
    ax.set_yticks(y)
    ax.set_yticklabels(labels, fontsize=style['fontsize']['y_ticks'])
    ax.set_xlabel('Energy (kcal/mol)', fontsize=style['fontsize']['x_label'])
    ax.set_title('Per-Residue Binding Contribution', fontsize=style['fontsize']['title'])
    ax.axvline(x=0, color='black', lw=0.5)
    if style['y_inverted']:
        ax.invert_yaxis()
    
    plt.tight_layout()
    plt.savefig(output_path, dpi=dpi, bbox_inches='tight')
    plt.close()


def plot_decomposition_heatmap(df_delta, output_path, top_n=30, dpi=300):
    if df_delta.empty:
        return
    
    style = GMX_MMPBSA_STYLE['Heatmap Plot']
    
    try:
        pivot = df_delta.pivot_table(index='frame', columns='residue', values='total', aggfunc='mean')
    except Exception:
        print("Warning: Could not create decomposition heatmap")
        return
    
    if pivot.empty or pivot.shape[1] == 0:
        print("Warning: No data for decomposition heatmap")
        return
    
    top_cols = pivot.var().sort_values(ascending=False).head(top_n).index
    pivot_top = pivot[top_cols]
    
    short_cols = {}
    for col in pivot_top.columns:
        parts = col.split(':')
        if len(parts) >= 4:
            short_cols[col] = f"{parts[2][0]}{parts[3]}"
        else:
            short_cols[col] = col
    pivot_top = pivot_top.rename(columns=short_cols)
    
    fig_width = style['figure_width']
    fig_height = style['figure_height']
    fig, ax = plt.subplots(figsize=(fig_width, fig_height))
    
    vmax = max(abs(pivot_top.values.min()), abs(pivot_top.values.max()))
    im = ax.imshow(pivot_top.values, aspect='auto', cmap=style['palette'],
                   vmin=-vmax, vmax=vmax, origin='lower')
    
    ax.set_xlabel('Residue', fontsize=style['fontsize']['x_label'])
    ax.set_ylabel('Frame', fontsize=style['fontsize']['y_label'])
    ax.set_title('Per-Residue Energy Contribution (kcal/mol)', fontsize=style['fontsize']['title'])
    
    ax.yaxis.tick_right()
    ax.yaxis.set_label_position('right')
    
    cbar = plt.colorbar(im, ax=ax, shrink=0.8, label='Energy (kcal/mol)')
    cbar.ax.tick_params(labelsize=style['fontsize']['y_ticks'])
    
    n_xticks = min(10, len(pivot_top.columns))
    if n_xticks > 0:
        xtick_positions = np.linspace(0, len(pivot_top.columns)-1, n_xticks).astype(int)
        ax.set_xticks(xtick_positions)
        ax.set_xticklabels([list(short_cols.values())[i] for i in xtick_positions], 
                          fontsize=style['fontsize']['x_ticks'], rotation=90, ha='center')
    
    n_yticks = min(10, len(pivot_top.index))
    if n_yticks > 0:
        ytick_positions = np.linspace(0, len(pivot_top.index)-1, n_yticks).astype(int)
        ax.set_yticks(ytick_positions)
        ax.set_yticklabels([str(pivot_top.index[i]) for i in ytick_positions], fontsize=style['fontsize']['y_ticks'])
    
    plt.tight_layout()
    plt.savefig(output_path, dpi=dpi, bbox_inches='tight')
    plt.close()


def plot_mmpbsa_results(mmpbsa_dir='.', output_dir='output/mmpbsa', dpi=300):
    mmpbsa_path = Path(mmpbsa_dir)
    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)
    
    results_csv = mmpbsa_path / 'FINAL_RESULTS_MMPBSA.csv'
    decomp_csv = mmpbsa_path / 'FINAL_DECOMP_MMPBSA.csv'
    
    if results_csv.exists():
        print(f"Parsing {results_csv}...")
        results = parse_results_csv(str(results_csv))
        
        if 'delta' in results:
            plot_energy_timecourse(results['delta'], 
                                   str(output_path / 'energy_timecourse.png'), dpi)
            plot_energy_summary(results['delta'],
                                str(output_path / 'energy_summary.png'), dpi)
    
    if decomp_csv.exists():
        print(f"Parsing {decomp_csv}...")
        decomp = parse_decomposition_csv(str(decomp_csv))
        
        if 'delta' in decomp:
            plot_decomposition_bar(decomp['delta'],
                                   str(output_path / 'decomposition_bar.png'), dpi)
            plot_decomposition_heatmap(decomp['delta'],
                                       str(output_path / 'decomposition_heatmap.png'), dpi)
    
    print(f"MMPBSA plots saved to {output_path}")


def detect_ligands(workspace_path, trials_file='trials.txt'):
    """Detect ligand directories containing trials.txt file."""
    workspace = Path(workspace_path)
    ligands = []
    
    for item in workspace.iterdir():
        if item.is_dir():
            trials_file_path = item / trials_file
            if trials_file_path.exists():
                ligands.append(item.name)
    
    return sorted(ligands)


def parse_trials_file(ligand_path, trials_file='trials.txt'):
    """Parse trials.txt to get list of trial directories."""
    trials_file_path = Path(ligand_path) / trials_file
    trials = []
    
    with open(trials_file_path, 'r') as f:
        for line in f:
            line = line.strip()
            if line and not line.startswith('#'):
                trial_path = Path(ligand_path) / line
                if trial_path.exists():
                    trials.append(line)
    
    return trials


def parse_results_dat(dat_path):
    """Parse FINAL_RESULTS_MMPBSA.dat to extract summary statistics with SD/SEM."""
    with open(dat_path, 'r') as f:
        content = f.read()
    
    results = {'complex': {}, 'receptor': {}, 'ligand': {}, 'delta': {}}
    current_section = None
    
    for line in content.split('\n'):
        if 'Complex:' in line and 'Energy Component' not in line:
            current_section = 'complex'
        elif 'Receptor:' in line and 'Energy Component' not in line:
            current_section = 'receptor'
        elif 'Ligand:' in line and 'Energy Component' not in line:
            current_section = 'ligand'
        elif 'Delta' in line and ':' in line:
            current_section = 'delta'
        elif current_section and line.strip() and not line.startswith('-') and not line.startswith('|'):
            parts = line.split()
            if len(parts) >= 6:
                try:
                    term = parts[0].replace('Δ', '')
                    avg = float(parts[1])
                    sd_prop = float(parts[2])
                    sd = float(parts[3])
                    sem_prop = float(parts[4])
                    sem = float(parts[5])
                    results[current_section][term] = {
                        'average': avg,
                        'sd_prop': sd_prop,
                        'sd': sd,
                        'sem_prop': sem_prop,
                        'sem': sem
                    }
                except (ValueError, IndexError):
                    continue
    
    return results


def collect_trial_data(trial_path, receptor_backbone_sel, receptor_sel, ligand_sel, 
                       contact_cutoff, hbond_dist, hbond_angle):
    """Collect all data from a single trial."""
    trial_path = Path(trial_path)
    
    tpr_file = trial_path / 'com.tpr'
    xtc_file = trial_path / 'com_traj.xtc'
    mmpbsa_results = trial_path / 'FINAL_RESULTS_MMPBSA.csv'
    mmpbsa_results_dat = trial_path / 'FINAL_RESULTS_MMPBSA.dat'
    mmpbsa_decomp = trial_path / 'FINAL_DECOMP_MMPBSA.csv'
    
    data = {
        'rmsd_receptor': None,
        'rmsd_ligand': None,
        'contacts': None,
        'hbonds': None,
        'mmpbsa': None,
        'mmpbsa_stats': None,
        'decomp': None
    }
    
    if tpr_file.exists() and xtc_file.exists():
        try:
            universe, ref_atoms, receptor, ligand = load_and_align_trajectory(
                str(tpr_file), str(xtc_file),
                align_selection=receptor_backbone_sel,
                receptor_sel=receptor_sel,
                ligand_sel=ligand_sel
            )
            
            df_rec = calculate_receptor_rmsd(universe, receptor_backbone_sel)
            df_lig = calculate_ligand_rmsd(universe, ligand_sel + ' and not name H*')
            data['rmsd_receptor'] = df_rec
            data['rmsd_ligand'] = df_lig
            
            contact_ts, times = calculate_per_residue_contacts(
                universe, receptor_sel, ligand_sel, contact_cutoff
            )
            df_contacts = pd.DataFrame(contact_ts, index=times)
            data['contacts'] = df_contacts
            
            df_time, df_residue = calculate_hydrogen_bonds(
                universe, receptor_sel, ligand_sel, hbond_dist, hbond_angle
            )
            data['hbonds'] = {'time': df_time, 'residue': df_residue}
            
        except Exception as e:
            print(f"Error processing trial {trial_path}: {e}")
    
    if mmpbsa_results.exists():
        try:
            results = parse_results_csv(str(mmpbsa_results))
            data['mmpbsa'] = results
        except Exception as e:
            print(f"Error parsing MMPBSA results {mmpbsa_results}: {e}")
    
    if mmpbsa_results_dat.exists():
        try:
            stats = parse_results_dat(str(mmpbsa_results_dat))
            data['mmpbsa_stats'] = stats
        except Exception as e:
            print(f"Error parsing MMPBSA dat {mmpbsa_results_dat}: {e}")
    
    if mmpbsa_decomp.exists():
        try:
            decomp = parse_decomposition_csv(str(mmpbsa_decomp))
            data['decomp'] = decomp
        except Exception as e:
            print(f"Error parsing decomposition {mmpbsa_decomp}: {e}")
    
    return data


def plot_multi_trial_rmsd(all_trial_data, trial_names, output_path, dpi=300):
    """Plot RMSD for multiple trials with average."""
    style = GMX_MMPBSA_STYLE['Line Plot']
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(style['figure_width'] * 2, style['figure_height']))
    
    palette = sns.color_palette('husl', n_colors=len(trial_names))
    
    rec_data = []
    lig_data = []
    
    for i, (name, data) in enumerate(zip(trial_names, all_trial_data)):
        if data['rmsd_receptor'] is not None:
            df = data['rmsd_receptor'].copy()
            df['trial'] = name
            rec_data.append(df)
            sns.lineplot(data=data['rmsd_receptor'], x='time_ps', y='rmsd',
                        color=palette[i], linewidth=style['line_width'] * 0.5,
                        alpha=0.6, label=name, ax=ax1)
        
        if data['rmsd_ligand'] is not None:
            df = data['rmsd_ligand'].copy()
            df['trial'] = name
            lig_data.append(df)
            sns.lineplot(data=data['rmsd_ligand'], x='time_ps', y='rmsd',
                        color=palette[i], linewidth=style['line_width'] * 0.5,
                        alpha=0.6, label=name, ax=ax2)
    
    if rec_data:
        df_rec_all = pd.concat(rec_data, ignore_index=True)
        rec_avg = df_rec_all.groupby('time_ps')['rmsd'].mean().reset_index()
        rec_std = df_rec_all.groupby('time_ps')['rmsd'].std().reset_index()
        ax1.plot(rec_avg['time_ps'], rec_avg['rmsd'], 'k-', linewidth=style['rolling_avg_width'],
                label='Average')
        ax1.fill_between(rec_avg['time_ps'], 
                        rec_avg['rmsd'] - rec_std['rmsd'],
                        rec_avg['rmsd'] + rec_std['rmsd'],
                        alpha=0.2, color='black')
    
    if lig_data:
        df_lig_all = pd.concat(lig_data, ignore_index=True)
        lig_avg = df_lig_all.groupby('time_ps')['rmsd'].mean().reset_index()
        lig_std = df_lig_all.groupby('time_ps')['rmsd'].std().reset_index()
        ax2.plot(lig_avg['time_ps'], lig_avg['rmsd'], 'k-', linewidth=style['rolling_avg_width'],
                label='Average')
        ax2.fill_between(lig_avg['time_ps'],
                        lig_avg['rmsd'] - lig_std['rmsd'],
                        lig_avg['rmsd'] + lig_std['rmsd'],
                        alpha=0.2, color='black')
    
    ax1.set_xlabel('Time (ps)', fontsize=style['fontsize']['x_label'])
    ax1.set_ylabel('RMSD (A)', fontsize=style['fontsize']['y_label'])
    ax1.set_title('Receptor Backbone RMSD', fontsize=style['fontsize']['title'])
    ax1.legend(loc='upper left', fontsize=style['fontsize']['legend'], 
               bbox_to_anchor=(1.02, 1), borderaxespad=0)
    
    ax2.set_xlabel('Time (ps)', fontsize=style['fontsize']['x_label'])
    ax2.set_ylabel('RMSD (A)', fontsize=style['fontsize']['y_label'])
    ax2.set_title('Ligand Heavy Atoms RMSD', fontsize=style['fontsize']['title'])
    ax2.legend(loc='upper left', fontsize=style['fontsize']['legend'],
               bbox_to_anchor=(1.02, 1), borderaxespad=0)
    
    plt.tight_layout()
    plt.savefig(output_path, dpi=dpi, bbox_inches='tight')
    plt.close()


def plot_multi_trial_mmpbsa_summary(all_trial_data, trial_names, output_path, dpi=300):
    """Plot MMPBSA energy summary for multiple trials with average and SD.
    
    Uses SD from .dat files and proper error propagation:
    - Within-trial variance: SD_prop² from each trial's .dat file
    - Between-trial variance: variance of trial means
    - Combined SD = sqrt(mean(within_trial_var) + between_trial_var)
    """
    style = GMX_MMPBSA_STYLE['Bar Plot']
    
    terms = ['VDWAALS', 'EEL', 'EGB', 'ESURF', 'GGAS', 'GSOLV', 'TOTAL']
    summary_data = []
    
    for term in terms:
        trial_means = []
        trial_sds = []
        trial_sems = []
        
        for name, data in zip(trial_names, all_trial_data):
            if data['mmpbsa_stats'] and 'delta' in data['mmpbsa_stats']:
                delta_stats = data['mmpbsa_stats']['delta']
                if term in delta_stats:
                    trial_means.append(delta_stats[term]['average'])
                    trial_sds.append(delta_stats[term]['sd_prop'])
                    trial_sems.append(delta_stats[term]['sem_prop'])
        
        if not trial_means:
            continue
        
        trial_means = np.array(trial_means)
        trial_sds = np.array(trial_sds)
        trial_sems = np.array(trial_sems)
        
        mean_of_means = np.mean(trial_means)
        
        within_trial_var = np.mean(trial_sds ** 2)
        between_trial_var = np.var(trial_means, ddof=1) if len(trial_means) > 1 else 0
        combined_var = within_trial_var + between_trial_var
        combined_sd = np.sqrt(combined_var)
        
        n_trials = len(trial_means)
        sem = combined_sd / np.sqrt(n_trials) if n_trials > 0 else 0
        
        summary_data.append({
            'term': term,
            'mean': mean_of_means,
            'sd': combined_sd,
            'sem': sem,
            'within_trial_sd': np.sqrt(within_trial_var),
            'between_trial_sd': np.sqrt(between_trial_var),
            'trial_means': trial_means,
            'trial_sds': trial_sds,
            'trial_sems': trial_sems,
            'n_trials': n_trials
        })
    
    if not summary_data:
        print("No MMPBSA data to plot")
        return None
    
    fig, ax = plt.subplots(figsize=(style['figure_width'], style['figure_height']))
    
    x = np.arange(len(summary_data))
    means = [d['mean'] for d in summary_data]
    sds = [d['sd'] for d in summary_data]
    labels = [d['term'] for d in summary_data]
    
    palette = sns.color_palette('husl', n_colors=len(labels))
    
    bars = ax.bar(x, means, yerr=sds, capsize=style['error_cap_size'],
                  color=palette, alpha=0.85,
                  error_kw=dict(elinewidth=style['error_line_width'],
                               ecolor=style['error_line_color'],
                               capsize=style['error_cap_size']))
    
    ax.set_xticks(x)
    ax.set_xticklabels(labels, rotation=style['x_rotation'], ha='right', fontsize=style['fontsize']['x_ticks'])
    ax.set_ylabel('Energy (kcal/mol)', fontsize=style['fontsize']['y_label'])
    ax.set_title('Binding Energy Summary (Average ± SD across trials)', fontsize=style['fontsize']['title'])
    ax.axhline(y=0, color='black', lw=0.5)
    
    plt.tight_layout()
    plt.savefig(output_path, dpi=dpi, bbox_inches='tight')
    plt.close()
    
    return summary_data


def save_mmpbsa_averaged_energies(summary_data, trial_names, all_trial_data, output_path):
    """Save averaged MMPBSA energies to text file with proper error propagation."""
    output_path = Path(output_path)
    
    with open(output_path, 'w') as f:
        f.write("=" * 110 + "\n")
        f.write("MMPBSA AVERAGED ENERGIES ACROSS TRIALS\n")
        f.write("Error propagation: Combined SD = sqrt(mean(within_trial_SD²) + var(trial_means))\n")
        f.write("=" * 110 + "\n\n")
        
        f.write(f"Number of trials: {len(trial_names)}\n")
        f.write(f"Trial names: {', '.join(trial_names)}\n\n")
        
        f.write("-" * 110 + "\n")
        f.write(f"{'Energy Term':<15} {'Average':>12} {'Combined SD':>12} {'SEM':>12} {'Within-trial SD':>15} {'Between-trial SD':>15}\n")
        f.write("-" * 110 + "\n")
        
        for d in summary_data:
            f.write(f"{d['term']:<15} {d['mean']:>12.2f} {d['sd']:>12.2f} {d['sem']:>12.2f} "
                   f"{d['within_trial_sd']:>15.2f} {d['between_trial_sd']:>15.2f}\n")
        
        f.write("-" * 110 + "\n")
        f.write("Combined SD = sqrt(mean(within_trial_SD_prop²) + var(trial_means))\n")
        f.write("Within-trial SD_prop from FINAL_RESULTS_MMPBSA.dat\n")
        f.write("SEM = Combined SD / sqrt(N_trials)\n\n")
        
        f.write("Per-trial statistics:\n")
        f.write("-" * 110 + "\n")
        
        header = f"{'Trial':<15}"
        for d in summary_data:
            header += f" {d['term']:>12}"
        f.write(header + "\n")
        f.write("-" * 110 + "\n")
        
        for i, name in enumerate(trial_names):
            row = f"{name:<15}"
            for d in summary_data:
                if i < len(d['trial_means']):
                    row += f" {d['trial_means'][i]:>12.2f}"
                else:
                    row += f" {'N/A':>12}"
            f.write(row + "\n")
        
        f.write("\n" + "-" * 110 + "\n")
        f.write("Per-trial SD(Prop.) from FINAL_RESULTS_MMPBSA.dat:\n")
        f.write("-" * 110 + "\n")
        
        header = f"{'Trial':<15}"
        for d in summary_data:
            header += f" {d['term']:>12}"
        f.write(header + " (SD)\n")
        f.write("-" * 110 + "\n")
        
        for i, name in enumerate(trial_names):
            row = f"{name:<15}"
            for d in summary_data:
                if i < len(d['trial_sds']):
                    row += f" {d['trial_sds'][i]:>12.2f}"
                else:
                    row += f" {'N/A':>12}"
            f.write(row + "\n")
        
        f.write("=" * 110 + "\n")


def plot_multi_ligand_comparison(ligand_data, ligand_names, output_path, dpi=300):
    """Plot comparison across multiple ligands."""
    style = GMX_MMPBSA_STYLE['Bar Plot']
    
    all_summary = []
    
    for ligand_name, trials_data in ligand_data.items():
        trial_summaries = []
        for trial_data in trials_data:
            if trial_data['mmpbsa'] and 'delta' in trial_data['mmpbsa']:
                df = trial_data['mmpbsa']['delta']
                if 'TOTAL' in df.columns:
                    trial_summaries.append(df['TOTAL'].mean())
        
        if trial_summaries:
            all_summary.append({
                'ligand': ligand_name,
                'mean': np.mean(trial_summaries),
                'std': np.std(trial_summaries),
                'n_trials': len(trial_summaries)
            })
    
    if not all_summary:
        print("No data for multi-ligand comparison")
        return
    
    df_summary = pd.DataFrame(all_summary)
    df_summary = df_summary.sort_values('mean')
    
    fig, ax = plt.subplots(figsize=(max(style['figure_width'], len(df_summary) * 0.8), 
                                         style['figure_height']))
    
    y = np.arange(len(df_summary))
    colors = [PALETTE['favorable'] if v < 0 else PALETTE['unfavorable'] 
              for v in df_summary['mean'].values]
    
    ax.barh(y, df_summary['mean'], xerr=df_summary['std'],
            color=colors, alpha=0.85,
            error_kw=dict(elinewidth=style['error_line_width'],
                         ecolor=style['error_line_color'],
                         capsize=style['error_cap_size']))
    
    ax.set_yticks(y)
    ax.set_yticklabels(df_summary['ligand'], fontsize=style['fontsize']['y_ticks'])
    ax.set_xlabel('ΔG (kcal/mol)', fontsize=style['fontsize']['x_label'])
    ax.set_title('Binding Free Energy Comparison', fontsize=style['fontsize']['title'])
    ax.axvline(x=0, color='black', lw=0.5)
    
    for i, (mean, n) in enumerate(zip(df_summary['mean'], df_summary['n_trials'])):
        ax.text(mean + 0.5, i, f'n={n}', va='center', fontsize=8)
    
    plt.tight_layout()
    plt.savefig(output_path, dpi=dpi, bbox_inches='tight')
    plt.close()


def run_single_trial_analysis(trial_path, output_dir, args):
    """Run complete analysis for a single trial (per-trial output)."""
    trial_path = Path(trial_path)
    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)
    
    tpr_file = trial_path / 'com.tpr'
    xtc_file = trial_path / 'com_traj.xtc'
    
    universe = None
    if tpr_file.exists() and xtc_file.exists() and args.ana in ['all', 'load', 'rmsd', 'contacts', 'hbonds']:
        universe, ref_atoms, receptor, ligand = load_and_align_trajectory(
            str(tpr_file), str(xtc_file),
            align_selection=args.receptor_backbone_sel,
            receptor_sel=args.receptor_sel,
            ligand_sel=args.ligand_sel
        )
        print(f"  Loaded trajectory: {universe.trajectory.n_frames} frames")
    
    if args.ana in ['all', 'rmsd'] and universe:
        calculate_rmsd(
            universe=universe,
            receptor_sel=args.receptor_backbone_sel,
            ligand_sel=args.ligand_sel + ' and not name H*',
            output_dir=str(output_path / 'rmsd')
        )
        print(f"  RMSD analysis complete")
    
    if args.ana in ['all', 'contacts'] and universe:
        calculate_contacts(
            universe=universe,
            receptor_sel=args.receptor_sel,
            ligand_sel=args.ligand_sel,
            cutoff=args.contact_cutoff,
            output_dir=str(output_path / 'contacts')
        )
        print(f"  Contact analysis complete")
    
    if args.ana in ['all', 'hbonds'] and universe:
        calculate_hydrogen_bonds_wrapper(
            universe=universe,
            receptor_sel=args.receptor_sel,
            ligand_sel=args.ligand_sel,
            distance_cutoff=args.hbond_distance_cutoff,
            angle_cutoff=args.hbond_angle_cutoff,
            output_dir=str(output_path / 'hydrogen_bonds')
        )
        print(f"  Hydrogen bond analysis complete")
    
    if args.ana in ['all', 'mmpbsa']:
        plot_mmpbsa_results(
            mmpbsa_dir=str(trial_path),
            output_dir=str(output_path / 'mmpbsa')
        )
        print(f"  MMPBSA plots complete")


def plot_multi_trial_contacts(all_trial_data, trial_names, output_path, dpi=300):
    """Plot contact analysis across multiple trials."""
    style = GMX_MMPBSA_STYLE['Bar Plot']
    
    all_contacts = []
    for name, data in zip(trial_names, all_trial_data):
        if data['contacts'] is not None:
            df = data['contacts'].copy()
            df['trial'] = name
            all_contacts.append(df)
    
    if not all_contacts:
        print("No contact data to plot")
        return
    
    df_all = pd.concat(all_contacts, ignore_index=False)
    df_avg = df_all.drop(columns='trial', errors='ignore').groupby(level=0).mean()
    
    avg_by_residue = df_avg.mean().sort_values(ascending=False).head(20)
    
    fig, ax = plt.subplots(figsize=(style['figure_width'], style['figure_height']))
    
    y = np.arange(len(avg_by_residue))
    palette = sns.color_palette('husl', n_colors=len(avg_by_residue))
    
    sns.barplot(x=avg_by_residue.values, y=y, palette=palette, ax=ax, orient='h')
    
    ax.set_yticks(y)
    ax.set_yticklabels(avg_by_residue.index, fontsize=style['fontsize']['y_ticks'])
    ax.set_xlabel('Average Contact Number', fontsize=style['fontsize']['x_label'])
    ax.set_title('Top 20 Residues by Average Contacts (Multi-Trial)', fontsize=style['fontsize']['title'])
    if style['y_inverted']:
        ax.invert_yaxis()
    
    plt.tight_layout()
    plt.savefig(output_path, dpi=dpi, bbox_inches='tight')
    plt.close()
    print(f"  Contact plot saved")


def plot_multi_trial_hbonds(all_trial_data, trial_names, output_path, dpi=300):
    """Plot hydrogen bond analysis across multiple trials."""
    style = GMX_MMPBSA_STYLE['Bar Plot']
    
    all_hbonds = []
    for name, data in zip(trial_names, all_trial_data):
        if data['hbonds'] is not None and not data['hbonds']['residue'].empty:
            df = data['hbonds']['residue'].copy()
            df['trial'] = name
            all_hbonds.append(df)
    
    if not all_hbonds:
        print("No hydrogen bond data to plot")
        return
    
    df_all = pd.concat(all_hbonds, ignore_index=True)
    
    avg_by_residue = df_all.groupby('residue')['hbond_count'].sum().sort_values(ascending=False).head(15)
    
    fig, ax = plt.subplots(figsize=(style['figure_width'], style['figure_height']))
    
    y = np.arange(len(avg_by_residue))
    palette = sns.color_palette('husl', n_colors=len(avg_by_residue))
    
    sns.barplot(x=avg_by_residue.values, y=y, palette=palette, ax=ax, orient='h')
    
    ax.set_yticks(y)
    ax.set_yticklabels(avg_by_residue.index, fontsize=style['fontsize']['y_ticks'])
    ax.set_xlabel('Total Hydrogen Bond Count', fontsize=style['fontsize']['x_label'])
    ax.set_title('Top 15 Residues by H-bond Frequency (Multi-Trial)', fontsize=style['fontsize']['title'])
    if style['y_inverted']:
        ax.invert_yaxis()
    
    plt.tight_layout()
    plt.savefig(output_path, dpi=dpi, bbox_inches='tight')
    plt.close()
    print(f"  Hydrogen bond plot saved")


def run_multi_trial_analysis(workspace_path, trials_file, output_dir, args):
    """Run analysis for multiple trials of a single ligand."""
    workspace = Path(workspace_path)
    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)
    
    trials = parse_trials_file(workspace, trials_file)
    if not trials:
        print(f"No trials found in {workspace}")
        return
    
    print(f"Found {len(trials)} trials: {', '.join(trials)}")
    
    all_trial_data = []
    for trial in trials:
        trial_path = workspace / trial
        print(f"\nProcessing trial: {trial}")
        
        trial_output = output_path / 'per_trial' / trial
        run_single_trial_analysis(trial_path, str(trial_output), args)
        
        data = collect_trial_data(
            trial_path,
            args.receptor_backbone_sel,
            args.receptor_sel,
            args.ligand_sel,
            args.contact_cutoff,
            args.hbond_distance_cutoff,
            args.hbond_angle_cutoff
        )
        all_trial_data.append(data)
    
    print("\nGenerating multi-trial combined plots...")
    
    if args.ana in ['all', 'rmsd']:
        rmsd_path = output_path / 'rmsd'
        rmsd_path.mkdir(exist_ok=True)
        plot_multi_trial_rmsd(all_trial_data, trials, str(rmsd_path / 'rmsd_multi_trial.png'),
                             dpi=args.plot_dpi)
        print(f"  RMSD combined plot saved")
    
    if args.ana in ['all', 'contacts']:
        contacts_path = output_path / 'contacts'
        contacts_path.mkdir(exist_ok=True)
        plot_multi_trial_contacts(all_trial_data, trials, 
                                  str(contacts_path / 'contacts_multi_trial.png'),
                                  dpi=args.plot_dpi)
    
    if args.ana in ['all', 'hbonds']:
        hbonds_path = output_path / 'hydrogen_bonds'
        hbonds_path.mkdir(exist_ok=True)
        plot_multi_trial_hbonds(all_trial_data, trials,
                                str(hbonds_path / 'hbonds_multi_trial.png'),
                                dpi=args.plot_dpi)
    
    if args.ana in ['all', 'mmpbsa']:
        mmpbsa_path = output_path / 'mmpbsa'
        mmpbsa_path.mkdir(exist_ok=True)
        
        summary_data = plot_multi_trial_mmpbsa_summary(
            all_trial_data, trials,
            str(mmpbsa_path / 'energy_summary_multi_trial.png'),
            dpi=args.plot_dpi
        )
        
        if summary_data:
            save_mmpbsa_averaged_energies(
                summary_data, trials, all_trial_data,
                str(mmpbsa_path / 'averaged_energies.txt')
            )
            print(f"  MMPBSA combined plots and averaged energies saved")
    
    return all_trial_data


def run_multi_ligand_analysis(workspace_path, trials_file, output_dir, args):
    """Run analysis for multiple ligands."""
    workspace = Path(workspace_path)
    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)
    
    ligands = detect_ligands(workspace, trials_file)
    if not ligands:
        print(f"No ligands found in {workspace}")
        return
    
    print(f"Found {len(ligands)} ligands: {', '.join(ligands)}")
    
    ligand_data = {}
    
    for ligand in ligands:
        ligand_path = workspace / ligand
        print(f"\n{'='*60}")
        print(f"Processing ligand: {ligand}")
        print(f"{'='*60}")
        
        ligand_output = output_path / 'per_ligand' / ligand
        trial_data = run_multi_trial_analysis(
            ligand_path, trials_file, str(ligand_output), args
        )
        
        if trial_data:
            ligand_data[ligand] = trial_data
    
    if len(ligands) > 1:
        print("\nGenerating multi-ligand comparison plots...")
        plot_multi_ligand_comparison(
            ligand_data, ligands,
            str(output_path / 'ligand_comparison.png'),
            dpi=args.plot_dpi
        )
        print(f"Ligand comparison saved to {output_path}")
    
    print(f"\nAll analyses complete. Output saved to: {output_path}")


def main():
    args = parse_arguments()
    
    output_dir = setup_output_directory(args.output_dir)
    
    if args.workspace:
        workspace_path = Path(args.workspace)
        trials_file_path = workspace_path / args.trials_file
        
        if trials_file_path.exists():
            print(f"Detected single ligand with multiple trials in {args.workspace}")
            run_multi_trial_analysis(args.workspace, args.trials_file, args.output_dir, args)
        else:
            print(f"Detected multi-ligand workspace: {args.workspace}")
            run_multi_ligand_analysis(args.workspace, args.trials_file, args.output_dir, args)
        return
    
    universe = None
    need_trajectory = args.ana in ['all', 'load', 'rmsd', 'contacts', 'hbonds']
    
    if need_trajectory:
        universe, reference_atoms, receptor, ligand = load_and_align_trajectory(
            topology=args.tpr,
            trajectory=args.xtc,
            index_file=args.ndx,
            align_selection=args.receptor_backbone_sel,
            receptor_sel=args.receptor_sel,
            ligand_sel=args.ligand_sel
        )
        print(f"Loaded and aligned trajectory: {universe.trajectory.n_frames} frames")
    
    if args.ana in ['all', 'load']:
        print("Load and align complete")
    
    if args.ana in ['all', 'rmsd']:
        rmsd_results = calculate_rmsd(
            universe=universe,
            receptor_sel=args.receptor_backbone_sel,
            ligand_sel=args.ligand_sel + ' and not name H*',
            output_dir=str(output_dir / 'rmsd')
        )
        print(f"RMSD analysis complete: {len(rmsd_results[0])} frames")
    
    if args.ana in ['all', 'contacts']:
        contact_results = calculate_contacts(
            universe=universe,
            receptor_sel=args.receptor_sel,
            ligand_sel=args.ligand_sel,
            cutoff=args.contact_cutoff,
            output_dir=str(output_dir / 'contacts')
        )
        print(f"Contact analysis complete")
    
    if args.ana in ['all', 'hbonds']:
        hbond_results = calculate_hydrogen_bonds_wrapper(
            universe=universe,
            receptor_sel=args.receptor_sel,
            ligand_sel=args.ligand_sel,
            distance_cutoff=args.hbond_distance_cutoff,
            angle_cutoff=args.hbond_angle_cutoff,
            output_dir=str(output_dir / 'hydrogen_bonds')
        )
        print(f"Hydrogen bond analysis complete")
    
    if args.ana in ['all', 'mmpbsa']:
        plot_mmpbsa_results(
            mmpbsa_dir=args.mmpbsa_dir,
            output_dir=str(output_dir / 'mmpbsa')
        )
        print(f"MMPBSA plots generated")
    
    print(f"\nAll analyses complete. Output saved to: {output_dir}")


if __name__ == '__main__':
    main()
