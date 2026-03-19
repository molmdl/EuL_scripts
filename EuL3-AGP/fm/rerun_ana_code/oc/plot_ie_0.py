import argparse
import os
import glob
import re
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm
import warnings

def get_resid(filename):
    match = re.search(r'_(\d+)\.xvg$', os.path.basename(filename))
    if match:
        return int(match.group(1))
    return -1

def process_replica(replica_dir, kt, N_bins, d1_min, d1_max):
    print(f"Processing replica: {replica_dir}")
    colvar_path = os.path.join(replica_dir, 'trj.COLVAR.o')
    if not os.path.exists(colvar_path):
        print(f"Warning: {colvar_path} not found. Skipping {replica_dir}.")
        return None, None, None

    with open(colvar_path, 'r') as f:
        fields = f.readline().strip().split()[2:]
    df = pd.read_csv(colvar_path, sep=r'\s+', comment='#', names=fields)
    
    d1 = df['d1'].values
    metad_bias = df['metad.bias'].values
    
    # Calculate weight from REWEIGHT_METAD logic: w = exp(bias / kT)
    as_bias = metad_bias / kt
    w = np.exp(as_bias - np.max(as_bias))
    n_colvar = len(d1)

    xvg_files = glob.glob(os.path.join(replica_dir, 'rerun_*.xvg'))
    if not xvg_files:
        print(f"Warning: No rerun_*.xvg files found in {replica_dir}.")
        return None, None, None

    xvg_files.sort(key=get_resid)
    resids = np.array([get_resid(f) for f in xvg_files])

    # Find common minimum frames to align arrays
    n_frames_list = []
    data_list = []
    for f in xvg_files:
        data = np.loadtxt(f, comments=['#', '@'])
        data_list.append(data)
        n_frames_list.append(data.shape[0])

    min_frames = min([n_colvar] + n_frames_list)
    print(f"  -> Aligning all data to {min_frames} frames.")

    IE_all = np.zeros((len(resids), min_frames))
    for i, data in enumerate(data_list):
        # IE = Coul-SR + LJ-SR (columns 1 and 2 of gmx energy default output)
        IE_all[i, :] = data[:min_frames, 1] + data[:min_frames, 2]

    d1 = d1[:min_frames]
    w = w[:min_frames]

    # Save details per frame for reference
    df_frame = pd.DataFrame({
        'frame': np.arange(min_frames),
        'd1': d1,
        'weight': w
    })
    df_frame.to_csv(os.path.join(replica_dir, 'frame_weights.csv'), index=False)
    print(f"  -> Saved frame weights to {replica_dir}/frame_weights.csv")

    bins = np.linspace(d1_min, d1_max, N_bins + 1)
    
    mean_IE = np.full((len(resids), N_bins), np.nan)
    std_IE = np.full((len(resids), N_bins), np.nan)

    for b in range(N_bins):
        bin_start = bins[b]
        bin_end = bins[b+1]
        
        if b == N_bins - 1:
            in_bin = (d1 >= bin_start) & (d1 <= bin_end)
        else:
            in_bin = (d1 >= bin_start) & (d1 < bin_end)
            
        if np.any(in_bin):
            w_bin = w[in_bin]
            ie_bin = IE_all[:, in_bin]
            
            sum_w = np.sum(w_bin)
            if sum_w > 0:
                mean = np.sum(ie_bin * w_bin, axis=1) / sum_w
                mean_IE[:, b] = mean
                variance = np.sum(w_bin * (ie_bin - mean[:, None])**2, axis=1) / sum_w
                std_IE[:, b] = np.sqrt(variance)

    # Save bin-wise stats
    np.savez(os.path.join(replica_dir, 'binned_IE.npz'), resids=resids, mean_IE=mean_IE, std_IE=std_IE, bins=bins)
    
    return resids, mean_IE, std_IE

def plot_heatmap(mean_IE, resids, N_bins, d1_min, d1_max, out_path, title_prefix=""):
    plt.figure(figsize=(10, 8))
    
    min_res = resids.min() if len(resids) > 0 else 0
    max_res = resids.max() if len(resids) > 0 else 100
    full_res_range = np.arange(min_res, max_res + 1)
    
    mapped_mean_IE = np.full((len(full_res_range), mean_IE.shape[1]), np.nan)
    for idx, res in enumerate(resids):
        row_idx = res - min_res
        mapped_mean_IE[row_idx, :] = mean_IE[idx, :]
    
    extent = [d1_min, d1_max, min_res - 0.5, max_res + 0.5]
    
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", category=RuntimeWarning)
        vmin = np.nanmin(mapped_mean_IE)
        vmax = np.nanmax(mapped_mean_IE)
        
    if np.isnan(vmin) or np.isnan(vmax):
        print(f"Warning: No valid data to plot for {out_path}.")
        return

    abs_max = max(abs(vmin), abs(vmax))
    # Prevent TwoSlopeNorm error if abs_max is 0
    if abs_max == 0:
        abs_max = 1e-3
        
    norm = TwoSlopeNorm(vcenter=0, vmin=-abs_max, vmax=abs_max)
    
    im = plt.imshow(mapped_mean_IE, aspect='auto', origin='lower', extent=extent,
                    cmap='bwr', norm=norm, interpolation='nearest')
    
    cbar = plt.colorbar(im)
    cbar.set_label("Average IE (kJ/mol)", fontsize=14)
    cbar.ax.tick_params(labelsize=12)
    
    plt.xlabel("D1 (nm)", fontsize=14)
    plt.ylabel("Residue ID", fontsize=14)
    
    plt.gca().set_xlim([d1_min, d1_max])
    plt.gca().set_ylim([min_res - 0.5, max_res + 0.5])
    
    x_major_ticks = np.arange(d1_min, d1_max + 0.001, 0.5)
    x_minor_ticks = np.arange(d1_min, d1_max + 0.001, 0.1)
    plt.gca().set_xticks(x_major_ticks)
    plt.gca().set_xticks(x_minor_ticks, minor=True)
    
    y_major_start = (min_res // 20) * 20
    y_major_ticks = np.arange(y_major_start, max_res + 20, 20)
    y_minor_ticks = np.arange(min_res, max_res + 1, 1)
    plt.gca().set_yticks(y_major_ticks)
    plt.gca().set_yticks(y_minor_ticks, minor=True)
    
    title = f"{title_prefix}Residue Interaction Energy"
    if title_prefix:
        title = f"{title_prefix.strip()} Residue Interaction Energy"
    plt.title(title, fontsize=16)
    plt.tick_params(axis='both', which='major', labelsize=12)
    plt.tight_layout()
    plt.savefig(out_path, dpi=300)
    print(f"Saved plot to {out_path}")
    plt.close()

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Process data from multiple replicas of funnel metadynamics")
    parser.add_argument('replicas', nargs='+', help='Directories of replicas to process (e.g., t1 t1_1)')
    parser.add_argument('--kt', type=float, default=2.494339, help='kT value for reweighting (kJ/mol)')
    parser.add_argument('--d1_min', type=float, default=0.0, help='Minimum D1 value')
    parser.add_argument('--d1_max', type=float, default=4.0, help='Maximum D1 value')
    parser.add_argument('--bins', type=int, default=40, help='Number of bins for D1 (default 40 for 0.1 spacing)')
    args = parser.parse_args()

    all_resids = None
    all_means = []
    
    for rep in args.replicas:
        resids, mean_IE, std_IE = process_replica(rep, args.kt, args.bins, args.d1_min, args.d1_max)
        if mean_IE is not None:
            if all_resids is None:
                all_resids = resids
            else:
                # Ensure all replicas have the same resids
                if not np.array_equal(all_resids, resids):
                    print(f"Warning: Residues in {rep} do not match the first replica. Intersection will not be handled perfectly here.")
            
            all_means.append(mean_IE)
            
            out_png = os.path.join(rep, 'res-ie.png')
            plot_heatmap(mean_IE, resids, args.bins, args.d1_min, args.d1_max, out_png, title_prefix=f"{rep} ")

    if len(all_means) > 1:
        print("Averaging across replicas...")
        # Average ignoring nans
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", category=RuntimeWarning)
            merged_mean = np.nanmean(np.array(all_means), axis=0)
        
        # Save merged to current directory
        plot_heatmap(merged_mean, all_resids, args.bins, args.d1_min, args.d1_max, 'merged_res-ie.png', title_prefix="Merged ")
