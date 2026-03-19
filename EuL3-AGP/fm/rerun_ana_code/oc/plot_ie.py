import argparse
import os
import glob
import re
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm
import warnings

# --- DEFAULT PARAMETERS ---
DEFAULT_KT = 2.494339
DEFAULT_D1_MIN = 0.0
DEFAULT_D1_MAX = 4.0
DEFAULT_BINS = 40
DEFAULT_LABEL_NEG = -40.0
DEFAULT_LABEL_POS = 10.0
# --------------------------

def get_resid_resname(filename):
    """Extract residue ID and 3-letter code from filename like 'rerun_Protein_ALA_128.xvg'"""
    match = re.search(r'_([A-Za-z0-9]+)_(\d+)\.xvg$', os.path.basename(filename))
    if match:
        return int(match.group(2)), match.group(1)
    return -1, "UNK"

def find_d1_ranges(mean_IE_row, bins, threshold, is_negative=True):
    ranges = []
    in_range = False
    start_d1 = None
    
    for i, val in enumerate(mean_IE_row):
        if np.isnan(val):
            condition = False
        else:
            condition = (val < threshold) if is_negative else (val > threshold)
            
        if condition:
            if not in_range:
                start_d1 = bins[i]
                in_range = True
        else:
            if in_range:
                end_d1 = bins[i]
                ranges.append(f"{start_d1:.1f}-{end_d1:.1f}")
                in_range = False
                
    if in_range:
        end_d1 = bins[-1]
        ranges.append(f"{start_d1:.1f}-{end_d1:.1f}")
        
    return ", ".join(ranges) if ranges else "N/A"

def process_replica(replica_dir, kt, N_bins, d1_min, d1_max):
    print(f"Processing replica: {replica_dir}")
    colvar_path = os.path.join(replica_dir, 'trj.COLVAR.o')
    if not os.path.exists(colvar_path):
        print(f"Warning: {colvar_path} not found. Skipping {replica_dir}.")
        return None, None, None, None

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
        return None, None, None, None

    xvg_files.sort(key=lambda x: get_resid_resname(x)[0])
    resids = []
    resnames = []
    for f in xvg_files:
        rid, rname = get_resid_resname(f)
        resids.append(rid)
        resnames.append(rname)
    resids = np.array(resids)
    resnames = np.array(resnames)

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
    np.savez(os.path.join(replica_dir, 'binned_IE.npz'), resids=resids, resnames=resnames, mean_IE=mean_IE, std_IE=std_IE, bins=bins)
    
    return resids, resnames, mean_IE, std_IE

def plot_heatmap(mean_IE, resids, resnames, bins, d1_min, d1_max, out_path, title_prefix="", 
                 neg_cutoff=DEFAULT_LABEL_NEG, pos_cutoff=DEFAULT_LABEL_POS, plot_labels=True, 
                 cbar_min=None, cbar_max=None):
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
        data_vmin = np.nanmin(mapped_mean_IE)
        data_vmax = np.nanmax(mapped_mean_IE)
        
    if np.isnan(data_vmin) or np.isnan(data_vmax):
        print(f"Warning: No valid data to plot for {out_path}.")
        return

    # User defined limits, or auto-detect min/max
    vmin = data_vmin if cbar_min is None else cbar_min
    vmax = data_vmax if cbar_max is None else cbar_max
    
    # Provide tiny fallback so norm scale always exists correctly
    if vmax <= 0:
        vmax = 0.1
    if vmin >= 0:
        vmin = -0.1

    norm = TwoSlopeNorm(vcenter=0, vmin=vmin, vmax=vmax)
    
    im = plt.imshow(mapped_mean_IE, aspect='auto', origin='lower', extent=extent,
                    cmap='bwr', norm=norm, interpolation='nearest')
    
    # Custom visually proportional colorbar ticks
    neg_ticks = np.linspace(vmin, 0, 5)
    pos_ticks = np.linspace(0, vmax, 5)
    cbar_ticks = np.concatenate([neg_ticks[:-1], pos_ticks])
    
    cbar_ticks_rounded = []
    for t in cbar_ticks:
        if t == 0:
            cbar_ticks_rounded.append(0)
        elif abs(t) >= 10:
            cbar_ticks_rounded.append(int(round(t)))
        else:
            cbar_ticks_rounded.append(round(t, 1))
    cbar_ticks_rounded = sorted(list(set(cbar_ticks_rounded)))
    
    cbar = plt.colorbar(im, ticks=cbar_ticks_rounded)
    cbar.set_label("Average IE (kJ/mol)", fontsize=18)
    cbar.ax.tick_params(labelsize=16)
    
    plt.xlabel("D1 (nm)", fontsize=18)
    plt.ylabel("Residue ID", fontsize=18)
    
    plt.gca().set_xlim([d1_min, d1_max])
    plt.gca().set_ylim([min_res - 0.5, max_res + 0.5])
    
    x_major_ticks = np.arange(d1_min, d1_max + 0.001, 0.5)
    x_minor_ticks = np.arange(d1_min, d1_max + 0.001, 0.1)
    plt.gca().set_xticks(x_major_ticks)
    plt.gca().set_xticks(x_minor_ticks, minor=True)
    
    y_major_start = (min_res // 20) * 20
    y_major_ticks = np.arange(y_major_start, max_res + 20, 20)
    # Reduce minor ticks: every 5 residues instead of every 1
    y_minor_ticks = np.arange(y_major_start, max_res + 1, 5)
    plt.gca().set_yticks(y_major_ticks)
    plt.gca().set_yticks(y_minor_ticks, minor=True)
    
    # Process major contributions
    labels_to_plot = []
    csv_rows = []
    
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", category=RuntimeWarning)
        for idx, res in enumerate(resids):
            rname = resnames[idx]
            row_IE = mean_IE[idx, :]
            
            min_val = np.nanmin(row_IE)
            if not np.isnan(min_val) and min_val < neg_cutoff:
                min_idx = np.nanargmin(row_IE)
                peak_d1 = (bins[min_idx] + bins[min_idx+1]) / 2.0
                d1_range = find_d1_ranges(row_IE, bins, neg_cutoff, True)
                csv_rows.append({'Residue': res, 'Code': rname, 'Type': 'Negative', 
                                 'Peak_IE': min_val, 'D1_Range': d1_range})
                labels_to_plot.append((peak_d1, res, f"{rname}{res} ({min_val:.2f} kJ/mol)"))
                
            max_val = np.nanmax(row_IE)
            if not np.isnan(max_val) and max_val > pos_cutoff:
                max_idx = np.nanargmax(row_IE)
                peak_d1 = (bins[max_idx] + bins[max_idx+1]) / 2.0
                d1_range = find_d1_ranges(row_IE, bins, pos_cutoff, False)
                csv_rows.append({'Residue': res, 'Code': rname, 'Type': 'Positive', 
                                 'Peak_IE': max_val, 'D1_Range': d1_range})
                labels_to_plot.append((peak_d1, res, f"{rname}{res} (+{max_val:.2f} kJ/mol)"))
    
    if csv_rows:
        df_labels = pd.DataFrame(csv_rows)
        csv_path = out_path.replace('.png', '_major_contributions.csv')
        df_labels.to_csv(csv_path, index=False)
        print(f"  -> Saved major contributions to {csv_path}")
        
        if plot_labels:
            labels_to_plot.sort(key=lambda x: x[1])
            occupied_y = []
            
            for peak_d1, res, text in labels_to_plot:
                # Plot circle (no fill, black edge)
                plt.plot(peak_d1, res, 'o', markerfacecolor='none', markeredgecolor='black', 
                         markersize=5, markeredgewidth=1.5)
                
                text_y = res + 2
                for oy in occupied_y:
                    if abs(text_y - oy) < 3.0:
                        text_y = oy + 3.5
                
                if text_y > max_res:
                    text_y = res - 3
                    
                occupied_y.append(text_y)
                
                plt.annotate(text, xy=(peak_d1, res), xytext=(peak_d1 + 0.1, text_y),
                             color='black', fontsize=10, fontweight='bold',
                             bbox=dict(boxstyle="round,pad=0.2", fc="white", alpha=0.5, ec="none"),
                             arrowprops=dict(arrowstyle="->", connectionstyle="arc3", color="black", alpha=0.7))

    title = f"{title_prefix}Residue Interaction Energy"
    if title_prefix:
        title = f"{title_prefix.strip()} Residue Interaction Energy"
    plt.title(title, fontsize=20)
    plt.tick_params(axis='both', which='major', labelsize=16)
    plt.tight_layout()
    plt.savefig(out_path, dpi=300)
    print(f"Saved plot to {out_path}")
    plt.close()

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Process data from multiple replicas of funnel metadynamics")
    parser.add_argument('replicas', nargs='+', help='Directories of replicas to process (e.g., t1 t1_1)')
    parser.add_argument('--kt', type=float, default=DEFAULT_KT, help=f'kT value for reweighting (default: {DEFAULT_KT})')
    parser.add_argument('--d1_min', type=float, default=DEFAULT_D1_MIN, help=f'Minimum D1 value (default: {DEFAULT_D1_MIN})')
    parser.add_argument('--d1_max', type=float, default=DEFAULT_D1_MAX, help=f'Maximum D1 value (default: {DEFAULT_D1_MAX})')
    parser.add_argument('--bins', type=int, default=DEFAULT_BINS, help=f'Number of bins for D1 (default: {DEFAULT_BINS})')
    parser.add_argument('--label_neg', type=float, default=DEFAULT_LABEL_NEG, help=f'Negative IE cutoff to label major contribution (default: {DEFAULT_LABEL_NEG})')
    parser.add_argument('--label_pos', type=float, default=DEFAULT_LABEL_POS, help=f'Positive IE cutoff to label major contribution (default: {DEFAULT_LABEL_POS})')
    parser.add_argument('--no_plot_labels', action='store_true', help='Disable overlaying labels on the plot itself')
    parser.add_argument('--cmin', type=float, default=None, help='Minimum value for colorbar (default: auto from data)')
    parser.add_argument('--cmax', type=float, default=None, help='Maximum value for colorbar (default: auto from data)')
    args = parser.parse_args()

    all_resids = None
    all_resnames = None
    all_means = []
    
    for rep in args.replicas:
        resids, resnames, mean_IE, std_IE = process_replica(rep, args.kt, args.bins, args.d1_min, args.d1_max)
        if mean_IE is not None:
            if all_resids is None:
                all_resids = resids
                all_resnames = resnames
            else:
                if not np.array_equal(all_resids, resids):
                    print(f"Warning: Residues in {rep} do not match the first replica. Intersection will not be handled perfectly here.")
            
            all_means.append(mean_IE)
            
            out_png = os.path.join(rep, 'res-ie.png')
            bins_arr = np.linspace(args.d1_min, args.d1_max, args.bins + 1)
            plot_heatmap(mean_IE, resids, resnames, bins_arr, args.d1_min, args.d1_max, out_png, 
                         title_prefix=f"{rep} ", neg_cutoff=args.label_neg, pos_cutoff=args.label_pos,
                         plot_labels=not args.no_plot_labels, cbar_min=args.cmin, cbar_max=args.cmax)

    if len(all_means) > 1:
        print("Averaging across replicas...")
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", category=RuntimeWarning)
            merged_mean = np.nanmean(np.array(all_means), axis=0)
        
        bins_arr = np.linspace(args.d1_min, args.d1_max, args.bins + 1)
        plot_heatmap(merged_mean, all_resids, all_resnames, bins_arr, args.d1_min, args.d1_max, 'merged_res-ie.png', 
                     title_prefix="Merged ", neg_cutoff=args.label_neg, pos_cutoff=args.label_pos,
                     plot_labels=not args.no_plot_labels, cbar_min=args.cmin, cbar_max=args.cmax)
