import logging
import os
import argparse
import tracemalloc
import numpy as np
import MDAnalysis as mda
from MDAnalysis.analysis import align
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.ticker import MaxNLocator
from scipy.ndimage import gaussian_filter

# Default variables
DEFAULT_TOPOLOGY = "../apo_100ns_0_pbc_0.pdb"
#DEFAULT_TRAJECTORIES = "apo_100ns_0_pbc.xtc,apo_100ns_1_pbc.xtc,apo_100ns_2_pbc.xtc,apo_100ns_3_pbc.xtc"
DEFAULT_TRAJECTORIES = "apo_100ns_merged.xtc"
DEFAULT_RECEPTOR_TYPE = "protein"
DEFAULT_SELECTION = "backbone" #"not name H*"
DEFAULT_WEIGHT_FILE = None
DEFAULT_TEMPERATURE = 300.0  # K, for free energy
DEFAULT_SMOOTH = False  # Smoothing off by default

# Constants
R = 0.008314  # kJ/mol/K
BINS = 100  # for 2D histograms
EPS = 1e-10  # to avoid log(0)
SMOOTH_SIGMA = 1.0  # for Gaussian smoothing

# Set up logging
logging.basicConfig(filename="analysis.log", level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s")

def main():
    logging.info("Starting PCA analysis script.")
    tracemalloc.start()

    try:
        # Parse command-line arguments
        parser = argparse.ArgumentParser(description="Perform PCA on biomolecular trajectories.")
        parser.add_argument("topology", type=str, nargs="?", default=DEFAULT_TOPOLOGY, help="Topology file (e.g., v1.pdb)")
        parser.add_argument("trajectories", type=str, nargs="?", default=DEFAULT_TRAJECTORIES, help='Comma-separated trajectories in quotes (e.g., "md_aq.xtc,md_org.xtc")')
        parser.add_argument("receptor_type", type=str, nargs="?", default=DEFAULT_RECEPTOR_TYPE, choices=["protein", "nucleic acid"], help="Receptor type")
        parser.add_argument("selection", type=str, nargs="?", default=DEFAULT_SELECTION, help="Selection for PCA (e.g., 'not name H*')")
        parser.add_argument("--weight", type=str, default=DEFAULT_WEIGHT_FILE, help="Weight file (e.g., weight.xvg)")
        parser.add_argument("--temperature", type=float, default=DEFAULT_TEMPERATURE, help="Temperature in K for free energy (e.g., 300)")
        parser.add_argument("--smooth", action="store_true", default=DEFAULT_SMOOTH, help="Enable Gaussian smoothing for 2D histograms")
        args = parser.parse_args()

        traj_str = args.trajectories.strip('"')
        traj_list = [t.strip() for t in traj_str.split(",")]
        sel_str = args.selection
        weight_file = args.weight
        T = args.temperature
        smooth = args.smooth

        logging.info(f"Inputs: topology={args.topology}, trajectories={traj_list}, receptor_type={args.receptor_type}, selection={sel_str}, weight={weight_file}, T={T}, smooth={smooth}")

        # Get individual trajectory lengths
        frame_starts = [0]
        for traj in traj_list:
            temp_u = mda.Universe(args.topology, traj)
            frame_starts.append(frame_starts[-1] + len(temp_u.trajectory))
        total_frames = frame_starts[-1]
        frame_ranges = list(zip(frame_starts[:-1], frame_starts[1:]))

        # Load combined universe
        logging.info("Loading combined trajectory.")
        u = mda.Universe(args.topology, traj_list)

        # Set reference from first frame
        ref_u = mda.Universe(args.topology)
        u.trajectory[0]
        align.alignto(u, ref_u, select=sel_str)
        ref = ref_u.select_atoms(sel_str)
        ref_coords = ref.positions.copy()

        # Load or generate weights
        if weight_file:
            logging.info(f"Loading weights from {weight_file}.")
            weights = np.loadtxt(weight_file)
        else:
            logging.info("Generating uniform weights.")
            weights = np.ones(total_frames)
        if len(weights) != total_frames:
            raise ValueError("Weights length does not match total frames.")

        # Collect aligned coordinates
        logging.info("Aligning frames and collecting coordinates.")
        atoms = u.select_atoms(sel_str)
        natoms = len(atoms)
        data = np.empty((total_frames, 3 * natoms))
        for i, ts in enumerate(u.trajectory):
            align.alignto(u, ref_u, select=sel_str)
            data[i] = atoms.positions.flatten()

        # Weighted PCA
        logging.info("Computing weighted PCA.")
        mean_coords = np.average(data, axis=0, weights=weights)
        centered = data - mean_coords
        cov = np.dot((weights[:, None] * centered).T, centered) / weights.sum()
        eigvals, eigvecs = np.linalg.eigh(cov)
        idx = np.argsort(eigvals)[::-1]
        eigvals = eigvals[idx]
        eigvecs = eigvecs[:, idx]
        fractions = eigvals / eigvals.sum()

        # Projections
        proj = np.dot(centered, eigvecs)

        # Save projections
        logging.info("Saving projections to pca.csv.")
        pca_data = np.column_stack((np.arange(total_frames), proj[:, :3]))
        np.savetxt("pca.csv", pca_data, delimiter=" ", fmt="%d %.6f %.6f %.6f")

        # Save eigenvalues (top 3 for simplicity, as per outputs)
        logging.info("Saving eigenvalues to eigenvalue.csv.")
        eig_data = []
        for i in range(3):
            eig_data.append([i+1, eigvals[i], fractions[i]])
        np.savetxt("eigenvalue.csv", eig_data, delimiter=" ", fmt="%d %.6f %.6f")

        # Save nmd files
        logging.info("Saving normal modes to pcN.nmd.")
        for i in range(3):
            mode_file = f"pc{i+1}.nmd"
            with open(mode_file, "w") as f:
                f.write(f"name PC{i+1}\n")
                f.write("coordinates " + " ".join(f"{c:.6f}" for c in ref_coords.flatten()) + "\n")
                scale = np.sqrt(eigvals[i])
                mode_vec = eigvecs[:, i].reshape(natoms, 3).flatten()
                f.write(f"mode {i+1} {scale:.6f} " + " ".join(f"{v:.6f}" for v in mode_vec) + "\n")

        # Free energy surfaces per trajectory
        if T is not None:
            logging.info("Computing free energy surfaces.")
            kT = R * T

            # Compute global min/max for PCs
            min_pcs = np.min(proj[:, :3], axis=0)
            max_pcs = np.max(proj[:, :3], axis=0)

            # Pairs: (0,1,"12"), (0,2,"13"), (1,2,"23")
            pairs = [(0,1, "12"), (0,2, "13"), (1,2, "23")]

            # Dict to store data for plotting: suffix -> list of (trj_name, X, Y, fe)
            suffix_to_data = {suffix: [] for _, _, suffix in pairs}

            # First loop: compute FE and save CSV
            for traj_idx, (start, end) in enumerate(frame_ranges):
                trj_name = os.path.splitext(os.path.basename(traj_list[traj_idx]))[0]
                proj_traj = proj[start:end, :3]
                weights_traj = weights[start:end]
                sum_w = weights_traj.sum()

                for px, py, suffix in pairs:
                    logging.info(f"Processing PC{suffix} for {trj_name}.")
                    x_min, x_max = min_pcs[px], max_pcs[px]
                    y_min, y_max = min_pcs[py], max_pcs[py]
                    xedges = np.linspace(x_min, x_max, BINS + 1)
                    yedges = np.linspace(y_min, y_max, BINS + 1)
                    hist, _, _ = np.histogram2d(proj_traj[:, px], proj_traj[:, py], bins=(xedges, yedges), weights=weights_traj)
                    if smooth:
                        logging.info(f"Applying Gaussian smoothing to PC{suffix} for {trj_name}.")
                        hist = gaussian_filter(hist, sigma=SMOOTH_SIGMA)
                    prob = hist / sum_w
                    fe = -kT * np.log(prob + EPS)
                    fe[hist == 0] = np.nan  # Mask empty bins

                    # Save values
                    X, Y = np.meshgrid((xedges[:-1] + xedges[1:]) / 2, (yedges[:-1] + yedges[1:]) / 2)
                    csv_data = np.column_stack((X.ravel(), Y.ravel(), fe.ravel()))
                    np.savetxt(f"pc{suffix}_{trj_name}.csv", csv_data, delimiter=" ", fmt="%.6f %.6f %.6f")

                    # Store for plotting
                    suffix_to_data[suffix].append((trj_name, X, Y, fe))

            # Compute global min/max FE per suffix
            suffix_to_ranges = {}
            for suffix, data_list in suffix_to_data.items():
                if data_list:
                    min_fe = min(np.nanmin(fe) for _, _, _, fe in data_list)
                    max_fe = max(np.nanmax(fe) for _, _, _, fe in data_list)
                    suffix_to_ranges[suffix] = (min_fe, max_fe)

            # Second loop: generate plots with consistent ranges
            for suffix, data_list in suffix_to_data.items():
                vmin, vmax = suffix_to_ranges.get(suffix, (0, 0))
                for trj_name, X, Y, fe in data_list:
                    logging.info(f"Plotting PC{suffix} for {trj_name} with global ranges.")
                    fig, ax = plt.subplots()
                    cmap = plt.get_cmap("cubehelix")
                    im = ax.pcolormesh(X, Y, fe, cmap=cmap, shading="auto", vmin=vmin, vmax=vmax)
                    cbar = fig.colorbar(im, ax=ax)
                    px, py = int(suffix[0]) - 1, int(suffix[1]) - 1  # 1-based to 0-based
                    ax.set_xlabel(f"PC{px+1}", fontsize=14)
                    ax.set_ylabel(f"PC{py+1}", fontsize=14)
                    ax.tick_params(labelsize=14)
                    ax.yaxis.set_major_locator(MaxNLocator(nbins=10))
                    cbar.set_label("Free Energy (kJ/mol)", fontsize=14)
                    cbar.ax.tick_params(labelsize=14)
                    plt.savefig(f"pc{suffix}_{trj_name}.png")
                    plt.close()

        logging.info("All tasks completed.")

    except Exception as e:
        logging.error(f"Error occurred: {str(e)}")
        raise

    finally:
        # Memory usage
        current, peak = tracemalloc.get_traced_memory()
        logging.info(f"Peak memory usage: {peak / 10**6:.2f} MB")
        tracemalloc.stop()

if __name__ == "__main__":
    main()
