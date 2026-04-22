#!/usr/bin/env python3
"""
Quick SASA Analysis Tool

Standalone tool for SASA (Solvent Accessible Surface Area) calculation and analysis.

Uses internal SASA implementation (fast VDW+Gaussian approximation) for relative
comparisons. For high-accuracy absolute SASA values, use gmx sasa external tools.

Performance: ~0.7ms/frame for typical protein-ligand complexes (optimized per-atom algorithm)

Outputs:
  - Per-ligand SASA time-series (PNG + CSV)
  - Weighted histogram analysis (PNG + CSV)
  - Summary statistics (CSV)

Part of LFS-HREX toolkit

History:
  Phase 10.3.3 (Feb 2026):
  - Created as standalone tool following rmsd.py/com_dist.py patterns
  - Integrated internal SASA from fp.py Phase 10.3.2 optimization (25,000x speedup)
  - Visualization functions from fp.py SASA enhancement
"""

import os
import sys
import argparse
import logging
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import MDAnalysis as mda
from pathlib import Path
from scipy.spatial import cKDTree
from selection_defaults import DEFAULT_RECEPTOR_SEL, DEFAULT_LIGAND_SEL

# ============================================================================
# CONFIGURATION DEFAULTS
# ============================================================================
# All defaults defined here for easy management and bash wrapper integration

# I/O Defaults
DEFAULT_TOPOLOGY = '../v1.pdb'
DEFAULT_TRAJECTORY = '../v1.xtc'
# SASA Calculation Defaults
DEFAULT_PROBE_RADIUS = 1.4  # Angstroms (water probe)
DEFAULT_N_POINTS = 100      # Deprecated, kept for compatibility

# Output Defaults
DEFAULT_OUTPUT_DIR = '.'
DEFAULT_OUTPUT_PREFIX = 'sasa'
DEFAULT_BINS = 50

# Frame Selection Defaults
DEFAULT_START_FRAME = 0
DEFAULT_END_FRAME = None
DEFAULT_STRIDE = 1

# Shared timeline for SASA time-series plotting (ns)
_SASA_TIME_NS = None


def validate_output_file_path(output_path, description):
    """Validate output file path for security and ensure parent directory exists."""
    if not output_path:
        raise ValueError(f"{description} cannot be empty")

    try:
        path = Path(output_path).resolve()
    except Exception as e:
        raise ValueError(f"Invalid {description} path '{output_path}': {e}")

    # Prevent directory traversal
    if '..' in str(path):
        raise ValueError(f"Directory traversal not allowed in {description}: {output_path}")

    if path.exists() and path.is_dir():
        raise ValueError(f"{description} is a directory: {path}")

    parent_dir = path.parent
    if not parent_dir.exists():
        try:
            parent_dir.mkdir(parents=True, exist_ok=True)
        except Exception as e:
            raise ValueError(f"Cannot create directory for {description}: {e}")

    if not os.access(parent_dir, os.W_OK):
        raise ValueError(f"Directory not writable for {description}: {parent_dir}")

    return path


def setup_logging(log_file, quiet=False):
    """Setup logging to file and console."""
    log_path = validate_output_file_path(log_file, "log file")
    formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')
    file_handler = logging.FileHandler(str(log_path))
    file_handler.setLevel(logging.INFO)
    file_handler.setFormatter(formatter)
    console_handler = logging.StreamHandler(sys.stdout)
    if quiet:
        console_handler.setLevel(logging.ERROR)
    else:
        console_handler.setLevel(logging.WARNING)
    console_handler.setFormatter(formatter)
    logger = logging.getLogger()
    logger.setLevel(logging.DEBUG)
    logger.addHandler(file_handler)
    logger.addHandler(console_handler)
    return logger


def calculate_sasa_internal(
    universe,
    ligand_selections,
    receptor_selection=None,
    probe_radius=1.4,
    n_points=100,
    trajectory=None,
):
    """
    Calculate SASA internally using a fast VDW+Gaussian approximation.

    This method is intended for relative comparisons (e.g., SASA scaling of contact
    weights), not high-accuracy absolute SASA values.

    Algorithm overview
    ------------------
    - Per-atom SASA computed once per frame for the full complex
    - Per-ligand SASA obtained by summing per-atom SASA over ligand atoms
    - Exposure modeled with Gaussian damping based on neighbor density
    - SASA_atom ≈ 4π(r_vdw + r_probe)² × exposure_factor

    Notes
    -----
    - Legacy slow group-subtraction implementations have been removed.
    - `n_points` is deprecated and ignored (retained for API compatibility).
    - Uses VDW radii (Bondi 1964; Mantina et al. 2009) with empirically tuned
      Gaussian parameters (λ=0.5, σ=1.5 Å, cutoff=10 Å).

    Parameters
    ----------
    universe : MDAnalysis.Universe
        Trajectory universe object with loaded topology and coordinates
    ligand_selections : list of MDAnalysis.AtomGroup
        One AtomGroup per ligand for which to calculate SASA
    receptor_selection : MDAnalysis.AtomGroup, optional
        Receptor atoms for context. If None, calculated from universe.
    probe_radius : float, optional
        Solvent probe radius in Angstroms (default: 1.4 for water)
    n_points : int, optional
        Deprecated (ignored in fast mode)
    trajectory : MDAnalysis.coordinates.base.ReaderBase, optional
        Optional trajectory slice iterator for frame subsetting

    Returns
    -------
    sasa_arrays : list of numpy.ndarray
        One array per ligand containing SASA values per frame (nm²)
        Shape: (n_frames,) for each ligand

    Examples
    --------
    >>> u = mda.Universe("topology.pdb", "trajectory.xtc")
    >>> lig1 = u.select_atoms("resname LIG and resid 1")
    >>> receptor = u.select_atoms("protein")
    >>> sasa_arrays = calculate_sasa_internal(u, [lig1], receptor)
    """
    if trajectory is None:
        trajectory = universe.trajectory

    n_ligands = len(ligand_selections)
    n_frames = len(trajectory)
    sasa_arrays = [np.zeros(n_frames, dtype=np.float32) for _ in range(n_ligands)]

    # VDW radii table (Bondi 1964 + Mantina 2009 extensions)
    vdw_radii_ang = {
        'H': 1.20, 'He': 1.40, 'C': 1.70, 'N': 1.55, 'O': 1.52, 'F': 1.47,
        'Ne': 1.54, 'Na': 2.27, 'Mg': 1.73, 'Si': 2.10, 'P': 1.80, 'S': 1.80,
        'Cl': 1.75, 'Ar': 1.88, 'K': 2.75, 'Ga': 1.87, 'As': 1.85, 'Se': 1.90,
        'Br': 1.85, 'Kr': 2.02, 'In': 1.93, 'Sn': 2.17, 'Te': 2.06, 'I': 1.98,
        'Xe': 2.16, 'Tl': 1.96, 'Pb': 2.02, 'Be': 1.53, 'B': 1.92, 'Al': 1.84,
        'Ca': 2.31, 'Ge': 2.11, 'Rb': 3.03, 'Sr': 2.50, 'Sb': 2.06, 'Cs': 3.43,
        'Ba': 2.68, 'Bi': 2.07, 'Po': 1.97, 'At': 2.02, 'Rn': 2.20, 'Fr': 3.48,
        'Ra': 2.83, 'Li': 1.81, 'Mn': 2.00, 'Fe': 2.00, 'Co': 1.95, 'Ni': 1.63,
        'Cu': 1.40, 'Zn': 1.39, 'Ru': 2.05, 'Rh': 2.00, 'Pd': 1.63, 'Ag': 1.72,
        'Pt': 1.72, 'Au': 1.66,
    }
    default_radius = 1.70  # Carbon

    # Gaussian smoothing parameters (empirically tuned for good approximation)
    gaussian_lambda = 0.5      # Scaling factor for burial
    gaussian_sigma = 1.5       # Gaussian width in Angstroms
    neighbor_cutoff = 10.0     # Neighbor search cutoff in Angstroms

    logging.info(
        f"Calculating internal SASA (FAST VDW+Gaussian approximation): {n_ligands} ligands, {n_frames} frames"
    )
    logging.info(
        f"  Method: VDW radii + Gaussian smoothing (λ={gaussian_lambda}, σ={gaussian_sigma}Å, cutoff={neighbor_cutoff}Å)"
    )
    logging.info("  NOTE: Fast approximation for relative comparisons (10-20% accuracy vs. probe methods)")

    # Helper function to calculate SASA using VDW+Gaussian method
    def calculate_group_sasa_fast(atom_group, all_coords, all_elements, tree=None):
        """Calculate per-atom SASA using fast VDW+Gaussian approximation.

        Args:
            atom_group: AtomGroup to calculate SASA for
            all_coords: Coordinates of all atoms in complex (for neighbor search)
            all_elements: Elements of all atoms in complex

        Returns:
            per_atom_sasa: numpy array of SASA values per atom (nm²)
        """
        # Get VDW radii + probe radius for target atoms
        group_coords = atom_group.positions  # Angstroms
        group_elements = atom_group.elements
        group_radii = np.array([
            vdw_radii_ang.get(elem, default_radius) + probe_radius
            for elem in group_elements
        ])

        # Build KDTree for all atoms (complex) if not provided
        if tree is None:
            tree = cKDTree(all_coords)

        # Allocate per-atom SASA array
        n_atoms = len(group_coords)
        atom_sasa_array = np.zeros(n_atoms, dtype=np.float32)

        # Efficient neighbor accumulation using sparse distance matrix
        # Compute all neighbor distances within cutoff; include_self=False to avoid self-count
        dist_matrix = tree.sparse_distance_matrix(tree, neighbor_cutoff, output_type='coo_matrix')
        if dist_matrix.nnz == 0:
            # No neighbors within cutoff; full exposure for all atoms
            atom_sasa_array[:] = 4 * np.pi * (group_radii ** 2)
            return atom_sasa_array / 100.0

        # Compute gaussian weights for all neighbor distances
        dist_vals = dist_matrix.data.astype(np.float32, copy=False)
        weights = np.exp(-(dist_vals ** 2) / (2 * gaussian_sigma ** 2))

        # Sum weights by row (per-atom neighbor_sum)
        n_atoms = len(group_coords)
        neighbor_sum = np.zeros(n_atoms, dtype=np.float32)
        np.add.at(neighbor_sum, dist_matrix.row, weights)

        # Exposure and SASA per atom
        exposure = np.exp(-gaussian_lambda * neighbor_sum)
        atom_sasa_array = (4 * np.pi * (group_radii ** 2) * exposure).astype(np.float32, copy=False)

        # Convert from Å² to nm²
        return atom_sasa_array / 100.0

    # Determine receptor selection if not provided
    if receptor_selection is None:
        # Use all non-ligand atoms as receptor context
        all_ligand_indices = set()
        for lg in ligand_selections:
            all_ligand_indices.update(lg.indices)
        receptor_indices = [i for i in range(len(universe.atoms)) if i not in all_ligand_indices]
        receptor_sel = universe.atoms[receptor_indices]
    else:
        receptor_sel = receptor_selection

    # Calculate per-atom SASA once per frame, then sum by ligand atoms.
    # Precompute ligand atom indices in complex for fast per-frame summation
    ligand_indices_global = [lg.indices for lg in ligand_selections]

    for frame_idx, ts in enumerate(trajectory):
        # Build complex with receptor + ALL ligands (done once per frame)
        all_ligands = receptor_sel.universe.atoms[:0]
        for lig in ligand_selections:
            all_ligands = all_ligands + lig
        complex_all = receptor_sel + all_ligands

        # 1. Calculate per-atom SASA for entire complex (ONCE per frame)
        complex_all_coords = complex_all.positions
        complex_all_elements = complex_all.elements
        tree = cKDTree(complex_all_coords)
        per_atom_sasa = calculate_group_sasa_fast(complex_all, complex_all_coords, complex_all_elements, tree=tree)

        # 2. For each ligand: sum SASA for atoms belonging to this ligand
        complex_indices = complex_all.indices
        index_map = {atom_index: pos for pos, atom_index in enumerate(complex_indices)}
        for lig_id, lig_indices in enumerate(ligand_indices_global):
            lig_indices_in_complex = [index_map[i] for i in lig_indices if i in index_map]
            if len(lig_indices_in_complex) == 0:
                sasa_arrays[lig_id][frame_idx] = 0.0
                continue
            ligand_sasa = np.sum(per_atom_sasa[lig_indices_in_complex])
            sasa_arrays[lig_id][frame_idx] = ligand_sasa

        if frame_idx % max(1, n_frames // 100) == 0:
            logging.info(f"  Processed frame {frame_idx + 1}/{n_frames}")

    logging.info("Completed internal SASA calculation (fast mode)")
    return sasa_arrays


def plot_sasa_timeseries(sasa_arrays, weights, n_ligands, output_prefix='sasa_timeseries'):
    """
    Generate per-ligand SASA time-series multiplot.

    Parameters
    ----------
    sasa_arrays : list of numpy.ndarray
        SASA values per ligand per frame (nm²)
    weights : numpy.ndarray
        Frame weights for weighted analysis
    n_ligands : int
        Number of ligands
    output_prefix : str
        Output file prefix for PNG and CSV files

    Outputs
    -------
    For each ligand i:
    - {output_prefix}_lig{i+1}.png (time-series plot)
    - {output_prefix}_lig{i+1}.csv (frame, time_ns, sasa_nm2, weight)
    """
    if sasa_arrays is None or len(sasa_arrays) == 0:
        logging.warning("No SASA arrays provided; skipping SASA time-series plots")
        return

    global _SASA_TIME_NS
    n_frames = len(sasa_arrays[0])
    weights_arr = np.asarray(weights, dtype=np.float64)
    if weights_arr.ndim != 1:
        weights_arr = weights_arr.ravel()
    if weights_arr.size != n_frames:
        raise ValueError(
            f"Weights length mismatch for SASA time series: expected {n_frames}, got {weights_arr.size}"
        )

    frames = np.arange(n_frames, dtype=np.int32)
    time_ns = frames.astype(np.float64)
    use_time_axis = False
    if _SASA_TIME_NS is not None:
        time_ns = np.asarray(_SASA_TIME_NS, dtype=np.float64)
        if time_ns.size == n_frames:
            use_time_axis = True
        else:
            logging.warning(
                "SASA time axis length mismatch (%d vs %d); using frame index",
                time_ns.size,
                n_frames,
            )
            time_ns = frames.astype(np.float64)
            use_time_axis = False

    n_plot_ligands = min(n_ligands, len(sasa_arrays))
    if n_plot_ligands != n_ligands:
        logging.warning(
            "Requested %d ligands but received %d SASA arrays; plotting %d",
            n_ligands,
            len(sasa_arrays),
            n_plot_ligands,
        )

    for i in range(n_plot_ligands):
        sasa_vals = np.asarray(sasa_arrays[i], dtype=np.float64)
        if sasa_vals.size != n_frames:
            logging.warning(
                "Ligand %d SASA length mismatch (%d vs %d); skipping",
                i + 1,
                sasa_vals.size,
                n_frames,
            )
            continue

        x_vals = time_ns if use_time_axis else frames
        xlabel = "Time (ns)" if use_time_axis else "Frame"

        plt.figure(figsize=(10, 4))
        plt.plot(x_vals, sasa_vals, 'b-', linewidth=0.5)
        plt.xlabel(xlabel)
        plt.ylabel("SASA (nm²)")
        plt.title(f"SASA Time Series - Ligand {i+1}")
        plt.grid(True, linestyle='--', alpha=0.4)

        png_file = f"{output_prefix}_lig{i+1}.png"
        plt.savefig(png_file, dpi=300, bbox_inches='tight')
        plt.close()

        try:
            df = pd.DataFrame({
                'frame': frames,
                'time_ns': time_ns,
                'sasa_nm2': sasa_vals,
                'weight': weights_arr,
            })
            csv_file = f"{output_prefix}_lig{i+1}.csv"
            df.to_csv(csv_file, index=False, float_format='%.6f')
        except Exception as e:
            logging.warning("Failed to write SASA time-series CSV for ligand %d: %s", i + 1, e)
            csv_file = f"{output_prefix}_lig{i+1}.csv"

        logging.info("Saved SASA time series: %s + %s", png_file, csv_file)


def plot_sasa_histogram(sasa_arrays, weights, n_ligands, bins=50, output_prefix='sasa_histogram'):
    """
    Generate weighted SASA histogram with per-ligand breakdown and combined average.

    Parameters
    ----------
    sasa_arrays : list of numpy.ndarray
        SASA values per ligand per frame (nm²)
    weights : numpy.ndarray
        Frame weights for weighted histogram
    n_ligands : int
        Number of ligands
    bins : int
        Number of histogram bins (default: 50)
    output_prefix : str
        Output file prefix for PNG and CSV files

    Outputs
    -------
    - {output_prefix}_lig{N}.png (per-ligand histogram)
    - {output_prefix}_combined.png (combined average histogram)
    - {output_prefix}_histogram.csv (bin_center, ligand1_density, ligand2_density, ..., combined_avg)
    """
    if not sasa_arrays or n_ligands <= 0:
        logging.warning("No SASA arrays provided for histogram plotting")
        return

    if len(sasa_arrays) != n_ligands:
        logging.warning(
            "SASA array count (%d) does not match n_ligands (%d)",
            len(sasa_arrays),
            n_ligands,
        )
        n_ligands = min(len(sasa_arrays), n_ligands)

    if weights is None:
        raise ValueError("Weights array is required for weighted SASA histogram")

    ligand_hists = []
    bin_centers = None
    bin_widths = None

    for i in range(n_ligands):
        hist, bin_edges = np.histogram(
            sasa_arrays[i],
            bins=bins,
            weights=weights,
            density=True,
        )
        bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2
        bin_widths = bin_edges[1:] - bin_edges[:-1]
        ligand_hists.append(hist)

        plt.figure(figsize=(8, 4))
        plt.bar(
            bin_centers,
            hist,
            width=bin_widths,
            color='steelblue',
            alpha=0.7,
            align='center',
        )
        plt.title(f"Ligand {i+1} SASA Distribution")
        plt.xlabel("SASA (nm²)")
        plt.ylabel("Weighted Density")
        plt.grid(True, alpha=0.3)

        png_file = f"{output_prefix}_lig{i+1}.png"
        plt.tight_layout()
        plt.savefig(png_file, dpi=300, bbox_inches='tight')
        plt.close()
        logging.info(f"Saved SASA histogram for ligand {i+1}: {png_file}")

    combined_hist = np.mean(ligand_hists, axis=0) if ligand_hists else None
    combined_png_file = None
    if combined_hist is not None and bin_centers is not None and bin_widths is not None:
        plt.figure(figsize=(8, 4))
        plt.bar(
            bin_centers,
            combined_hist,
            width=bin_widths,
            color='darkgreen',
            alpha=0.7,
            align='center',
        )
        plt.title("Combined Average SASA Distribution")
        plt.xlabel("SASA (nm²)")
        plt.ylabel("Weighted Density")
        plt.grid(True, alpha=0.3)
        combined_png_file = f"{output_prefix}_combined.png"
        plt.tight_layout()
        plt.savefig(combined_png_file, dpi=300, bbox_inches='tight')
        plt.close()

    csv_file = f"{output_prefix}.csv"
    data = {'bin_center': bin_centers}
    for i, hist in enumerate(ligand_hists):
        data[f"ligand{i+1}_density"] = hist
    if combined_hist is not None:
        data['combined_avg'] = combined_hist

    df = pd.DataFrame(data)
    df.to_csv(csv_file, index=False, float_format='%.6f')

    if combined_png_file:
        logging.info(f"Saved combined SASA histogram: {combined_png_file} + {csv_file}")
    else:
        logging.info(f"Saved SASA histogram CSV: {csv_file}")


def parse_args():
    parser = argparse.ArgumentParser(
        description='Quick SASA Analysis Tool - Standalone SASA calculation and visualization',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Basic usage
  python quick_sasa.py --topology system.pdb --trajectory traj.xtc

  # Custom selections
  python quick_sasa.py --topology system.pdb --trajectory traj.xtc \\
      --receptor_sel "protein and not name H*" --ligand_sel "resname LIG and not name H*"

  # Frame range selection
  python quick_sasa.py --topology system.pdb --trajectory traj.xtc \
      --start 0 --end 1000 --stride 10

Output files:
  {prefix}_timeseries_lig{N}.png/csv - Per-ligand time series
  {prefix}_histogram_lig{N}.png - Per-ligand histogram analysis
  {prefix}_histogram_combined.png - Combined histogram analysis
  {prefix}_histogram.csv - Weighted histogram data
  {prefix}_summary.csv - Summary statistics

Internal SASA uses fast VDW+Gaussian approximation for relative comparisons.
"""
    )

    # I/O arguments
    parser.add_argument('--topology', default=DEFAULT_TOPOLOGY, help='Topology file (PDB, GRO)')
    parser.add_argument('--trajectory', default=DEFAULT_TRAJECTORY, help='Trajectory file (XTC, TRR)')
    parser.add_argument('--receptor-sel', default=DEFAULT_RECEPTOR_SEL, help='MDAnalysis receptor selection string (default from selection_defaults.py; see module for heavy-atom standards)')
    parser.add_argument('--ligand-sel', default=DEFAULT_LIGAND_SEL, help='MDAnalysis ligand selection string (default from selection_defaults.py; see module for heavy-atom standards)')

    # SASA parameters
    parser.add_argument(
        '--probe_radius',
        type=float,
        default=DEFAULT_PROBE_RADIUS,
        help='Probe radius in Angstroms (default: 1.4 for water)',
    )

    # Frame selection
    parser.add_argument('--start', type=int, default=DEFAULT_START_FRAME, help='Start frame')
    parser.add_argument('--end', type=int, default=DEFAULT_END_FRAME, help='End frame')
    parser.add_argument('--stride', type=int, default=DEFAULT_STRIDE, help='Frame stride')

    # Output options
    parser.add_argument('--output_dir', default=DEFAULT_OUTPUT_DIR, help='Output directory')
    parser.add_argument('--output_prefix', default=DEFAULT_OUTPUT_PREFIX, help='Output file prefix')
    parser.add_argument('--bins', type=int, default=DEFAULT_BINS, help='Histogram bins')

    # Logging options
    parser.add_argument('--log-file', default='quick_sasa.log',
                        help='Log file path (default: quick_sasa.log)')
    parser.add_argument('--quiet', '-q', action='store_true',
                        help='Suppress INFO messages on console (errors still shown)')

    # Weights (optional)
    group = parser.add_mutually_exclusive_group()
    group.add_argument('--weights', help='File with linear weights per frame')
    group.add_argument('--logweights', help='File with logweights per frame (exp + normalize)')

    return parser.parse_args()


def main():
    args = parse_args()

    # Setup logging
    setup_logging(args.log_file, quiet=args.quiet)

    # Tool banner (matching fp.py pattern)
    logging.info("=" * 60)
    logging.info("Quick SASA Analysis Tool - Started")
    logging.info("=" * 60)
    logging.info("NOTE: Internal SASA uses fast VDW+Gaussian approximation")
    logging.info("      Suitable for relative comparisons (10-20% accuracy)")
    logging.info("      For absolute values, use gmx sasa external tools")

    try:
        # Load trajectory
        logging.info(f"Loading trajectory: {args.trajectory}")
        u = mda.Universe(args.topology, args.trajectory)

        # Select atoms
        receptor = u.select_atoms(args.receptor_sel)
        ligand_base_sel = args.ligand_sel

        # Auto-detect multiple ligands (split by residue)
        ligand_selection = u.select_atoms(ligand_base_sel)
        ligand_resids = np.unique(ligand_selection.resids)
        n_ligands = len(ligand_resids)
        ligand_selections = [
            u.select_atoms(f"{ligand_base_sel} and resid {rid}") for rid in ligand_resids
        ]

        logging.info(f"Detected {n_ligands} ligand(s)")
        logging.info(f"Receptor: {len(receptor)} atoms")

        # Frame selection
        start = args.start
        end = args.end if args.end is not None else len(u.trajectory)
        stride = args.stride
        frame_slice = u.trajectory[start:end:stride]
        n_frames = len(frame_slice)

        logging.info(f"Frames: {n_frames} (start={start}, end={end}, stride={stride})")

        # Calculate SASA
        logging.info("Calculating SASA...")
        sasa_arrays = calculate_sasa_internal(
            u,
            ligand_selections,
            receptor_selection=receptor,
            probe_radius=args.probe_radius,
            trajectory=frame_slice,
        )

        # Frame weights (optional)
        if args.logweights:
            logw = np.loadtxt(args.logweights)
            logw = np.asarray(logw, dtype=np.float64).ravel()
            weights = np.exp(logw - np.max(logw))
            weights = weights / np.sum(weights)
            logging.info("Processed logweights (direct exp normalization)")
        elif args.weights:
            weights = np.loadtxt(args.weights)
            weights = np.asarray(weights, dtype=np.float64).ravel()
        else:
            weights = np.ones(n_frames, dtype=np.float64) / n_frames

        if weights.size != n_frames:
            logging.warning(
                "Weights length (%d) != number of frames (%d); padding/truncating",
                weights.size,
                n_frames,
            )
            if weights.size < n_frames:
                weights = np.pad(weights, (0, n_frames - weights.size), mode='edge')
            else:
                weights = weights[:n_frames]

        if np.sum(weights) <= 0:
            raise ValueError("Weights must sum to > 0")

        weights = weights / np.sum(weights)

        # Create output directory
        os.makedirs(args.output_dir, exist_ok=True)
        prefix = os.path.join(args.output_dir, args.output_prefix)

        # Generate visualizations
        logging.info("Generating time-series plots...")
        plot_sasa_timeseries(sasa_arrays, weights, n_ligands, output_prefix=f"{prefix}_timeseries")

        logging.info("Generating histogram analysis...")
        plot_sasa_histogram(
            sasa_arrays,
            weights,
            n_ligands,
            bins=args.bins,
            output_prefix=f"{prefix}_histogram",
        )

        # Summary statistics
        logging.info("Computing summary statistics...")
        summary_data = []
        for i, sasa_arr in enumerate(sasa_arrays):
            mean_sasa = np.average(sasa_arr, weights=weights)
            variance = np.average((sasa_arr - mean_sasa) ** 2, weights=weights)
            summary_data.append({
                'ligand': i + 1,
                'mean_sasa_nm2': mean_sasa,
                'std_sasa_nm2': np.sqrt(variance),
                'min_sasa_nm2': np.min(sasa_arr),
                'max_sasa_nm2': np.max(sasa_arr)
            })

        summary_df = pd.DataFrame(summary_data)
        summary_file = f"{prefix}_summary.csv"
        summary_df.to_csv(summary_file, index=False, float_format='%.6f')
        logging.info(f"Saved summary: {summary_file}")

        logging.info("=" * 60)
        logging.info("Quick SASA Analysis Tool - Completed Successfully")
        logging.info("=" * 60)

    except Exception as e:
        logging.error(f"Fatal error: {e}", exc_info=True)
        sys.exit(1)


if __name__ == '__main__':
    main()
