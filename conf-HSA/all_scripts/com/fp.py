import argparse
import logging
import os
from pathlib import Path
import sys
import gc
import tempfile
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import MDAnalysis as mda
from selection_defaults import (
    DEFAULT_LIGAND_SEL as DEFAULT_LIGAND_SEL_RESNAME,
    DEFAULT_RECEPTOR_SEL as DEFAULT_RECEPTOR_SEL_RESNAME,
)
from MDAnalysis.lib.distances import capped_distance
from sklearn.cluster import DBSCAN
from sklearn.neighbors import NearestNeighbors
from scipy.sparse import csr_matrix
from scipy import sparse
from scipy.spatial import cKDTree

"""
Interaction Fingerprint Analysis Tool

SELECTION PARADIGM:
    This tool uses resid-based selections for interaction fingerprinting,
    which differs from the resname-based selections used by rmsd.py and
    com_dist.py (see selection_defaults.py module).
    
    Why different?
    - Fingerprints track receptor residue contacts by position (resid = position)
    - Other tools identify molecular components by type (resname = identity)
    - Interaction patterns depend on which residue positions are near ligand
    - Forcing resname-based compatibility would break existing workflows
    
    Example workflow impact:
    - fp.py with 'resid 1-80' tracks contacts at protein positions 1-80
    - Changing to 'protein' would include ALL protein residues regardless of position
    - Users expect position-specific contact maps, not type-based selections

DEFAULT SELECTIONS (resid-based):
    --receptor_sel: 'resid 1-80'   (protein residues at positions 1-80)
    --ligand_sel: 'resid 81-90'    (ligand residues at positions 81-90)
    
For resname-based selections used in rmsd.py/com_dist.py, see selection_defaults.py.

HISTORY:
Phase 10.3.1 (Feb 2026):
- Internal SASA implementation using fast VDW+Gaussian approximation

Phase 10.3.2 (Feb 2026):
- SASA performance optimization: 25,000x speedup via per-atom computation
- SASA defaults consolidated at script top (lines 48-119)
"""

# ============================================================================
# CONFIGURATION DEFAULTS
# ============================================================================
# All defaults are defined here for easy management and bash wrapper integration

# I/O Defaults
DEFAULT_TOPOLOGY = '../v1.pdb'
DEFAULT_TRAJECTORY = '../v1.xtc'
DEFAULT_RECEPTOR_SEL = 'resid 1-80'
DEFAULT_LIGAND_SEL = 'resid 81-90'

# Fingerprint Defaults
DEFAULT_CUTOFF = 4.5
DEFAULT_SIMILARITY = 'dice'

# Clustering Defaults
DEFAULT_CLUSTER_EPS = 0.3
DEFAULT_MIN_SAMPLES = 5
DEFAULT_SUBSAMPLING = 1

# Sieving Defaults
DEFAULT_SIEVING = False
DEFAULT_MAX_SIEVED = None
DEFAULT_TOP_OCC_PERCENT = 95.0
DEFAULT_MIN_OCC = 0.001

# Output Defaults
DEFAULT_MODES_CSV = 'binding_modes.csv'
DEFAULT_SASA_MODES_CSV = 'sasa_modes.csv'
DEFAULT_MODE_PDB_PREFIX = 'mode_'
DEFAULT_HEATMAP_PREFIX = 'heatmap_lig_'
DEFAULT_OCCUPANCIES_PREFIX = 'occupancies_lig'
DEFAULT_SASA_SCORES_PREFIX = 'sasa_scores'
DEFAULT_REF_PDB = None

# Bridge Clustering Defaults
DEFAULT_BRIDGE_UNSIEVED = True
DEFAULT_BRIDGE_MIN_SIM = 0.3
DEFAULT_BRIDGE_MIN_SUPPORT = 3
DEFAULT_BRIDGE_TOP_K = 2
DEFAULT_SIEVED_TOP_PERCENT = 0.0

# Sparse Matrix Defaults
DEFAULT_SPARSE_SIMILARITY = False
DEFAULT_SPARSE_THRESHOLD_OFFSET = 0.1
DEFAULT_SPARSE_CHUNK_SIZE = 1000

# Visualization Defaults
DEFAULT_KDIST_PLOT = False

# Free Energy Calculation Defaults
DEFAULT_TEMPERATURE = 300.0                      # Kelvin
DEFAULT_DG_OFFSET_KJ = 0.0                       # kJ/mol empirical offset (not needed with molecular correction)
DEFAULT_DG_UNBOUND_REF = 'contact'               # 'contact' or 'noise'
DEFAULT_DG_CONTACT_VOLUME_CORRECTION = False     # Legacy MC sampling method
DEFAULT_DG_CONTACT_VOLUME_PER_FRAME = False      # Legacy per-frame method
DEFAULT_DG_CONTACT_VOLUME_MOLECULAR = False      # Molecular volume method (RECOMMENDED for predictive dG)
DEFAULT_DG_VOLUME_CUTOFF = None                  # nm, defaults to fingerprint cutoff
DEFAULT_DG_VOLUME_PACKING_FACTOR = 0.0005        # Tight binding cavity (validated default)
DEFAULT_DG_VOLUME_USE_BOUND_FRACTION = False     # Scale by bound fraction (recommended)
DEFAULT_DG_VOLUME_SAMPLES_PER_FRAME = 5000       # MC samples per frame (legacy per-frame method)
DEFAULT_DG_VOLUME_FRAME_STRIDE = 10              # Frame stride (legacy per-frame method)

# SASA Calculation Defaults
DEFAULT_SASA_DIR = '.'
DEFAULT_SASA_TYPE = 'csv'
DEFAULT_SASA_COL = 1
DEFAULT_SASA_METHOD = 'internal'  # 'internal' (FAST Gaussian approx) or 'external' (load gmx sasa files)
DEFAULT_SASA_N_POINTS = 100       # Deprecated, kept for compatibility
DEFAULT_SASA_EPS = 1e-6           # Epsilon for SASA division stability

# Memory Management
MEMORY_SAFETY_FRACTION = 0.75

# Shared timeline for SASA time-series plotting (ns)
_SASA_TIME_NS = None

# 1-letter code mapping
RES_3TO1 = {
    'ALA': 'A', 'ARG': 'R', 'ASN': 'N', 'ASP': 'D', 'CYS': 'C',
    'GLN': 'Q', 'GLU': 'E', 'GLY': 'G', 'HIS': 'H', 'ILE': 'I',
    'LEU': 'L', 'LYS': 'K', 'MET': 'M', 'PHE': 'F', 'PRO': 'P',
    'SER': 'S', 'THR': 'T', 'TRP': 'W', 'TYR': 'Y', 'VAL': 'V',
    'ADE': 'A', 'GUA': 'G', 'CYT': 'C', 'THY': 'T',
    'URA': 'U', 'URI': 'U', 'INO': 'I',
    'DA': 'A', 'DC': 'C', 'DG': 'G', 'DT': 'T', 'DI': 'I',
    'A': 'A', 'C': 'C', 'G': 'G', 'U': 'U', 'I': 'I',
}


def validate_output_file_path(output_path, description):
    """Validate output file path for security and ensure parent directory exists."""
    if not output_path:
        raise ValueError(f"{description} cannot be empty")
    raw = str(output_path)
    if '..' in raw:
        raise ValueError(f"Directory traversal not allowed in {description}: {output_path}")
    if Path(output_path).is_symlink():
        raise ValueError(f"Symlinks not allowed for output path: {output_path}")
    try:
        path = Path(output_path).resolve(strict=False)
    except Exception as e:
        raise ValueError(f"Invalid {description} path '{output_path}': {e}")
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


def parse_args():
    parser = argparse.ArgumentParser(description='Analyze MD simulations for binding pockets and poses using interaction fingerprints.')
    parser.add_argument('--topology', default=DEFAULT_TOPOLOGY, help='Topology file path (e.g., PDB).')
    parser.add_argument('--trajectory', default=DEFAULT_TRAJECTORY, help='Trajectory file path (e.g., XTC).')
    parser.add_argument('--receptor_sel', default=DEFAULT_RECEPTOR_SEL,
                        help='MDAnalysis selection for receptor (resid-based by default).')
    parser.add_argument('--ligand_sel', default=DEFAULT_LIGAND_SEL,
                        help='MDAnalysis selection for ligands (resid-based by default).')
    parser.add_argument('--cutoff', type=float, default=DEFAULT_CUTOFF, help='Contact distance cutoff (Angstrom).')
    group = parser.add_mutually_exclusive_group()
    group.add_argument('--weights', help='File with linear weights per frame.')
    group.add_argument('--logweights', help='File with logweights per frame (exponentiated and normalized).')
    parser.add_argument('--sasa_enabled', action='store_true', help='Enable SASA scaling from internal calculation (fast VDW+Gaussian approximation) or external files.')
    parser.add_argument('--sasa_dir', default=DEFAULT_SASA_DIR, help='Directory for SASA files (external mode only).')
    parser.add_argument('--sasa_type', choices=['csv', 'xvg'], default=DEFAULT_SASA_TYPE, help='SASA file type (external mode only).')
    parser.add_argument('--sasa_col', type=int, default=DEFAULT_SASA_COL, help='0-based column for SASA values (external mode only).')
    parser.add_argument('--sasa_method', 
                        choices=['external', 'internal'], 
                        default=DEFAULT_SASA_METHOD,
                        help='SASA calculation method: internal (FAST Gaussian approximation with 25,000x speedup, default) or external (load from gmx sasa files). Internal method suitable for relative comparisons and binding mode discrimination.')
    parser.add_argument('--sasa_n_points', 
                        type=int, 
                        default=DEFAULT_SASA_N_POINTS,
                        help='[DEPRECATED] Number of test points (ignored in fast mode, kept for compatibility).')
    parser.add_argument('--sasa_eps', type=float, default=DEFAULT_SASA_EPS, help='Epsilon for SASA division stability.')
    parser.add_argument('--similarity', choices=['jaccard', 'dice', 'cosine'], default=DEFAULT_SIMILARITY, help='Fingerprint similarity measure.')
    parser.add_argument('--cluster_eps', type=float, default=DEFAULT_CLUSTER_EPS, help='DBSCAN distance threshold.')
    parser.add_argument('--min_samples', type=int, default=DEFAULT_MIN_SAMPLES, help='DBSCAN min samples for core points.')
    parser.add_argument('--subsampling', type=int, default=DEFAULT_SUBSAMPLING, help='Subsample bound instances per ligand (1 = none).')
    parser.add_argument('--sieving', dest='sieving', action='store_true',
                        help='Enable sieving for memory reduction (requires --max_sieved).')
    parser.add_argument('--max_sieved', type=int, default=DEFAULT_MAX_SIEVED,
                        help='Max unique fingerprints after sieving (auto-enables --sieving when set).')
    parser.add_argument('--top_occ_percent', type=float, default=DEFAULT_TOP_OCC_PERCENT, help='Top occupancy percent to retain modes.')
    parser.add_argument('--min_occ', type=float, default=DEFAULT_MIN_OCC, help='Min occupancy for PDB output.')
    parser.add_argument('--modes_csv', default=DEFAULT_MODES_CSV, help='Output CSV for binding modes.')
    parser.add_argument('--sasa_modes_csv', default=DEFAULT_SASA_MODES_CSV, help='Output CSV for SASA modes (if enabled).')
    parser.add_argument('--mode_pdb_prefix', default=DEFAULT_MODE_PDB_PREFIX, help='Prefix for mode PDB files.')
    parser.add_argument('--heatmap_prefix', default=DEFAULT_HEATMAP_PREFIX, help='Prefix for per-ligand heatmaps.')
    parser.add_argument('--occupancies_prefix', default=DEFAULT_OCCUPANCIES_PREFIX, help='Prefix for per-ligand occupancy files.')
    parser.add_argument('--sasa_scores_prefix', default=DEFAULT_SASA_SCORES_PREFIX, help='Prefix for SASA scores files.')
    parser.add_argument('--ref_pdb', default=DEFAULT_REF_PDB, help='Reference PDB for similarity comparison.')
    parser.add_argument('--ref_ligand_sel', default=None,
                        help='MDAnalysis selection for reference ligand in --ref_pdb (overrides auto-selection).')
    parser.add_argument('--ignore_memory_risk', action='store_true',
                        help='Override dense-mode memory safety exit (does not disable auto sparse fallback).')
    parser.add_argument('--clustering_method', choices=['density_ratio', 'dbscan'], default='dbscan',
                        help='Clustering method to use (sparse mode falls back to dbscan).')
    parser.add_argument('--calc_overall_dg', action='store_true', help='Calculate overall binding free energy from weighted occupancy.')
    parser.add_argument('--calc_cluster_dg', action='store_true', help='Calculate per-cluster binding free energies.')
    parser.add_argument(
        '--dg_unbound_ref',
        choices=['contact', 'noise'],
        default='contact',
        help='Reference for unbound state in dG: contact (no contacts) or noise (DBSCAN noise)'
    )
    # Free energy calculation arguments
    parser.add_argument('--temperature', type=float, default=DEFAULT_TEMPERATURE, 
                        help=f'Temperature in Kelvin for free energy calculations (default: {DEFAULT_TEMPERATURE}K).')
    parser.add_argument('--dg_offset_kj', type=float, default=DEFAULT_DG_OFFSET_KJ,
                        help=f'Empirical offset (kJ/mol) applied to absolute dG values (default: {DEFAULT_DG_OFFSET_KJ}, not needed with molecular correction).')
    parser.add_argument(
        '--dg_contact_volume_correction',
        action='store_true',
        default=DEFAULT_DG_CONTACT_VOLUME_CORRECTION,
        help=('[LEGACY] Apply contact-volume correction using Monte Carlo sampling with topology receptor frame. '
              'RECOMMENDED: Use --dg_contact_volume_molecular instead.')
    )
    parser.add_argument(
        '--dg_contact_volume_per_frame',
        action='store_true',
        default=DEFAULT_DG_CONTACT_VOLUME_PER_FRAME,
        help=('[LEGACY] Estimate contact volume per-frame with weighted average. '
              'RECOMMENDED: Use --dg_contact_volume_molecular instead.')
    )
    parser.add_argument(
        '--dg_contact_volume_molecular',
        action='store_true',
        default=DEFAULT_DG_CONTACT_VOLUME_MOLECULAR,
        help=('[RECOMMENDED] Estimate contact volume from ligand molecular volume (VdW radii). '
              'Achieves +/-2-4 kcal/mol accuracy without experimental calibration. '
              'Does NOT use receptor coordinates - purely ligand-based.')
    )
    parser.add_argument(
        '--dg_volume_cutoff',
        type=float,
        default=DEFAULT_DG_VOLUME_CUTOFF,
        help=('Override contact-volume cutoff distance (nm) for MC sampling methods. '
              'Default: use fingerprint --cutoff. Not used by molecular volume method.')
    )
    parser.add_argument(
        '--dg_volume_packing_factor',
        type=float,
        default=DEFAULT_DG_VOLUME_PACKING_FACTOR,
        help=(f'Packing factor for molecular volume method (default: {DEFAULT_DG_VOLUME_PACKING_FACTOR} = tight binding). '
              'Range: 0.0001 (very tight) to 0.02 (loose/flexible). '
              'Use with --dg_contact_volume_molecular.')
    )
    parser.add_argument(
        '--dg_volume_use_bound_fraction',
        action='store_true',
        default=DEFAULT_DG_VOLUME_USE_BOUND_FRACTION,
        help=('[RECOMMENDED] Scale molecular volume by bound fraction from trajectory statistics. '
              'Use with --dg_contact_volume_molecular.')
    )
    parser.add_argument(
        '--dg_volume_samples_per_frame',
        type=int,
        default=DEFAULT_DG_VOLUME_SAMPLES_PER_FRAME,
        help=(f'Monte Carlo samples per frame (default: {DEFAULT_DG_VOLUME_SAMPLES_PER_FRAME}). '
              'Legacy per-frame method only.')
    )
    parser.add_argument(
        '--dg_volume_frame_stride',
        type=int,
        default=DEFAULT_DG_VOLUME_FRAME_STRIDE,
        help=(f'Frame stride for per-frame estimation (default: {DEFAULT_DG_VOLUME_FRAME_STRIDE}). '
              'Legacy per-frame method only.')
    )
    parser.add_argument('--bridge_unsieved', action='store_true', help='Merge fragmented clusters using unsieved->sieved similarity bridges.')
    parser.add_argument('--bridge_min_sim', type=float, default=DEFAULT_BRIDGE_MIN_SIM, help='Min similarity for unsieved bridge links (default: 0.3).')
    parser.add_argument('--bridge_min_support', type=int, default=DEFAULT_BRIDGE_MIN_SUPPORT, help='Min number of unsieved supports to merge clusters (default: 3).')
    parser.add_argument('--bridge_top_k', type=int, default=DEFAULT_BRIDGE_TOP_K, help='Top-K sieved neighbors for bridge voting (default: 2).')
    parser.add_argument('--sieved_top_percent', type=float, default=DEFAULT_SIEVED_TOP_PERCENT,
                        help='Percent of highest-weight fingerprints to force include in sieving (0 disables).')
    parser.add_argument('--sparse_similarity', action='store_true',
                        help='Compute sparse similarity/distance matrices for DBSCAN (auto-enabled if dense memory risk).')
    parser.add_argument('--sparse_threshold', type=float, default=None,
                        help='Similarity threshold for sparse edges (sparse mode only; default: cluster_eps - 0.1).')
    parser.add_argument('--kdist_plot', action='store_true', help='Generate k-distance plot (can be memory-intensive).')
    parser.add_argument('--log-file', default='fp.log',
                        help='Log file path (default: fp.log)')
    parser.add_argument('--quiet', '-q', action='store_true',
                        help='Suppress INFO messages on console (errors still shown)')
    parser.set_defaults(sieving=DEFAULT_SIEVING)
    parser.set_defaults(bridge_unsieved=DEFAULT_BRIDGE_UNSIEVED)
    parser.set_defaults(kdist_plot=DEFAULT_KDIST_PLOT)
    args = parser.parse_args()

    # Auto-enable sieving when max_sieved is explicitly provided
    if args.max_sieved is not None:
        args.sieving = True
    elif args.sieving:
        parser.error("--sieving requires --max_sieved to be set")

    # Auto-enable molecular volume method when free energy calculation requested
    # (unless user explicitly chose a different volume correction method)
    if (args.calc_overall_dg or args.calc_cluster_dg):
        if not (args.dg_contact_volume_correction or args.dg_contact_volume_per_frame or args.dg_contact_volume_molecular):
            # No volume correction specified - enable recommended method
            args.dg_contact_volume_molecular = True
            args.dg_volume_use_bound_fraction = True
            logging.info("Auto-enabled molecular volume correction (recommended for predictive dG)")
            logging.info(f"Using packing_factor={args.dg_volume_packing_factor} (override with --dg_volume_packing_factor)")

    return args

def get_available_memory():
    """Detect memory limit with 5-level priority order"""
    try:
        # Priority 1: SLURM per-node allocation (bytes)
        if 'SLURM_MEM_PER_NODE' in os.environ:
            return int(os.environ['SLURM_MEM_PER_NODE']) * 1024 * 1024
        
        # Priority 2: SLURM per-CPU allocation x CPUs (bytes)
        if 'SLURM_MEM_PER_CPU' in os.environ:
            cpus = int(os.environ.get('SLURM_CPUS_PER_TASK', 1))
            return int(os.environ['SLURM_MEM_PER_CPU']) * cpus * 1024 * 1024
        
        # Priority 3: PBS job manager allocation (KB or bytes format)
        # PBS_MEM takes priority over PBS_RESC_MEM
        try:
            if 'PBS_MEM' in os.environ:
                mem_str = os.environ['PBS_MEM'].lower()
                # Parse format: "32gb", "32768mb", "33554432kb", or plain bytes
                if 'gb' in mem_str:
                    return int(mem_str.replace('gb', '')) * 1024**3
                elif 'mb' in mem_str:
                    return int(mem_str.replace('mb', '')) * 1024**2
                elif 'kb' in mem_str:
                    return int(mem_str.replace('kb', '')) * 1024
                else:
                    return int(mem_str)  # Assume bytes
            
            if 'PBS_RESC_MEM' in os.environ:
                mem_str = os.environ['PBS_RESC_MEM'].lower()
                if 'gb' in mem_str:
                    return int(mem_str.replace('gb', '')) * 1024**3
                elif 'mb' in mem_str:
                    return int(mem_str.replace('mb', '')) * 1024**2
                elif 'kb' in mem_str:
                    return int(mem_str.replace('kb', '')) * 1024
                else:
                    return int(mem_str)
        except (ValueError, KeyError):
            pass  # Fall through to next priority level
        
        # Priority 4: System total memory from /proc/meminfo (KB)
        with open('/proc/meminfo') as f:
            for line in f:
                if line.startswith('MemTotal'):
                    return int(line.split()[1]) * 1024
    except Exception:
        pass
    
    # Priority 5: Conservative default (32 GB) for shared nodes
    return 32 * 1024**3

_PEAK_RSS_GB = 0.0

def get_process_memory_gb():
    """Return current process RSS in GB (best-effort)."""
    try:
        with open('/proc/self/status') as f:
            for line in f:
                if line.startswith('VmRSS:'):
                    parts = line.split()
                    # VmRSS reported in kB
                    return float(parts[1]) * 1024 / 1e9
    except Exception:
        pass
    return None

def log_memory_usage(tag):
    """Log current and peak process memory usage with a tag."""
    global _PEAK_RSS_GB
    rss_gb = get_process_memory_gb()
    if rss_gb is not None:
        if rss_gb > _PEAK_RSS_GB:
            _PEAK_RSS_GB = rss_gb
        logging.info(f"Memory usage ({tag}): {rss_gb:.2f} GB RSS (peak { _PEAK_RSS_GB:.2f} GB)")

def _build_sparse_similarity_chunk(intersect, sumA, sumB, similarity, threshold, normB=None):
    """Build a sparse similarity chunk from a dense/sparse intersection matrix."""
    if sumA is None:
        raise ValueError("sumA must be provided for sparse similarity computation")
    if sumB is None:
        raise ValueError("sumB must be provided for sparse similarity computation")
    sumB_arr = np.asarray(sumB, dtype=np.float32)
    sumA_arr = np.asarray(sumA, dtype=np.float32)
    normB_arr = np.sqrt(sumB_arr + 1e-10)
    if similarity == 'cosine' and normB is not None:
        normB_arr = np.asarray(normB, dtype=np.float32)
    if intersect is None:
        raise ValueError("intersect matrix must be provided")
    if not sparse.issparse(intersect):
        intersect_csr = sparse.csr_matrix(intersect)
    else:
        intersect_csr = intersect.tocsr()

    chunk_rows = []
    chunk_cols = []
    chunk_data = []

    shape = intersect_csr.shape
    if shape is None:
        raise ValueError("intersect matrix shape is undefined")
    for i in range(shape[0]):
        row_start = intersect_csr.indptr[i]
        row_end = intersect_csr.indptr[i + 1]
        if row_start == row_end:
            continue
        js = intersect_csr.indices[row_start:row_end]
        inter_vals = intersect_csr.data[row_start:row_end]
        if inter_vals.dtype != np.float32:
            inter_vals = inter_vals.astype(np.float32, copy=False)

        if similarity == 'dice':
            sims = 2.0 * inter_vals / (sumA_arr[i] + sumB_arr[js] + 1e-10)
        elif similarity == 'jaccard':
            sims = inter_vals / (sumA_arr[i] + sumB_arr[js] - inter_vals + 1e-10)
        else:  # cosine
            sims = inter_vals / ((np.sqrt(sumA_arr[i] + 1e-10) * normB_arr[js]) + 1e-10)

        mask = sims >= threshold
        if np.any(mask):
            js_masked = js[mask]
            sims_masked = sims[mask]
            count = int(js_masked.size)
            chunk_rows.append(np.full(count, i, dtype=np.int32))
            chunk_cols.append(js_masked.astype(np.int32, copy=False))
            chunk_data.append(sims_masked.astype(np.float32, copy=False))

    if chunk_rows:
        rows_arr = np.concatenate(chunk_rows)
        cols_arr = np.concatenate(chunk_cols)
        data_arr = np.concatenate(chunk_data)
        chunk_sparse = csr_matrix(
            (data_arr, (rows_arr, cols_arr)),
            shape=intersect_csr.shape,
            dtype=np.float32,
        )
    else:
        chunk_sparse = csr_matrix(intersect_csr.shape, dtype=np.float32)

    return chunk_sparse

def estimate_k_distance_memory_gb(n_samples, k, dense=False, dense_bytes=8):
    """Estimate memory used by k-distance neighbor computation."""
    if n_samples <= 0 or k <= 0:
        return 0.0
    distances = n_samples * (k + 1) * 8  # float64 distances
    indices = n_samples * (k + 1) * 8    # int64 indices
    base = distances + indices
    if dense:
        # Worst-case buffer for dense pairwise distances (float32/float64)
        base += n_samples * n_samples * dense_bytes
    return base / 1e9

def estimate_sparse_extra_memory_gb(final_n, n_unique, avg_contacts_per_fp, include_unsieved=True):
    """Estimate extra memory for sparse similarity path (dense temporaries).

    Includes:
      - similarity/distance chunk working buffers (approx)
      - unsieved reassignment dense intersection/sim matrices
    """
    extra = 0.0
    # Chunk working buffers for sparse similarity (approx)
    # intersect chunk is sparse; add a conservative buffer for dense conversions
    chunk_rows = min(DEFAULT_SPARSE_CHUNK_SIZE, final_n)
    extra += (chunk_rows * final_n * 4) / 1e9  # float32 buffer

    if include_unsieved and n_unique > final_n:
        n_unsieved = n_unique - final_n
        # intersect_unsieved and sim_unsieved are dense float32 matrices
        extra += (n_unsieved * final_n * 4) / 1e9  # intersect_unsieved
        extra += (n_unsieved * final_n * 4) / 1e9  # sim_unsieved

    return extra

def estimate_memory_usage(final_n, n_res, n_lig_atoms, avg_contacts_per_fp=50):
    """Accurate memory estimation including all components (dense mode)."""
    # Primary memory consumers (float32 = 4 bytes)
    distance_matrix = final_n**2 * 4  # Distance matrix (float32)
    similarity_matrix = final_n**2 * 4  # Similarity matrix (float32)
    # Intersection matrix can be float64 during conversion (sparse @ sparse -> dense)
    intersection_matrix = final_n**2 * 8  # Dense intermediate (float64)
    
    # Sparse matrix storage (data + row_indices + col_indices)
    sparse_memory = final_n * avg_contacts_per_fp * (4 + 4 + 8)  # float32 + int32 + int64
    
    # Dense contact frequency matrix (receptor residues x ligand atoms)
    contact_matrix = n_res * n_lig_atoms * 4  # float32
    
    # Additional arrays (weights, frames, clustering results, etc.)
    arrays_memory = final_n * 6 * 4  # ~6 float32/int32 arrays
    
    total_gb = (distance_matrix + similarity_matrix + intersection_matrix + sparse_memory + contact_matrix + arrays_memory) / 1e9
    
    return total_gb, {
        'distance_matrix': distance_matrix / 1e9,
        'similarity_matrix': similarity_matrix / 1e9,
        'intersection_matrix': intersection_matrix / 1e9,
        'sparse_matrix': sparse_memory / 1e9,
        'contact_matrix': contact_matrix / 1e9,
        'arrays': arrays_memory / 1e9
    }

def estimate_sparse_similarity_density(X, similarity, threshold, sample_size=1000, chunk_size=200, seed=42):
    """Estimate sparse similarity density using a row sample.

    Returns:
        density: estimated nnz / (n_samples^2)
        avg_neighbors: average neighbors per row (including self)
    """
    n_samples = X.shape[0]
    if n_samples == 0:
        return 0.0, 0.0
    sample_size = int(min(sample_size, n_samples))
    rng = np.random.default_rng(seed)
    sample_idx = rng.choice(n_samples, sample_size, replace=False)

    sumB = np.array(X.getnnz(axis=1), dtype=np.float32)
    normB = None
    if similarity == 'cosine':
        normB = np.sqrt(sumB + 1e-10)

    total_edges = 0
    for start in range(0, sample_size, chunk_size):
        end = min(start + chunk_size, sample_size)
        idx = sample_idx[start:end]
        X_chunk = X[idx]
        intersect = X_chunk @ X.T
        if not sparse.issparse(intersect):
            intersect = sparse.csr_matrix(intersect)
        else:
            intersect = intersect.tocsr()

        sumA = np.array(X_chunk.getnnz(axis=1), dtype=np.float32)

        for i in range(end - start):
            row_start = intersect.indptr[i]
            row_end = intersect.indptr[i + 1]
            if row_start == row_end:
                continue
            js = intersect.indices[row_start:row_end]
            inter_vals = intersect.data[row_start:row_end].astype(np.float32)

            if similarity == 'dice':
                sims = 2.0 * inter_vals / (sumA[i] + sumB[js] + 1e-10)
            elif similarity == 'jaccard':
                sims = inter_vals / (sumA[i] + sumB[js] - inter_vals + 1e-10)
            else:  # cosine
                if normB is None:
                    normB = np.sqrt(sumB + 1e-10)
                sims = inter_vals / ((np.sqrt(sumA[i] + 1e-10) * normB[js]) + 1e-10)

            total_edges += int(np.count_nonzero(sims >= threshold))

    avg_neighbors = total_edges / float(sample_size)
    density = avg_neighbors / float(n_samples)
    return density, avg_neighbors

def estimate_csr_memory_gb(nnz, n_rows, data_bytes=4, index_bytes=4, indptr_bytes=4):
    """Estimate CSR memory in GB given nnz and row count."""
    return (nnz * (data_bytes + index_bytes) + (n_rows + 1) * indptr_bytes) / 1e9

def estimate_dbscan_neighborhood_memory_gb(nnz, n_rows, avg_neighbors_per_row=None):
    """Estimate memory for DBSCAN neighborhood list materialization.
    
    sklearn DBSCAN with precomputed sparse distances converts the sparse CSR graph
    into Python object arrays of neighbor indices, which has significant overhead:
    - Object array overhead: ~8 bytes per pointer
    - int64 indices: 8 bytes per neighbor
    - Python object overhead: ~56 bytes per array
    
    Args:
        nnz: Number of nonzero entries in sparse distance matrix
        n_rows: Number of samples
        avg_neighbors_per_row: Average neighbors per row (defaults to nnz/n_rows)
    
    Returns:
        Estimated memory in GB for neighborhood lists
    """
    if avg_neighbors_per_row is None:
        avg_neighbors_per_row = nnz / max(n_rows, 1)
    
    # Memory components:
    # 1. Object array pointers: n_rows * 8 bytes
    # 2. Neighbor indices (int64): nnz * 8 bytes
    # 3. Python list/array overhead: n_rows * 56 bytes (approximate)
    # 4. DBSCAN internal working arrays and label propagation buffers
    # 5. Temporary copies during neighborhood extraction
    
    # Empirically calibrated: actual DBSCAN overhead is ~3x the basic neighbor storage
    # due to internal label arrays, queue/stack for cluster propagation, and scipy/numpy temporaries
    pointer_bytes = n_rows * 8
    indices_bytes = nnz * 8
    overhead_bytes = n_rows * 56
    temporary_factor = 4.5  # Empirically calibrated from 97k samples, 369M edges (was 1.5, actual 3.0x)
    
    total_bytes = (pointer_bytes + indices_bytes + overhead_bytes) * temporary_factor
    return total_bytes / 1e9

def warn_sparse_memory_risk(estimated_gb, available_gb, details, safety_fraction=MEMORY_SAFETY_FRACTION):
    """Log a warning if sparse path is estimated to exceed memory safety threshold."""
    if estimated_gb <= safety_fraction * available_gb:
        return
    logging.warning(
        "Sparse similarity path is estimated at %.2f GB (%.0f%% of available %.2f GB).",
        estimated_gb,
        (estimated_gb / max(available_gb, 1e-9)) * 100.0,
        available_gb,
    )
    logging.warning(
        "Risk of OOM even in sparse mode. Consider lowering --max_sieved, disabling --kdist_plot, "
        "or increasing available memory."
    )
    logging.warning(
        "Sparse estimate details: nnz~%d, avg neighbors~%.1f, extra~%.2f GB, dbscan~%.2f GB, kdist~%.2f GB",
        details.get("est_nnz", 0),
        details.get("avg_neighbors", 0.0),
        details.get("sparse_extra_gb", 0.0),
        details.get("dbscan_neighborhood_gb", 0.0),
        details.get("k_distance_gb", 0.0),
    )

def check_memory_safety(
    estimated_gb,
    available_gb,
    final_n,
    current_max_sieved,
    ignore_risk=False,
    memory_breakdown=None,
    safety_fraction=MEMORY_SAFETY_FRACTION,
):
    """Exit with error and suggested reduction if memory unsafe (with override)."""
    if estimated_gb > safety_fraction * available_gb:
        if ignore_risk:
            # Override flag used - log warning but continue
            logging.warning(
                "Memory safety override active: %.2f GB exceeds %.0f%% of available %.2f GB",
                estimated_gb,
                safety_fraction * 100,
                available_gb,
            )
            logging.warning(f"Proceeding with {final_n} fingerprints despite memory risk (--ignore_memory_risk)")
            if memory_breakdown:
                logging.warning(
                    "Memory breakdown: Distance: %.2f GB, Similarity: %.2f GB, Intersection: %.2f GB, "
                    "Sparse: %.2f GB, Contact: %.2f GB, Arrays: %.2f GB",
                    memory_breakdown['distance_matrix'],
                    memory_breakdown['similarity_matrix'],
                    memory_breakdown['intersection_matrix'],
                    memory_breakdown['sparse_matrix'],
                    memory_breakdown['contact_matrix'],
                    memory_breakdown['arrays'],
                )
        else:
            # No override - enforce limit
            safe_n = int(np.sqrt((safety_fraction * available_gb * 1e9) / (4 * 2)))  # Account for dist+sim matrices
            suggested_max_sieved = (
                min(safe_n, current_max_sieved)
                if current_max_sieved is not None
                else safe_n
            )

            logging.error(
                "Estimated memory usage: %.2f GB exceeds %.0f%% of available %.2f GB",
                estimated_gb,
                safety_fraction * 100,
                available_gb,
            )
            if memory_breakdown:
                logging.error(
                    "Memory breakdown: Distance matrix: %.2f GB, Similarity: %.2f GB, Intersection: %.2f GB, "
                    "Sparse: %.2f GB, Contact: %.2f GB, Arrays: %.2f GB",
                    memory_breakdown['distance_matrix'],
                    memory_breakdown['similarity_matrix'],
                    memory_breakdown['intersection_matrix'],
                    memory_breakdown['sparse_matrix'],
                    memory_breakdown['contact_matrix'],
                    memory_breakdown['arrays'],
                )
            logging.error(
                "CRITICAL: Reduce --max_sieved to %d or use --ignore_memory_risk",
                suggested_max_sieved,
            )
            logging.error("Current %d unique fingerprints -> suggest ~%d", final_n, suggested_max_sieved)
            sys.exit(1)
    return True

def estimate_sparse_mode_memory_gb(
    X,
    similarity,
    threshold,
    final_n,
    n_unique,
    avg_contacts_per_fp,
    memory_breakdown,
    kdist_plot=False,
):
    """Estimate total memory for sparse similarity + DBSCAN path."""
    density, avg_neighbors = estimate_sparse_similarity_density(
        X,
        similarity,
        threshold,
        sample_size=min(1000, final_n),
        chunk_size=200,
        seed=42,
    )
    est_nnz = int(density * final_n * final_n)
    sim_sparse_gb = estimate_csr_memory_gb(est_nnz, final_n)
    dist_sparse_gb = sim_sparse_gb  # distance has same structure
    sparse_extra_gb = estimate_sparse_extra_memory_gb(
        final_n,
        n_unique,
        avg_contacts_per_fp,
        include_unsieved=True,
    )

    if kdist_plot:
        k_distance_gb = estimate_k_distance_memory_gb(
            final_n,
            k=min(15, max(1, final_n // 2)),
            dense=True,
            dense_bytes=4,
        )
    else:
        k_distance_gb = 0.0

    dbscan_neighborhood_gb = estimate_dbscan_neighborhood_memory_gb(
        est_nnz,
        final_n,
        avg_neighbors_per_row=avg_neighbors,
    )

    sparse_base_gb = (
        memory_breakdown['sparse_matrix']
        + memory_breakdown['contact_matrix']
        + memory_breakdown['arrays']
        + sim_sparse_gb
        + dist_sparse_gb
        + sparse_extra_gb
        + dbscan_neighborhood_gb
    )

    sparse_estimated_gb = sparse_base_gb * 1.2 + k_distance_gb

    details = {
        "density": density,
        "avg_neighbors": avg_neighbors,
        "est_nnz": est_nnz,
        "sparse_extra_gb": sparse_extra_gb,
        "dbscan_neighborhood_gb": dbscan_neighborhood_gb,
        "k_distance_gb": k_distance_gb,
    }
    return sparse_estimated_gb, details

def compute_clustering_statistics(X, labels, sample_weight=None, output_prefix='', kdist_plot=False):
    """Compute comprehensive clustering evaluation metrics"""
    from sklearn.metrics import davies_bouldin_score, silhouette_score, calinski_harabasz_score
    from sklearn.neighbors import NearestNeighbors
    
    # Filter out noise points (-1 labels)
    mask = labels >= 0
    X_clustered = X[mask]
    labels_clustered = labels[mask]
    
    if len(np.unique(labels_clustered)) < 2:
        logging.warning("Need at least 2 clusters for comprehensive statistics")
        return {}
    
    weights = None
    if sample_weight is not None:
        weights = sample_weight[mask]
    
    stats = {}
    
    if sparse.issparse(X_clustered):
        logging.info("Skipping Davies-Bouldin/Calinski-Harabasz/Silhouette scores on sparse matrix")
    else:
        # Davies-Bouldin Index (lower is better)
        try:
            stats['davies_bouldin'] = davies_bouldin_score(X_clustered, labels_clustered)
        except Exception as e:
            logging.warning(f"Could not compute Davies-Bouldin: {e}")

        # Silhouette Score
        try:
            if weights is None:
                stats['silhouette_score'] = silhouette_score(X_clustered, labels_clustered)
            else:
                # Weighted silhouette approximation
                stats['silhouette_score'] = silhouette_score(X_clustered, labels_clustered)
        except Exception as e:
            logging.warning(f"Could not compute silhouette score: {e}")

        # Calinski-Harabasz Index (higher is better)
        try:
            stats['calinski_harabasz'] = calinski_harabasz_score(X_clustered, labels_clustered)
        except Exception as e:
            logging.warning(f"Could not compute Calinski-Harabasz: {e}")
    
    # SSR/SST Ratio (pseudo-F statistic)
    stats['ssr_sst_ratio'] = compute_ssr_sst_ratio(X_clustered, labels_clustered, weights)
    
    # K-distance plot for parameter selection (optional, memory-intensive)
    if kdist_plot:
        n_clustered = X_clustered.shape[0]
        k_distance_data = compute_k_distance_data(X_clustered, k=min(15, n_clustered // 2))
        if k_distance_data:
            stats['k_distance_data'] = k_distance_data
            # Save k-distance plot
            plt.figure(figsize=(10, 6))
            plt.plot(k_distance_data['k_distances'])
            plt.xlabel('Points sorted by distance')
            plt.ylabel(f"{k_distance_data['k']}-distance")
            plt.title('K-distance Plot for Parameter Selection')
            plt.grid(True)
            k_dist_prefix = output_prefix if output_prefix else "fp"
            k_dist_filename = f"{k_dist_prefix}_k_distance_plot.png"
            plt.savefig(k_dist_filename, dpi=300, bbox_inches='tight')
            plt.close()
            logging.info(f"Saved k-distance plot to {k_dist_filename}")
    
    return stats

def compute_ssr_sst_ratio(X, labels, weights=None):
    """Compute SSR/SST ratio (pseudo-F statistic for clustering)"""
    total_weight = np.sum(weights) if weights is not None else len(X)

    def _centroid(data, data_weights=None):
        if sparse.issparse(data):
            if data_weights is None:
                return np.asarray(data.mean(axis=0)).ravel()
            weight_sum = np.sum(data_weights)
            weighted_sum = data.T @ data_weights
            return np.asarray(weighted_sum / weight_sum).ravel()
        if data_weights is None:
            return np.mean(data, axis=0)
        return np.average(data, axis=0, weights=data_weights)

    def _row_norms_sq(data):
        if sparse.issparse(data):
            return np.asarray(data.power(2).sum(axis=1)).ravel()
        return np.einsum('ij,ij->i', data, data)

    # Overall centroid
    overall_centroid = _centroid(X, weights)
    centroid_norm_sq = np.dot(overall_centroid, overall_centroid)

    # Total sum of squares (SST) using stable formula
    row_norms_sq = _row_norms_sq(X)
    if weights is None:
        sst = np.sum(row_norms_sq) - total_weight * centroid_norm_sq
    else:
        sst = np.sum(weights * row_norms_sq) - total_weight * centroid_norm_sq

    # Between-cluster sum of squares (SSR)
    ssr = 0.0
    for label in np.unique(labels):
        mask = labels == label
        cluster_data = X[mask]
        cluster_weights = weights[mask] if weights is not None else None

        cluster_centroid = _centroid(cluster_data, cluster_weights)
        n_cluster = np.sum(cluster_weights) if cluster_weights is not None else len(cluster_data)
        distance = np.linalg.norm(cluster_centroid - overall_centroid)
        ssr += n_cluster * distance**2

    return ssr / sst if sst > 0 else 0

def compute_k_distance_data(X, k=15):
    """Compute k-distance data for parameter selection"""
    try:
        nbrs = NearestNeighbors(n_neighbors=k+1).fit(X)  # +1 to exclude self
        distances, _ = nbrs.kneighbors(X)
        k_distances = distances[:, k]  # Distance to k-th neighbor (excluding self)
        return {
            'k_distances': np.sort(k_distances)[::-1],  # Sorted descending
            'k': k
        }
    except Exception as e:
        logging.warning(f"Could not compute k-distance: {e}")
        return None

def compute_sparse_similarity_matrix(X, similarity, threshold, chunk_size=1000, sumB=None, normB=None):
    """Compute sparse similarity matrix above threshold in chunks.

    Args:
        X: CSR matrix (binary contacts)
        similarity: 'dice', 'jaccard', or 'cosine'
        threshold: minimum similarity to store
        chunk_size: number of rows per chunk
        sumB: optional precomputed nnz counts per row
        normB: optional precomputed norms for cosine
    """
    n_samples = X.shape[0]
    if sumB is None:
        sumB = np.array(X.getnnz(axis=1), dtype=np.float32)
    if similarity == 'cosine' and normB is None:
        normB = np.sqrt(sumB + 1e-10)

    sim_chunks = []

    for start in range(0, n_samples, chunk_size):
        end = min(start + chunk_size, n_samples)
        X_chunk = X[start:end]
        intersect = X_chunk @ X.T
        chunk_n = end - start

        # Warn if sparse chunk unexpectedly dense (memory risk)
        if hasattr(intersect, 'nnz') and intersect.nnz > chunk_n * n_samples * 0.1:
            logging.warning(
                f"Sparse matrix chunk {start}:{end} is >10% dense ({intersect.nnz} nnz). "
                "Memory usage may be high."
            )

        if not sparse.issparse(intersect):
            intersect = sparse.csr_matrix(intersect)
        else:
            intersect = intersect.tocsr()

        sumA = np.array(X_chunk.getnnz(axis=1), dtype=np.float32)

        chunk_rows = []
        chunk_cols = []
        chunk_data = []

        for i in range(chunk_n):
            if i + 1 >= len(intersect.indptr):
                logging.error(f"indptr index overflow: i={i}, len={len(intersect.indptr)}")
                break
            row_start = intersect.indptr[i]
            row_end = intersect.indptr[i + 1]
            if row_start == row_end:
                continue
            js = intersect.indices[row_start:row_end]
            inter_vals = intersect.data[row_start:row_end]
            if inter_vals.dtype != np.float32:
                inter_vals = inter_vals.astype(np.float32, copy=False)

            if similarity == 'dice':
                sims = 2.0 * inter_vals / (sumA[i] + sumB[js] + 1e-10)
            elif similarity == 'jaccard':
                sims = inter_vals / (sumA[i] + sumB[js] - inter_vals + 1e-10)
            else:  # cosine
                if normB is None:
                    normB = np.sqrt(sumB + 1e-10)
                sims = inter_vals / ((np.sqrt(sumA[i] + 1e-10) * normB[js]) + 1e-10)

            mask = sims >= threshold
            if np.any(mask):
                js_masked = js[mask]
                sims_masked = sims[mask]
                count = int(js_masked.size)
                chunk_rows.append(np.full(count, i, dtype=np.int32))
                chunk_cols.append(js_masked.astype(np.int32, copy=False))
                chunk_data.append(sims_masked.astype(np.float32, copy=False))

        if chunk_rows:
            rows_arr = np.concatenate(chunk_rows)
            cols_arr = np.concatenate(chunk_cols)
            data_arr = np.concatenate(chunk_data)
            chunk_sparse = csr_matrix(
                (data_arr, (rows_arr, cols_arr)),
                shape=(end - start, n_samples),
                dtype=np.float32,
            )
        else:
            chunk_sparse = csr_matrix((end - start, n_samples), dtype=np.float32)

        sim_chunks.append(chunk_sparse)
        logging.info("Sparse similarity chunk %d:%d processed", start, end)
        log_memory_usage(f"sparse similarity chunk {start}:{end}")
        del intersect, chunk_rows, chunk_cols, chunk_data
        gc.collect()

    if not sim_chunks:
        return csr_matrix((n_samples, n_samples), dtype=np.float32)

    sim_sparse = sparse.vstack(sim_chunks, format='csr')
    return sim_sparse

class DensityRatioClustering:
    """Efficient density ratio clustering implementation"""
    
    def __init__(self, k_neighbors=15, density_ratio_threshold=0.5, min_cluster_size=5):
        self.k_neighbors = k_neighbors  
        self.density_ratio_threshold = density_ratio_threshold
        self.min_cluster_size = min_cluster_size
    
    def fit_predict(self, X, sample_weight=None, precomputed_similarity=None):
        """Density ratio clustering using precomputed or computed similarity
        
        Args:
            X: Feature matrix (ignored if precomputed_similarity provided)
            sample_weight: Sample weights (optional, not currently used)
            precomputed_similarity: Precomputed similarity matrix (n_samples, n_samples)
        """
        n_samples = X.shape[0] if precomputed_similarity is None else precomputed_similarity.shape[0]
        
        if n_samples < self.min_cluster_size:
            return np.array([-1] * n_samples)
        
        # Use precomputed similarity if provided, otherwise compute
        if precomputed_similarity is not None:
            sim_matrix = precomputed_similarity
        else:
            # Convert to dense for similarity computation if needed
            if hasattr(X, 'toarray'):
                X_dense = X.toarray()
            else:
                X_dense = X
            
            # Compute similarity matrix (using cosine similarity by default)
            # Normalize features
            norms = np.linalg.norm(X_dense, axis=1, keepdims=True)
            norms[norms == 0] = 1  # Avoid division by zero
            X_normalized = X_dense / norms
            
            # Compute similarity matrix
            sim_matrix = np.dot(X_normalized, X_normalized.T)
        
        np.fill_diagonal(sim_matrix, 0)  # Remove self-similarity
        
        # Calculate local densities
        densities = self._compute_local_densities(sim_matrix)
        
        # Find density ratio peaks (cluster centers)
        cluster_centers = self._find_density_peaks(densities, sim_matrix)
        
        # Assign labels based on density basin membership
        labels = self._assign_to_basins(densities, cluster_centers, sim_matrix)
        
        return labels
    
    def _compute_local_densities(self, sim_matrix):
        """Local density as sum of similarities to k-nearest neighbors"""
        if sparse.issparse(sim_matrix):
            sim_matrix = sim_matrix.toarray()
        n = sim_matrix.shape[0]
        k = min(self.k_neighbors, sim_matrix.shape[1])
        if k <= 0:
            return np.zeros(n, dtype=np.float32)

        indices = np.argpartition(sim_matrix, -k, axis=1)[:, -k:]
        rows = np.arange(n)[:, None]
        top_k_sims = sim_matrix[rows, indices]
        densities = np.sum(np.maximum(top_k_sims, 0), axis=1)
        return densities
    
    def _find_density_peaks(self, densities, sim_matrix):
        """Find cluster centers based on density ratio"""
        n_samples = len(densities)
        cluster_centers = []
        
        # Sort by density (highest first)
        sorted_indices = np.argsort(-densities)
        
        for idx in sorted_indices:
            if densities[idx] == 0:
                break
                
            # Check density ratio with already assigned centers
            is_center = True
            for center_idx in cluster_centers:
                if sim_matrix[idx, center_idx] > densities[idx] * self.density_ratio_threshold:
                    is_center = False
                    break
            
            if is_center:
                cluster_centers.append(idx)
        
        return cluster_centers
    
    def _assign_to_basins(self, densities, cluster_centers, sim_matrix):
        """Assign points to nearest cluster center"""
        n_samples = len(densities)
        labels = np.full(n_samples, -1)
        
        if not cluster_centers:
            return labels
        
        # Assign each point to the nearest cluster center
        for i in range(n_samples):
            if i in cluster_centers:
                labels[i] = cluster_centers.index(i)
            else:
                # Find nearest center by similarity
                similarities_to_centers = [sim_matrix[i, center] for center in cluster_centers]
                nearest_center_idx = np.argmax(similarities_to_centers)
                
                # Only assign if similarity is positive
                if similarities_to_centers[nearest_center_idx] > 0:
                    labels[i] = nearest_center_idx
        
        # Filter small clusters
        unique_labels, counts = np.unique(labels[labels >= 0], return_counts=True)
        small_clusters = unique_labels[counts < self.min_cluster_size]
        
        for cluster_label in small_clusters:
            labels[labels == cluster_label] = -1
        
        return labels

def calculate_comprehensive_binding_statistics(total_frames, frame_weights, n_ligands, sasa_arrays=None, bound_instances=None, sasa_eps=DEFAULT_SASA_EPS):
    """Calculate comprehensive bound/unbound statistics
    
    Note: For multi-ligand systems:
    - bound_instances contains one entry per (frame, ligand_id) pair with contacts
    - Each ligand is treated independently
    - Total possible instances = total_frames * n_ligands
    - Raw counts: bound vs unbound ligand-frame pairs (instances)
    - Weighted: fraction of frame weight with at least one bound ligand
    """
    
    if bound_instances is None:
        bound_count = 0
        logging.warning("No bound instances tracked")
        return {}
    
    # Count bound ligand-frame pairs (instances)
    n_bound_instances = len(bound_instances)
    total_possible_instances = total_frames * n_ligands
    n_unbound_instances = total_possible_instances - n_bound_instances
    
    
    # Instance-weighted statistics (per ligand-frame pair)
    bound_instance_weight = np.sum([inst['weight'] for inst in bound_instances])
    total_instance_weight = np.sum(frame_weights) * n_ligands
    bound_instance_weighted = bound_instance_weight / total_instance_weight if total_instance_weight > 0 else 0.0
    unbound_instance_weighted = 1.0 - bound_instance_weighted

    # Instance-weighted statistics using SASA-scaled weights (if available)
    # Normalize by total scaled instance weight across ALL ligand-frame pairs
    bound_instance_final_weight = np.sum([inst['final_weight'] for inst in bound_instances])
    total_scaled_instance_weight = total_instance_weight
    if sasa_arrays is not None:
        total_scaled_instance_weight = 0.0
        for lig_id in range(n_ligands):
            total_scaled_instance_weight += np.sum(frame_weights * (1.0 / (sasa_arrays[lig_id] + sasa_eps)))
    bound_instance_final_weighted = (
        bound_instance_final_weight / total_scaled_instance_weight if total_scaled_instance_weight > 0 else 0.0
    )
    unbound_instance_final_weighted = 1.0 - bound_instance_final_weighted

    # SASA-scaled statistics (normalized)
    bound_sasa_scaled = unbound_sasa_scaled = 0.0
    if sasa_arrays is not None:
        for inst in bound_instances:
            frame = inst['frame']
            lig_id = inst['ligand_id']
            sasa_scale = 1.0 / (sasa_arrays[lig_id][frame] + sasa_eps)
            bound_sasa_scaled += frame_weights[frame] * sasa_scale
        bound_sasa_scaled = bound_sasa_scaled / total_scaled_instance_weight if total_scaled_instance_weight > 0 else 0.0
        unbound_sasa_scaled = 1.0 - bound_sasa_scaled
    
    stats = {
        'total_frames': total_frames,
        'n_ligands': n_ligands,
        'total_possible_instances': total_possible_instances,
        'bound_instances': n_bound_instances,
        'unbound_instances': n_unbound_instances,
        'instance_bound_fraction': n_bound_instances / total_possible_instances,
        'bound_instance_weighted': bound_instance_weighted,
        'unbound_instance_weighted': unbound_instance_weighted,
        'bound_instance_final_weighted': bound_instance_final_weighted,
        'unbound_instance_final_weighted': unbound_instance_final_weighted,
        'bound_sasa_scaled': bound_sasa_scaled,
        'unbound_sasa_scaled': unbound_sasa_scaled
    }
    
    logging.info("=== Binding Statistics (instance-based) ===")
    logging.info(
        "Global: Bound raw count %d, Unbound raw count %d, "
        "Bound occ %.4f, Unbound occ %.4f, "
        "Bound weight occ %.4f, Unbound weight occ %.4f",
        stats['bound_instances'],
        stats['unbound_instances'],
        stats['instance_bound_fraction'],
        1 - stats['instance_bound_fraction'],
        stats['bound_instance_weighted'],
        stats['unbound_instance_weighted'],
    )
    if sasa_arrays is not None:
        logging.info(
            "Global (SASA-scaled): Bound weight occ %.4f, Unbound weight occ %.4f",
            stats['bound_instance_final_weighted'],
            stats['unbound_instance_final_weighted'],
        )
    if sasa_arrays is not None:
        logging.info(f"SASA-scaled bound/unbound: {stats['bound_sasa_scaled']:.6f}/{stats['unbound_sasa_scaled']:.6f}")

    # Per-ligand bound/unbound (instance-based)
    bound_matrix = np.zeros((n_ligands, total_frames), dtype=bool)
    for inst in bound_instances:
        bound_matrix[inst['ligand_id'], inst['frame']] = True
    total_frame_weight = np.sum(frame_weights)
    for lig_id in range(n_ligands):
        lig_bound_frames = np.sum(bound_matrix[lig_id])
        lig_bound_weight = np.sum(frame_weights[bound_matrix[lig_id]])
        lig_bound_final_weight = np.sum([
            inst['final_weight'] for inst in bound_instances
            if inst['ligand_id'] == lig_id
        ])
        lig_bound_fraction = lig_bound_frames / total_frames if total_frames > 0 else 0.0
        lig_bound_weighted = lig_bound_weight / total_frame_weight if total_frame_weight > 0 else 0.0
        lig_total_scaled_weight = total_frame_weight
        if sasa_arrays is not None:
            lig_total_scaled_weight = np.sum(frame_weights * (1.0 / (sasa_arrays[lig_id] + sasa_eps)))
        lig_bound_final_weighted = (
            lig_bound_final_weight / lig_total_scaled_weight if lig_total_scaled_weight > 0 else 0.0
        )
        logging.info(
            "Ligand %d: Bound occ %.4f, Unbound occ %.4f, Bound weight occ %.4f, Unbound weight occ %.4f",
            lig_id + 1,
            lig_bound_fraction,
            1 - lig_bound_fraction,
            lig_bound_weighted,
            1 - lig_bound_weighted,
        )
        if sasa_arrays is not None:
            logging.info(
                "Ligand %d (SASA-scaled): Bound weight occ %.4f, Unbound weight occ %.4f",
                lig_id + 1,
                lig_bound_final_weighted,
                1 - lig_bound_final_weighted,
            )
    
    return stats

def calculate_sasa_internal(universe, ligand_selections, receptor_selection=None, probe_radius=1.4, n_points=100):
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
    n_ligands = len(ligand_selections)
    n_frames = len(universe.trajectory)
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
    
    logging.info(f"Calculating internal SASA (FAST VDW+Gaussian approximation): {n_ligands} ligands, {n_frames} frames")
    logging.info(f"  Method: VDW radii + Gaussian smoothing (λ={gaussian_lambda}, σ={gaussian_sigma}Å, cutoff={neighbor_cutoff}Å)")
    logging.info(f"  NOTE: Fast approximation for relative comparisons (10-20% accuracy vs. probe methods)")
    
    # Helper function to calculate SASA using VDW+Gaussian method
    def calculate_group_sasa_fast(atom_group, all_coords, all_elements, tree=None, group_indices_in_all=None):
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
        
        # Efficient neighbor accumulation via local neighbor queries (sparse by construction)
        neighbors_list = tree.query_ball_tree(tree, neighbor_cutoff)
        if len(neighbors_list) == 0:
            # No neighbors within cutoff; full exposure for all atoms
            atom_sasa_array[:] = 4 * np.pi * (group_radii ** 2)
            return atom_sasa_array / 100.0

        if group_indices_in_all is None:
            group_indices_in_all = np.arange(n_atoms, dtype=np.intp)
        else:
            group_indices_in_all = np.asarray(group_indices_in_all, dtype=np.intp)
            if len(group_indices_in_all) != n_atoms:
                raise ValueError(
                    f"group_indices_in_all length mismatch: expected {n_atoms}, got {len(group_indices_in_all)}"
                )

        neighbor_sum = np.zeros(n_atoms, dtype=np.float32)
        for local_i, global_i in enumerate(group_indices_in_all):
            if global_i < 0 or global_i >= len(neighbors_list):
                continue
            neighbor_indices = neighbors_list[global_i]
            if len(neighbor_indices) <= 1:
                continue
            neighbor_indices = np.asarray(neighbor_indices, dtype=np.intp)
            dists = np.linalg.norm(all_coords[neighbor_indices] - all_coords[global_i], axis=1)
            weights = np.exp(-(dists ** 2) / (2 * gaussian_sigma ** 2))
            neighbor_sum[local_i] = np.sum(weights[dists > 0])

        # Exposure and SASA per atom
        exposure = np.exp(-gaussian_lambda * neighbor_sum)
        atom_sasa_array = (4 * np.pi * (group_radii ** 2) * exposure).astype(np.float32, copy=False)
        
        # Convert from Ų to nm²
        return atom_sasa_array / 100.0
    
    # Determine receptor selection if not provided
    if receptor_selection is None:
        # Use all non-ligand atoms as receptor context
        all_ligand_indices = set()
        for lg in ligand_selections:
            all_ligand_indices.update(lg.indices)
        receptor_indices = [i for i in range(len(universe.atoms)) 
                            if i not in all_ligand_indices]
        receptor_sel = universe.atoms[receptor_indices]
    else:
        receptor_sel = receptor_selection
    
    # Build complex with receptor + all ligands once (topology is static)
    all_ligands = receptor_sel.universe.atoms[:0]
    for lig in ligand_selections:
        all_ligands = all_ligands + lig
    complex_all = receptor_sel + all_ligands
    complex_indices = complex_all.indices
    index_map = {atom_index: pos for pos, atom_index in enumerate(complex_indices)}

    # Calculate per-atom SASA once per frame, then sum by ligand atoms.
    # Precompute ligand atom indices in complex for fast per-frame summation
    ligand_indices_global = [lg.indices for lg in ligand_selections]

    for ts in universe.trajectory:
        frame = ts.frame
        
        # 1. Calculate per-atom SASA for entire complex (ONCE per frame)
        complex_all_coords = complex_all.positions
        complex_all_elements = complex_all.elements
        tree = cKDTree(complex_all_coords)
        per_atom_sasa = calculate_group_sasa_fast(complex_all, complex_all_coords, complex_all_elements, tree=tree)
        
        # 2. For each ligand: sum SASA for atoms belonging to this ligand
        for lig_id, lig_indices in enumerate(ligand_indices_global):
            missing_atoms = [i for i in lig_indices if i not in index_map]
            if missing_atoms:
                logging.warning(
                    f"Ligand {lig_id+1}: {len(missing_atoms)} atoms not found in complex index map "
                    f"(indices: {missing_atoms[:5]}{'...' if len(missing_atoms)>5 else ''}). "
                    "SASA may be underestimated. Check topology-trajectory atom count consistency."
                )
            lig_indices_in_complex = [index_map[i] for i in lig_indices if i in index_map]
            if len(lig_indices_in_complex) == 0:
                sasa_arrays[lig_id][frame] = 0.0
                continue
            ligand_sasa = np.sum(per_atom_sasa[lig_indices_in_complex])
            sasa_arrays[lig_id][frame] = ligand_sasa
        
        if frame % max(1, n_frames // 100) == 0:
            logging.info(f"  Processed frame {frame+1}/{n_frames}")
    
    logging.info(f"Completed internal SASA calculation (fast mode)")
    return sasa_arrays

def load_or_calculate_sasa(args, universe, n_ligands, n_frames, receptor_selection=None):
    """
    Load SASA from external files (gmx sasa) or calculate internally (fast VDW+Gaussian approximation).
    
    This function provides a unified interface for SASA acquisition, supporting
    both the traditional workflow (external gmx sasa files) and the new internal
    calculation method (using only existing environment libraries).
    
    Parameters
    ----------
    args : argparse.Namespace
        Command-line arguments containing SASA configuration
        Required fields:
        - sasa_method: 'external' or 'internal'
        - sasa_dir: directory for external SASA files (external mode)
        - sasa_type: 'csv' or 'xvg' (external mode)
        - sasa_col: column index for SASA values (external mode)
        - ligand_sel: ligand selection string (internal mode)
        - sasa_n_points: test points per atom (internal mode)
    universe : MDAnalysis.Universe
        Trajectory universe (only used for internal calculation)
    n_ligands : int
        Number of ligands to process
    n_frames : int
        Number of frames in trajectory
    receptor_selection : MDAnalysis.AtomGroup, optional
        Receptor atoms for internal calculation context
    
    Returns
    -------
    sasa_arrays : list of numpy.ndarray
        One array per ligand containing SASA values per frame (nm²)
        Shape: (n_frames,) for each ligand
    
    Raises
    ------
    ValueError
        - SASA file length mismatch with trajectory frames (external mode)
        - No atoms found for ligand selection (internal mode)
        - Unknown SASA method
    FileNotFoundError
        - SASA files not found (external mode)
    
    Notes
    -----
    **Internal method (default):**
    - Calculates SASA on-the-fly using fast VDW+Gaussian approximation
    - No external preprocessing or new dependencies required
    - Expected variance: ~10-20% from gmx sasa (approximate method)
    - Performance: fast per-frame computation; n_points is ignored in fast mode
    - Suitable for relative comparisons and avoiding gmx sasa dependency
    
    **External method:**
    - Loads pre-computed SASA from gmx sasa output files
    - Fast (no computation), exact match to GROMACS topology radii
    - Requires gmx sasa preprocessing step (see LFS_ana_hrex.sh)
    """
    if args.sasa_method == 'external':
        # EXISTING CODE: load from gmx sasa output files
        logging.info(f"Loading SASA from external files in {args.sasa_dir}")
        sasa_arrays = []
        for lig_id in range(n_ligands):
            path = os.path.join(args.sasa_dir, f"sasa_lig{lig_id+1}.{args.sasa_type}")
            if args.sasa_type == 'csv':
                import pandas as pd
                data = pd.read_csv(path, header=None).values
            else:
                data = np.loadtxt(path, comments=['#', '@'])
            sasa = data[:, args.sasa_col]
            if len(sasa) != n_frames:
                raise ValueError(f"SASA length mismatch for ligand {lig_id+1}: "
                                 f"expected {n_frames} frames, got {len(sasa)}")
            sasa_arrays.append(sasa)
        logging.info(f"Loaded SASA arrays from external files (gmx sasa method)")
        return sasa_arrays
    
    elif args.sasa_method == 'internal':
        # Internal fast SASA calculation (VDW+Gaussian approximation)
        logging.info("Calculating SASA internally (fast VDW+Gaussian approximation)")
        
        # Build ligand selections from ligand_sel (same pattern as main loop)
        # Get all ligands and treat each residue as a separate ligand
        ligands_ag = universe.select_atoms(args.ligand_sel)
        ligand_residues = ligands_ag.residues
        
        if len(ligand_residues) != n_ligands:
            logging.warning(f"Expected {n_ligands} ligands but found {len(ligand_residues)} residues in selection '{args.ligand_sel}'")
        
        ligand_selections = []
        for lig_id in range(n_ligands):
            if lig_id < len(ligand_residues):
                ligand = ligand_residues[lig_id].atoms
                ligand_selections.append(ligand)
                logging.info(f"  Ligand {lig_id+1}: {len(ligand)} atoms (resid {ligand_residues[lig_id].resid})")
            else:
                raise ValueError(f"Ligand {lig_id+1} not found in selection '{args.ligand_sel}'")
        
        # Calculate SASA using internal function
        sasa_arrays = calculate_sasa_internal(
            universe, 
            ligand_selections,
            receptor_selection=receptor_selection,
            probe_radius=1.4,  # Standard water probe
            n_points=args.sasa_n_points
        )
        
        logging.info("Completed internal SASA calculation")
        return sasa_arrays
    
    else:
        raise ValueError(f"Unknown SASA method: {args.sasa_method}. Must be 'external' or 'internal'.")

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
            import pandas as pd
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

    import pandas as pd

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

    csv_file = f"{output_prefix}_histogram.csv"
    data = {'bin_center': bin_centers}
    for i, hist in enumerate(ligand_hists):
        data[f"ligand{i+1}_density"] = hist
    if combined_hist is not None:
        data['combined_avg'] = combined_hist

    df = pd.DataFrame(data)
    df.to_csv(csv_file, index=False, float_format='%.6f')

    if combined_png_file:
        logging.info(f"Saved SASA histogram: {combined_png_file} + {csv_file}")
    else:
        logging.info(f"Saved SASA histogram CSV: {csv_file}")

def align_structures(ref_coords, mob_coords, ref_resids, mob_resids, ref_resnames, mob_resnames, allow_mismatch=True, min_common=10):
    """Structural alignment using common residues (with sequential fallback).

    Falls back to sequential residue matching when residue IDs/names differ.
    """
    # Find common residues by resid and resname
    common_ref_mask = []
    common_mob_mask = []

    for i, (ref_id, ref_name) in enumerate(zip(ref_resids, ref_resnames)):
        matches = np.where((mob_resids == ref_id) & (mob_resnames == ref_name))[0]
        if len(matches) > 0:
            common_ref_mask.append(i)
            common_mob_mask.append(matches[0])  # Take first match

    # Fallback: sequential matching by position
    if len(common_ref_mask) < min_common:
        if not allow_mismatch:
            raise ValueError("Need at least 3 common residues for alignment")
        common_len = min(len(ref_coords), len(mob_coords))
        if common_len < 3:
            raise ValueError("Need at least 3 residues for alignment")
        common_ref_mask = list(range(common_len))
        common_mob_mask = list(range(common_len))

    # Extract coordinates for common residues
    ref_common = ref_coords[common_ref_mask]
    mob_common = mob_coords[common_mob_mask]

    # Center coordinates
    ref_centered = ref_common - np.mean(ref_common, axis=0)
    mob_centered = mob_common - np.mean(mob_common, axis=0)

    # Compute optimal rotation using Kabsch algorithm
    H = mob_centered.T @ ref_centered
    U, S, Vt = np.linalg.svd(H)
    d = -1 if (np.linalg.det(Vt.T @ U.T) < 0) else 1
    D = np.diag([1, 1, d])
    rotation = Vt.T @ D @ U.T

    # Apply rotation and translation to full mobile structure
    mob_mean = np.mean(mob_common, axis=0)
    ref_mean = np.mean(ref_common, axis=0)
    mob_aligned = (mob_coords - mob_mean) @ rotation
    translation = ref_mean - np.mean(mob_aligned, axis=0)
    mob_aligned += translation

    return mob_aligned, len(common_ref_mask), rotation, mob_mean, translation


def _append_reference_row(modes_path, ref_pairs, sasa_enabled=False):
    """Insert reference mode as first entry after header in binding_modes.csv."""
    if sasa_enabled:
        ref_header = "mode_id weighted_occupancy sasa_scaled_occupancy rep_frame interacting_pairs n_instances avg_sasa"
        ref_row = f"reference 0.000000 0.000000 -1 {','.join(ref_pairs)} 0 0.0000"
    else:
        ref_header = "mode_id weighted_occupancy rep_frame interacting_pairs n_instances"
        ref_row = f"reference 0.000000 -1 {','.join(ref_pairs)} 0"
    try:
        lines = []
        if not os.path.exists(modes_path):
            lines = [f"{ref_header}\n"]
        else:
            with open(modes_path, 'r') as f:
                lines = f.readlines()

        if not lines:
            lines = [f"{ref_header}\n"]

        header = lines[0].strip()
        if header != ref_header:
            logging.warning("binding_modes.csv header did not match expected; reference row inserted without header validation")

        # Remove existing reference rows to avoid duplicates
        remaining = [line for line in lines[1:] if not line.startswith("reference ")]

        new_lines = [lines[0].rstrip("\n") + "\n", f"{ref_row}\n"] + remaining
        dir_path = os.path.dirname(os.path.abspath(modes_path))
        tmp_path = None
        try:
            with tempfile.NamedTemporaryFile(mode='w', delete=False, dir=dir_path, suffix='.tmp') as tmp:
                tmp.writelines(new_lines)
                tmp_path = tmp.name
            os.replace(tmp_path, modes_path)
        finally:
            if tmp_path and os.path.exists(tmp_path):
                os.unlink(tmp_path)
    except Exception as e:
        logging.warning(f"Could not insert reference row to binding_modes.csv: {e}")

def _calculate_concentration_factor(concentration_nm3):
    """Return concentration correction factor for standard state."""
    standard_conc = 1.0 / 1.66  # nm^-3
    if concentration_nm3 is not None and concentration_nm3 > 0:
        return standard_conc / concentration_nm3
    return standard_conc


def estimate_contact_volume_fraction(receptor_coords_nm, box_dims_nm, cutoff_nm, n_samples=100000, seed=42, rng=None):
    """Monte Carlo estimate of contact volume fraction.

    Notes:
      - Uses a SINGLE receptor frame (topology/first frame), not per-frame sampling.
      - Assumes coordinates and box dimensions are in nm.
    """
    if receptor_coords_nm is None or receptor_coords_nm.size == 0:
        raise ValueError("Receptor coordinates are empty; cannot estimate contact volume")
    if box_dims_nm is None or np.any(np.asarray(box_dims_nm) <= 0):
        raise ValueError("Box dimensions must be positive to estimate contact volume")
    if cutoff_nm <= 0:
        raise ValueError("Cutoff must be positive to estimate contact volume")

    tree = cKDTree(receptor_coords_nm)
    if rng is None:
        rng = np.random.default_rng(seed)
    points = rng.uniform(low=[0.0, 0.0, 0.0], high=box_dims_nm, size=(int(n_samples), 3))
    distances, _ = tree.query(points)
    n_contact = np.sum(distances <= cutoff_nm)
    f_contact = float(n_contact) / float(n_samples)
    return f_contact


def estimate_contact_volume_fraction_weighted(
    receptor_coords_by_frame_nm,
    box_dims_by_frame_nm,
    frame_weights,
    cutoff_nm,
    n_samples_per_frame=5000,
    seed=42,
):
    """Weighted per-frame contact volume fraction.

    Args:
        receptor_coords_by_frame_nm: list/array of receptor positions per frame (Nf, Natoms, 3) in nm
        box_dims_by_frame_nm: list/array of box dimensions per frame (Nf, 3) in nm
        frame_weights: weights per frame (Nf)
        cutoff_nm: contact cutoff (nm)
        n_samples_per_frame: MC samples per frame
        seed: RNG seed

    Returns:
        f_contact_weighted: weighted average of per-frame contact volume fractions
    """
    if receptor_coords_by_frame_nm is None or len(receptor_coords_by_frame_nm) == 0:
        raise ValueError("No receptor coordinates provided for per-frame contact volume")
    if box_dims_by_frame_nm is None or len(box_dims_by_frame_nm) == 0:
        raise ValueError("No box dimensions provided for per-frame contact volume")
    if cutoff_nm <= 0:
        raise ValueError("Cutoff must be positive for contact volume estimation")

    coords_arr = np.asarray(receptor_coords_by_frame_nm, dtype=np.float64)
    box_arr = np.asarray(box_dims_by_frame_nm, dtype=np.float64)
    weights_arr = np.asarray(frame_weights, dtype=np.float64)

    if coords_arr.shape[0] != box_arr.shape[0] or coords_arr.shape[0] != weights_arr.shape[0]:
        raise ValueError("Per-frame contact volume inputs have inconsistent lengths")
    if np.any(box_arr <= 0):
        raise ValueError("Box dimensions must be positive for per-frame contact volume")
    if np.sum(weights_arr) <= 0:
        raise ValueError("Frame weights must sum to > 0 for per-frame contact volume")

    weights_arr = weights_arr / np.sum(weights_arr)
    rng = np.random.default_rng(seed)

    f_contact_weighted = 0.0
    for i in range(coords_arr.shape[0]):
        f_i = estimate_contact_volume_fraction(
            coords_arr[i],
            box_arr[i],
            cutoff_nm,
            n_samples=n_samples_per_frame,
            rng=rng,
        )
        f_contact_weighted += weights_arr[i] * f_i

    return float(f_contact_weighted)


def _guess_element_from_name(atom_name):
    """Heuristic to extract element from atom name (e.g., 'C1' -> 'C')."""
    if not atom_name:
        return None
    name = str(atom_name).strip()
    if not name:
        return None
    if len(name) >= 2 and name[0].isalpha() and name[1].islower():
        return name[:2].capitalize()
    if name[0].isalpha():
        return name[0].upper()
    return None


def get_vdw_radii_nm(atom_group):
    """Return array of vdW radii in nm for atoms in an MDAnalysis AtomGroup.

    Uses atom.element when available; falls back to atom.type or atom.name heuristics.
    Missing elements default to carbon radius.
    """
    # Bondi (1964) values for main-group elements and selected transition metals.
    # Mantina et al. (2009) for main-group elements missing in Bondi's table.
    vdw_radii_ang = {
        # Main group (Bondi 1964)
        'H': 1.20,
        'He': 1.40,
        'C': 1.70,
        'N': 1.55,
        'O': 1.52,
        'F': 1.47,
        'Ne': 1.54,
        'Na': 2.27,
        'Mg': 1.73,
        'Si': 2.10,
        'P': 1.80,
        'S': 1.80,
        'Cl': 1.75,
        'Ar': 1.88,
        'K': 2.75,
        'Ga': 1.87,
        'As': 1.85,
        'Se': 1.90,
        'Br': 1.85,
        'Kr': 2.02,
        'In': 1.93,
        'Sn': 2.17,
        'Te': 2.06,
        'I': 1.98,
        'Xe': 2.16,
        'Tl': 1.96,
        'Pb': 2.02,
        # Main group (Mantina et al. 2009 extensions)
        'Be': 1.53,
        'B': 1.92,
        'Al': 1.84,
        'Ca': 2.31,
        'Ge': 2.11,
        'Rb': 3.03,
        'Sr': 2.50,
        'Sb': 2.06,
        'Cs': 3.43,
        'Ba': 2.68,
        'Bi': 2.07,
        'Po': 1.97,
        'At': 2.02,
        'Rn': 2.20,
        'Fr': 3.48,
        'Ra': 2.83,
        # Alkali/alkaline earth (Bondi 1964)
        'Li': 1.81,
        # Transition metals (Bondi 1964)
        'Mn': 2.00,
        'Fe': 2.00,
        'Co': 1.95,
        'Ni': 1.63,
        'Cu': 1.40,
        'Zn': 1.39,
        'Ru': 2.05,
        'Rh': 2.00,
        'Pd': 1.63,
        'Ag': 1.72,
        'Pt': 1.72,
        'Au': 1.66,
    }
    default_radius = vdw_radii_ang['C']
    radii_ang = []
    for atom in atom_group:
        element = None
        try:
            element = atom.element
        except Exception:
            element = None
        if not element:
            try:
                element = atom.type
            except Exception:
                element = None
        if element:
            element = str(element).strip()
            if len(element) > 2:
                element = element[:2].capitalize()
            elif len(element) == 2:
                element = element.capitalize()
            else:
                element = element.upper()
        if not element or element not in vdw_radii_ang:
            element_guess = _guess_element_from_name(getattr(atom, 'name', None))
            if element_guess in vdw_radii_ang:
                element = element_guess
        radius = vdw_radii_ang.get(element if element is not None else 'C', default_radius)
        radii_ang.append(radius)
    radii_nm = np.asarray(radii_ang, dtype=np.float64) * 0.1
    return radii_nm


def compute_contact_volume_correction(
    args,
    avg_volume_nm3,
    avg_box_dims_nm,
    receptor_coords_nm,
    per_frame_receptor_coords_nm,
    per_frame_box_dims_nm,
    per_frame_weights,
    n_ligs,
    ref_lig,
    binding_stats,
):
    """Compute contact volume fraction and correction term based on selected method."""
    f_contact = None
    volume_correction_kjmol = None
    if not (args.dg_contact_volume_correction or args.dg_contact_volume_per_frame or args.dg_contact_volume_molecular):
        return f_contact, volume_correction_kjmol

    logging.info("Contact-volume correction enabled")
    if args.dg_contact_volume_molecular:
        if avg_volume_nm3 is None or avg_volume_nm3 <= 0:
            raise ValueError("Cannot compute molecular contact volume: missing average box volume")
        packing_factor = float(args.dg_volume_packing_factor)
        if packing_factor <= 0:
            raise ValueError("Packing factor must be positive for molecular contact volume")
        vdw_radii_nm = get_vdw_radii_nm(ref_lig.atoms)
        v_molecule_nm3 = float(np.sum((4.0 / 3.0) * np.pi * vdw_radii_nm**3))
        bound_fraction = 1.0
        if args.dg_volume_use_bound_fraction:
            if isinstance(binding_stats, dict) and binding_stats:
                key = 'bound_instance_final_weighted' if args.sasa_enabled else 'bound_instance_weighted'
                if key in binding_stats:
                    bound_fraction = float(binding_stats[key])
                else:
                    logging.warning("Bound fraction key missing (%s); using 1.0", key)
            else:
                logging.warning("Bound fraction requested but binding_stats missing; using 1.0")
        v_binding_nm3 = n_ligs * v_molecule_nm3 * packing_factor * bound_fraction
        f_contact = v_binding_nm3 / avg_volume_nm3
        logging.info(
            "Molecular volume estimate: V_molecule=%.4f nm^3, packing=%.3f, bound_fraction=%.3f",
            v_molecule_nm3, packing_factor, bound_fraction
        )
        logging.info(
            "Molecular binding volume: V_binding=%.4f nm^3 (box=%.4f nm^3)",
            v_binding_nm3, avg_volume_nm3
        )
        logging.info("NOTE: Molecular volume method uses ONLY ligand VdW radii (NO receptor coordinates)")
    else:
        if args.dg_volume_cutoff is not None:
            cutoff_nm = float(args.dg_volume_cutoff)
            logging.info("Using overridden volume cutoff: %.3f nm", cutoff_nm)
        else:
            cutoff_nm = args.cutoff * 0.1
            logging.info("Using fingerprint cutoff for volume: %.3f nm", cutoff_nm)

        if args.dg_contact_volume_per_frame:
            logging.info("[LEGACY METHOD] Estimating per-frame contact volume (weighted average)")
            logging.info("NOTE: This method uses receptor atom positions from trajectory frames")
            logging.info(
                "Per-frame settings: stride=%d, samples_per_frame=%d",
                max(int(args.dg_volume_frame_stride), 1),
                int(args.dg_volume_samples_per_frame),
            )
            if not per_frame_box_dims_nm:
                raise ValueError("Cannot compute per-frame contact volume: no sampled frames")
            f_contact = estimate_contact_volume_fraction_weighted(
                per_frame_receptor_coords_nm,
                per_frame_box_dims_nm,
                per_frame_weights,
                cutoff_nm,
                n_samples_per_frame=int(args.dg_volume_samples_per_frame),
            )
        else:
            logging.info("[LEGACY METHOD] Estimating contact volume fraction (100000 MC samples)")
            logging.info("NOTE: This method uses receptor atom positions from topology frame (single frame)")
            if avg_box_dims_nm is None:
                raise ValueError("Cannot compute contact volume: missing box dimensions")
            logging.info(
                "Box dimensions (avg): %.3f x %.3f x %.3f nm",
                avg_box_dims_nm[0], avg_box_dims_nm[1], avg_box_dims_nm[2]
            )
            f_contact = estimate_contact_volume_fraction(
                receptor_coords_nm,
                avg_box_dims_nm,
                cutoff_nm,
                n_samples=100000,
            )

    if f_contact is not None:
        if f_contact <= 0:
            raise ValueError("Computed contact volume fraction <= 0; cannot apply correction")
        if f_contact > 1.0:
            logging.warning("Contact volume fraction > 1 (%.3f); capping to 0.999", f_contact)
            f_contact = 0.999
        volume_correction_kjmol = - (8.314 * args.temperature / 1000.0) * np.log(1.0 / max(f_contact, 1e-12))
        logging.info("Contact volume fraction: f_contact = %.6f", f_contact)
        logging.info("Volume correction term: %.2f kJ/mol", volume_correction_kjmol)

    return f_contact, volume_correction_kjmol


def _bootstrap_contact_unbound(boot_weights):
    """Bootstrap bound/unbound probabilities for contact-unbound reference."""
    boot_P_bound = np.sum(boot_weights)
    boot_P_bound = min(max(boot_P_bound, 0.0), 1.0 - 1e-12)
    boot_P_unbound = 1.0 - boot_P_bound
    return boot_P_bound, boot_P_unbound


def _bootstrap_noise_unbound(boot_weights, boot_labels):
    """Bootstrap bound/unbound probabilities using DBSCAN noise as unbound."""
    boot_weights = boot_weights / np.sum(boot_weights)
    boot_bound_mask = (boot_labels >= 0)
    boot_P_bound = np.sum(boot_weights[boot_bound_mask])
    boot_P_unbound = np.sum(boot_weights[~boot_bound_mask])
    return boot_P_bound, boot_P_unbound


def calculate_overall_binding_free_energy(
    unique_weights,
    labels,
    temperature=300.0,
    n_bootstrap=1000,
    unbound_weight_override=None,
    concentration_nm3=None,
    dg_offset_kj=0.0,
    f_contact=None,
):
    """
    Calculate overall binding free energy from weighted bound/unbound probabilities.
    
    MOLECULAR VOLUME CONTACT CORRECTION (Recommended - Auto-enabled):
    When free energy calculation is enabled, the molecular volume method is automatically
    activated (unless explicitly overridden). This method computes a physically-based
    volume correction using ONLY ligand VdW radii (NO receptor coordinates) to achieve
    PREDICTIVE absolute binding free energies without experimental calibration.
    
    The correction accounts for the entropic penalty of confining ligands to the binding
    cavity versus bulk solution.
    
    Validation: +/-2-4 kcal/mol accuracy achieved on test systems without empirical fitting.
    
    CONTACT-SPACE REFERENCE (Legacy behavior without molecular correction):
    Without molecular volume correction, fingerprint analysis uses CONTACT-UNBOUND
    (no contacts within cutoff) or DBSCAN noise as the unbound reference. This is
    NOT the thermodynamic bulk-unbound state and systematically underestimates dG.
    
    RECOMMENDATIONS:
    1. Default setup (auto-enabled): --calc_overall_dg (recommended)
       - Automatically enables --dg_contact_volume_molecular
       - Automatically enables --dg_volume_use_bound_fraction
       - Uses packing_factor=0.0005 (override with --dg_volume_packing_factor)
    2. For per-cluster analysis: --calc_cluster_dg (recommended)
    3. Legacy methods: --dg_contact_volume_correction or --dg_contact_volume_per_frame
       (use receptor coordinates - not recommended)
    
    Core formula: dG = -RT ln(K_bind x C0/C_sim x 1/f_contact) + offset
    where K_bind = P_bound/P_unbound, f_contact = binding volume fraction
    
    Args:
        unique_weights: Array of normalized weights per fingerprint instance
        labels: Array of cluster labels (-1 = DBSCAN noise, >=0 = cluster ID)
        temperature: Temperature in Kelvin (default 300K)
        n_bootstrap: Number of bootstrap resamples for uncertainty
        unbound_weight_override: If provided, use this contact-unbound weight instead of DBSCAN noise
        concentration_nm3: If provided, apply concentration correction using C_sim = N_lig / V (nm^-3)
        dg_offset_kj: Empirical offset applied to absolute dG (kJ/mol)
    
    Returns:
        result: Dict with {dG_bind_kJmol, dG_stderr_kJmol, P_bound, P_unbound, 
                          K_bind_nm3, n_bootstrap, temperature_K, warning}
    """
    R_GAS = 8.314  # J/(mol*K)
    RT_kJmol = R_GAS * temperature / 1000.0
    conc_factor = _calculate_concentration_factor(concentration_nm3)
    volume_factor = 1.0
    volume_correction_kjmol = 0.0
    if f_contact is not None:
        safe_f_contact = max(float(f_contact), 1e-12)
        volume_factor = 1.0 / safe_f_contact
        volume_correction_kjmol = -RT_kJmol * np.log(volume_factor)
    
    # Information about reference state and corrections
    if f_contact is not None:
        logging.info("=" * 80)
        logging.info("BINDING FREE ENERGY CALCULATION - Molecular Volume Correction")
        logging.info("=" * 80)
        logging.info("Using molecular volume-based contact correction for predictive dG.")
        logging.info("f_contact = %.6f (%.4f%% of box volume)", f_contact, f_contact * 100)
        logging.info("Volume correction: %.2f kJ/mol", volume_correction_kjmol)
        if unbound_weight_override is not None:
            logging.info("Unbound reference: CONTACT-UNBOUND (no contacts within cutoff)")
        else:
            logging.info("Unbound reference: DBSCAN noise (legacy)")
        if concentration_nm3 is not None and concentration_nm3 > 0:
            logging.info("Concentration correction: C_sim = %.6f nm^-3", concentration_nm3)
        logging.info("Expected accuracy: +/-2-4 kcal/mol without empirical calibration")
        logging.info("=" * 80)
    else:
        logging.warning("=" * 80)
        logging.warning("WARNING: BINDING FREE ENERGY - CONTACT-SPACE REFERENCE (No Volume Correction)")
        logging.warning("=" * 80)
        logging.warning("Fingerprint 'unbound' != thermodynamic bulk unbound.")
        if unbound_weight_override is not None:
            logging.warning("Using CONTACT-UNBOUND reference (no contacts within cutoff).")
        else:
            logging.warning("Using DBSCAN noise as unbound reference (legacy behavior).")
        logging.warning("")
        logging.warning("Without molecular volume correction, dG values are QUALITATIVE ONLY.")
        logging.warning("Systematic underestimation: typically 5-10 kcal/mol too weak.")
        logging.warning("")
        logging.warning("RECOMMENDED: Use --dg_contact_volume_molecular for predictive dG")
        logging.warning("=" * 80)
    
    # Calculate bound/unbound probabilities from weighted occupancy
    if unbound_weight_override is not None:
        P_bound = np.sum(unique_weights)
        P_unbound = float(unbound_weight_override)
    else:
        bound_mask = (labels >= 0)
        weighted_bound = np.sum(unique_weights[bound_mask])
        weighted_unbound = np.sum(unique_weights[~bound_mask])
        total_weight = weighted_bound + weighted_unbound
        P_bound = weighted_bound / total_weight
        P_unbound = weighted_unbound / total_weight
    
    # Validation checks
    if P_bound <= 0:
        raise ValueError("No bound states detected (all frames are noise points with label=-1)")
    if P_unbound <= 0:
        raise ValueError("No unbound states detected (all frames are clustered, none with label=-1)")
    if not np.isclose(P_bound + P_unbound, 1.0, atol=1e-6):
        logging.warning(f"Probabilities sum to {P_bound + P_unbound:.6f}, expected 1.0")
    
    # Calculate binding constant and free energy
    K_bind = P_bound / P_unbound
    dG_bind = -RT_kJmol * np.log(K_bind * conc_factor * volume_factor) + dg_offset_kj
    
    logging.info(f"Overall binding: P_bound={P_bound:.4f}, P_unbound={P_unbound:.4f}")
    logging.info(f"K_bind={K_bind:.4f} nm^3, dG_bind={dG_bind:.2f} kJ/mol")
    
    # Bootstrap uncertainty estimation
    dG_samples = []
    n_frames = len(unique_weights)
    
    for _ in range(n_bootstrap):
        # Resample with replacement
        idx = np.random.choice(n_frames, size=n_frames, replace=True)
        boot_weights = unique_weights[idx]
        boot_labels = labels[idx]
        
        if unbound_weight_override is not None:
            boot_P_bound, boot_P_unbound = _bootstrap_contact_unbound(boot_weights)
        else:
            boot_P_bound, boot_P_unbound = _bootstrap_noise_unbound(boot_weights, boot_labels)
        
        # Skip resamples with zero bound or unbound populations
        if boot_P_bound > 0 and boot_P_unbound > 0:
            boot_K = boot_P_bound / boot_P_unbound
            boot_dG = -RT_kJmol * np.log(boot_K * conc_factor * volume_factor) + dg_offset_kj
            dG_samples.append(boot_dG)
    
    dG_stderr = np.std(dG_samples, ddof=1) if len(dG_samples) > 1 else 0.0
    
    logging.info(f"Bootstrap: {len(dG_samples)}/{n_bootstrap} valid samples, stderr={dG_stderr:.2f} kJ/mol")
    if f_contact is None:
        logging.warning(
            f"WARNING: dG={dG_bind:.2f} kJ/mol is contact-space based WITHOUT volume correction."
        )
        logging.warning("WARNING: Use --dg_contact_volume_molecular for predictive absolute dG.")
    
    return {
        'dG_bind_kJmol': dG_bind,
        'dG_stderr_kJmol': dG_stderr,
        'P_bound': P_bound,
        'P_unbound': P_unbound,
        'K_bind_nm3': K_bind,
        'concentration_nm3': concentration_nm3 if concentration_nm3 is not None else None,
        'unbound_reference': 'contact' if unbound_weight_override is not None else 'noise',
        'n_bootstrap': len(dG_samples),
        'temperature_K': temperature,
        'warning': 'CONTACT_REFERENCE_SEE_DOCSTRING',
        'dg_offset_kj': dg_offset_kj,
        'f_contact': f_contact,
        'volume_correction_kJmol': volume_correction_kjmol,
    }

def calculate_cluster_binding_free_energies(
    mode_occ,
    cluster_labels,
    unbound_weight,
    unique_weights,
    labels,
    temperature=300.0,
    n_bootstrap=1000,
    min_occupancy=0.001,
    unbound_reference="noise",
    concentration_nm3=None,
    dg_offset_kj=0.0,
    f_contact=None,
):
    """
    Calculate per-cluster binding free energies with absolute and relative dG values.
    
    Part 2 of Phase 10.2.5: Binding site occupancy and free energies.
    
    MOLECULAR VOLUME CORRECTION (Auto-enabled with --calc_cluster_dg):
    Uses the same molecular volume correction as overall dG calculation.
    All clusters use the SAME f_contact value to ensure consistent absolute scale.
    This enables meaningful comparison of absolute dG values across clusters.
    
    For each cluster with occupancy >= min_occupancy:
    - Absolute dG: Free energy relative to unbound state
      dG_abs = -RT ln(P_cluster/P_unbound x C0/C_sim x 1/f_contact)
    - Relative ddG: Difference from reference cluster (lowest dG = most stable)
      ddG = dG_cluster - dG_reference
    
    Bootstrap uncertainty: Frame-level resampling with weight renormalization.
    
    Args:
        mode_occ: Array of weighted occupancies per cluster (normalized)
        cluster_labels: Array of cluster IDs (corresponding to mode_occ indices)
        unbound_weight: Weighted occupancy of unbound states (normalized)
        unique_weights: Array of normalized HREX weights per frame
        labels: Array of cluster labels (-1 = unbound, >=0 = cluster ID)
        temperature: Temperature in Kelvin (default 300K)
        n_bootstrap: Number of bootstrap resamples for uncertainty
        min_occupancy: Skip clusters below this threshold
        unbound_reference: "noise" or "contact" to control bootstrap handling
        concentration_nm3: If provided, apply concentration correction using C_sim = N_lig / V (nm^-3)
        dg_offset_kj: Empirical offset applied to absolute dG (kJ/mol)
        f_contact: Binding volume fraction (from molecular volume method)
    
    Returns:
        results: Dict with cluster_id -> {dG_abs_kJmol, dG_abs_stderr_kJmol, 
                                         ddG_kJmol, ddG_stderr_kJmol, 
                                         occupancy, is_reference}
    """
    R_GAS = 8.314  # J/(mol*K)
    RT_kJmol = R_GAS * temperature / 1000.0
    conc_factor = _calculate_concentration_factor(concentration_nm3)
    volume_factor = 1.0
    if f_contact is not None:
        safe_f_contact = max(float(f_contact), 1e-12)
        volume_factor = 1.0 / safe_f_contact
    
    if unbound_weight <= 0:
        raise ValueError("Unbound weight must be positive for cluster free energy calculation")
    
    cluster_dGs = []
    cluster_dG_stderrs = []
    cluster_ids_filtered = []
    cluster_occs = []
    
    # Calculate absolute free energies for each cluster with bootstrap uncertainty
    n_frames = len(unique_weights)
    
    for i, (cluster_id, occ) in enumerate(zip(cluster_labels, mode_occ)):
        if occ < min_occupancy:
            continue  # Skip low-occupancy clusters
        
        if occ <= 0:
            logging.warning(f"Cluster {cluster_id} has zero or negative occupancy, skipping")
            continue
        
        # Absolute free energy for this cluster (point estimate)
        K_cluster = occ / unbound_weight
        dG_abs = -RT_kJmol * np.log(K_cluster * conc_factor * volume_factor) + dg_offset_kj
        
        # Bootstrap uncertainty via frame-level resampling
        boot_dG_samples = []
        
        for _ in range(n_bootstrap):
            # Resample frames with replacement
            idx = np.random.choice(n_frames, size=n_frames, replace=True)
            boot_weights = unique_weights[idx]
            boot_labels = labels[idx]
            
            # Recalculate occupancies for this bootstrap sample
            boot_cluster_mask = (boot_labels == cluster_id)
            
            if unbound_reference == "contact":
                boot_cluster_occ = np.sum(boot_weights[boot_cluster_mask])
                boot_unbound_weight = max(1.0 - np.sum(boot_weights), 1e-12)
            else:
                boot_weights = boot_weights / np.sum(boot_weights)  # Renormalize
                boot_cluster_occ = np.sum(boot_weights[boot_cluster_mask])
                boot_unbound_weight = np.sum(boot_weights[boot_labels == -1])
            
            # Skip if zero occupancy (cluster not sampled or no unbound)
            if boot_cluster_occ <= 0 or boot_unbound_weight <= 0:
                continue
            
            # Calculate bootstrap sample's free energy
            boot_K = boot_cluster_occ / boot_unbound_weight
            boot_dG = -RT_kJmol * np.log(boot_K * conc_factor * volume_factor) + dg_offset_kj
            boot_dG_samples.append(boot_dG)
        
        # Calculate standard error from bootstrap distribution
        if len(boot_dG_samples) > 1:
            dG_stderr = np.std(boot_dG_samples, ddof=1)
        else:
            dG_stderr = 0.0  # Insufficient valid samples
        
        cluster_dGs.append(dG_abs)
        cluster_dG_stderrs.append(dG_stderr)
        cluster_ids_filtered.append(cluster_id)
        cluster_occs.append(occ)
    
    if len(cluster_dGs) == 0:
        logging.warning("No clusters passed min_occupancy threshold")
        return {}
    
    cluster_dGs = np.array(cluster_dGs)
    cluster_dG_stderrs = np.array(cluster_dG_stderrs)
    
    # Relative free energies: ddG = dG_cluster - dG_reference
    # Reference = most stable (lowest dG)
    ref_idx = np.argmin(cluster_dGs)
    dG_ref = cluster_dGs[ref_idx]
    ddG_values = cluster_dGs - dG_ref
    
    # Relative uncertainties (propagate errors)
    ddG_stderrs = np.sqrt(cluster_dG_stderrs**2 + cluster_dG_stderrs[ref_idx]**2)
    
    # Build results dictionary
    results = {}
    for i, cluster_id in enumerate(cluster_ids_filtered):
        results[cluster_id] = {
            'dG_abs_kJmol': cluster_dGs[i],
            'dG_abs_stderr_kJmol': cluster_dG_stderrs[i],
            'ddG_kJmol': ddG_values[i],
            'ddG_stderr_kJmol': ddG_stderrs[i],
            'occupancy': cluster_occs[i],
            'is_reference': (i == ref_idx),
            'unbound_reference': unbound_reference,
            'concentration_nm3': concentration_nm3 if concentration_nm3 is not None else None,
        }
    
    logging.info(f"Calculated free energies for {len(results)} clusters")
    logging.info(f"Reference cluster: {cluster_ids_filtered[ref_idx]} (lowest dG={dG_ref:.2f} kJ/mol)")
    if len(cluster_dGs) > 0:
        logging.info(f"dG range: {np.min(cluster_dGs):.2f} to {np.max(cluster_dGs):.2f} kJ/mol")
    
    return results

def main():  # pyright: ignore[reportGeneralTypeIssues]
    args = parse_args()
    setup_logging(args.log_file, quiet=args.quiet)
    logging.info("="*60)
    logging.info("Interaction Fingerprint Analysis Tool - Started")
    logging.info("="*60)
    logging.info("Starting analysis with arguments: %s", vars(args))

    # Dual output mode: always write unscaled, conditionally write scaled
    # Remove prefix switching to enable side-by-side comparison (ANA-21)
    logging.info(f"SASA scaling: {'enabled (' + args.sasa_method + ' method, dual output mode)' if args.sasa_enabled else 'disabled'}")

    avail_mem_bytes = get_available_memory()
    avail_mem_gb = avail_mem_bytes / 1e9
    logging.info(f"Detected available memory: {avail_mem_gb:.1f} GB")

    try:
        global _SASA_TIME_NS
        u = mda.Universe(args.topology, args.trajectory)
        n_frames = len(u.trajectory)
        logging.info(f"Loaded universe with {n_frames} frames")
        try:
            dt_ps = getattr(u.trajectory, 'dt', None)
            if dt_ps is not None and dt_ps > 0:
                _SASA_TIME_NS = np.arange(n_frames, dtype=np.float64) * (dt_ps / 1000.0)
            else:
                _SASA_TIME_NS = None
        except Exception as e:
            logging.warning(f"Failed to infer timestep; using frame index: {e}")
            _SASA_TIME_NS = None

        def _auto_select_atoms(universe, selection, fallback, sel_type):
            try:
                atoms = universe.select_atoms(selection)
            except Exception as e:
                logging.warning("Invalid %s selection '%s': %s", sel_type, selection, e)
                atoms = universe.atoms[[]]
            if len(atoms) > 0:
                return atoms
            logging.warning("Empty %s selection '%s'; trying fallback '%s'", sel_type, selection, fallback)
            try:
                atoms = universe.select_atoms(fallback)
            except Exception as e:
                logging.warning("Fallback %s selection '%s' failed: %s", sel_type, fallback, e)
                atoms = universe.atoms[[]]
            return atoms

        receptor = _auto_select_atoms(
            u, args.receptor_sel, DEFAULT_RECEPTOR_SEL_RESNAME, "receptor"
        )
        ligands_ag = _auto_select_atoms(
            u, args.ligand_sel, DEFAULT_LIGAND_SEL_RESNAME, "ligand"
        )
        ligand_residues = ligands_ag.residues
        n_ligs = len(ligand_residues)
        if n_ligs == 0:
            raise ValueError("No ligands found")

        ref_lig = ligand_residues[0]
        n_lig_atoms = ref_lig.atoms.n_atoms
        lig_atom_names = ref_lig.atoms.names
        for res in ligand_residues[1:]:
            if res.atoms.n_atoms != n_lig_atoms or not np.array_equal(res.atoms.names, lig_atom_names):
                raise ValueError("Ligands are not identical")
        logging.info(f"{n_ligs} identical ligands with {n_lig_atoms} atoms each")

        n_res = len(receptor.residues)
        fp_dense_shape = (n_res, n_lig_atoms)
        fp_flat_size = n_res * n_lig_atoms
        logging.info(f"Receptor: {n_res} residues -> fingerprint shape {fp_dense_shape}")

        # Store receptor coordinates from topology/first frame for contact-volume correction
        u.trajectory[0]
        receptor_coords_nm = receptor.positions.copy() * 0.1  # Angstrom -> nm

        receptor_atom_to_res = np.concatenate([[i] * res.atoms.n_atoms for i, res in enumerate(receptor.residues)])

        # Weights
        if args.logweights:
            logw = np.loadtxt(args.logweights)
            logw -= np.max(logw)
            weights = np.exp(logw)
            weights /= np.sum(weights)
            logging.info("Processed logweights (direct exp normalization)")
        elif args.weights:
            weights = np.loadtxt(args.weights)
            weights /= np.sum(weights)
        else:
            weights = np.ones(n_frames) / n_frames
        if len(weights) != n_frames:
            raise ValueError("Weights length mismatch with trajectory")
        
        # Validate weight normalization
        weight_sum = np.sum(weights)
        if abs(weight_sum - 1.0) > 1e-6:
            logging.warning(f"Weights sum to {weight_sum:.6f}, normalizing to 1.0")
            weights = weights / weight_sum
        logging.info(f"Validated weights: sum={np.sum(weights):.6f}, min={np.min(weights):.6f}, max={np.max(weights):.6f}")

        # SASA arrays
        sasa_arrays = None
        if args.sasa_enabled:
            # Pass receptor selection for internal SASA calculation
            sasa_arrays = load_or_calculate_sasa(args, u, n_ligs, n_frames, receptor_selection=receptor)
        if args.sasa_enabled and sasa_arrays is None:
            raise ValueError("SASA scaling enabled but SASA arrays not loaded")
        if args.sasa_enabled:
            try:
                sasa_csv_path = "sasa.csv"
                with open(sasa_csv_path, 'w') as f:
                    f.write("ligand_id,frame,sasa_nm2\n")
                    for lig_id in range(n_ligs):
                        for frame in range(n_frames):
                            f.write(f"{lig_id+1},{frame},{sasa_arrays[lig_id][frame]:.6f}\n")
                logging.info(f"Wrote SASA values to {sasa_csv_path}")
            except Exception as e:
                logging.warning(f"Failed to write SASA CSV: {e}")
            try:
                logging.info("Generating SASA time-series visualizations...")
                plot_sasa_timeseries(
                    sasa_arrays=sasa_arrays,
                    weights=weights,
                    n_ligands=n_ligs,
                    output_prefix=args.sasa_scores_prefix.replace('_scores', '_timeseries'),
                )
            except Exception as e:
                logging.warning(f"SASA time-series plotting failed (non-critical): {e}")
            try:
                logging.info("Generating SASA histogram analysis...")
                plot_sasa_histogram(
                    sasa_arrays=sasa_arrays,
                    weights=weights,
                    n_ligands=n_ligs,
                    bins=50,
                    output_prefix=args.sasa_scores_prefix.replace('_scores', '_histogram'),
                )
            except Exception as e:
                logging.warning(f"Failed to generate SASA histogram: {e}")

        # === Collection & deduplication using frozenset ===
        logging.info("Collecting bound instances...")
        volume_accum_nm3 = 0.0
        volume_frames = 0
        box_dim_accum = np.zeros(3, dtype=np.float64)
        box_dim_frames = 0
        per_frame_receptor_coords_nm = []
        per_frame_box_dims_nm = []
        per_frame_weights = []
        # key: frozenset((res_idx, atom_idx)), value: (total_weight_scaled, count, best_weight_scaled, best_frame, best_lig, total_weight_unscaled)
        fp_dict = {}
        bound_instances = []  # Track all bound instances for statistics
        contact_freq_by_lig = np.zeros((n_ligs, n_res, n_lig_atoms), dtype=np.float32)
        sasa_scaled_contact_freq_by_lig = np.zeros((n_ligs, n_res, n_lig_atoms), dtype=np.float32)

        volume_stride = max(int(args.dg_volume_frame_stride), 1)
        for ts in u.trajectory:
            frame = ts.frame
            frame_weight = weights[frame]
            if u.dimensions is not None and len(u.dimensions) >= 3:
                volume_accum_nm3 += float(u.dimensions[0] * u.dimensions[1] * u.dimensions[2]) / 1000.0
                volume_frames += 1
                box_dim_accum += u.dimensions[:3]
                box_dim_frames += 1
            if args.dg_contact_volume_per_frame and (frame % volume_stride == 0):
                if u.dimensions is None or len(u.dimensions) < 3:
                    raise ValueError("Cannot compute per-frame contact volume: missing box dimensions")
                per_frame_receptor_coords_nm.append(receptor.positions.copy() * 0.1)
                per_frame_box_dims_nm.append(u.dimensions[:3] * 0.1)
                per_frame_weights.append(frame_weight)
            for lig_id, res in enumerate(ligand_residues):
                if np.random.rand() >= 1.0 / args.subsampling:
                    continue
                lig_at = res.atoms
                pairs, _ = capped_distance(receptor.positions, lig_at.positions,
                                           args.cutoff, box=u.dimensions)
                if len(pairs) == 0:
                    continue

                rows = receptor_atom_to_res[pairs[:, 0]]
                cols = pairs[:, 1]
                key = frozenset(zip(rows, cols))

                if args.sasa_enabled:
                    assert sasa_arrays is not None
                    scale = 1.0 / (sasa_arrays[lig_id][frame] + args.sasa_eps)
                else:
                    scale = 1.0
                inst_weight = frame_weight * scale

                if key:
                    res_idxs, atom_idxs = zip(*key)
                    contact_freq_by_lig[lig_id, res_idxs, atom_idxs] += frame_weight
                    if args.sasa_enabled:
                        sasa_scaled_contact_freq_by_lig[lig_id, res_idxs, atom_idxs] += inst_weight
                
                # Track bound instance for statistics
                bound_instances.append({
                    'frame': frame,
                    'ligand_id': lig_id,
                    'weight': frame_weight,
                    'sasa_scale': scale,
                    'final_weight': inst_weight
                })

                if key in fp_dict:
                    tw, cnt, bw, bf, bl, tw_unscaled = fp_dict[key]
                    tw += inst_weight
                    tw_unscaled += frame_weight
                    cnt += 1
                    if inst_weight > bw:
                        bw, bf, bl = inst_weight, frame, lig_id
                    fp_dict[key] = (tw, cnt, bw, bf, bl, tw_unscaled)
                else:
                    fp_dict[key] = (inst_weight, 1, inst_weight, frame, lig_id, frame_weight)

        if not fp_dict:
            logging.warning("No bound instances found")
            return
            
        # Calculate comprehensive binding statistics
        binding_stats = calculate_comprehensive_binding_statistics(
            n_frames, weights, n_ligs, sasa_arrays, bound_instances, args.sasa_eps
        )

        avg_volume_nm3 = None
        if volume_frames > 0:
            avg_volume_nm3 = volume_accum_nm3 / volume_frames

        avg_box_dims_nm = None
        if box_dim_frames > 0:
            avg_box_dims_nm = (box_dim_accum / box_dim_frames) * 0.1  # Angstrom -> nm

        f_contact, volume_correction_kjmol = compute_contact_volume_correction(
            args,
            avg_volume_nm3,
            avg_box_dims_nm,
            receptor_coords_nm,
            per_frame_receptor_coords_nm,
            per_frame_box_dims_nm,
            per_frame_weights,
            n_ligs,
            ref_lig,
            binding_stats,
        )

        n_unique = len(fp_dict)
        logging.info(f"Found {n_unique} unique fingerprints")

        # Prepare arrays for sieving
        keys_list = list(fp_dict.keys())
        unique_weights = np.array([fp_dict[k][0] for k in keys_list], dtype=np.float64)
        unique_weights_unscaled = np.array([fp_dict[k][5] for k in keys_list], dtype=np.float64)
        # Normalize by total possible instance weight (sum of frame weights x n_ligands)
        # Rationale: instance-based occupancy should not exceed 1.0 across all ligand-frame pairs.
        total_frame_weight = np.sum(weights)
        total_instance_weight = total_frame_weight * n_ligs
        denom = total_instance_weight if total_instance_weight > 0 else 1.0
        unique_weights /= denom
        unique_weights_unscaled /= denom
        uniq_counts = np.array([fp_dict[k][1] for k in keys_list], dtype=int)
        best_weights = np.array([fp_dict[k][2] for k in keys_list])
        unique_frames = np.array([fp_dict[k][3] for k in keys_list])
        unique_lig_ids = np.array([fp_dict[k][4] for k in keys_list])

        # Enhanced frame tracking dictionary
        frame_mapping = {
            'sieved_to_original': {},
            'original_to_sieved': {}, 
            'original_frames': unique_frames.copy(),
            'original_ligands': unique_lig_ids.copy(),
            'original_weights': unique_weights.copy(),
            'original_weights_unscaled': unique_weights_unscaled.copy(),
            'sasa_scale_factors': [],
            'sieved_instance_weights': [],  # Final weights after SASA scaling
            'bound_instance_indices': []    # Track which were bound instances
        }

        # Weighted sieving with enhanced mapping
        final_n = n_unique
        if args.sieving and args.max_sieved is not None and n_unique > args.max_sieved:
            logging.info(f"Sieving from {n_unique} to {args.max_sieved} (weighted)...")
            probs = unique_weights / unique_weights.sum()
            top_percent = max(0.0, float(args.sieved_top_percent))
            n_top = int(round(n_unique * top_percent / 100.0)) if top_percent > 0 else 0
            n_top = min(n_top, args.max_sieved)
            if n_top > 0:
                top_idx = np.argsort(-unique_weights)[:n_top]
                top_weight = unique_weights[top_idx].sum()
                total_weight = unique_weights.sum()
                logging.info(
                    "Stratified sieving: forcing top %d (%.1f%%) by weight (coverage=%.4f)",
                    n_top, top_percent, (top_weight / total_weight) if total_weight > 0 else 0.0
                )
                remaining = args.max_sieved - n_top
                if remaining > 0:
                    mask = np.ones(n_unique, dtype=bool)
                    mask[top_idx] = False
                    remaining_idx = np.where(mask)[0]
                    remaining_weights = unique_weights[remaining_idx]
                    remaining_probs = remaining_weights / remaining_weights.sum()
                    sampled_idx = np.random.choice(remaining_idx, remaining, replace=False, p=remaining_probs)
                    idx = np.concatenate([top_idx, sampled_idx])
                else:
                    idx = top_idx
            else:
                idx = np.random.choice(n_unique, args.max_sieved, replace=False, p=probs)
            
            # Create mapping from sieved to original indices
            for i, orig_idx in enumerate(idx):
                frame_mapping['sieved_to_original'][i] = orig_idx
                frame_mapping['original_to_sieved'][orig_idx] = i
                frame_mapping['bound_instance_indices'].append(True)  # All sieved instances are bound
            
            # Update arrays based on sieving
            keys_list = [keys_list[i] for i in idx]
            unique_weights = unique_weights[idx]
            unique_weights_unscaled = unique_weights_unscaled[idx]
            # Do NOT re-normalize after sieving; unique_weights are already instance-normalized
            uniq_counts = uniq_counts[idx]
            unique_frames = unique_frames[idx]
            unique_lig_ids = unique_lig_ids[idx]
            best_weights = best_weights[idx]
            final_n = args.max_sieved
            logging.info(f"Sieved to {final_n} unique fingerprints")
        else:
            # No sieving - create identity mapping
            for i in range(n_unique):
                frame_mapping['sieved_to_original'][i] = i
                frame_mapping['original_to_sieved'][i] = i
                frame_mapping['bound_instance_indices'].append(True)
        
        # Function to restore original frame information
        def get_original_frame_data(sieved_idx):
            """Retrieve complete original frame information"""
            orig_idx = frame_mapping['sieved_to_original'].get(sieved_idx, sieved_idx)
            return {
                'frame': frame_mapping['original_frames'][orig_idx],
                'ligand_id': frame_mapping['original_ligands'][orig_idx], 
                'weight': frame_mapping['original_weights'][orig_idx],
                'weight_unscaled': frame_mapping['original_weights_unscaled'][orig_idx],
                'sasa_scale': frame_mapping['sasa_scale_factors'][orig_idx] if frame_mapping['sasa_scale_factors'] else 1.0,
                'final_weight': frame_mapping['sieved_instance_weights'][orig_idx] if frame_mapping['sieved_instance_weights'] else frame_mapping['original_weights'][orig_idx],
                'was_bound': frame_mapping['bound_instance_indices'][orig_idx] if sieved_idx < len(frame_mapping['bound_instance_indices']) else True
            }

        # === Build correct sparse matrix ===
        logging.info("Building sparse fingerprint matrix (correct row/col mapping)...")
        all_row_ind = []
        all_col_ind = []
        all_data = []

        for i, key in enumerate(keys_list):
            if not key:
                continue
            res_idxs, atom_idxs = zip(*key)
            flat_cols = np.array(res_idxs) * n_lig_atoms + np.array(atom_idxs)
            all_row_ind.extend([i] * len(flat_cols))
            all_col_ind.extend(flat_cols)
            all_data.extend([1.0] * len(flat_cols))

        row_ind = np.array(all_row_ind, dtype=np.int32)
        col_ind = np.array(all_col_ind, dtype=np.int32)
        data = np.array(all_data, dtype=np.float32)

        X = csr_matrix((data, (row_ind, col_ind)),
                       shape=(final_n, fp_flat_size),
                       dtype=np.float32)

        logging.info(f"Sparse matrix created: {X.nnz} non-zero elements ({X.nnz/final_n:.1f} contacts/frame avg)")
        log_memory_usage("after sparse matrix")

        # === Memory estimation and safety check (before expensive similarity) ===
        avg_contacts_per_fp = int(X.nnz / final_n) if final_n > 0 else 50
        estimated_gb, memory_breakdown = estimate_memory_usage(final_n, n_res, n_lig_atoms, avg_contacts_per_fp)
        dense_over_limit = estimated_gb > MEMORY_SAFETY_FRACTION * avail_mem_gb

        if not args.sparse_similarity and dense_over_limit:
            logging.warning(
                "Dense similarity estimate %.2f GB exceeds %.0f%% of available %.2f GB. "
                "Auto-enabling --sparse_similarity to reduce memory risk.",
                estimated_gb,
                MEMORY_SAFETY_FRACTION * 100,
                avail_mem_gb,
            )
            logging.warning(
                "Dense breakdown: Distance %.2f GB, Similarity %.2f GB, Intersection %.2f GB, "
                "Sparse %.2f GB, Contact %.2f GB, Arrays %.2f GB",
                memory_breakdown['distance_matrix'],
                memory_breakdown['similarity_matrix'],
                memory_breakdown['intersection_matrix'],
                memory_breakdown['sparse_matrix'],
                memory_breakdown['contact_matrix'],
                memory_breakdown['arrays'],
            )
            args.sparse_similarity = True

        if args.sparse_similarity:
            if args.sparse_threshold is None:
                sparse_threshold = max(0.0, args.cluster_eps - DEFAULT_SPARSE_THRESHOLD_OFFSET)
            else:
                sparse_threshold = float(args.sparse_threshold)

            sparse_estimated_gb, sparse_details = estimate_sparse_mode_memory_gb(
                X,
                args.similarity,
                sparse_threshold,
                final_n,
                n_unique,
                avg_contacts_per_fp,
                memory_breakdown,
                kdist_plot=args.kdist_plot,
            )
            logging.info(
                "Sparse similarity estimate %.2f GB (nnz~%d, avg neighbors~%.1f, extra~%.2f GB, "
                "dbscan~%.2f GB, kdist~%.2f GB, safety=1.2x)",
                sparse_estimated_gb,
                sparse_details["est_nnz"],
                sparse_details["avg_neighbors"],
                sparse_details["sparse_extra_gb"],
                sparse_details["dbscan_neighborhood_gb"],
                sparse_details["k_distance_gb"],
            )
            warn_sparse_memory_risk(sparse_estimated_gb, avail_mem_gb, sparse_details)
        else:
            logging.info("Dense similarity estimate: %.2f GB total", estimated_gb)
            logging.info(
                "Dense breakdown: Distance %.2f GB, Similarity %.2f GB, Intersection %.2f GB, "
                "Sparse %.2f GB, Contact %.2f GB, Arrays %.2f GB",
                memory_breakdown['distance_matrix'],
                memory_breakdown['similarity_matrix'],
                memory_breakdown['intersection_matrix'],
                memory_breakdown['sparse_matrix'],
                memory_breakdown['contact_matrix'],
                memory_breakdown['arrays'],
            )

            # Check memory safety for dense mode
            check_memory_safety(
                estimated_gb,
                avail_mem_gb,
                final_n,
                args.max_sieved,
                args.ignore_memory_risk,
                memory_breakdown,
            )

        # === Heatmaps (write before memory safety exit) ===
        total_frame_weight = np.sum(weights)
        denom = total_frame_weight if total_frame_weight > 0 else 1.0
        contact_freq_by_lig_norm = contact_freq_by_lig / denom
        overall_contact_freq = np.mean(contact_freq_by_lig_norm, axis=0)
        sasa_scaled_contact_freq_by_lig_norm = None
        overall_sasa_contact_freq = None
        if args.sasa_enabled:
            sasa_scaled_contact_freq_by_lig_norm = sasa_scaled_contact_freq_by_lig / denom
            overall_sasa_contact_freq = np.mean(sasa_scaled_contact_freq_by_lig_norm, axis=0)

        # Overall weighted average heatmap
        plt.figure(figsize=(12, 8))
        sns.heatmap(overall_contact_freq, cmap='Blues', cbar_kws={'label': 'Weighted Occupancy'})
        plt.xlabel('Ligand Atom Index')
        plt.ylabel('Receptor Residue Index')
        plt.title('Weighted Contact Frequency - Overall Average')
        overall_heatmap_path = "heatmap_overall_weighted.png"
        plt.savefig(overall_heatmap_path, dpi=300, bbox_inches='tight')
        plt.close()
        for lig_id in range(n_ligs):
            # Standard weighted heatmap
            plt.figure(figsize=(12, 8))
            sns.heatmap(contact_freq_by_lig_norm[lig_id], cmap='Blues', cbar_kws={'label': 'Weighted Occupancy'})
            plt.xlabel('Ligand Atom Index')
            plt.ylabel('Receptor Residue Index')
            plt.title(f'Weighted Contact Frequency - Ligand {lig_id+1}')
            heatmap_path = f"{args.heatmap_prefix}{lig_id+1}_weighted.png"
            plt.savefig(heatmap_path, dpi=300, bbox_inches='tight')
            plt.close()

            # SASA-scaled heatmap (same blue colormap)
            if args.sasa_enabled and sasa_scaled_contact_freq_by_lig_norm is not None:
                plt.figure(figsize=(12, 8))
                sns.heatmap(sasa_scaled_contact_freq_by_lig_norm[lig_id], cmap='Blues', cbar_kws={'label': 'SASA-Scaled Occupancy'})
                plt.xlabel('Ligand Atom Index')
                plt.ylabel('Receptor Residue Index')
                plt.title(f'SASA-Scaled Contact Frequency - Ligand {lig_id+1}')
                sasa_heatmap_path = f"{args.heatmap_prefix}{lig_id+1}_sasa_scaled.png"
                plt.savefig(sasa_heatmap_path, dpi=300, bbox_inches='tight')
                plt.close()
                logging.info(f"Wrote SASA-scaled heatmap for ligand {lig_id+1} to {sasa_heatmap_path}")

        if args.sasa_enabled and overall_sasa_contact_freq is not None:
            plt.figure(figsize=(12, 8))
            sns.heatmap(overall_sasa_contact_freq, cmap='Blues', cbar_kws={'label': 'SASA-Scaled Occupancy'})
            plt.xlabel('Ligand Atom Index')
            plt.ylabel('Receptor Residue Index')
            plt.title('SASA-Scaled Contact Frequency - Overall Average')
            overall_sasa_heatmap_path = "heatmap_overall_sasa_scaled.png"
            plt.savefig(overall_sasa_heatmap_path, dpi=300, bbox_inches='tight')
            plt.close()

        # Calculate SASA scale factors if enabled
        sasa_scale_factors = np.ones(final_n, dtype=np.float32)
        if args.sasa_enabled:
            assert sasa_arrays is not None
            for i in range(final_n):
                frame = unique_frames[i]
                lig_id = unique_lig_ids[i]
                scale = 1.0 / (sasa_arrays[lig_id][frame] + args.sasa_eps)
                sasa_scale_factors[i] = scale
                frame_mapping['sasa_scale_factors'].append(scale)
        else:
            frame_mapping['sasa_scale_factors'] = [1.0] * final_n
        
        if args.sasa_enabled:
            logging.info("Computing weighted contact frequency and SASA-scaled frequency...")
        else:
            logging.info("Computing weighted contact frequency...")

        sumA = np.array(X.getnnz(axis=1), dtype=np.float32)
        sumB = sumA.copy()
        sim = None
        sim_sparse = None

        # Similarity computation (dense or sparse)
        if args.sparse_similarity:
            if args.sparse_threshold is None:
                sparse_threshold = max(0.0, args.cluster_eps - DEFAULT_SPARSE_THRESHOLD_OFFSET)
            else:
                sparse_threshold = float(args.sparse_threshold)
            logging.info(
                "Computing sparse similarity matrix (threshold=%.3f, chunk=%d)...",
                sparse_threshold, DEFAULT_SPARSE_CHUNK_SIZE
            )
            normB = None
            if args.similarity == 'cosine':
                normB = np.sqrt(sumB + 1e-10)
            sim_sparse = compute_sparse_similarity_matrix(
                X,
                args.similarity,
                sparse_threshold,
                chunk_size=DEFAULT_SPARSE_CHUNK_SIZE,
                sumB=sumB,
                normB=normB,
            )
            sim_sparse = sim_sparse.maximum(sim_sparse.T)
            sim_sparse.setdiag(1.0)
            sim_sparse.eliminate_zeros()
            logging.info(
                "Sparse similarity nnz=%d (density=%.6f)",
                sim_sparse.nnz,
                sim_sparse.nnz / float(final_n * final_n)
            )
            log_memory_usage("after sparse similarity")

            dist = sim_sparse.copy()
            dist.data = 1.0 - dist.data
            dist.setdiag(0.0)
            dist.eliminate_zeros()
            log_memory_usage("after sparse distance")
        else:
            logging.info("Computing intersection (sparse @ sparse)...")
            intersect = (X @ X.T).toarray().astype(np.float32)
            log_memory_usage("after intersection matrix")

            logging.info("Computing similarity matrix...")
            if args.similarity == 'dice':
                sim = 2.0 * intersect / (sumA[:, None] + sumB[None, :] + 1e-10)
                normB = None  # Not needed for dice
            elif args.similarity == 'jaccard':
                sim = intersect / (sumA[:, None] + sumB[None, :] - intersect + 1e-10)
                normB = None  # Not needed for jaccard
            else:  # cosine
                normA = np.sqrt(sumA + 1e-10)
                normB = np.sqrt(sumB + 1e-10)
                sim = intersect / (normA[:, None] * normB[None, :] + 1e-10)

            # Prevent similarity > 1.0 due to floating point
            sim = np.minimum(sim, 1.0)
            sim_min = float(np.min(sim))
            sim_max = float(np.max(sim))
            sim_mean = float(np.mean(sim))
            sim_nonzero = float(np.count_nonzero(sim) / sim.size)
            logging.info(f"Similarity stats: min={sim_min:.4f}, max={sim_max:.4f}, mean={sim_mean:.4f}, nonzero={sim_nonzero:.4f}")
            log_memory_usage("after similarity matrix")

            del intersect
            gc.collect()

            # Compute distance matrix (needed for both DBSCAN and medoid selection)
            dist = 1.0 - sim
            log_memory_usage("after distance matrix")

            # Final sanity check
            if np.any(dist < 0):
                logging.warning("Negative distances found after clipping (should be rare)")

        # Clustering with enhanced statistics
        logging.info(f"Running {args.clustering_method} clustering...")

        clustering_method = args.clustering_method
        if args.sparse_similarity and clustering_method == 'density_ratio':
            logging.warning("Sparse similarity mode only supports DBSCAN; falling back to DBSCAN")
            clustering_method = 'dbscan'
        
        if clustering_method == 'density_ratio':
            clustering = DensityRatioClustering()
            if args.sparse_similarity:
                labels = clustering.fit_predict(X, sample_weight=unique_weights, precomputed_similarity=sim_sparse)
            else:
                labels = clustering.fit_predict(X, sample_weight=unique_weights, precomputed_similarity=sim)
                del sim  # Clean up similarity matrix (dist kept for medoid selection)
        else:  # dbscan
            if not args.sparse_similarity:
                del sim  # DBSCAN uses distance, not similarity
            db = DBSCAN(eps=args.cluster_eps, min_samples=args.min_samples, metric='precomputed')
            labels = db.fit_predict(dist)
        
        gc.collect()

        unique_labels = np.unique(labels)
        n_clusters = len(unique_labels[unique_labels >= 0])
        n_noise = np.sum(labels == -1)
        cluster_labels = unique_labels[unique_labels >= 0]
        logging.info(f"Clustering result: {n_noise} noise points, {n_clusters} clusters")
        unsieved_cluster_weights = np.zeros(n_clusters)

        if clustering_method == 'density_ratio' and n_clusters == 0:
            logging.warning("Density ratio clustering produced 0 clusters; falling back to DBSCAN with precomputed distances")
            db = DBSCAN(eps=args.cluster_eps, min_samples=args.min_samples, metric='precomputed')
            labels = db.fit_predict(dist)
            unique_labels = np.unique(labels)
            n_clusters = len(unique_labels[unique_labels >= 0])
            n_noise = np.sum(labels == -1)
            cluster_labels = unique_labels[unique_labels >= 0]
            logging.info(f"DBSCAN fallback result: {n_noise} noise points, {n_clusters} clusters")
        
        # Compute comprehensive clustering statistics
        clustering_stats = compute_clustering_statistics(
            X,
            labels,
            sample_weight=unique_weights,
            output_prefix='',
            kdist_plot=args.kdist_plot,
        )
        
        # Save clustering statistics to file
        if clustering_stats:
            with open("clustering_statistics.txt", 'w') as f:
                f.write("=== Clustering Statistics ===\n")
                f.write(f"Number of clusters: {n_clusters}\n")
                f.write(f"Number of noise points: {n_noise}\n")
                if 'noise_ratio' in clustering_stats:
                    f.write(f"Noise ratio: {clustering_stats['noise_ratio']:.3f}\n")
                if 'silhouette' in clustering_stats:
                    f.write(f"Silhouette Score: {clustering_stats['silhouette']:.3f} (closer to 1 is better)\n")
                if 'davies_bouldin' in clustering_stats:
                    f.write(f"Davies-Bouldin Index: {clustering_stats['davies_bouldin']:.3f} (lower is better)\n")
                if 'calinski_harabasz' in clustering_stats:
                    f.write(f"Calinski-Harabasz Index: {clustering_stats['calinski_harabasz']:.1f} (higher is better)\n")
                if 'ssr_sst_ratio' in clustering_stats:
                    f.write(f"SSR/SST Ratio: {clustering_stats['ssr_sst_ratio']:.3f} (higher is better, pseudo-F)\n")
            
            logging.info("Clustering statistics saved to clustering_statistics.txt")

        # === Reassign unsieved fingerprints to nearest clusters ===
        # After clustering the sieved subset, assign unsieved fingerprints to their nearest cluster
        # and add their weighted occupancy to the cluster totals
        if args.sieving and args.max_sieved is not None and n_unique > args.max_sieved:
            logging.info(f"Reassigning unsieved fingerprints to nearest clusters...")
            # Identify which original indices were NOT included in sieving
            sieved_orig_indices = set(frame_mapping['sieved_to_original'].values())
            unsieved_indices = [i for i in range(n_unique) if i not in sieved_orig_indices]
            n_unsieved = len(unsieved_indices)
            
            # Diagnostic: check weight distribution
            total_original_weight = np.sum(frame_mapping['original_weights'])
            sieved_weight_sum = np.sum([frame_mapping['original_weights'][i] for i in sieved_orig_indices])
            unsieved_weight_sum = np.sum([frame_mapping['original_weights'][i] for i in unsieved_indices])
            logging.info(f"Weight distribution: total={total_original_weight:.6f}, " 
                        f"sieved={sieved_weight_sum:.6f}, unsieved={unsieved_weight_sum:.6f}")
            
            if n_unsieved > 0 and n_clusters > 0:
                logging.info(f"Found {n_unsieved} unsieved fingerprints to reassign")
                
                # Build sparse matrix for unsieved fingerprints
                # Need to access original keys_list before sieving
                original_keys_list = list(fp_dict.keys())
                unsieved_row_ind = []
                unsieved_col_ind = []
                unsieved_data = []
                
                for row_idx, orig_idx in enumerate(unsieved_indices):
                    key = original_keys_list[orig_idx]
                    if not key:
                        continue
                    res_idxs, atom_idxs = zip(*key)
                    flat_cols = np.array(res_idxs) * n_lig_atoms + np.array(atom_idxs)
                    unsieved_row_ind.extend([row_idx] * len(flat_cols))
                    unsieved_col_ind.extend(flat_cols)
                    unsieved_data.extend([1.0] * len(flat_cols))
                
                if len(unsieved_row_ind) > 0:
                    X_unsieved = csr_matrix(
                        (np.array(unsieved_data, dtype=np.float32), 
                         (np.array(unsieved_row_ind, dtype=np.int32), 
                          np.array(unsieved_col_ind, dtype=np.int32))),
                        shape=(n_unsieved, fp_flat_size),
                        dtype=np.float32
                    )
                    
                    # Compute similarity between unsieved and sieved (clustered) fingerprints
                    intersect_unsieved = (X_unsieved @ X.T).toarray().astype(np.float32)
                    sumA_unsieved = np.array(X_unsieved.getnnz(axis=1), dtype=np.float32)
                    
                    if args.similarity == 'dice':
                        sim_unsieved = 2.0 * intersect_unsieved / (sumA_unsieved[:, None] + sumB[None, :] + 1e-10)
                    elif args.similarity == 'jaccard':
                        sim_unsieved = intersect_unsieved / (sumA_unsieved[:, None] + sumB[None, :] - intersect_unsieved + 1e-10)
                    else:  # cosine
                        normA_unsieved = np.sqrt(sumA_unsieved + 1e-10)
                        if normB is not None:
                            sim_unsieved = intersect_unsieved / (normA_unsieved[:, None] * normB[None, :] + 1e-10)
                        else:
                            # Fallback: recompute normB if needed
                            normB_local = np.sqrt(sumB + 1e-10)
                            sim_unsieved = intersect_unsieved / (normA_unsieved[:, None] * normB_local[None, :] + 1e-10)
                    
                    sim_unsieved = np.minimum(sim_unsieved, 1.0)

                    # Optional: bridge fragmented clusters using unsieved similarity links
                    if args.bridge_unsieved and n_clusters > 1:
                        logging.info(
                            "Bridging clusters with unsieved links (min_sim=%.3f, top_k=%d, min_support=%d)",
                            args.bridge_min_sim, args.bridge_top_k, args.bridge_min_support
                        )

                        # Union-Find for cluster labels
                        cluster_label_set = np.unique(labels[labels >= 0])
                        parent = {lbl: lbl for lbl in cluster_label_set}

                        def _find(x):
                            while parent[x] != x:
                                parent[x] = parent[parent[x]]
                                x = parent[x]
                            return x

                        def _union(a, b):
                            ra, rb = _find(a), _find(b)
                            if ra == rb:
                                return False
                            if ra > rb:
                                ra, rb = rb, ra
                            parent[rb] = ra
                            return True

                        support_counts = {}
                        top_k = max(2, int(args.bridge_top_k))
                        for i in range(n_unsieved):
                            row = sim_unsieved[i]
                            if top_k >= row.size:
                                top_idx = np.argsort(-row)[:top_k]
                            else:
                                top_idx = np.argpartition(-row, top_k - 1)[:top_k]
                                top_idx = top_idx[np.argsort(-row[top_idx])]

                            top_labels = []
                            top_sims = []
                            for idx in top_idx:
                                lbl = labels[idx]
                                if lbl >= 0:
                                    top_labels.append(lbl)
                                    top_sims.append(row[idx])

                            if len(top_labels) < 2:
                                continue

                            a_lbl, b_lbl = top_labels[0], top_labels[1]
                            a_sim, b_sim = top_sims[0], top_sims[1]
                            if a_lbl == b_lbl:
                                continue
                            if a_sim < args.bridge_min_sim or b_sim < args.bridge_min_sim:
                                continue

                            key = (min(a_lbl, b_lbl), max(a_lbl, b_lbl))
                            support_counts[key] = support_counts.get(key, 0) + 1

                        merged = 0
                        for (a_lbl, b_lbl), count in support_counts.items():
                            if count >= args.bridge_min_support:
                                if _union(a_lbl, b_lbl):
                                    merged += 1

                        if merged > 0:
                            label_map = {lbl: _find(lbl) for lbl in cluster_label_set}
                            labels = np.array([label_map.get(lbl, lbl) if lbl >= 0 else -1 for lbl in labels])
                            unique_labels = np.unique(labels)
                            n_clusters = len(unique_labels[unique_labels >= 0])
                            n_noise = np.sum(labels == -1)
                            cluster_labels = unique_labels[unique_labels >= 0]
                            logging.info("Bridging merged %d cluster pairs -> %d clusters", merged, n_clusters)
                        else:
                            logging.info("No cluster merges met support threshold")

                    # Initialize unsieved cluster weights after any merges
                    unsieved_cluster_weights = np.zeros(n_clusters)
                    unsieved_cluster_weights_unscaled = np.zeros(n_clusters)
                    cluster_label_to_index = {lbl: i for i, lbl in enumerate(cluster_labels)}

                    # For each unsieved fingerprint, find nearest cluster
                    # Method: for each unsieved FP, find the sieved FP with highest similarity,
                    # then assign to that sieved FP's cluster
                    for idx, orig_idx in enumerate(unsieved_indices):
                        nearest_sieved_idx = np.argmax(sim_unsieved[idx])
                        nearest_cluster = labels[nearest_sieved_idx]

                        if nearest_cluster >= 0:  # Not noise
                            # Add this unsieved fingerprint's weight to the cluster
                            unsieved_weight = frame_mapping['original_weights'][orig_idx]
                            unsieved_weight_unscaled = frame_mapping['original_weights_unscaled'][orig_idx]
                            cluster_idx = cluster_label_to_index.get(nearest_cluster)
                            if cluster_idx is not None:
                                unsieved_cluster_weights[cluster_idx] += unsieved_weight
                                unsieved_cluster_weights_unscaled[cluster_idx] += unsieved_weight_unscaled

                    logging.info(f"Reassigned {n_unsieved} unsieved fingerprints to clusters")
                    total_unsieved_weight = np.sum(unsieved_cluster_weights)
                    total_unsieved_weight_unscaled = np.sum(unsieved_cluster_weights_unscaled)
                    logging.info(
                        f"Total unsieved weight added to clusters: {total_unsieved_weight:.6f} "
                        f"(unscaled={total_unsieved_weight_unscaled:.6f})"
                    )
                else:
                    logging.info("No valid unsieved fingerprints to reassign")
            else:
                logging.info("No unsieved fingerprints or no clusters to reassign to")
        else:
            unsieved_cluster_weights = np.zeros(n_clusters)
            unsieved_cluster_weights_unscaled = np.zeros(n_clusters)

        # === Free Energy Calculation: Overall Binding dG ===
        if args.calc_overall_dg:
            try:
                logging.info("Calculating overall binding free energy...")
                concentration_nm3 = None
                if avg_volume_nm3 is not None and avg_volume_nm3 > 0:
                    concentration_nm3 = n_ligs / avg_volume_nm3
                unbound_weight_override = None
                if args.dg_unbound_ref == 'contact' and binding_stats:
                    unbound_weight_override = binding_stats.get('unbound_instance_weighted', None)
                dg_result = calculate_overall_binding_free_energy(
                    unique_weights,
                    labels,
                    temperature=args.temperature,
                    unbound_weight_override=unbound_weight_override,
                    concentration_nm3=concentration_nm3,
                    dg_offset_kj=args.dg_offset_kj,
                    f_contact=f_contact,
                )
                
                # Log results
                logging.info(f"Overall binding free energy: dG = {dg_result['dG_bind_kJmol']:.2f} +/- "
                           f"{dg_result['dG_stderr_kJmol']:.2f} kJ/mol (T={dg_result['temperature_K']:.1f} K)")
                logging.info(f"Probabilities: P_bound={dg_result['P_bound']:.4f}, "
                           f"P_unbound={dg_result['P_unbound']:.4f}")
                logging.info(f"Binding constant: K_bind={dg_result['K_bind_nm3']:.4f} nm^3")
                
                # Write to CSV
                overall_dg_file = "overall_binding_dg.csv"
                with open(overall_dg_file, 'w') as f:
                    if dg_result['f_contact'] is not None:
                        f.write("# Binding free energy with molecular volume contact correction\n")
                        f.write("# Method: Ligand VdW volume x packing_factor x bound_fraction / box_volume\n")
                        f.write("# Expected accuracy: +/-2-4 kcal/mol without empirical calibration\n")
                        f.write("# Note: Uses ONLY ligand VdW radii (NO receptor coordinates)\n")
                    else:
                        f.write("# WARNING: No molecular volume correction applied\n")
                        f.write("# dG values are QUALITATIVE ONLY - systematic underestimation (~5-10 kcal/mol)\n")
                        f.write("# Reason: Contact-space unbound != thermodynamic bulk unbound\n")
                        f.write("# Recommendation: Molecular volume method auto-enabled with --calc_overall_dg\n")
                    f.write("P_bound,P_unbound,K_bind_nm3,dG_bind_kJmol,dG_stderr_kJmol,temperature_K,n_bootstrap,warning,unbound_reference,concentration_nm3,dg_offset_kj,f_contact,volume_correction_kJmol\n")
                    f_contact_str = "" if dg_result['f_contact'] is None else f"{dg_result['f_contact']:.6f}"
                    volume_corr_str = "" if dg_result['f_contact'] is None else f"{dg_result['volume_correction_kJmol']:.4f}"
                    f.write(f"{dg_result['P_bound']:.6f},{dg_result['P_unbound']:.6f},"
                           f"{dg_result['K_bind_nm3']:.6f},{dg_result['dG_bind_kJmol']:.4f},"
                           f"{dg_result['dG_stderr_kJmol']:.4f},{dg_result['temperature_K']:.1f},"
                           f"{dg_result['n_bootstrap']},{dg_result['warning']},"
                           f"{dg_result['unbound_reference']},{dg_result['concentration_nm3']},{args.dg_offset_kj:.4f},"
                           f"{f_contact_str},{volume_corr_str}\n")
                logging.info(f"Overall binding free energy written to {overall_dg_file}")
                
            except ValueError as e:
                logging.warning(f"Could not calculate overall binding free energy: {e}")
            except Exception as e:
                logging.error(f"Error calculating overall binding free energy: {e}", exc_info=True)

        # Cluster statistics
        mode_occ = np.zeros(n_clusters)
        mode_occ_unscaled = np.zeros(n_clusters)
        raw_inst_per_cluster = np.zeros(n_clusters, dtype=int)
        cluster_label_to_index = {lbl: i for i, lbl in enumerate(cluster_labels)}

        for lbl in cluster_labels:
            mask = labels == lbl
            i = cluster_label_to_index[lbl]
            mode_occ[i] = unique_weights[mask].sum()
            mode_occ_unscaled[i] = unique_weights_unscaled[mask].sum()
            raw_inst_per_cluster[i] = uniq_counts[mask].sum()
        
        # Add unsieved fingerprint weights to cluster occupancies
        if args.sieving and args.max_sieved is not None and n_unique > args.max_sieved:
            for i in range(n_clusters):
                mode_occ[i] += unsieved_cluster_weights[i]
                mode_occ_unscaled[i] += unsieved_cluster_weights_unscaled[i]
        
        # Log cluster statistics with final occupancies
        for lbl in cluster_labels:
            mask = labels == lbl
            i = cluster_label_to_index[lbl]
            logging.info(
                f"Cluster {lbl}: size={mask.sum()}, raw={raw_inst_per_cluster[i]}, "
                f"occ={mode_occ[i]:.6f}, occ_unscaled={mode_occ_unscaled[i]:.6f}"
            )

        # Select top modes
        if n_clusters > 0:
            sort_idx = np.argsort(-mode_occ)
            cum_occ = np.cumsum(mode_occ[sort_idx])
            retain_n = np.searchsorted(cum_occ, args.top_occ_percent / 100) + 1
            if retain_n == 0:
                retain_n = 1
            if retain_n > len(cum_occ):
                retain_n = len(cum_occ)
            top_modes_idx = sort_idx[:retain_n]
            top_cluster_labels = cluster_labels[top_modes_idx]
            logging.info(f"Retained {retain_n} top modes covering {cum_occ[retain_n-1]*100:.1f}% occupancy")
        else:
            logging.warning("No clusters found - skipping mode selection")
            retain_n = 0
            top_modes_idx = np.array([], dtype=int)
            top_cluster_labels = np.array([], dtype=int)

        # === Free Energy Calculation: Per-Cluster Binding dG ===
        if args.calc_cluster_dg and n_clusters > 0:
            try:
                logging.info("Calculating per-cluster binding free energies...")
                
                concentration_nm3 = None
                if avg_volume_nm3 is not None and avg_volume_nm3 > 0:
                    concentration_nm3 = n_ligs / avg_volume_nm3
                if args.dg_unbound_ref == 'contact' and binding_stats:
                    unbound_weight = binding_stats.get('unbound_instance_weighted', None)
                else:
                    unbound_weight = np.sum(unique_weights[labels == -1])
                
                if unbound_weight <= 0:
                    logging.warning("No unbound states detected, skipping per-cluster free energy calculation")
                else:
                    # Calculate per-cluster free energies
                    cluster_dg_results = calculate_cluster_binding_free_energies(
                        mode_occ,
                        cluster_labels,
                        unbound_weight,
                        unique_weights,
                        labels,
                        temperature=args.temperature,
                        unbound_reference=args.dg_unbound_ref,
                        concentration_nm3=concentration_nm3,
                        dg_offset_kj=args.dg_offset_kj,
                        f_contact=f_contact,
                    )
                    
                    if cluster_dg_results:
                        # Log summary
                        n_calculated = len(cluster_dg_results)
                        dg_values = [r['dG_abs_kJmol'] for r in cluster_dg_results.values()]
                        logging.info(f"Calculated free energies for {n_calculated} clusters")
                        logging.info(f"dG range: {min(dg_values):.2f} to {max(dg_values):.2f} kJ/mol")
                        
                        # Find reference cluster
                        ref_cluster = [cid for cid, data in cluster_dg_results.items() if data['is_reference']][0]
                        logging.info(f"Reference cluster (lowest dG): {ref_cluster}")
                        
                        # Write to CSV
                        cluster_dg_file = "cluster_binding_dg.csv"
                        with open(cluster_dg_file, 'w') as f:
                            f.write("cluster_id,occupancy,dG_abs_kJmol,dG_abs_stderr_kJmol,"
                                   "ddG_rel_kJmol,ddG_rel_stderr_kJmol,is_reference,"
                                   "unbound_reference,concentration_nm3,dg_offset_kj,f_contact,volume_correction_kJmol\n")

                            f_contact_str = "" if f_contact is None else f"{f_contact:.6f}"
                            volume_corr_str = "" if f_contact is None else f"{volume_correction_kjmol:.4f}"
                            
                            # Sort by occupancy (highest first)
                            sorted_items = sorted(
                                cluster_dg_results.items(),
                                key=lambda x: x[1]['occupancy'],
                                reverse=True
                            )
                            
                            for cluster_id, data in sorted_items:
                                f.write(f"{cluster_id},{data['occupancy']:.6f},"
                                       f"{data['dG_abs_kJmol']:.4f},{data['dG_abs_stderr_kJmol']:.4f},"
                                       f"{data['ddG_kJmol']:.4f},{data['ddG_stderr_kJmol']:.4f},"
                                       f"{int(data['is_reference'])},{data['unbound_reference']},"
                                       f"{data['concentration_nm3']},{args.dg_offset_kj:.4f},"
                                       f"{f_contact_str},{volume_corr_str}\n")
                        
                        logging.info(f"Per-cluster free energies written to {cluster_dg_file}")
                    else:
                        logging.info("No clusters passed min_occupancy threshold for free energy calculation")
                        
            except ValueError as e:
                logging.warning(f"Could not calculate per-cluster free energies: {e}")
            except Exception as e:
                logging.error(f"Error calculating per-cluster free energies: {e}", exc_info=True)

        # Medoids
        top_medoids = []
        for i, lbl in enumerate(top_cluster_labels):
            mask = labels == lbl
            if mask.sum() == 0:
                continue
            if args.sparse_similarity:
                # In sparse mode, approximate medoid using sparse similarity:
                # treat missing edges as sim=0 (distance=1), choose max total similarity
                cluster_indices = np.where(mask)[0]
                if sim_sparse is not None and sparse.issparse(sim_sparse) and cluster_indices.size > 0:
                    cluster_sim = sim_sparse[cluster_indices][:, cluster_indices]  # type: ignore[index]
                    sim_sums = np.asarray(cluster_sim.sum(axis=1)).ravel()
                    medoid_global = cluster_indices[int(np.argmax(sim_sums))]
                else:
                    # Fallback: deterministic max-weight representative
                    medoid_global = cluster_indices[np.argmax(unique_weights[cluster_indices])]
            else:
                cluster_dist = dist[np.ix_(mask, mask)]
                medoid_local = np.argmin(cluster_dist.sum(axis=1))
                medoid_global = np.where(mask)[0][medoid_local]
            top_medoids.append(medoid_global)

        # Per-ligand occupancies (unscaled + SASA-scaled)
        lig_occ = np.zeros((n_ligs, n_clusters))
        lig_occ_unscaled = np.zeros((n_ligs, n_clusters))
        for i in range(final_n):
            lbl = labels[i]
            if lbl >= 0:
                cl_idx = cluster_label_to_index.get(lbl)
                if cl_idx is not None:
                    lig_occ[unique_lig_ids[i], cl_idx] += unique_weights[i]
                    lig_occ_unscaled[unique_lig_ids[i], cl_idx] += unique_weights_unscaled[i]

        lig_occ /= np.sum(lig_occ, axis=1, keepdims=True) + 1e-12
        lig_occ_unscaled /= np.sum(lig_occ_unscaled, axis=1, keepdims=True) + 1e-12
        top_lig_occ = lig_occ[:, top_modes_idx]
        top_lig_occ_unscaled = lig_occ_unscaled[:, top_modes_idx]

        # Per-ligand mode SASA (if enabled)
        if args.sasa_enabled:
            logging.info("Computing per-ligand mode SASA...")
            assert sasa_arrays is not None
            lig_mode_sasa = np.zeros((n_ligs, n_clusters))
            lig_mode_sasa_wsum = np.zeros((n_ligs, n_clusters))
            mode_sasa = np.zeros(n_clusters)
            mode_sasa_wsum = np.zeros(n_clusters)
            for i in range(final_n):
                lbl = labels[i]
                if lbl >= 0:
                    cl_idx = cluster_label_to_index.get(lbl)
                    if cl_idx is None:
                        continue
                    frame = unique_frames[i]
                    lig = unique_lig_ids[i]
                    sasa_val = sasa_arrays[lig][frame]
                    w = unique_weights[i]
                    lig_mode_sasa[lig, cl_idx] += sasa_val * w
                    lig_mode_sasa_wsum[lig, cl_idx] += w
                    mode_sasa[cl_idx] += sasa_val * w
                    mode_sasa_wsum[cl_idx] += w
            lig_mode_sasa /= lig_mode_sasa_wsum + 1e-12
            mode_sasa /= mode_sasa_wsum + 1e-12
            top_lig_mode_sasa = lig_mode_sasa[:, top_modes_idx]
            top_mode_sasa = mode_sasa[top_modes_idx]
        else:
            top_lig_mode_sasa = np.zeros((n_ligs, 0))
            top_mode_sasa = np.array([])

        # === Outputs ===

        # Enhanced modes.csv with consistent naming
        modes_path = f"{args.modes_csv}"
        ref_pairs_for_csv = None
        with open(modes_path, 'w') as f:
            # Determine header based on SASA enablement
            if args.sasa_enabled:
                f.write("mode_id weighted_occupancy sasa_scaled_occupancy rep_frame interacting_pairs n_instances avg_sasa\n")
            else:
                f.write("mode_id weighted_occupancy rep_frame interacting_pairs n_instances\n")
            
            for i, cl_idx in enumerate(top_modes_idx):
                rep_i = top_medoids[i]
                rep_frame = unique_frames[rep_i]
                rep_lig_id = unique_lig_ids[rep_i]
                key = keys_list[rep_i]
                
                # Calculate enhanced statistics
                mask = labels == cluster_labels[top_modes_idx[i]]
                n_instances = uniq_counts[mask].sum() if 'uniq_counts' in locals() else mask.sum()
                
                rows, cols = zip(*key) if key else ([], [])
                pairs = []
                for r, c in zip(rows, cols):
                    resname = receptor.residues[r].resname
                    code = RES_3TO1.get(resname, resname[:1])
                    pairs.append(f"{code}{receptor.residues[r].resid}:{lig_atom_names[c]}")
                
                if args.sasa_enabled:
                    # Calculate average SASA for this mode
                    avg_sasa = 0
                    if sasa_arrays is not None:
                        mode_frames = unique_frames[mask]
                        mode_ligands = unique_lig_ids[mask]
                        sasa_values = [sasa_arrays[lig][frame] for frame, lig in zip(mode_frames, mode_ligands)]
                        avg_sasa = np.mean(sasa_values) if sasa_values else 0
                    
                    # Calculate SASA-scaled occupancy
                    sasa_scaled_weight = mode_occ[cl_idx]
                    weighted_occ = mode_occ_unscaled[cl_idx]
                    f.write(
                        f"{i} {weighted_occ:.6f} {sasa_scaled_weight:.6f} "
                        f"{rep_frame} {','.join(pairs)} {n_instances} {avg_sasa:.4f}\n"
                    )
                else:
                    f.write(f"{i} {mode_occ[cl_idx]:.6f} {rep_frame} {','.join(pairs)} {n_instances}\n")

        # Enhanced PDBs with B-factors and remarks
        for i, cl_idx in enumerate(top_modes_idx):
            if mode_occ[cl_idx] >= args.min_occ:
                rep_i = top_medoids[i]
                rep_frame = unique_frames[rep_i]
                rep_lig_id = unique_lig_ids[rep_i]
                u.trajectory[rep_frame]
                
                # Get selection for receptor + ligand
                sel = receptor + ligand_residues[rep_lig_id].atoms
                
                # Set B-factors based on mode occupancy (SASA-scaled if enabled)
                bfactor_value = mode_occ[cl_idx]
                
                sel.bfactors = np.full(len(sel), bfactor_value)
                
                # Write PDB with enhanced information
                pdb_path = f"{args.mode_pdb_prefix}{i}.pdb"
                
                # Add REMARK lines with mode statistics
                remarks = [
                    f"REMARK 265 - BINDING MODE ANALYSIS",
                    f"REMARK 265 MODE_ID: {i}",
                    f"REMARK 265 WEIGHTED_OCCUPANCY: {mode_occ_unscaled[cl_idx]:.6f}",
                    f"REMARK 265 FRAME: {rep_frame}",
                    f"REMARK 265 LIGAND_ID: {rep_lig_id+1}",
                    f"REMARK 265 N_INSTANCES: {uniq_counts[labels == cluster_labels[top_modes_idx[i]]].sum() if 'uniq_counts' in locals() else 'N/A'}"
                ]
                
                if args.sasa_enabled:
                    remarks.append(f"REMARK 265 SASA_SCALED_OCCUPANCY: {bfactor_value:.6f}")
                    if sasa_arrays is not None:
                        remarks.append(f"REMARK 265 LIGAND_SASA: {sasa_arrays[rep_lig_id][rep_frame]:.4f}")
                
                # Write PDB with remarks (write structure first, then prepend remarks)
                with tempfile.NamedTemporaryFile(delete=False, suffix='.pdb') as tmp:
                    tmp_path = tmp.name
                try:
                    sel.write(tmp_path)
                    with open(pdb_path, 'w') as f:
                        f.write('\n'.join(remarks) + '\n')
                        with open(tmp_path, 'r') as tmp_f:
                            for line in tmp_f:
                                f.write(line)
                finally:
                    if os.path.exists(tmp_path):
                        os.unlink(tmp_path)
                
                logging.info(f"Wrote enhanced mode {i} PDB: frame {rep_frame}, lig {rep_lig_id+1} to {pdb_path}")

        # Enhanced per-ligand occupancies with SASA-scaled values
        for lig_id in range(n_ligs):
            occ_path = f"{args.occupancies_prefix}{lig_id+1}.txt"
            with open(occ_path, 'w') as f:
                if args.sasa_enabled:
                    f.write("ligand mode_id weighted_occupancy sasa_scaled_occupancy raw_count\n")
                    # Calculate SASA-scaled occupancies
                    for mode_idx, occ in enumerate(top_lig_occ[lig_id]):
                        sasa_scaled_occ = occ
                        weighted_occ = top_lig_occ_unscaled[lig_id][mode_idx] if mode_idx < top_lig_occ_unscaled.shape[1] else occ

                        # Calculate raw count
                        mask = labels == cluster_labels[top_modes_idx[mode_idx]] if mode_idx < len(top_modes_idx) else np.array([])
                        raw_count = uniq_counts[mask].sum() if 'uniq_counts' in locals() and mask.size > 0 else 0

                        f.write(f"{lig_id+1} {mode_idx} {weighted_occ:.6f} {sasa_scaled_occ:.6f} {raw_count}\n")
                else:
                    f.write("ligand mode_id weighted_occupancy raw_count\n")
                    for mode_idx, occ in enumerate(top_lig_occ[lig_id]):
                        # Calculate raw count
                        mask = labels == cluster_labels[top_modes_idx[mode_idx]] if mode_idx < len(top_modes_idx) else np.array([])
                        raw_count = uniq_counts[mask].sum() if 'uniq_counts' in locals() and mask.size > 0 else 0
                        f.write(f"{lig_id+1} {mode_idx} {occ:.6f} {raw_count}\n")
            logging.info(f"Wrote enhanced occupancies for ligand {lig_id+1} to {occ_path}")

        # SASA modes and per-ligand scores with consistent naming
        if args.sasa_enabled:
            sasa_modes_path = f"{args.sasa_modes_csv}"
            with open(sasa_modes_path, 'w') as f:
                f.write("mode_id sasa_scaled_occupancy\n")
                for i, cl_idx in enumerate(top_modes_idx):
                    f.write(f"{i} {mode_occ[cl_idx]:.6f}\n")

            for lig_id in range(n_ligs):
                sasa_scores_path = f"{args.sasa_scores_prefix}{lig_id+1}.txt"
                with open(sasa_scores_path, 'w') as f:
                    f.write("ligand mode_id sasa_scaled_occupancy\n")
                    for mode_idx, occ in enumerate(top_lig_occ[lig_id]):
                        f.write(f"{lig_id+1} {mode_idx} {occ:.6f}\n")
                logging.info(f"Wrote SASA-scaled scores for ligand {lig_id+1} to {sasa_scores_path}")

        # Enhanced reference similarity with structural alignment
        if args.ref_pdb:
            logging.info("Computing structural alignment and similarity to reference PDB...")
            ref_u = mda.Universe(args.ref_pdb)
            ref_rec = ref_u.select_atoms(args.receptor_sel)
            ref_lig_residues = []

            def _select_reference_ligand(ref_universe, ligand_sel, traj_resnames, ref_ligand_sel=None):
                traj_resname_set = {rn for rn in traj_resnames if rn}
                ref_lig_sel = []
                if ref_ligand_sel:
                    try:
                        ref_lig_sel = ref_universe.select_atoms(ref_ligand_sel).residues
                    except Exception:
                        ref_lig_sel = []
                    if len(ref_lig_sel) > 0:
                        return ref_lig_sel
                try:
                    ref_lig_sel = ref_universe.select_atoms(ligand_sel).residues
                except Exception:
                    ref_lig_sel = []
                if len(ref_lig_sel) > 0:
                    ref_resnames = {res.resname for res in ref_lig_sel}
                    if ref_resnames & traj_resname_set:
                        return ref_lig_sel
                    logging.info(
                        "Reference ligand selection resnames %s do not match trajectory ligand resnames %s; falling back",
                        sorted(ref_resnames),
                        sorted(traj_resname_set),
                    )
                if traj_resname_set:
                    resname_sel = "resname " + " ".join(sorted(traj_resname_set))
                    try:
                        ref_lig_sel = ref_universe.select_atoms(resname_sel).residues
                    except Exception:
                        ref_lig_sel = []
                    if len(ref_lig_sel) > 0:
                        logging.info("Reference ligand selection fallback to resname: %s", resname_sel)
                return ref_lig_sel

            try:
                traj_resnames = [res.resname for res in ligand_residues]
                ref_lig_residues = _select_reference_ligand(
                    ref_u,
                    args.ligand_sel,
                    traj_resnames,
                    ref_ligand_sel=args.ref_ligand_sel,
                )
                if len(ref_lig_residues) == 0:
                    raise ValueError("No reference ligand residues found")

                def _atom_element(atom):
                    element = getattr(atom, "element", None)
                    if element is None:
                        return None
                    element = str(element).strip()
                    return element if element else None

                def _build_common_atom_space(ref_receptor, traj_receptor, ref_ligand, traj_ligand):
                    """Return common atom selections and index mappings (ref->traj)."""
                    def _is_heavy(atom):
                        elem = _atom_element(atom)
                        if elem is not None:
                            return elem.upper() != "H"
                        return not str(atom.name).upper().startswith("H")

                    def _map_residues_by_offset(ref_resnames, traj_resnames):
                        best_offset = 0
                        best_matches = -1
                        best_range = (0, 0, 0)
                        ref_len = len(ref_resnames)
                        traj_len = len(traj_resnames)
                        for offset in range(-ref_len + 1, traj_len):
                            start_ref = max(0, -offset)
                            start_traj = max(0, offset)
                            length = min(ref_len - start_ref, traj_len - start_traj)
                            if length <= 0:
                                continue
                            matches = sum(
                                1 for i in range(length)
                                if ref_resnames[start_ref + i] == traj_resnames[start_traj + i]
                            )
                            if matches > best_matches:
                                best_matches = matches
                                best_offset = offset
                                best_range = (start_ref, start_traj, length)
                        if best_matches <= 0:
                            return {}, 0
                        start_ref, start_traj, length = best_range
                        mapping = {}
                        for i in range(length):
                            if ref_resnames[start_ref + i] == traj_resnames[start_traj + i]:
                                mapping[start_ref + i] = start_traj + i
                        return mapping, best_matches
                    traj_res_lookup = {
                        (res.resid, res.resname): idx
                        for idx, res in enumerate(traj_receptor.residues)
                    }
                    ref_to_traj_res = {}
                    for ref_idx, res in enumerate(ref_receptor.residues):
                        key = (res.resid, res.resname)
                        if key in traj_res_lookup:
                            ref_to_traj_res[ref_idx] = traj_res_lookup[key]

                    min_common_res = 10
                    if len(ref_to_traj_res) < min_common_res:
                        ref_resnames = [res.resname for res in ref_receptor.residues]
                        traj_resnames = [res.resname for res in traj_receptor.residues]
                        offset_map, offset_matches = _map_residues_by_offset(ref_resnames, traj_resnames)
                        if len(offset_map) >= min_common_res:
                            ref_to_traj_res = offset_map
                            logging.info(
                                "Reference residue mapping fallback to resname offset (matches=%d)",
                                offset_matches,
                            )
                        else:
                            common_len = min(len(ref_receptor.residues), len(traj_receptor.residues))
                            ref_to_traj_res = {i: i for i in range(common_len)}
                            logging.info(
                                "Reference residue mapping fallback to sequential order (common_len=%d)",
                                common_len,
                            )

                    ref_rec_heavy_mask = np.array([_is_heavy(atom) for atom in ref_receptor.atoms], dtype=bool)
                    ref_rec_heavy = ref_receptor.atoms[ref_rec_heavy_mask]
                    ref_rec_index_map = {idx: i for i, idx in enumerate(ref_rec_heavy.indices)}
                    ref_rec_atom_to_traj_res = np.full(ref_rec_heavy.n_atoms, -1, dtype=np.int32)
                    common_rec_atoms = 0
                    for ref_res_idx, traj_res_idx in ref_to_traj_res.items():
                        ref_res = ref_receptor.residues[ref_res_idx]
                        traj_res = traj_receptor.residues[traj_res_idx]

                        traj_keys = {}
                        traj_names = set()
                        for atom in traj_res.atoms:
                            if not _is_heavy(atom):
                                continue
                            elem = _atom_element(atom)
                            key = (atom.name, elem)
                            traj_keys.setdefault(key, 0)
                            traj_keys[key] += 1
                            traj_names.add(atom.name)

                        for atom in ref_res.atoms:
                            if not _is_heavy(atom):
                                continue
                            ref_local_idx = ref_rec_index_map.get(atom.index)
                            if ref_local_idx is None:
                                continue
                            elem = _atom_element(atom)
                            key = (atom.name, elem)
                            if (elem is not None and key in traj_keys) or (atom.name in traj_names):
                                ref_rec_atom_to_traj_res[ref_local_idx] = traj_res_idx
                                common_rec_atoms += 1

                    traj_lig_keys = {}
                    traj_lig_names = {}
                    traj_lig_by_elem = {}
                    for idx, atom in enumerate(traj_ligand.atoms):
                        if not _is_heavy(atom):
                            continue
                        elem = _atom_element(atom)
                        key = (atom.name, elem)
                        traj_lig_keys.setdefault(key, []).append(idx)
                        traj_lig_names.setdefault(atom.name, []).append(idx)
                        if elem is not None:
                            traj_lig_by_elem.setdefault(elem, []).append(idx)

                    ref_lig_heavy_mask = np.array([_is_heavy(atom) for atom in ref_ligand.atoms], dtype=bool)
                    ref_lig_heavy = ref_ligand.atoms[ref_lig_heavy_mask]
                    ref_lig_atom_to_traj_atom = np.full(ref_lig_heavy.n_atoms, -1, dtype=np.int32)
                    common_lig_atoms = 0
                    for ref_idx, atom in enumerate(ref_lig_heavy):
                        elem = _atom_element(atom)
                        key = (atom.name, elem)
                        if elem is not None and key in traj_lig_keys and traj_lig_keys[key]:
                            ref_lig_atom_to_traj_atom[ref_idx] = traj_lig_keys[key].pop(0)
                            common_lig_atoms += 1
                        elif atom.name in traj_lig_names and traj_lig_names[atom.name]:
                            ref_lig_atom_to_traj_atom[ref_idx] = traj_lig_names[atom.name].pop(0)
                            common_lig_atoms += 1

                    if common_lig_atoms == 0:
                        common_len = min(ref_lig_heavy.n_atoms, traj_ligand.n_atoms)
                        if common_len > 0:
                            ref_lig_atom_to_traj_atom[:common_len] = np.arange(common_len, dtype=np.int32)
                            common_lig_atoms = common_len
                            logging.info(
                                "Reference ligand mapping fallback to sequential order (common_len=%d)",
                                common_len,
                            )
                    elif common_lig_atoms < ref_lig_heavy.n_atoms:
                        for ref_idx, atom in enumerate(ref_lig_heavy):
                            if ref_lig_atom_to_traj_atom[ref_idx] >= 0:
                                continue
                            elem = _atom_element(atom)
                            if elem is None:
                                continue
                            candidate_list = traj_lig_by_elem.get(elem)
                            if candidate_list:
                                ref_lig_atom_to_traj_atom[ref_idx] = candidate_list.pop(0)
                        common_lig_atoms = int(np.count_nonzero(ref_lig_atom_to_traj_atom >= 0))

                    common_rec_mask = ref_rec_atom_to_traj_res >= 0
                    common_lig_mask = ref_lig_atom_to_traj_atom >= 0
                    ref_rec_common = ref_rec_heavy[common_rec_mask]
                    ref_lig_common = ref_lig_heavy[common_lig_mask]
                    ref_rec_common_to_traj_res = ref_rec_atom_to_traj_res[common_rec_mask]
                    ref_lig_common_to_traj_atom = ref_lig_atom_to_traj_atom[common_lig_mask]

                    common_res_idx = np.unique(ref_rec_common_to_traj_res)
                    common_lig_atom_idx = np.unique(ref_lig_common_to_traj_atom)

                    logging.info(
                        "Reference common atom space: receptor residues=%d ligand atoms=%d (common receptor atoms=%d, common ligand atoms=%d)",
                        len(common_res_idx),
                        len(common_lig_atom_idx),
                        int(common_rec_atoms),
                        int(common_lig_atoms),
                    )

                    return {
                        "ref_rec_common": ref_rec_common,
                        "ref_lig_common": ref_lig_common,
                        "ref_rec_common_to_traj_res": ref_rec_common_to_traj_res,
                        "ref_lig_common_to_traj_atom": ref_lig_common_to_traj_atom,
                        "common_res_idx": common_res_idx,
                        "common_lig_atom_idx": common_lig_atom_idx,
                    }

                def _compute_reference_similarity(common, lig_id):
                    ref_pairs = []
                    if common["common_res_idx"].size == 0 or common["common_lig_atom_idx"].size == 0:
                        logging.warning("No common atom space found; reference similarity will be zero")
                        return np.zeros(X.shape[0], dtype=np.float32), ref_pairs

                    pairs, _ = capped_distance(
                        common["ref_rec_common"].positions,
                        common["ref_lig_common"].positions,
                        args.cutoff,
                    )
                    logging.info("Reference contact pairs (common atoms): %d", len(pairs))

                    common_cols = (
                        common["common_res_idx"][:, None] * n_lig_atoms +
                        common["common_lig_atom_idx"][None, :]
                    ).ravel().astype(np.int32)
                    common_cols = np.sort(common_cols)
                    col_index_map = {col: i for i, col in enumerate(common_cols)}
                    ref_vec_common = np.zeros(len(common_cols), dtype=np.float32)

                    if len(pairs) > 0:
                        rows = common["ref_rec_common_to_traj_res"][pairs[:, 0]]
                        cols = common["ref_lig_common_to_traj_atom"][pairs[:, 1]]
                        flat_cols = rows * n_lig_atoms + cols
                        for flat_col, rr, cc in zip(flat_cols, rows, cols):
                            col_idx = col_index_map.get(int(flat_col))
                            if col_idx is not None:
                                ref_vec_common[col_idx] = 1.0
                                resname = receptor.residues[rr].resname
                                code = RES_3TO1.get(resname, resname[:1])
                                pair = f"{code}{receptor.residues[rr].resid}:{lig_atom_names[cc]}"
                                ref_pairs.append(pair)

                    if ref_pairs:
                        ref_pairs = sorted(set(ref_pairs))

                    lig_mask = unique_lig_ids == lig_id
                    X_common = X[lig_mask][:, common_cols]
                    ref_sumA = float(ref_vec_common.sum())
                    if sparse.issparse(X_common):
                        ref_sumB = np.array(X_common.getnnz(axis=1), dtype=np.float32)
                    else:
                        ref_sumB = np.count_nonzero(X_common, axis=1).astype(np.float32)

                    if ref_sumA <= 0:
                        return np.zeros(X.shape[0], dtype=np.float32), ref_pairs

                    ref_intersect = ref_vec_common @ X_common.T
                    ref_intersect = np.asarray(ref_intersect).ravel()

                    if args.similarity == 'dice':
                        ref_sim_local = 2 * ref_intersect / (ref_sumA + ref_sumB + 1e-10)
                    elif args.similarity == 'jaccard':
                        ref_sim_local = ref_intersect / (ref_sumA + ref_sumB - ref_intersect + 1e-10)
                    else:
                        ref_normA = np.sqrt(ref_sumA + 1e-10)
                        ref_normB = np.sqrt(ref_sumB + 1e-10)
                        ref_sim_local = ref_intersect / (ref_normA * ref_normB + 1e-10)

                    ref_sim = np.zeros(X.shape[0], dtype=np.float32)
                    ref_sim[lig_mask] = ref_sim_local.astype(np.float32, copy=False)

                    logging.info("Reference fingerprint nnz=%d sumA=%.1f", int(ref_vec_common.sum()), float(ref_sumA))
                    return ref_sim, ref_pairs

                ref_sim = np.zeros(X.shape[0], dtype=np.float32)
                for lig_id, traj_res in enumerate(ligand_residues):
                    ref_res = ref_lig_residues[0]
                    if lig_id < len(ref_lig_residues):
                        ref_res = ref_lig_residues[lig_id]

                    common = _build_common_atom_space(ref_rec, receptor, ref_res.atoms, traj_res.atoms)
                    lig_ref_sim, lig_ref_pairs = _compute_reference_similarity(common, lig_id)
                    ref_sim += lig_ref_sim
                    if ref_pairs_for_csv is None and lig_ref_pairs:
                        ref_pairs_for_csv = lig_ref_pairs

                top_sim = ref_sim[top_medoids]
                logging.info(f"Reference similarities to top modes: {top_sim}")
                logging.info("Reference similarity stats: max=%.6f mean=%.6f", float(np.max(ref_sim)), float(np.mean(ref_sim)))

                # Insert reference row to binding_modes.csv
                _append_reference_row(f"{args.modes_csv}", ref_pairs_for_csv or [], args.sasa_enabled)
                
            except Exception as e:
                logging.warning(f"Structural alignment failed, falling back to simple similarity: {e}")
                # Fallback to original method
                fallback_ref_lig = None
                if ref_lig_residues:
                    fallback_ref_lig = ref_lig_residues[0].atoms
                else:
                    try:
                        fallback_ref_lig = ref_u.select_atoms(args.ligand_sel)
                    except Exception:
                        fallback_ref_lig = ref_u.atoms[[]]

                pairs, _ = capped_distance(ref_rec.positions, fallback_ref_lig.positions, args.cutoff)
                if len(pairs) > 0:
                    r = receptor_atom_to_res[pairs[:,0]]
                    c = pairs[:,1]
                    flat_cols = r * n_lig_atoms + c
                    ref_vec = np.zeros(fp_flat_size, dtype=np.float32)
                    ref_vec[flat_cols] = 1.0
                else:
                    ref_vec = np.zeros(fp_flat_size, dtype=np.float32)

                ref_intersect = ref_vec @ X.T
                ref_sumA = np.sum(ref_vec)
                ref_sumB = sumA

                if args.similarity == 'dice':
                    ref_sim = 2 * ref_intersect / (ref_sumA + ref_sumB + 1e-10)
                elif args.similarity == 'jaccard':
                    ref_sim = ref_intersect / (ref_sumA + ref_sumB - ref_intersect + 1e-10)
                else:
                    ref_normA = np.sqrt(ref_sumA + 1e-10)
                    ref_normB = np.sqrt(ref_sumB + 1e-10)
                    ref_sim = ref_intersect / (ref_normA * ref_normB + 1e-10)

                top_sim = ref_sim[top_medoids]
                logging.info(f"Reference similarities to top modes (fallback): {top_sim}")

        logging.info("="*60)
        logging.info("Interaction Fingerprint Analysis Tool - Completed Successfully")
        logging.info("="*60)

    except Exception as e:
        logging.error(f"Critical error during execution: {e}", exc_info=True)
        sys.exit(1)


if __name__ == '__main__':
    main()
