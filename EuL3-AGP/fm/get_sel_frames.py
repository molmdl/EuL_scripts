#!/usr/bin/env python

import argparse
import logging
import numpy as np
import pandas as pd
import scipy.interpolate
import scipy.ndimage
from scipy.signal import argrelextrema
import MDAnalysis as mda

# Defaults defined here (overridable by arguments)
DEFAULT_FES = 'fes.dat'
DEFAULT_TOP = 'v1.pdb'
DEFAULT_TRAJ = 'v1.xtc'
DEFAULT_COLVAR = 'trj.COLVAR'
DEFAULT_FE_MIN = 0.0
DEFAULT_FE_MAX = 5.0
DEFAULT_DELTA_FE = 1.0
HIGH_VALUE = 80.0  # kJ/mol, for inf/NaN/excluded regions
DEFAULT_FILTER_CVS = ['dm0.min', 'dm1.min', 'dm2.min']
DEFAULT_FILTER_MIN = 0.0
DEFAULT_FILTER_MAX = 0.6

def main():
    parser = argparse.ArgumentParser(description='Extract frames at FES minima basins with optional per-CV filtering')
    parser.add_argument('--fes', default=DEFAULT_FES, help='FES file (default: %(default)s)')
    parser.add_argument('--top', default=DEFAULT_TOP, help='Topology (default: %(default)s)')
    parser.add_argument('--traj', default=DEFAULT_TRAJ, help='Trajectory (default: %(default)s)')
    parser.add_argument('--colvar', default=DEFAULT_COLVAR, help='COLVAR file (default: %(default)s)')
    parser.add_argument('--fe_min', type=float, default=DEFAULT_FE_MIN, help='Min FE (default: %(default)s kJ/mol)')
    parser.add_argument('--fe_max', type=float, default=DEFAULT_FE_MAX, help='Max FE for minima (default: %(default)s kJ/mol)')
    parser.add_argument('--delta_fe', type=float, default=DEFAULT_DELTA_FE, help='Basin width in FE (default: %(default)s kJ/mol)')
    parser.add_argument('--cv_cols', nargs='+', required=True, help='CV column names for FES in order (e.g. d1 fps.ld)')
    parser.add_argument('--filter', action='store_true', help='Enable additional per-CV filtering with defaults')
    parser.add_argument('--filter_cvs', nargs='*', help='Override filter CV column names (default if --filter: dm0.min dm1.min dm2.min)')
    parser.add_argument('--filter_ranges', nargs='*', help='Override filter ranges as min-max strings (default: 0-0.7 each)')
    args = parser.parse_args()

    logging.basicConfig(level=logging.INFO,
                        format='%(asctime)s %(levelname)s: %(message)s')
    logger = logging.getLogger(__name__)

    try:
        logger.info('Starting FES minima extraction with optional filtering')
        logger.info(f'Inputs: FES={args.fes}, TOP={args.top}, TRAJ={args.traj}, COLVAR={args.colvar}')
        logger.info(f'CVs for FES: {" ".join(args.cv_cols)}')
        logger.info(f'FE range: {args.fe_min:.1f} to {args.fe_max:.1f} kJ/mol, delta={args.delta_fe:.1f} kJ/mol')
        logger.info('Non-periodic CVs assumed for simplicity.')

        # Filtering setup
        filter_cvs = None
        filter_mins = None
        filter_maxs = None
        if args.filter or args.filter_cvs is not None:
            if args.filter_cvs is None:
                filter_cvs = DEFAULT_FILTER_CVS[:]
            else:
                filter_cvs = args.filter_cvs
            n_filter = len(filter_cvs)
            if n_filter == 0:
                logger.info('Empty filter_cvs list – disabling additional filtering')
                filter_cvs = None
            else:
                if args.filter_ranges is not None and len(args.filter_ranges) != n_filter:
                    raise ValueError('Number of --filter_ranges must match number of filter CVs')
                if args.filter_ranges is None:
                    filter_mins = np.full(n_filter, DEFAULT_FILTER_MIN, dtype=float)
                    filter_maxs = np.full(n_filter, DEFAULT_FILTER_MAX, dtype=float)
                    range_display = '0-0.7'
                else:
                    filter_mins = np.empty(n_filter, dtype=float)
                    filter_maxs = np.empty(n_filter, dtype=float)
                    for j, rg in enumerate(args.filter_ranges):
                        mn, mx = map(float, rg.split('-'))
                        filter_mins[j] = mn
                        filter_maxs[j] = mx
                    range_display = ' '.join(args.filter_ranges)
                logger.info(f'Additional per-CV filtering enabled: CVs {" ".join(filter_cvs)} | ranges {range_display}')

        # Load FES and parse fields
        fields = None
        with open(args.fes) as f:
            for line in f:
                if line.startswith('#! FIELDS'):
                    fields = line.split()[2:]
                    break
        if fields is None:
            raise ValueError('No #! FIELDS line found in FES file')
        
        df_fes = pd.read_csv(args.fes, delim_whitespace=True, comment='#', names=fields)
        dim = len(args.cv_cols)
        if dim not in (1, 2):
            raise ValueError('Only 1D or 2D FES supported')
        
        fe_col = fields[dim]  # FE column immediately after CVs
        logger.info(f'Detected dimension: {dim}D, FE column: {fe_col}')

        # Handle inf/NaN and clip > fe_max
        fe = df_fes[fe_col].to_numpy(dtype=float)
        fe[np.isinf(fe) | np.isnan(fe)] = HIGH_VALUE
        fe[fe > args.fe_max] = HIGH_VALUE

        # Grids (first CV in cv_cols varies fastest per PLUMED convention)
        grids = [np.sort(df_fes[cv].unique()) for cv in args.cv_cols]
        n_bins = [len(g) for g in grids]
        if np.prod(n_bins) != len(df_fes):
            raise ValueError('Incomplete or irregular grid in FES')

        # Prepare for RegularGridInterpolator (slowest dimension first)
        grids_interp = tuple(reversed(grids))
        n_bins_interp = tuple(reversed(n_bins))
        fe_grid = fe.reshape(n_bins_interp)

        # Minima detection
        if dim == 1:
            min_idx = argrelextrema(fe_grid, np.less, order=1)[0]
            min_fe_raw = fe_grid[min_idx]
            min_pos_raw = grids[0][min_idx]
            min_pos_raw = min_pos_raw[:, np.newaxis]  # (n,1) for consistency
        else:
            min_filter = scipy.ndimage.minimum_filter(fe_grid, size=3,
                                                      mode='constant', cval=HIGH_VALUE)
            minima_mask = (min_filter == fe_grid)
            slow_i, fast_i = np.where(minima_mask)
            min_fe_raw = fe_grid[slow_i, fast_i]
            min_pos_raw = np.column_stack((grids[0][fast_i], grids[1][slow_i]))

        # Filter valid minima
        valid = min_fe_raw < HIGH_VALUE
        min_pos = min_pos_raw[valid]
        min_fe = min_fe_raw[valid]
        sort_idx = np.argsort(min_fe)
        min_pos = min_pos[sort_idx]
        min_fe = min_fe[sort_idx]

        logger.info(f'Found {len(min_fe)} local minima within FE < {args.fe_max} kJ/mol')
        for i, (pos, energy) in enumerate(zip(min_pos, min_fe), 1):
            pos_str = ' '.join(f'{p:.6f}' for p in pos)
            logger.info(f'Minimum {i}: CV(s) {pos_str}  FE {energy:.2f} kJ/mol')

        if len(min_fe) == 0:
            logger.warning('No valid minima found – exiting')
            return

        # Load COLVAR
        cv_fields = None
        with open(args.colvar) as f:
            for line in f:
                if line.startswith('#! FIELDS'):
                    cv_fields = line.split()[2:]
                    break
        df_cv = pd.read_csv(args.colvar, delim_whitespace=True, comment='#', names=cv_fields)
        if 'time' not in df_cv.columns:
            logger.warning('No "time" column detected in COLVAR')
        cv_frames = df_cv[args.cv_cols].to_numpy()  # (n_frames, dim)

        # Validate/filter columns if filtering enabled
        filter_frames = None
        if filter_cvs is not None:
            missing = [c for c in filter_cvs if c not in df_cv.columns]
            if missing:
                raise ValueError(f'Missing filter CV columns in COLVAR: {missing}')
            filter_frames = df_cv[filter_cvs].to_numpy()

        # Load trajectory once
        u = mda.Universe(args.top, args.traj)
        if len(u.trajectory) != len(df_cv):
            raise ValueError(f'Frame count mismatch: traj {len(u.trajectory)} vs COLVAR {len(df_cv)}')

        # Interpolator
        interpolator = scipy.interpolate.RegularGridInterpolator(
            grids_interp, fe_grid,
            method='linear', bounds_error=False, fill_value=HIGH_VALUE)

        query_points = cv_frames[:, ::-1] if dim == 2 else cv_frames
        fe_frames = interpolator(query_points)

        # Voronoi assignment in CV space
        diffs = cv_frames[:, np.newaxis, :] - min_pos[np.newaxis, :, :]
        dists = np.linalg.norm(diffs, axis=2)
        closest_min_idx = np.argmin(dists, axis=1)

        # Extraction per minimum
        for k in range(len(min_fe)):
            basin_mask = (closest_min_idx == k) & (fe_frames <= min_fe[k] + args.delta_fe)
            selected_frames = np.flatnonzero(basin_mask)
            if len(selected_frames) == 0:
                logger.warning(f'No frames assigned to minimum {k+1}')
                continue
            selected_frames = np.sort(selected_frames)
            logger.info(f'Minimum {k+1}: {len(selected_frames)} frames assigned (main basin)')
            out_xtc = f'minima_{k+1}.xtc'
            with mda.Writer(out_xtc, u.atoms.n_atoms) as W:
                for fr in selected_frames:
                    u.trajectory[fr]
                    W.write(u.atoms)
            logger.info(f'Written {out_xtc}')

            # Additional per-CV filtered trajectories
            if filter_frames is not None:
                for i, cv_name in enumerate(filter_cvs):
                    filter_mask = (filter_frames[:, i] >= filter_mins[i]) & (filter_frames[:, i] <= filter_maxs[i])
                    filtered_mask = basin_mask & filter_mask
                    filt_selected = np.flatnonzero(filtered_mask)
                    if len(filt_selected) == 0:
                        logger.warning(f'No frames for minimum {k+1} with filter {cv_name}')
                        continue
                    filt_selected = np.sort(filt_selected)
                    cv_safe = cv_name.replace('.', '_')
                    out_filt = f'minima_{k+1}_filter_{cv_safe}.xtc'
                    with mda.Writer(out_filt, u.atoms.n_atoms) as W:
                        for fr in filt_selected:
                            u.trajectory[fr]
                            W.write(u.atoms)
                    logger.info(f'Written filtered {out_filt} ({len(filt_selected)} frames)')

        logger.info('Analysis complete')

    except Exception as e:
        logger.error(f'Fatal error: {e}')
        raise

if __name__ == '__main__':
    main()
