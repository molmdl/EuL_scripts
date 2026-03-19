import sys
import os
import argparse
import logging
import numpy as np

# Defaults
DEFAULT_TEMP = 300.0
KB = 0.008314462618  # kJ/mol/K
DEFAULT_KBT = KB * DEFAULT_TEMP
DEFAULT_PREFIX = "merged"
DEFAULT_MIN_TO_ZERO = True

def setup_logging():
    logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

def parse_arguments():
    parser = argparse.ArgumentParser(description="Average PLUMED histograms and compute free energy surface.")
    parser.add_argument("list_file", type=str, help="Text file containing list of histogram file paths (one per line).")
    parser.add_argument("--temp", type=float, default=None, help="Temperature in K (default: 300.0).")
    parser.add_argument("--kbt", type=float, default=None, help="kBT value in kJ/mol (overrides temp).")
    parser.add_argument("--prefix", type=str, default=DEFAULT_PREFIX, help="Output file prefix (default: 'merged').")
    parser.add_argument("--mintozero", action="store_true", default=DEFAULT_MIN_TO_ZERO,
                        help="Shift FES so minimum is 0 kJ/mol (default: True)")
    parser.add_argument("--no-mintozero", dest="mintozero", action="store_false",
                        help="Do NOT shift FES minimum to zero")
    return parser.parse_args()

def parse_header(headers):
    info = {}
    for line in headers:
        if line.startswith("#! FIELDS"):
            info['fields'] = line.split()[2:]
        elif line.startswith("#! SET"):
            parts = line.split()
            key = parts[2]
            value = parts[3]
            if key.startswith('min_') or key.startswith('max_'):
                info[key] = float(value)
            elif key.startswith('nbins_'):
                info[key] = int(value)
            elif key == 'normalisation':
                info[key] = float(value)
            else:
                info[key] = value
    return info

def validate_file_consistency(ref_info, info, file_path):
    if len(ref_info['fields']) != len(info['fields']):
        raise ValueError(f"Dimension mismatch in {file_path}")
    for key in ref_info:
        if key.startswith('min_') or key.startswith('max_') or key.startswith('nbins_') or key.startswith('periodic_'):
            if ref_info[key] != info[key]:
                raise ValueError(f"Grid mismatch in {file_path} for {key}")

def main():
    setup_logging()
    args = parse_arguments()

    try:
        # Determine kbt
        if args.kbt is not None:
            kbt = args.kbt
        elif args.temp is not None:
            kbt = KB * args.temp
        else:
            kbt = DEFAULT_KBT
        logging.info(f"Using kBT = {kbt:.6f} kJ/mol")

        # Read list of files
        with open(args.list_file, 'r') as f:
            hist_files = [line.strip() for line in f if line.strip()]

        if not hist_files:
            raise ValueError("No histogram files listed.")

        # Check file existence
        for file_path in hist_files:
            if not os.path.exists(file_path):
                raise FileNotFoundError(f"File not found: {file_path}")

        logging.info(f"Processing {len(hist_files)} histogram files.")

        # Process first file to get reference
        first_file = hist_files[0]
        headers = []
        with open(first_file, 'r') as f:
            line = f.readline()
            while line.startswith('#!'):
                headers.append(line.strip())
                line = f.readline()
        skiprows = len(headers)

        ref_info = parse_header(headers)
        fields = ref_info['fields']
        num_fields = len(fields)

        if num_fields not in [3, 5]:
            raise ValueError("Only 1D (3 fields) or 2D (5 fields) histograms supported.")

        num_cvs = (num_fields - 1) // 2
        hist_col = num_cvs   # index of probability column (0-based)

        logging.info(f"Detected {num_cvs}D histogram with {num_fields} columns.")

        data = np.loadtxt(first_file, skiprows=skiprows)
        if data.shape[1] != num_fields:
            raise ValueError(f"First file has wrong number of columns: expected {num_fields}, got {data.shape[1]}")

        coords = data[:, :num_cvs]
        hists = [data[:, hist_col]]

        n_points = data.shape[0]
        logging.info(f"Reference grid has {n_points} bins.")

        # Process remaining files
        for file_path in hist_files[1:]:
            headers = []
            with open(file_path, 'r') as f:
                line = f.readline()
                while line.startswith('#!'):
                    headers.append(line.strip())
                    line = f.readline()

            info = parse_header(headers)
            validate_file_consistency(ref_info, info, file_path)

            data = np.loadtxt(file_path, skiprows=len(headers))

            if data.shape[0] != n_points or data.shape[1] != num_fields:
                raise ValueError(
                    f"Data shape mismatch in {file_path}: "
                    f"expected ({n_points}, {num_fields}), got {data.shape}"
                )

            hists.append(data[:, hist_col])

        # Average
        stack = np.stack(hists, axis=0)
        avg_hist = np.mean(stack, axis=0)

        # Gentle normalization
        total = np.sum(avg_hist)
        if total == 0:
            raise ValueError("Average histogram sum is zero.")
        avg_hist /= total
        logging.info("Averaged and normalized histogram.")

        # Compute FES
        fes = -kbt * np.log(avg_hist)

        # Cap infinities to maximum finite value
        finite_mask = np.isfinite(fes)
        if np.any(finite_mask):
            max_fes = np.max(fes[finite_mask])
            fes[~finite_mask] = max_fes
        else:
            fes[:] = 0.0  # extremely unlikely
        logging.info("Capped infinite FES values to maximum finite value.")

        # Optional: shift minimum to zero
        if args.mintozero:
            finite_mask = np.isfinite(fes)  # re-check after capping
            if np.any(finite_mask):
                min_fes = np.min(fes[finite_mask])
                fes -= min_fes
                logging.info(f"Shifted FES by -{min_fes:.3f} kJ/mol so minimum is 0.")
            else:
                logging.warning("No finite FES values → cannot shift to zero.")
        else:
            logging.info("FES minimum NOT shifted to zero (as requested).")

        # Output formats
        if num_cvs == 1:
            fmt = ['%.6f', '%.10e']
        else:
            fmt = ['%.6f', '%.6f', '%.10e']

        # Save averaged hist (no derivatives)
        hist_out = f"{args.prefix}_hist.dat"
        np.savetxt(hist_out, np.column_stack((coords, avg_hist)), fmt=fmt, delimiter=' ')
        logging.info(f"Saved averaged histogram to {hist_out}")

        # Save FES (no derivatives)
        fes_out = f"{args.prefix}_fes.dat"
        np.savetxt(fes_out, np.column_stack((coords, fes)), fmt=fmt, delimiter=' ')
        logging.info(f"Saved FES to {fes_out}")

    except Exception as e:
        logging.error(f"Error: {str(e)}")
        sys.exit(1)

if __name__ == "__main__":
    main()
