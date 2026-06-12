#!/usr/bin/env python3
"""
T4 — N-C-C-N Torsion Classification (16-System Summary)

Classifies SAP vs TSAP geometry in each frame of 16 trajectories using
N-C-C-N ring torsions and chromophore torsions.

Usage:
    python scripts/T4_classify_torsion.py
    python scripts/T4_classify_torsion.py --data-dir data --out-dir analysis
"""

import argparse
import csv
import json
import os
import sys
from pathlib import Path

import numpy as np
import MDAnalysis as mda

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt


# ═══════════════════════════════════════════════════════
#  DIHEDRAL COMPUTATION (from metal_geo_analysis.py)
# ═══════════════════════════════════════════════════════

def calc_dihedrals_batch(pos: np.ndarray, ijkl: np.ndarray) -> np.ndarray:
    """
    Compute T dihedral angles in one vectorised pass.

    pos  : (N_atoms, 3) positions
    ijkl : (N_torsions, 4) integer index array, each row = [i, j, k, l]

    Returns (N_torsions,) array of angles in degrees [-180, 180].
    NaN is returned for degenerate cases.
    """
    i, j, k, l = ijkl[:, 0], ijkl[:, 1], ijkl[:, 2], ijkl[:, 3]
    b1 = pos[j] - pos[i]          # (T, 3)
    b2 = pos[k] - pos[j]          # (T, 3)
    b3 = pos[l] - pos[k]          # (T, 3)
    n1 = np.cross(b1, b2)         # (T, 3)
    n2 = np.cross(b2, b3)         # (T, 3)
    nn1 = np.linalg.norm(n1, axis=1, keepdims=True)   # (T, 1)
    nn2 = np.linalg.norm(n2, axis=1, keepdims=True)
    b2n = b2 / np.linalg.norm(b2, axis=1, keepdims=True)
    # avoid division by zero
    valid = (nn1[:, 0] > 1e-9) & (nn2[:, 0] > 1e-9)
    nn1 = np.where(nn1 > 1e-9, nn1, 1.0)
    nn2 = np.where(nn2 > 1e-9, nn2, 1.0)
    n1 /= nn1
    n2 /= nn2
    # IUPAC/GROMACS convention: sin term = (n1 × n2) · b2_hat
    n1xn2 = np.cross(n1, n2)       # (T, 3)
    dot_sin = np.einsum('ti,ti->t', n1xn2, b2n)
    dot_n1_n2 = np.einsum('ti,ti->t', n1, n2)
    angles = np.degrees(np.arctan2(dot_sin, dot_n1_n2))
    angles[~valid] = np.nan
    return angles


# ═══════════════════════════════════════════════════════
#  SYSTEM DISCOVERY
# ═══════════════════════════════════════════════════════

def discover_systems(data_dir: str, indices: dict) -> list:
    """
    Discover all systems with traj.xyz files in data_dir.

    Returns list of dicts: [{"name": str, "traj": str, "indices": dict}, ...]
    """
    data_path = Path(data_dir)
    traj_files = sorted(data_path.glob("*/traj.xyz"))

    if len(traj_files) != 16:
        raise ValueError(
            f"Expected exactly 16 traj.xyz files, found {len(traj_files)}"
        )

    systems = []
    for traj_path in traj_files:
        name = traj_path.parent.name
        if "_sap" in name and "_tsap" not in name:
            idx_table = indices["sap"]
        elif "_tsap" in name:
            idx_table = indices["tsap"]
        else:
            raise ValueError(
                f"System '{name}' does not contain '_sap' or '_tsap' in its name"
            )
        systems.append({
            "name": name,
            "traj": str(traj_path),
            "indices": idx_table,
        })

    return systems


# ═══════════════════════════════════════════════════════
#  TORSION CLASSIFICATION PER SYSTEM
# ═══════════════════════════════════════════════════════

def classify_system(system: dict) -> list:
    """
    Classify every frame in a single trajectory.

    Returns list of dicts, one per frame:
        {"system": name, "frame": int,
         "t1": float|nan, "t2": float|nan, "t3": float|nan, "t4": float|nan,
         "tc": float|nan, "classification": str}
    """
    name = system["name"]
    traj_path = system["traj"]
    idx_table = system["indices"]

    # Build (5, 4) index array: rows 0-3 = ring torsions, row 4 = chrom torsion
    ring_torsions = idx_table["ring_torsions"]
    chrom_torsion = idx_table["chrom_torsion"]

    if len(ring_torsions) != 4:
        raise ValueError(
            f"System {name}: expected 4 ring_torsions, got {len(ring_torsions)}"
        )

    _TORSION_IJKL = np.zeros((5, 4), dtype=int)
    for ridx, rt in enumerate(ring_torsions):
        atoms = rt["atoms"]
        if len(atoms) != 4:
            raise ValueError(
                f"System {name}: ring torsion {ridx} has {len(atoms)} atoms, expected 4"
            )
        _TORSION_IJKL[ridx, :] = atoms

    chrom_atoms = chrom_torsion["atoms"]
    if len(chrom_atoms) != 4:
        raise ValueError(
            f"System {name}: chrom_torsion has {len(chrom_atoms)} atoms, expected 4"
        )
    _TORSION_IJKL[4, :] = chrom_atoms

    # Load trajectory
    u = mda.Universe(traj_path)
    records = []

    for ts in u.trajectory:
        pos = u.atoms.positions.copy()   # (N_atoms, 3)
        all_tors = calc_dihedrals_batch(pos, _TORSION_IJKL)  # (5,)

        t1, t2, t3, t4, tc = all_tors
        ring_arr = all_tors[0:4]
        valid = ring_arr[~np.isnan(ring_arr)]

        if len(valid) == 0 or np.isnan(tc):
            geom = 'UNK'
        else:
            mean_ring = float(valid.mean())
            tc_f = float(tc)
            if abs(mean_ring) < 10.0 or abs(tc_f) < 10.0:
                geom = 'UNK'
            else:
                geom = 'SAP' if np.sign(mean_ring) != np.sign(tc_f) else 'TSAP'

        records.append({
            "system": name,
            "frame": ts.frame,
            "t1": t1,
            "t2": t2,
            "t3": t3,
            "t4": t4,
            "tc": tc,
            "classification": geom,
        })

    return records


# ═══════════════════════════════════════════════════════
#  OUTPUT WRITERS
# ═══════════════════════════════════════════════════════

def write_csv(records: list, out_path: str):
    """Write frame records to CSV."""
    with open(out_path, 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(["system", "frame", "t1", "t2", "t3", "t4", "tc", "classification"])
        for rec in records:
            row = [
                rec["system"],
                rec["frame"],
                f"{rec['t1']:.2f}" if not np.isnan(rec["t1"]) else "",
                f"{rec['t2']:.2f}" if not np.isnan(rec["t2"]) else "",
                f"{rec['t3']:.2f}" if not np.isnan(rec["t3"]) else "",
                f"{rec['t4']:.2f}" if not np.isnan(rec["t4"]) else "",
                f"{rec['tc']:.2f}" if not np.isnan(rec["tc"]) else "",
                rec["classification"],
            ]
            writer.writerow(row)


def compute_statistics(records: list) -> list:
    """
    Compute per-system SAP/TSAP/UNK fractions.

    Returns list of dicts sorted alphabetically by system name.
    """
    from collections import defaultdict

    stats = defaultdict(lambda: {"total": 0, "SAP": 0, "TSAP": 0, "UNK": 0})
    for rec in records:
        sys_name = rec["system"]
        stats[sys_name]["total"] += 1
        stats[sys_name][rec["classification"]] += 1

    result = []
    for sys_name in sorted(stats.keys()):
        s = stats[sys_name]
        total = s["total"]
        result.append({
            "system": sys_name,
            "total": total,
            "n_sap": s["SAP"],
            "n_tsap": s["TSAP"],
            "n_unk": s["UNK"],
            "frac_sap": s["SAP"] / total,
            "frac_tsap": s["TSAP"] / total,
            "frac_unk": s["UNK"] / total,
        })
    return result


def plot_summary(stats: list, out_path: str):
    """Generate combined bar plot of SAP/TSAP/UNK fractions."""
    n_systems = len(stats)
    system_names = [s["system"] for s in stats]
    frac_sap = [s["frac_sap"] for s in stats]
    frac_tsap = [s["frac_tsap"] for s in stats]
    frac_unk = [s["frac_unk"] for s in stats]

    x = np.arange(n_systems)
    width = 0.25

    fig, ax = plt.subplots(figsize=(14, 5))
    bars_sap = ax.bar(x - width, frac_sap, width, label='SAP', color='#4682B4')
    bars_tsap = ax.bar(x, frac_tsap, width, label='TSAP', color='#BA55D3')
    bars_unk = ax.bar(x + width, frac_unk, width, label='UNK', color='#D3D3D3')

    ax.set_ylabel('Fraction of frames')
    ax.set_title('NCCN Torsion-Based Geometry Classification (16 Systems)')
    ax.set_xticks(x)
    ax.set_xticklabels(system_names, rotation=45, ha='right')
    ax.set_ylim(0, 1.05)
    ax.legend()

    # Add numeric labels on top of bars
    for bars in [bars_sap, bars_tsap, bars_unk]:
        for bar in bars:
            height = bar.get_height()
            ax.annotate(f'{height:.2f}',
                        xy=(bar.get_x() + bar.get_width() / 2, height),
                        xytext=(0, 3),
                        textcoords="offset points",
                        ha='center', va='bottom',
                        fontsize=6)

    fig.tight_layout()
    fig.savefig(out_path, dpi=150)
    plt.close(fig)


def print_summary(stats: list, csv_path: str, png_path: str):
    """Print terminal summary table."""
    print()
    print(f"{'System':<20} {'Frames':>8} {'SAP':>8} {'TSAP':>8} {'UNK':>8}")
    print("-" * 56)
    for s in stats:
        print(
            f"{s['system']:<20} "
            f"{s['total']:>8} "
            f"{s['frac_sap']:>8.3f} "
            f"{s['frac_tsap']:>8.3f} "
            f"{s['frac_unk']:>8.3f}"
        )
    print()
    print(f"CSV output: {os.path.abspath(csv_path)}")
    print(f"PNG output: {os.path.abspath(png_path)}")

    # Print mismatch warning
    print()
    print("Mismatch check (systems where classification doesn't match starting label):")
    mismatch_found = False
    for s in stats:
        name = s["system"]
        if "_sap" in name and "_tsap" not in name:
            expected = "SAP"
        elif "_tsap" in name:
            expected = "TSAP"
        else:
            continue
        dominant = max([("SAP", s["frac_sap"]), ("TSAP", s["frac_tsap"]), ("UNK", s["frac_unk"])], key=lambda x: x[1])[0]
        if dominant != expected and s["frac_unk"] < 0.5:
            mismatch_found = True
            print(f"  WARNING: {name} (label={expected}) has dominant class {dominant} "
                  f"(SAP={s['frac_sap']:.2f}, TSAP={s['frac_tsap']:.2f}, UNK={s['frac_unk']:.2f})")
    if not mismatch_found:
        print("  None detected (all SAP-labeled systems dominated by SAP, TSAP by TSAP).")
    print()


# ═══════════════════════════════════════════════════════
#  VALIDATION
# ═══════════════════════════════════════════════════════

def validate_outputs(records: list, stats: list, csv_path: str, png_path: str):
    """Run runtime validation checks."""
    # V5: For every system, n_sap + n_tsap + n_unk == n_total
    from collections import defaultdict
    counts = defaultdict(lambda: {"total": 0, "SAP": 0, "TSAP": 0, "UNK": 0})
    for rec in records:
        counts[rec["system"]]["total"] += 1
        counts[rec["system"]][rec["classification"]] += 1

    for sys_name, c in counts.items():
        total = c["total"]
        summed = c["SAP"] + c["TSAP"] + c["UNK"]
        assert summed == total, f"V5 failed for {sys_name}: {summed} != {total}"

    # V6: CSV file written and has header + data rows
    assert os.path.exists(csv_path), f"V6 failed: CSV not found at {csv_path}"
    with open(csv_path) as f:
        lines = f.readlines()
    assert len(lines) >= 2, f"V6 failed: CSV has only {len(lines)} lines"
    header = lines[0].strip()
    assert header == "system,frame,t1,t2,t3,t4,tc,classification", f"V6 failed: bad header: {header}"

    # V7: PNG file exists and size > 0
    assert os.path.exists(png_path), f"V7 failed: PNG not found at {png_path}"
    assert os.path.getsize(png_path) > 0, f"V7 failed: PNG is empty"

    # V8: Fractions per system sum to 1.0 (within 1e-6)
    for s in stats:
        total_frac = s["frac_sap"] + s["frac_tsap"] + s["frac_unk"]
        assert abs(total_frac - 1.0) < 1e-6, (
            f"V8 failed for {s['system']}: fractions sum to {total_frac}"
        )

    print("All validation checks passed (V5-V8).")


# ═══════════════════════════════════════════════════════
#  MAIN
# ═══════════════════════════════════════════════════════

def main():
    parser = argparse.ArgumentParser(
        description="T4 - NCCN Torsion Classification for 16 xTB MD trajectories"
    )
    parser.add_argument(
        "--data-dir", default="data",
        help="Directory containing system subdirectories with traj.xyz files"
    )
    parser.add_argument(
        "--out-dir", default="analysis",
        help="Directory to write output CSV and PNG"
    )
    parser.add_argument(
        "--indices", default="analysis/indices.json",
        help="Path to indices.json with SAP/TSAP atom index tables"
    )
    args = parser.parse_args()

    # Load indices
    if not os.path.exists(args.indices):
        raise FileNotFoundError(f"Indices file not found: {args.indices}")
    with open(args.indices) as f:
        indices = json.load(f)

    # V2: indices.json contains both sap and tsap keys
    assert "sap" in indices, "V2 failed: 'sap' key missing from indices.json"
    assert "tsap" in indices, "V2 failed: 'tsap' key missing from indices.json"

    # V3: Each index table has exactly 4 ring_torsions and 1 chrom_torsion
    for key in ("sap", "tsap"):
        table = indices[key]
        assert len(table.get("ring_torsions", [])) == 4, (
            f"V3 failed: {key} ring_torsions count != 4"
        )
        assert "chrom_torsion" in table, (
            f"V3 failed: {key} missing chrom_torsion"
        )

    # Discover systems
    systems = discover_systems(args.data_dir, indices)
    print(f"Discovered {len(systems)} systems.")

    # Process all systems
    all_records = []
    for sys in systems:
        print(f"Processing {sys['name']} ...", end=" ", flush=True)
        records = classify_system(sys)
        all_records.extend(records)
        print(f"{len(records)} frames")

    # Ensure output directory exists
    os.makedirs(args.out_dir, exist_ok=True)

    # Write CSV
    csv_path = os.path.join(args.out_dir, "torsion_classification.csv")
    write_csv(all_records, csv_path)

    # Compute statistics and plot
    stats = compute_statistics(all_records)
    png_path = os.path.join(args.out_dir, "plot_torsion_classification_summary.png")
    plot_summary(stats, png_path)

    # Print summary
    print_summary(stats, csv_path, png_path)

    # Validate
    validate_outputs(all_records, stats, csv_path, png_path)

    return 0


if __name__ == '__main__':
    sys.exit(main())
