#!/usr/bin/env python3
"""
HA3_rmsd_rmsf_allheavy.py — RMSD/RMSF with all-heavy-atom alignment
and core-aligned comparison for all 16 xTB trajectories.

Replaces MDAnalysis AlignTraj with manual Kabsch alignment on 68 scaffold
heavy atoms (same as HA1). Also computes core-aligned (Eu+8 = 9 atoms)
RMSD for comparison.

Outputs:
  analysis/rmsd_allheavy.csv        (32000 rows: system, frame, time_ps,
                                      rmsd_allheavy_A, rmsd_core_A,
                                      species, stereo, conformer)
  analysis/rmsf_allheavy.csv        (per-atom RMSF after all-heavy alignment)
  analysis/plot_rmsd_timeseries_allheavy.png
  analysis/plot_rmsd_distribution_allheavy.png  (KDE allheavy vs core)
  analysis/plot_rmsf_scaffold_allheavy.png      (per-atom RMSF bar chart)
  analysis/plot_rmsf_element_avg_allheavy.png

Usage:
    python scripts/HA3_rmsd_rmsf_allheavy.py \
        --data-dir data \
        --output-dir analysis \
        --mappings ../analysis/atom_mappings.json \
        --scaffold ../analysis/scaffold_mapping.json
"""

import argparse
import json
import os
import shutil
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

import MDAnalysis as mda

FRAME_INTERVAL_PS = 0.05

CORE_EU8_IDX = sorted([0, 1, 2, 3, 4, 5, 6, 7, 54])

SYSTEMS = [
    "me_rrrD_sap", "me_rrrL_sap", "me_sssD_sap", "me_sssL_sap",
    "me_rrrD_tsap", "me_rrrL_tsap", "me_sssD_tsap", "me_sssL_tsap",
    "phe_rrrD_sap", "phe_rrrL_sap", "phe_sssD_sap", "phe_sssL_sap",
    "phe_rrrD_tsap", "phe_rrrL_tsap", "phe_sssD_tsap", "phe_sssL_tsap",
]


def kabsch(P, Q):
    P_centroid = P.mean(axis=0)
    Q_centroid = Q.mean(axis=0)
    Pc = P - P_centroid
    Qc = Q - Q_centroid
    H = Pc.T @ Qc
    U, S, Vt = np.linalg.svd(H)
    d = np.sign(np.linalg.det(U @ Vt))
    D = np.diag([1, 1, d])
    R = U @ D @ Vt
    return R, P_centroid, Q_centroid


def backup_if_exists(path):
    if os.path.exists(path):
        backup_path = path + ".bak"
        shutil.copy2(path, backup_path)
        print(f"  Backed up existing file to {backup_path}")


def load_xyz_trajectory(traj_path):
    u = mda.Universe(traj_path, format="XYZ")
    n_frames = len(u.trajectory)
    n_atoms = u.atoms.n_atoms
    coords = np.empty((n_frames, n_atoms, 3), dtype=np.float64)
    for i, ts in enumerate(u.trajectory):
        coords[i] = ts.positions
    elements = [a.name for a in u.atoms]
    return coords, elements


def load_xyz_single(xyz_path):
    u = mda.Universe(xyz_path, format="XYZ")
    ts = u.trajectory[0]
    positions = ts.positions.copy()
    elements = [a.name for a in u.atoms]
    return positions, elements


def apply_tsap_mapping(coords, elements, mapping):
    mapped_coords = np.empty_like(coords)
    mapped_coords[:, mapping, :] = coords
    mapped_elements_arr = [None] * len(elements)
    for j, sap_idx in enumerate(mapping):
        mapped_elements_arr[sap_idx] = elements[j]
    return mapped_coords, mapped_elements_arr


def parse_system_info(system):
    parts = system.split("_")
    species = parts[0]
    stereo = parts[1]
    conformer = parts[2].upper()
    return species, stereo, conformer


def compute_rmsd(positions, ref_positions):
    return np.sqrt(np.mean((positions - ref_positions) ** 2, axis=(1, 2)))


def compute_rmsd_per_frame(aligned_coords, ref_frame, atom_indices):
    n_frames = aligned_coords.shape[0]
    ref_pos = ref_frame[np.array(atom_indices)]
    rmsds = np.empty(n_frames)
    for f in range(n_frames):
        diff = aligned_coords[f, atom_indices] - ref_pos
        rmsds[f] = np.sqrt(np.mean(np.sum(diff ** 2, axis=1)))
    return rmsds


def parse_args():
    parser = argparse.ArgumentParser(
        description="RMSD/RMSF with all-heavy alignment + core-aligned comparison")
    parser.add_argument("--data-dir", type=str, default="data",
                        help="Directory containing system subdirectories")
    parser.add_argument("--output-dir", type=str, default="analysis",
                        help="Output directory for CSVs and plots")
    parser.add_argument("--mappings", type=str, default="../analysis/atom_mappings.json",
                        help="Path to atom_mappings.json")
    parser.add_argument("--scaffold", type=str, default="../analysis/scaffold_mapping.json",
                        help="Path to scaffold_mapping.json")
    return parser.parse_args()


def main():
    args = parse_args()
    data_dir = Path(args.data_dir)
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    print("=" * 65)
    print("HA3: RMSD/RMSF with All-Heavy Alignment")
    print("=" * 65)

    # ── 1. Load mappings ──────────────────────────────────────
    print(f"\n[1/5] Loading mapping files ...")

    with open(args.mappings) as f:
        atom_mappings = json.load(f)
    with open(args.scaffold) as f:
        scaffold_mapping = json.load(f)

    n_scaffold = scaffold_mapping["n_scaffold"]
    print(f"  Scaffold atoms: {n_scaffold}")

    # ── 2. Load reference ──────────────────────────────────────
    print(f"\n[2/5] Loading reference (me_rrrD_sap xtbopt.xyz) ...")

    ref_path = data_dir / "me_rrrD_sap" / "xtbopt.xyz"
    ref_frame, ref_elements = load_xyz_single(str(ref_path))
    n_ref = len(ref_elements)
    print(f"  Reference: {n_ref} atoms")

    SCAFFOLD_HEAVY_IDX = [i for i in range(n_scaffold) if ref_elements[i] != "H"]
    print(f"  Scaffold heavy atoms: {len(SCAFFOLD_HEAVY_IDX)}")
    if len(SCAFFOLD_HEAVY_IDX) != 68:
        print(f"  [WARN] Expected 68 scaffold heavy atoms, got {len(SCAFFOLD_HEAVY_IDX)}")

    ref_allheavy_pos = ref_frame[np.array(SCAFFOLD_HEAVY_IDX)]
    ref_core_pos = ref_frame[np.array(CORE_EU8_IDX)]

    # ── 3. Process all systems ────────────────────────────────
    print(f"\n[3/5] Processing all 16 systems ...")

    all_rmsd_dfs = []
    all_rmsf_rows = []

    for system in SYSTEMS:
        species, stereo, conformer = parse_system_info(system)
        start_conf = "sap" if system.endswith("_sap") else "tsap"

        traj_path = data_dir / system / "traj.xyz"
        print(f"  Loading {system} ...", end=" ")
        coords, elements = load_xyz_trajectory(str(traj_path))
        n_frames, n_atoms, _ = coords.shape
        print(f"{n_frames} frames, {n_atoms} atoms")

        if start_conf == "tsap":
            mapping_key = f"{species}_tsap_to_sap"
            mapping = atom_mappings[mapping_key]
            coords, elements = apply_tsap_mapping(coords, elements, mapping)

        # ── All-heavy alignment (68 scaffold heavy atoms) ──
        aligned_allheavy = np.empty_like(coords)
        for f in range(n_frames):
            mobile_align = coords[f, np.array(SCAFFOLD_HEAVY_IDX)]
            R, P_cent, Q_cent = kabsch(mobile_align, ref_allheavy_pos)
            aligned_allheavy[f] = (coords[f] - P_cent) @ R + Q_cent

        # RMSD of 68 scaffold heavy atoms after all-heavy alignment
        rmsd_allheavy_A = compute_rmsd_per_frame(
            aligned_allheavy, ref_frame, SCAFFOLD_HEAVY_IDX)

        # ── Core alignment (Eu+8 = 9 atoms) ──
        aligned_core = np.empty_like(coords)
        for f in range(n_frames):
            mobile_align = coords[f, np.array(CORE_EU8_IDX)]
            R, P_cent, Q_cent = kabsch(mobile_align, ref_core_pos)
            aligned_core[f] = (coords[f] - P_cent) @ R + Q_cent

        # RMSD of 68 scaffold heavy atoms after core alignment
        rmsd_core_A = compute_rmsd_per_frame(
            aligned_core, ref_frame, SCAFFOLD_HEAVY_IDX)

        # ── Build RMSD DataFrame ──
        rmsd_df = pd.DataFrame({
            "system": [system] * n_frames,
            "frame": np.arange(n_frames),
            "time_ps": np.arange(n_frames) * FRAME_INTERVAL_PS,
            "rmsd_allheavy_A": rmsd_allheavy_A,
            "rmsd_core_A": rmsd_core_A,
            "species": [species] * n_frames,
            "stereo": [stereo] * n_frames,
            "conformer": [conformer] * n_frames,
        })
        all_rmsd_dfs.append(rmsd_df)

        # ── RMSF computation (after all-heavy alignment) ──
        mean_coords = aligned_allheavy[:, SCAFFOLD_HEAVY_IDX, :].mean(axis=0)
        deviations = aligned_allheavy[:, SCAFFOLD_HEAVY_IDX, :] - mean_coords
        per_atom_rmsf = np.sqrt(np.mean(np.sum(deviations ** 2, axis=2), axis=0))

        for i, aidx in enumerate(SCAFFOLD_HEAVY_IDX):
            all_rmsf_rows.append({
                "system": system,
                "atom_index": aidx,
                "element": elements[aidx],
                "species": species,
                "stereo": stereo,
                "conformer": conformer,
                "rmsf_A": per_atom_rmsf[i],
            })

        # Print per-system summary
        print(f"    allheavy RMSD: mean={rmsd_allheavy_A.mean():.4f}, "
              f"max={rmsd_allheavy_A.max():.4f} A")
        print(f"    core RMSD:     mean={rmsd_core_A.mean():.4f}, "
              f"max={rmsd_core_A.max():.4f} A")

        del coords, aligned_allheavy, aligned_core

    all_rmsd = pd.concat(all_rmsd_dfs, ignore_index=True)
    rmsf_df = pd.DataFrame(all_rmsf_rows)

    # ── Save CSVs ──────────────────────────────────────────────
    print(f"\n[4/5] Saving CSVs ...")

    rmsd_path = output_dir / "rmsd_allheavy.csv"
    backup_if_exists(str(rmsd_path))
    all_rmsd.to_csv(rmsd_path, index=False)
    print(f"  Saved {rmsd_path} ({len(all_rmsd)} rows)")

    rmsf_path = output_dir / "rmsf_allheavy.csv"
    backup_if_exists(str(rmsf_path))
    rmsf_df.to_csv(rmsf_path, index=False)
    print(f"  Saved {rmsf_path} ({len(rmsf_df)} rows)")

    nan_count = all_rmsd.isnull().sum().sum()
    print(f"  NaN values in RMSD: {nan_count}")

    # ── 5. Generate plots ─────────────────────────────────────
    print(f"\n[5/5] Generating plots ...")

    # --- Plot 1: RMSD timeseries ---
    fig, axes = plt.subplots(1, 2, figsize=(14, 5), sharey=True)
    for idx, conf in enumerate(["SAP", "TSAP"]):
        ax = axes[idx]
        sub = all_rmsd[all_rmsd["conformer"] == conf]
        for species, color in [("me", "#2a9d8f"), ("phe", "#e76f51")]:
            for stereo in ["rrrD", "rrrL", "sssD", "sssL"]:
                mask = (sub["species"] == species) & (sub["stereo"] == stereo)
                data = sub[mask]
                if len(data) > 0:
                    label = f"{species}_{stereo}"
                    alpha = 0.3 if stereo.endswith("D") else 0.3
                    ls = "-" if stereo.endswith("D") else "--"
                    ax.plot(data["time_ps"], data["rmsd_allheavy_A"],
                            color=color, alpha=alpha, linewidth=0.5, linestyle=ls,
                            label=label)
        ax.set_title(f"{conf}")
        ax.set_xlabel("Time (ps)")
        if idx == 0:
            ax.set_ylabel("RMSD (Å) — all-heavy alignment")
        handles, labels = ax.get_legend_handles_labels()
        by_label = dict(zip(labels, handles))
        ax.legend(by_label.values(), by_label.keys(), fontsize=7, framealpha=0.9)
    fig.suptitle("RMSD Timeseries — 68 scaffold heavy atoms after all-heavy alignment", fontsize=11)
    fig.tight_layout()
    p1_path = output_dir / "plot_rmsd_timeseries_allheavy.png"
    backup_if_exists(str(p1_path))
    fig.savefig(p1_path, dpi=300)
    plt.close(fig)
    print(f"  Saved {p1_path}")

    # --- Plot 2: RMSD distribution (KDE allheavy vs core) ---
    fig, axes = plt.subplots(1, 2, figsize=(12, 5), sharey=True)
    for idx, conf in enumerate(["SAP", "TSAP"]):
        ax = axes[idx]
        sub = all_rmsd[all_rmsd["conformer"] == conf]
        for species, color in [("me", "#2a9d8f"), ("phe", "#e76f51")]:
            data_ah = sub[sub["species"] == species]["rmsd_allheavy_A"]
            data_cr = sub[sub["species"] == species]["rmsd_core_A"]
            if len(data_ah) > 0:
                from scipy.stats import gaussian_kde
                try:
                    kde_ah = gaussian_kde(data_ah)
                    kde_cr = gaussian_kde(data_cr)
                    x_range = np.linspace(0, max(data_cr.max(), data_ah.max()) * 1.1, 200)
                    ax.plot(x_range, kde_ah(x_range), color=color, linestyle="-",
                            label=f"{species} allheavy", alpha=0.8)
                    ax.plot(x_range, kde_cr(x_range), color=color, linestyle="--",
                            label=f"{species} core", alpha=0.8)
                except Exception:
                    ax.hist(data_ah, bins=30, alpha=0.3, color=color,
                            label=f"{species} allheavy")
                    ax.hist(data_cr, bins=30, alpha=0.3, color=color,
                            label=f"{species} core", histtype="step")
        ax.set_title(f"{conf}")
        ax.set_xlabel("RMSD (Å)")
        if idx == 0:
            ax.set_ylabel("Density")
        handles, labels = ax.get_legend_handles_labels()
        by_label = dict(zip(labels, handles))
        ax.legend(by_label.values(), by_label.keys(), fontsize=7, framealpha=0.9)
    fig.suptitle("RMSD Distribution — all-heavy vs core (Eu+8) alignment\n"
                 "(solid = all-heavy, dashed = core)", fontsize=11)
    fig.tight_layout()
    p2_path = output_dir / "plot_rmsd_distribution_allheavy.png"
    backup_if_exists(str(p2_path))
    fig.savefig(p2_path, dpi=300)
    plt.close(fig)
    print(f"  Saved {p2_path}")

    # --- Plot 3: RMSF per scaffold atom ---
    element_colors = {
        "O": "#e63946", "N": "#457b9d", "C": "#a8a8a8", "Eu": "#d4a017"
    }

    # Average RMSF per atom across all systems
    avg_rmsf = rmsf_df.groupby("atom_index")["rmsf_A"].mean().reset_index()
    std_rmsf = rmsf_df.groupby("atom_index")["rmsf_A"].std(ddof=1).reset_index()
    avg_rmsf = avg_rmsf.merge(std_rmsf, on="atom_index", suffixes=("_mean", "_std"))
    avg_rmsf = avg_rmsf.merge(
        rmsf_df[["atom_index", "element"]].drop_duplicates(), on="atom_index"
    )
    avg_rmsf = avg_rmsf.sort_values("atom_index")

    fig, ax = plt.subplots(figsize=(14, 5))
    x = np.arange(len(avg_rmsf))
    colors = [element_colors.get(e, "#999999") for e in avg_rmsf["element"]]
    ax.bar(x, avg_rmsf["rmsf_A_mean"], yerr=avg_rmsf["rmsf_A_std"],
           color=colors, edgecolor="none", alpha=0.8, capsize=1)
    ax.set_xlabel("Scaffold heavy atom (sorted by atom_index)")
    ax.set_ylabel("Mean RMSF (Å)")
    ax.set_title("Per-atom RMSF after all-heavy alignment (mean ± std across 16 systems)")
    handles = [plt.Rectangle((0, 0), 1, 1, color=element_colors[e])
               for e in ["Eu", "O", "N", "C"]]
    ax.legend(handles, ["Eu", "O", "N", "C"], title="Element", loc="upper left")
    fig.tight_layout()
    p3_path = output_dir / "plot_rmsf_scaffold_allheavy.png"
    backup_if_exists(str(p3_path))
    fig.savefig(p3_path, dpi=300)
    plt.close(fig)
    print(f"  Saved {p3_path}")

    # --- Plot 4: Element-averaged RMSF ---
    elem_avg = rmsf_df.groupby(["species", "conformer", "element"])["rmsf_A"].mean().reset_index()
    elem_avg = elem_avg[elem_avg["element"].isin(["Eu", "O", "N", "C"])]

    fig, axes = plt.subplots(1, 2, figsize=(12, 5), sharey=True)
    for idx, sp in enumerate(["me", "phe"]):
        ax = axes[idx]
        sub = elem_avg[elem_avg["species"] == sp]
        conf_colors = {"SAP": "#2a9d8f", "TSAP": "#e76f51"}
        elements = ["Eu", "O", "N", "C"]
        x = np.arange(len(elements))
        width = 0.35
        for ci, conf in enumerate(["SAP", "TSAP"]):
            cdata = sub[sub["conformer"] == conf]
            vals = [cdata[cdata["element"] == e]["rmsf_A"].values[0]
                    if len(cdata[cdata["element"] == e]) > 0 else 0
                    for e in elements]
            ax.bar(x + ci * width, vals, width,
                   label=conf, color=conf_colors[conf], alpha=0.8, edgecolor="white")
        ax.set_xticks(x + width / 2)
        ax.set_xticklabels(elements)
        ax.set_xlabel("Element")
        if idx == 0:
            ax.set_ylabel("Mean RMSF (Å)")
        ax.set_title(f"Species: {sp}")
        ax.legend(title="Conformer", fontsize=8)
    fig.suptitle("Element-averaged RMSF after all-heavy alignment", fontsize=11)
    fig.tight_layout()
    p4_path = output_dir / "plot_rmsf_element_avg_allheavy.png"
    backup_if_exists(str(p4_path))
    fig.savefig(p4_path, dpi=300)
    plt.close(fig)
    print(f"  Saved {p4_path}")

    # ── Validation ────────────────────────────────────────────
    print("\n" + "=" * 65)
    print("VALIDATION")
    print("=" * 65)
    checks_passed = 0
    checks_total = 6

    correct_rows = len(all_rmsd) == 32000
    print(f"[{'PASS' if correct_rows else 'FAIL'}] rmsd_allheavy.csv has {len(all_rmsd)} rows (expected 32000)")
    checks_passed += correct_rows

    rmsf_exists = rmsf_path.exists()
    print(f"[{'PASS' if rmsf_exists else 'FAIL'}] rmsf_allheavy.csv exists")
    checks_passed += rmsf_exists

    plot_files = [
        output_dir / "plot_rmsd_timeseries_allheavy.png",
        output_dir / "plot_rmsd_distribution_allheavy.png",
        output_dir / "plot_rmsf_scaffold_allheavy.png",
        output_dir / "plot_rmsf_element_avg_allheavy.png",
    ]
    all_plots = all(p.exists() for p in plot_files)
    print(f"[{'PASS' if all_plots else 'FAIL'}] All 4 PNG plots created")
    checks_passed += all_plots

    # allheavy RMSD should be near-zero for alignment atoms
    frame0_rmsd = all_rmsd[all_rmsd["frame"] == 0]["rmsd_allheavy_A"]
    near_zero = (frame0_rmsd < 0.01).all()
    print(f"[{'PASS' if near_zero else 'WARN'}] Frame 0 allheavy RMSD near-zero: "
          f"max={frame0_rmsd.max():.6f} A")

    # core RMSD should be > 0 (core alignment doesn't optimize scaffold-wide fit)
    nonzero_core = (all_rmsd["rmsd_core_A"] > 0).all()
    print(f"[{'PASS' if nonzero_core else 'FAIL'}] Core-aligned RMSD > 0 for all frames")
    checks_passed += nonzero_core

    no_nan = nan_count == 0
    print(f"[{'PASS' if no_nan else 'FAIL'}] No NaN values in RMSD CSV")
    checks_passed += no_nan

    print(f"\n{'=' * 65}")
    print(f"Checks passed: {checks_passed}/{checks_total}")
    print(f"{'=' * 65}")

    # ── Summary statistics ─────────────────────────────────────
    print("\n" + "=" * 65)
    print("RMSD SUMMARY")
    print("=" * 65)
    for system in SYSTEMS:
        sub = all_rmsd[all_rmsd["system"] == system]
        ah_mean = sub["rmsd_allheavy_A"].mean()
        ah_max = sub["rmsd_allheavy_A"].max()
        cr_mean = sub["rmsd_core_A"].mean()
        cr_max = sub["rmsd_core_A"].max()
        ratio = cr_mean / ah_mean if ah_mean > 0 else float("inf")
        print(f"  {system:25s}  allheavy={ah_mean:.4f}/{ah_max:.4f}  "
              f"core={cr_mean:.4f}/{cr_max:.4f}  ratio={ratio:.2f}x")

    print("\n" + "=" * 65)
    print("RMSF SUMMARY (top 10 most flexible scaffold atoms)")
    print("=" * 65)
    top_atoms = avg_rmsf.nlargest(10, "rmsf_A_mean")
    for _, row in top_atoms.iterrows():
        print(f"  Atom {row['atom_index']:3d} ({row['element']:2s}): "
              f"RMSF = {row['rmsf_A_mean']:.4f} ± {row['rmsf_A_std']:.4f} Å")

    print("\n" + "=" * 65)
    print("OUTPUTS")
    print("=" * 65)
    print(f"  {rmsd_path}")
    print(f"  {rmsf_path}")
    for p in plot_files:
        print(f"  {p}")
    print("=" * 65)

    return 0 if checks_passed >= checks_total - 1 else 1


if __name__ == "__main__":
    sys.exit(main())
