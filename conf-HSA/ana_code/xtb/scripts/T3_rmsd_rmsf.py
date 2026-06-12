#!/usr/bin/env python3
"""
RMSD and RMSF Analysis Script

Computes RMSD vs xtbopt reference (heavy atoms and Eu+O+N core) and per-atom RMSF after core alignment.
Produces CSV tables and publication-quality plots.

Usage:
    python scripts/T3_rmsd_rmsf.py --data-dir data --output-dir analysis

Requires: MDAnalysis 2.7.0, numpy, pandas, matplotlib
"""

import argparse
import pathlib
import sys
from collections import defaultdict

import numpy as np
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
import seaborn as sns

import MDAnalysis as mda
from MDAnalysis.analysis import align, rms

matplotlib.use("Agg")

# ────────────────────────── Constants ──────────────────────────
# Core selection for RMSD alignment: all O/N donor atoms + Eu.
# This is a SUPERSET of the Eu+8 PCA alignment atoms [0,1,2,3,4,5,6,7,54]
# because it includes ALL O/N atoms in the molecule (not just the 8 scaffold
# donors). This is standard practice for core RMSD — a broader core definition
# that captures the full coordination sphere. The Eu+8 benchmark is the
# alignment used for joint PCA (see T3_joint_pca_eu8.py), not for RMSD.
CORE_SEL = "name Eu O N"
HEAVY_SEL = "not name H"
FRAME_INTERVAL_PS = 0.05  # 50 fs between saved frames (subsampled from 1 fs MD step, 100 ps / 2000 frames)

SYSTEMS = [
    "me_rrrD_sap", "me_rrrL_sap", "me_sssD_sap", "me_sssL_sap",
    "me_rrrD_tsap", "me_rrrL_tsap", "me_sssD_tsap", "me_sssL_tsap",
    "phe_rrrD_sap", "phe_rrrL_sap", "phe_sssD_sap", "phe_sssL_sap",
    "phe_rrrD_tsap", "phe_rrrL_tsap", "phe_sssD_tsap", "phe_sssL_tsap",
]

CONFORMER_TYPES = {"sap": "SAP", "tsap": "TSAP"}


# ────────────────────────── Helpers ──────────────────────────
def parse_system_info(system: str):
    """Parse system name into (species, stereochemistry, conformer)."""
    parts = system.split("_")
    species = parts[0]
    stereo = parts[1]
    conformer = CONFORMER_TYPES[parts[2]]
    return species, stereo, conformer


def compute_rmsd(positions, ref_positions):
    """Compute RMSD (Å) per frame: positions [n_frames, n_atoms, 3], ref_positions [n_atoms, 3]."""
    return np.sqrt(np.mean((positions - ref_positions) ** 2, axis=(1, 2)))


def process_system(data_dir: pathlib.Path, system: str):
    """
    Process a single system.
    Returns dict with rmsd_df, rmsf_df, rmsf_values, n_frames, atom_names, rmsd_stats.
    """
    traj_path = data_dir / system / "traj.xyz"
    ref_path = data_dir / system / "xtbopt.xyz"

    u = mda.Universe(str(traj_path))
    xref = mda.Universe(str(ref_path))

    n_frames = u.trajectory.n_frames
    n_atoms = len(u.atoms)

    print(f"  Processing {system}: {n_frames} frames, {n_atoms} atoms")

    if n_frames != 2000:
        print(f"    WARNING: {system} has {n_frames} frames (expected 2000)", file=sys.stderr)

    # ── Align on core ──
    align.AlignTraj(u, xref, select=CORE_SEL, in_memory=True).run()

    # ── Core RMSD (sanity check) ──
    core_ag = u.select_atoms(CORE_SEL)
    ref_core = xref.select_atoms(CORE_SEL)
    core_pos = np.array([core_ag.positions.copy() for _ in u.trajectory])
    ref_core_pos = ref_core.positions
    core_rmsd = compute_rmsd(core_pos, ref_core_pos)

    # sanity check – reset trajectory after array comprehension
    u.trajectory[0]

    # ── Heavy atom RMSD ──
    heavy_ag = u.select_atoms(HEAVY_SEL)
    ref_heavy = xref.select_atoms(HEAVY_SEL)
    heavy_pos = np.array([heavy_ag.positions.copy() for _ in u.trajectory])
    ref_heavy_pos = ref_heavy.positions
    heavy_rmsd = compute_rmsd(heavy_pos, ref_heavy_pos)

    # reset trajectory again
    u.trajectory[0]

    # ── Per-atom RMSF (ALL atoms) ──
    rmsf_all = rms.RMSF(u.atoms).run()
    rmsf_values = rmsf_all.results.rmsf

    # ── Build RMSD DataFrame ──
    rmsd_df = pd.DataFrame({
        "system": [system] * n_frames,
        "frame": np.arange(n_frames),
        "time_ps": np.arange(n_frames) * FRAME_INTERVAL_PS,
        "rmsd_heavy_A": heavy_rmsd,
        "rmsd_core_A": core_rmsd,
    })

    # Parse metadata
    species, stereo, conformer = parse_system_info(system)
    rmsd_df["species"] = species
    rmsd_df["stereo"] = stereo
    rmsd_df["conformer"] = conformer

    # ── Build RMSF DataFrame ──
    atom_names = u.atoms.names
    rmsf_df = pd.DataFrame({
        "system": [system] * n_atoms,
        "atom_index": np.arange(n_atoms),
        "element": atom_names,
        "rmsf_A": rmsf_values,
    })

    # per-system stats for stdout
    rmsd_stats = {
        "n_frames": n_frames,
        "rmsd_heavy_mean": float(heavy_rmsd.mean()),
        "rmsd_heavy_std": float(heavy_rmsd.std()),
        "rmsd_heavy_min": float(heavy_rmsd.min()),
        "rmsd_heavy_max": float(heavy_rmsd.max()),
        "rmsd_core_first": float(core_rmsd[0]),
        "rmsd_core_max": float(core_rmsd.max()),
    }

    return {
        "rmsd_df": rmsd_df,
        "rmsf_df": rmsf_df,
        "rmsf_values": rmsf_values,
        "n_frames": n_frames,
        "n_atoms": n_atoms,
        "atom_names": atom_names,
        "rmsd_stats": rmsd_stats,
        "species": species,
        "conformer": conformer,
    }


def build_scaffold_comparison(results, conformer_type):
    """
    Build scaffold comparison DataFrame for a given conformer type (SAP or TSAP).
    Common scaffold = first min_atoms indices within conformer type.
    Returns DataFrame with columns:
      [conformer_type, atom_index, element, <system>_rmsf for each system]
    """
    systems = [s for s in SYSTEMS if s.endswith(conformer_type.lower())]
    if not systems:
        return pd.DataFrame()

    min_atoms = min(results[s]["n_atoms"] for s in systems)

    df_data = {
        "conformer_type": [CONFORMER_TYPES[conformer_type.lower()]] * min_atoms,
        "atom_index": np.arange(min_atoms),
        "element": results[systems[0]]["atom_names"][:min_atoms],
    }

    for s in systems:
        df_data[f"{s}_rmsf"] = results[s]["rmsf_values"][:min_atoms]

    return pd.DataFrame(df_data)


def build_element_averaged(results):
    """Build element-averaged RMSF DataFrame for cross-conformer comparison."""
    rows = []
    for s in SYSTEMS:
        r = results[s]
        elements = np.unique(r["atom_names"])
        for elem in elements:
            mask = r["atom_names"] == elem
            avg_rmsf = r["rmsf_values"][mask].mean()
            rows.append({
                "system": s,
                "species": r["species"],
                "conformer": r["conformer"],
                "element": elem,
                "mean_rmsf_A": avg_rmsf,
            })
    return pd.DataFrame(rows)


# ────────────────────────── Plotting ──────────────────────────
def plot_rmsd_timeseries(all_rmsd: pd.DataFrame, output_dir: pathlib.Path):
    """Plot 1: RMSD time series, faceted by conformer type."""
    sns.set_theme(style="whitegrid", context="paper")

    g = sns.relplot(
        data=all_rmsd,
        x="time_ps",
        y="rmsd_heavy_A",
        col="conformer",
        hue="species",
        style="stereo",
        kind="line",
        ci=None,
        palette={"me": "#2a9d8f", "phe": "#e76f51"},
        height=4,
        aspect=1.2,
        facet_kws={"sharey": True, "sharex": True},
    )
    g.set_axis_labels("Time (ps)", "RMSDheavy (Å)")
    g.set_titles("{col_name}")
    g._legend.set_title("Species / Stereo")
    plt.tight_layout()
    g.savefig(output_dir / "plot_rmsd_timeseries.png", dpi=300)
    plt.close()
    print("  → Saved plot_rmsd_timeseries.png")


def plot_rmsd_distribution(all_rmsd: pd.DataFrame, output_dir: pathlib.Path):
    """Plot 2: RMSD distribution (KDE) grouped by species and conformer."""
    sns.set_theme(style="whitegrid", context="paper")

    fig, axes = plt.subplots(1, 2, figsize=(10, 4.5), sharey=True)
    for idx, conf in enumerate(["SAP", "TSAP"]):
        ax = axes[idx]
        sub = all_rmsd[all_rmsd["conformer"] == conf]
        for species, color in [("me", "#2a9d8f"), ("phe", "#e76f51")]:
            data = sub[sub["species"] == species]["rmsd_heavy_A"]
            if len(data) > 0:
                sns.kdeplot(data=data, ax=ax, color=color, label=species, fill=True, alpha=0.3)
        ax.set_title(conf)
        ax.set_xlabel("RMSDheavy (Å)")
        if idx == 0:
            ax.set_ylabel("Density")
        ax.legend(title="Species")
    plt.tight_layout()
    fig.savefig(output_dir / "plot_rmsd_distribution.png", dpi=300)
    plt.close()
    print("  → Saved plot_rmsd_distribution.png")


def plot_rmsf_scaffold(scaffold_sap, scaffold_tsap, output_dir: pathlib.Path):
    """Plot 3: RMSF bar chart for scaffold atoms."""
    sns.set_theme(style="whitegrid", context="paper")

    element_colors = {"O": "#e63946", "N": "#457b9d", "C": "#a8a8a8", "Eu": "#d4a017", "H": "#eeeeee"}

    fig, axes = plt.subplots(2, 1, figsize=(12, 10))
    for ax, conf, df in [(axes[0], "SAP", scaffold_sap), (axes[1], "TSAP", scaffold_tsap)]:
        if df.empty:
            ax.set_title(f"{conf} — no data")
            continue

        n_atoms = len(df)
        x = np.arange(n_atoms)
        systems = [c for c in df.columns if c.endswith("_rmsf")]
        mean_rmsf = df[systems].mean(axis=1).values
        std_rmsf = df[systems].std(axis=1).values

        colors = [element_colors.get(e, "#999999") for e in df["element"]]

        ax.bar(x, mean_rmsf, yerr=std_rmsf, color=colors, edgecolor="none", alpha=0.8)
        ax.set_xlabel("Scaffold atom index")
        ax.set_ylabel("Mean RMSF (Å)")
        ax.set_title(f"{conf} — scaffold RMSF (mean ± std across {len(systems)} trajectories)")

        # legend
        handles = [plt.Rectangle((0, 0), 1, 1, color=element_colors[e]) for e in ["O", "N", "C", "Eu", "H"]]
        ax.legend(handles, ["O", "N", "C", "Eu", "H"], title="Element", loc="upper left")

    plt.tight_layout()
    fig.savefig(output_dir / "plot_rmsf_scaffold.png", dpi=300)
    plt.close()
    print("  → Saved plot_rmsf_scaffold.png")


def plot_rmsf_element_avg(element_df: pd.DataFrame, output_dir: pathlib.Path):
    """Plot 4: Element-averaged RMSF comparison across conformers."""
    sns.set_theme(style="whitegrid", context="paper")

    # Filter to heavy elements only
    heavy_elements = ["Eu", "O", "N", "C"]
    df = element_df[element_df["element"].isin(heavy_elements)].copy()

    g = sns.catplot(
        data=df,
        x="element",
        y="mean_rmsf_A",
        hue="conformer",
        col="species",
        kind="bar",
        palette={"SAP": "#2a9d8f", "TSAP": "#e76f51"},
        height=4,
        aspect=1.0,
    )
    g.set_axis_labels("Element", "Mean RMSF (Å)")
    g.set_titles("Species: {col_name}")
    g._legend.set_title("Conformer")
    plt.tight_layout()
    g.savefig(output_dir / "plot_rmsf_element_avg.png", dpi=300)
    plt.close()
    print("  → Saved plot_rmsf_element_avg.png")


# ────────────────────────── Main ──────────────────────────
def main():
    parser = argparse.ArgumentParser(description="Compute RMSD/RMSF for xTB MD trajectories")
    parser.add_argument("--data-dir", type=pathlib.Path, default=pathlib.Path("data"),
                        help="Directory containing subdirectories for each system")
    parser.add_argument("--output-dir", type=pathlib.Path, default=pathlib.Path("analysis"),
                        help="Output directory for CSVs and plots")
    args = parser.parse_args()

    data_dir = args.data_dir
    output_dir = args.output_dir
    output_dir.mkdir(parents=True, exist_ok=True)

    print(f"Data directory:  {data_dir.resolve()}")
    print(f"Output directory: {output_dir.resolve()}")
    print(f"Systems to process: {len(SYSTEMS)}")
    print("-" * 50)

    # Process each system
    results = {}
    all_rmsd_dfs = []
    all_rmsf_dfs = []

    for system in SYSTEMS:
        result = process_system(data_dir, system)
        results[system] = result
        all_rmsd_dfs.append(result["rmsd_df"])
        all_rmsf_dfs.append(result["rmsf_df"])

    # Combine
    all_rmsd = pd.concat(all_rmsd_dfs, ignore_index=True)
    all_rmsf = pd.concat(all_rmsf_dfs, ignore_index=True)

    # ── Save combined RMSD ──
    rmsd_path = output_dir / "rmsd_all.csv"
    all_rmsd.to_csv(rmsd_path, index=False)
    print(f"\nSaved: {rmsd_path} ({len(all_rmsd)} rows)")

    # ── Save per-system RMSF CSVs ──
    rmsf_per_system_dir = output_dir / "rmsf_per_system"
    rmsf_per_system_dir.mkdir(exist_ok=True)
    for _, row in all_rmsf.groupby("system"):
        system_name = row["system"].iloc[0]
        fname = f"rmsf_{system_name}.csv"
        row.to_csv(rmsf_per_system_dir / fname, index=False)
    print(f"Saved: {rmsf_per_system_dir} (16 files)")

    # Also save all-in-one RMSF
    rmsf_all_path = output_dir / "rmsf_all.csv"
    all_rmsf.to_csv(rmsf_all_path, index=False)
    print(f"Saved: {rmsf_all_path} ({len(all_rmsf)} rows)")

    # ── Scaffold comparisons ──
    for conf in ["sap", "tsap"]:
        df = build_scaffold_comparison(results, conf)
        out_path = output_dir / f"rmsf_scaffold_comparison_{conf}.csv"
        df.to_csv(out_path, index=False)
        print(f"Saved: {out_path}")

    # ── Element-averaged ──
    element_avg = build_element_averaged(results)
    elem_path = output_dir / "rmsf_element_avg.csv"
    element_avg.to_csv(elem_path, index=False)
    print(f"Saved: {elem_path}")

    # ── Generate plots ──
    print("\nGenerating plots...")
    plot_rmsd_timeseries(all_rmsd, output_dir)
    plot_rmsd_distribution(all_rmsd, output_dir)
    plot_rmsf_scaffold(
        build_scaffold_comparison(results, "sap"),
        build_scaffold_comparison(results, "tsap"),
        output_dir,
    )
    plot_rmsf_element_avg(element_avg, output_dir)

    # ── Summary Statistics ──
    print("\n" + "=" * 60)
    print("SUMMARY STATISTICS")
    print("=" * 60)
    for system in SYSTEMS:
        st = results[system]["rmsd_stats"]
        status = "OK"
        if st["rmsd_core_first"] > 0.05:
            status = "ALIGNMENT SUSPECT"
        if st["rmsd_heavy_max"] > 5.0:
            status += " | POSSIBLE DISSOCIATION"
        print(
            f"{system:20s}  heavy={st['rmsd_heavy_mean']:.2f}±{st['rmsd_heavy_std']:.2f} Å  "
            f"[min/max={st['rmsd_heavy_min']:.2f}/{st['rmsd_heavy_max']:.2f}]  "
            f"core(first)={st['rmsd_core_first']:.4f} Å  [{status}]"
        )

    # ── Global checks ──
    print("\n" + "=" * 60)
    print("VALIDATION")
    print("=" * 60)
    checks_passed = 0
    checks_total = 6

    # 1. Frame count
    all_2000 = all(results[s]["n_frames"] == 2000 for s in SYSTEMS)
    print(f"[{'PASS' if all_2000 else 'FAIL'}] All systems have 2000 frames")
    checks_passed += all_2000

    # 2. RMSD sanity
    all_sane = all(results[s]["rmsd_stats"]["rmsd_core_first"] < 0.05 for s in SYSTEMS)
    print(f"[{'PASS' if all_sane else 'FAIL'}] First-frame core RMSD < 0.05 Å for all systems")
    checks_passed += all_sane

    # 3. RMSD range (exclude frame 0 where aligned reference RMSD is ~0)
    heavy_nonzero = all_rmsd[all_rmsd["frame"] > 0]["rmsd_heavy_A"]
    in_range = 0.05 <= heavy_nonzero.min() and heavy_nonzero.max() <= 5.0
    print(f"[{'PASS' if in_range else 'WARN'}] Heavy RMSD range (excl. frame 0): {heavy_nonzero.min():.3f}–{heavy_nonzero.max():.3f} Å (expected 0.05–5.0)")
    checks_passed += in_range

    # 4. CSV rows
    correct_rows = len(all_rmsd) == 32000
    print(f"[{'PASS' if correct_rows else 'FAIL'}] rmsd_all.csv has {len(all_rmsd)} rows (expected 32000)")
    checks_passed += correct_rows

    # 5. CSV columns
    expected_cols = {"system", "frame", "time_ps", "rmsd_heavy_A", "rmsd_core_A"}
    correct_cols = expected_cols <= set(all_rmsd.columns)
    print(f"[{'PASS' if correct_cols else 'FAIL'}] rmsd_all.csv has required columns")
    checks_passed += correct_cols

    # 6. Plots exist
    plot_files = [
        output_dir / "plot_rmsd_timeseries.png",
        output_dir / "plot_rmsd_distribution.png",
        output_dir / "plot_rmsf_scaffold.png",
        output_dir / "plot_rmsf_element_avg.png",
    ]
    all_plots = all(p.exists() for p in plot_files)
    print(f"[{'PASS' if all_plots else 'FAIL'}] All 4 PNG plots created")
    checks_passed += all_plots

    print(f"\n{'=' * 60}")
    print(f"Checks passed: {checks_passed}/{checks_total}")
    print(f"{'=' * 60}")

    if checks_passed < checks_total:
        print("\nNOTE: Some checks did not pass. Review warnings above.", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
