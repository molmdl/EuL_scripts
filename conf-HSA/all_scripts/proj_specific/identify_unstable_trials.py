#!/usr/bin/env python3
"""
identify_unstable_trials.py
───────────────────────────
Scan all per-trial contacts_total.csv (and optionally rmsd_ligand.csv) under
all_ana/per_ligand/ and flag trials whose time-series shows sustained
loss-of-binding events.

INSTABILITY CRITERIA (primary = contacts, secondary = RMSD)
────────────────────────────────────────────────────────────
1. CONTACT DROP  (primary, decisive)
   A frame is "contact-unstable" when:
       total_contacts < trial_mean − N_SIGMA × trial_std
   where N_SIGMA = 2.0 by default.

   Rationale: contacts vary with a coefficient of variation of ~7 % for
   well-bound ligands (std/mean ≈ 40/600).  Mean − 2σ sits well inside the
   Gaussian tail and marks a genuine drop, not thermal noise.  For trials
   where the ligand has partly dissociated, the distribution is bimodal
   (bound + unbound) with a high std; mean − 2σ can become negative, which
   means *no* individual frame ever crosses the threshold — but those trials
   are caught by the ABSOLUTE_CONTACT_FLOOR safety net (see below).

2. ABSOLUTE CONTACT FLOOR  (catches dissociated trials)
   A trial is always flagged if its mean total contacts < FLOOR_CONTACTS.
   FLOOR_CONTACTS = 200 by default.

   Rationale: well-bound trials across all 14 ligands range from 250 to 830
   mean contacts.  Trials with mean < 200 (e.g. 75, 89, 142) correspond to
   near-complete unbinding for the majority of the trajectory.

3. RMSD EXCEEDANCE  (secondary, corroborating)
   A frame is "rmsd-unstable" when:
       ligand_RMSD > RMSD_CUTOFF_A  (default 5 Å)

   NOTE: these ligands are large and the binding pocket is flexible; RMSD
   values of 5–20 Å are common even for bound states.  RMSD is therefore
   used as a *corroborating* signal only: a trial is reported as
   RMSD-unstable if it passes the RMSD window test, but RMSD alone does NOT
   flag a trial — only the contact criteria do.

4. SUSTAINED WINDOW  (filters transient fluctuations)
   A trial is contact-flagged if it contains ≥ WINDOW_FRAMES consecutive
   contact-unstable frames.
   WINDOW_FRAMES = 50 by default (= 5 ns at 100 ps/frame).

   Rationale: the user requested high tolerance because the pocket is
   flexible.  Short dips are normal; 5 ns of consecutive low contact is a
   genuine unbinding event.  The absolute-floor criterion is window-
   independent (it flags based on mean, not streaks).

All thresholds are tunable via command-line arguments.
"""

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd


# ─── helpers ──────────────────────────────────────────────────────────────────

def longest_consecutive_run(mask: np.ndarray) -> int:
    """Return the length of the longest True-run in a boolean array."""
    if not mask.any():
        return 0
    # pad with False so that changes at borders are detected
    padded = np.concatenate(([False], mask, [False]))
    diff = np.diff(padded.astype(int))
    starts = np.where(diff == 1)[0]
    ends   = np.where(diff == -1)[0]
    return int((ends - starts).max())


def fraction_unstable(mask: np.ndarray) -> float:
    return float(mask.mean())


def load_contacts(trial_dir: Path) -> pd.Series | None:
    p = trial_dir / "contacts" / "contacts_total.csv"
    if not p.exists():
        return None
    df = pd.read_csv(p)
    return df["total_contacts"]


def load_rmsd(trial_dir: Path) -> pd.Series | None:
    p = trial_dir / "rmsd" / "rmsd_ligand.csv"
    if not p.exists():
        return None
    df = pd.read_csv(p)
    return df["rmsd"]


# ─── main logic ───────────────────────────────────────────────────────────────

def analyse_trial(
    trial_dir: Path,
    n_sigma: float,
    floor_contacts: float,
    rmsd_cutoff: float,
    window_frames: int,
) -> dict:
    """Return a dict of metrics for one trial."""
    result = {
        "path": str(trial_dir),
        "contacts_mean": np.nan,
        "contacts_std": np.nan,
        "contacts_threshold": np.nan,
        "longest_contact_window": 0,
        "frac_contact_unstable": np.nan,
        "below_floor": False,
        "rmsd_longest_window": 0,
        "frac_rmsd_unstable": np.nan,
        "flagged_contacts": False,
        "flagged_rmsd": False,
        "flagged": False,
        "flag_reason": "",
    }

    # ── contacts ──────────────────────────────────────────────────────────────
    contacts = load_contacts(trial_dir)
    if contacts is not None:
        mean_c = contacts.mean()
        std_c  = contacts.std()
        result["contacts_mean"]      = round(mean_c, 1)
        result["contacts_std"]       = round(std_c, 1)
        threshold = mean_c - n_sigma * std_c
        result["contacts_threshold"] = round(threshold, 1)

        below_floor = mean_c < floor_contacts
        result["below_floor"] = below_floor

        mask_low = contacts.values < threshold
        longest  = longest_consecutive_run(mask_low)
        frac     = fraction_unstable(mask_low)
        result["longest_contact_window"] = longest
        result["frac_contact_unstable"]  = round(frac, 3)

        reasons = []
        if below_floor:
            reasons.append(
                f"mean contacts {mean_c:.0f} < floor {floor_contacts:.0f}"
            )
        if longest >= window_frames:
            reasons.append(
                f"contact drop sustained {longest} frames "
                f"(>= {window_frames}) @ threshold {threshold:.0f}"
            )
        if reasons:
            result["flagged_contacts"] = True
            result["flag_reason"] = "; ".join(reasons)

    # ── RMSD (corroborating only) ──────────────────────────────────────────
    rmsd = load_rmsd(trial_dir)
    if rmsd is not None:
        mask_rmsd = rmsd.values > rmsd_cutoff
        longest_r = longest_consecutive_run(mask_rmsd)
        frac_r    = fraction_unstable(mask_rmsd)
        result["rmsd_longest_window"] = longest_r
        result["frac_rmsd_unstable"]  = round(frac_r, 3)
        if longest_r >= window_frames:
            result["flagged_rmsd"] = True

    # ── overall flag ──────────────────────────────────────────────────────────
    result["flagged"] = result["flagged_contacts"]
    return result


def find_trials(base_dir: Path):
    """Yield (ligand_name, trial_name, trial_dir) for all trials."""
    per_ligand = base_dir / "per_ligand"
    for lig_dir in sorted(per_ligand.iterdir()):
        if not lig_dir.is_dir():
            continue
        trial_root = lig_dir / "per_trial"
        if not trial_root.exists():
            continue
        for trial_dir in sorted(trial_root.iterdir()):
            if trial_dir.is_dir():
                yield lig_dir.name, trial_dir.name, trial_dir


# ─── CLI ──────────────────────────────────────────────────────────────────────

def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    p.add_argument(
        "--base-dir",
        default=Path(__file__).parent / "all_ana",
        type=Path,
        help="Path to all_ana/ directory (default: ./all_ana)",
    )
    p.add_argument(
        "--n-sigma",
        type=float,
        default=2.0,
        metavar="N",
        help="Contact drop threshold: flag frames below mean − N×std  [default: 2.0]",
    )
    p.add_argument(
        "--floor-contacts",
        type=float,
        default=200.0,
        metavar="F",
        help="Always flag trials whose mean contacts < F (catches near-dissociation) "
             "[default: 200]",
    )
    p.add_argument(
        "--rmsd-cutoff",
        type=float,
        default=5.0,
        metavar="Å",
        help="RMSD threshold in Å for corroborating instability  [default: 5.0]",
    )
    p.add_argument(
        "--window",
        type=int,
        default=50,
        metavar="W",
        help="Minimum consecutive unstable frames to flag a trial  [default: 50 = 5 ns]",
    )
    p.add_argument(
        "--csv",
        type=Path,
        default=None,
        metavar="FILE",
        help="Optional: write full results table to a CSV file",
    )
    p.add_argument(
        "--verbose", "-v",
        action="store_true",
        help="Print all trials (not just flagged ones)",
    )
    return p


def main():
    args = build_parser().parse_args()
    base_dir: Path = args.base_dir

    if not base_dir.exists():
        sys.exit(f"ERROR: base directory not found: {base_dir}")

    print(
        f"\nInstability detection settings\n"
        f"  Contact threshold : mean − {args.n_sigma}×std per trial\n"
        f"  Absolute floor    : mean contacts < {args.floor_contacts:.0f}\n"
        f"  RMSD cutoff       : > {args.rmsd_cutoff} Å  (corroborating only)\n"
        f"  Sustained window  : ≥ {args.window} consecutive frames ({args.window * 0.1:.0f} ns)\n"
    )

    records = []
    flagged_paths = []

    for lig, trial, trial_dir in find_trials(base_dir):
        res = analyse_trial(
            trial_dir,
            n_sigma=args.n_sigma,
            floor_contacts=args.floor_contacts,
            rmsd_cutoff=args.rmsd_cutoff,
            window_frames=args.window,
        )
        res["ligand"] = lig
        res["trial"]  = trial
        records.append(res)

        if res["flagged"]:
            flagged_paths.append(res["path"])
            corr = ""
            if res["flagged_rmsd"]:
                corr = f"  [+RMSD: {res['rmsd_longest_window']} frames > {args.rmsd_cutoff} Å]"
            print(
                f"UNSTABLE  {lig}/{trial}\n"
                f"          {res['flag_reason']}{corr}\n"
                f"          path: {res['path']}"
            )
        elif args.verbose:
            print(
                f"stable    {lig}/{trial}  "
                f"mean={res['contacts_mean']:.0f}  "
                f"longest_drop={res['longest_contact_window']}"
            )

    n_total   = len(records)
    n_flagged = len(flagged_paths)

    print(f"\nSummary: {n_flagged} / {n_total} trials flagged as unstable\n")

    if n_flagged > 0:
        print("Unstable trial paths:")
        for p in flagged_paths:
            print(f"  {p}")

    if args.csv:
        cols = [
            "ligand", "trial", "path",
            "contacts_mean", "contacts_std", "contacts_threshold",
            "below_floor", "longest_contact_window", "frac_contact_unstable",
            "rmsd_longest_window", "frac_rmsd_unstable",
            "flagged_contacts", "flagged_rmsd", "flagged", "flag_reason",
        ]
        df = pd.DataFrame(records)[cols]
        df.to_csv(args.csv, index=False)
        print(f"\nFull results written to: {args.csv}")

    return 0


if __name__ == "__main__":
    sys.exit(main())
