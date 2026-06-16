import os
import glob
import math
from datetime import date

R = 8.314e-3  # kJ/(mol·K)
T = 298.15    # K
KCAL_TO_KJ = 4.184

def dG_from_logK(logK):
    return -2.303 * R * T * logK

def scale_factor(dG_exp_kJ, dG_comp_kJ):
    return dG_exp_kJ / dG_comp_kJ

def process_mmpbsa(mmpbsa_dir="mmpbsa",
                   ref_ligand="ibp_r",
                   logK_ref=5.7,
                   T=298.15):
    files = sorted(glob.glob(os.path.join(mmpbsa_dir, "*.dat")))
    raw = []

    for f in files:
        name = os.path.basename(f).replace("_avg_mmpbsa.dat", "")
        with open(f) as fh:
            for line in fh:
                if line.startswith("TOTAL"):
                    parts = line.split()
                    avg_kcal = float(parts[1])
                    sd_kcal = float(parts[2])
                    avg_kJ = avg_kcal * KCAL_TO_KJ
                    sd_kJ = sd_kcal * KCAL_TO_KJ
                    raw.append((name, avg_kcal, sd_kcal, avg_kJ, sd_kJ))
                    break

    dG_exp = dG_from_logK(logK_ref)

    ref_kJ = None
    for name, avg_kcal, sd_kcal, avg_kJ, sd_kJ in raw:
        if name == ref_ligand:
            ref_kJ = avg_kJ
            break

    if ref_kJ is None:
        raise ValueError(f"Reference ligand '{ref_ligand}' not found in data")

    f_scale = scale_factor(dG_exp, ref_kJ)

    today = date.today().isoformat()
    outname = f"mmpbsa_dG_scaled_{today}.csv"

    with open(outname, "w") as out:
        out.write("ligand,dG_kcal,SD_kcal,dG_kJ,SD_kJ,dG_scaled_kJ,SD_scaled_kJ\n")
        for name, avg_kcal, sd_kcal, avg_kJ, sd_kJ in raw:
            dG_scaled = f_scale * avg_kJ
            sd_scaled = abs(f_scale) * sd_kJ
            out.write(f"{name},{avg_kcal:.2f},{sd_kcal:.2f},{avg_kJ:.2f},{sd_kJ:.2f},{dG_scaled:.2f},{sd_scaled:.2f}\n")

    print(f"Scale factor f = {f_scale:.4f}")
    print(f"ΔG_exp = {dG_exp:.2f} kJ/mol")
    print(f"Reference {ref_ligand}: ΔG_comp = {ref_kJ:.2f} kJ/mol")
    print(f"Written: {outname}")

if __name__ == "__main__":
    process_mmpbsa()
