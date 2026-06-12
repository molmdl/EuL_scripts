import os
import glob

files = sorted(glob.glob("mmpbsa/*.dat"))

with open("mmpbsa_dG.csv", "w") as out:
    out.write("ligand,total_dG\n")
    for f in files:
        name = os.path.basename(f).replace("_avg_mmpbsa.dat", "")
        with open(f) as fh:
            for line in fh:
                if line.startswith("TOTAL"):
                    parts = line.split()
                    out.write(f"{name},{parts[1]}\n")
                    break
