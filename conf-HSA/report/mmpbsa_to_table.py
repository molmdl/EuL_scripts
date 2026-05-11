import os
import glob

files = sorted(glob.glob("mmpbsa/*.dat"))
data = []

for f in files:
    name = os.path.basename(f).replace("_avg_mmpbsa.dat", "")
    with open(f) as fh:
        lines = fh.readlines()[1:]  # skip header
    
    row = {"File": name}
    for line in lines:
        parts = line.split()
        if len(parts) >= 3:
            term = parts[0]
            avg = parts[1]
            sd = parts[2]
            row[term] = f"{avg} ± {sd}"
    data.append(row)

# Write CSV
with open("mmpbsa_summary.csv", "w") as f:
    header = ["File", "VDWAALS", "EEL", "GGAS", "GSOLV", "TOTAL"]
    f.write(",".join(header) + "\n")
    for row in data:
        f.write(",".join([row.get(h, "") for h in header]) + "\n")

# Write Markdown
with open("mmpbsa_summary.md", "w") as f:
    header = ["File", "VDWAALS", "EEL", "GGAS", "GSOLV", "TOTAL"]
    f.write("| " + " | ".join(header) + " |\n")
    f.write("|" + "|".join(["---"]*6) + "|\n")
    for row in data:
        f.write("| " + " | ".join([row.get(h, "") for h in header]) + " |\n")

print("Created mmpbsa_summary.csv and mmpbsa_summary.md")
