# Aim
Write a python script to analyze gromacs + gmx_MMPBSA results of a protein-ligand system. The following functions should be included:
0. load topology and  trajectory, select only atoms in trajectory from topology to load successfully, align frames on receptor backbone
1. rmsd: receptor backbone, ligand heavy atom, raw csv data and png plots (overlap both lines on same plot)
2. contact: distace-based per-residue contact number between receptor and ligand, time series, csv and png plot.
3. hbond: per-residue hbond count between receptor and ligand, time series, csv and png plot.:
4. plot mmpbsa results without the gui: generate plots as in gmx_mmpbsa_ana (e.g. energy, decomposition...) but without the need of gui

# Workflow
This is complicated multi-stage work. split into research stage, plan stage, execute stagw, check stage. Spawn subagents for that
1. research on how to do aims 0-4 + combine. write research context to file
2. plan: read research context, write multiple small plans, check plans
3. execute: execute plans to complete the script
4. check: test, debug the complete script.

# Explanation of test files
This directory has the proccessed trajectory and topology for mmpbsa, and results of mmpbsa calculation.

## system topology and trajectory
com.tpr - full topology of the solvated system
com_traj.xtc stripped trj with only receptor and ligand
dzp.itp
posre.itp
rec.itp
sys.top
index.ndx

## scripts for trj processing and running gmx_MMPBSA
trj4mmpbsa.sh
mmpbsa.in
mmpbsa.sh

## Results of mmgbsa
COMPACT_MMXSA_RESULTS.mmxsa
FINAL_DECOMP_MMPBSA.csv
FINAL_DECOMP_MMPBSA.dat
FINAL_RESULTS_MMPBSA.csv
FINAL_RESULTS_MMPBSA.dat
gmx_MMPBSA.log

# Behavior
* Do not use the rm command, do not edit files listed above
* be profesional, clear, concise. be critical to own work when checking. 
* generalize so the code work on any specified protein-ligand system
- take cmd flags, and define default variables at the top of the script
- do not hardcode when it can be variable
- no double code around variables that are paths
- the code should be efficient and consistent. using vectorized approach for speed, while maintain scientific accuracy. use only libraries in the current environment.
