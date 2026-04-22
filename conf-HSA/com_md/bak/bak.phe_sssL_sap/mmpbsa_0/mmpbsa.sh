#!/bin/bash
#SBATCH -J gmxMMPBSA
#SBATCH -p cpu
#SBATCH -n 16
#SBATCH -c 1

mpirun -np 16 gmx_MMPBSA -O -i mmpbsa.in -cs com.tpr -ct com_traj.xtc -ci ../index.ndx -cg 1 12 -cp ../bypass_sys.top -o FINAL_RESULTS_MMPBSA.dat -eo FINAL_RESULTS_MMPBSA.csv -do FINAL_DECOMP_MMPBSA.dat -deo FINAL_DECOMP_MMPBSA.csv -nogui
