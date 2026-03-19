#!/bin/bash
#SBATCH -J cluster
#SBATCH -p cpu
#SBATCH --mem=32GB
#SBATCH -n 1
#SBSTCH -c 1

source ~/data/scripts/gmx2023.5_plumed2.9.3_gpu_shared.env

echo 1 1 | gmx cluster -method gromos -f ../apo_100ns_merged_fit.xtc -s ../apo_100ns_0.tpr -rmsmin 0.2 -cutoff 0.4 -sz -clid -cl clusters.xtc
echo 1 | gmx trjconv -s ../apo_100ns_0.tpr -f clusters.xtc -sep -pbc mol -o aagp2_E64V.gro
