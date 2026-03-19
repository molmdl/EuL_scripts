#!/bin/bash
#SBATCH -J plumed_reweight
#SBATCH -p cpu
#SBATCH -n 1
#SBATCH -c 1

source /data/nglokwan/scripts/gmx2023.5_plumed2.9.3_gpu_shared.env

#plumed driver --plumed ../plumed_reweight_single.dat --kt 2.494339 --mf_xtc ./fm_1200ns_${1}.xtc --trajectory-stride 20
plumed driver --plumed ./plumed_reweight_merged.dat --kt 2.494339 --noatoms
