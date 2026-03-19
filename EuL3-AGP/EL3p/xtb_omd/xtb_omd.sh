#!/bin/bash
#SBATCH -J xtb_opt
#SBATCH -n 1
#SBATCH -c 16
#SBATCH -p cpu

source /share/apps/intel/oneapi/setvars.sh
source ~/scripts/xtb.env
export MKL_NUM_THREADS=16
export OMP_NUM_THREADS=16

xtb E3P_w.xyz -omd --input xtb_md.inp --gfn2 --chrg 1 --gbsa h2o  > xtb_omd.log # --gbsa h2o  
#xtb EL1_w.xyz -md --input xtb_md.inp --gfn0  > xtb_omd.log # --gbsa h2o  

