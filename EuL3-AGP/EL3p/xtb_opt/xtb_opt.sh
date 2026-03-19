#!/bin/bash
#SBATCH -J xtb_opt
#SBATCH -n 1
#SBATCH -c 16
#SBATCH -p cpu

source /share/apps/intel/oneapi/setvars.sh
source /data/nglokwan/scripts/xtb.env
export MKL_NUM_THREADS=16
export OMP_NUM_THREADS=16

xtb $1 -o verytight --gfn2 --gbsa h2o -c 1  > xtb_opt.log
