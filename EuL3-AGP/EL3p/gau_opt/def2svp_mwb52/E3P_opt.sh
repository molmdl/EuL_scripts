#!/bin/bash
#SBATCH -J EuL3_opt
#SBATCH -n 1
#SBATCH -c 48
#SBATCH -p cpu

source /share/apps/intel/oneapi/setvars.sh
export g16root=/data/zhaojy/g16/
export GAUSS_SCRDIR=/data/nglokwan/.scr/
export GAUSS_EXEDIR=/data/zhaojy/g16/g16/bsd:/share/home/zhaojy/g16/g16
source /data/zhaojy/g16/g16/bsd/g16.profile

#g16 < optfreq.gjf > optfreq.log
g16 < SP_stableopt.gjf > SP_stableopt.log
