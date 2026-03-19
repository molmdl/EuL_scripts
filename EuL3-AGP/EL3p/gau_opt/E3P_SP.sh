#!/bin/bash
#SBATCH -J EuL3_SP
#SBATCH -n 1
#SBATCH -c 48
#SBATCH -p cpu
#SBATCH -w node2
#SBATCH --mem-per-cpu=2500 #--mem=120GB

source /share/apps/intel/oneapi/setvars.sh
export g16root=/data/zhaojy/g16/
export GAUSS_SCRDIR=/data/nglokwan/.scr/
export GAUSS_EXEDIR=/data/zhaojy/g16/g16/bsd:/share/home/zhaojy/g16/g16
source /data/zhaojy/g16/g16/bsd/g16.profile

#g16 < optfreq_S0.gjf > optfreq_S0.log
#g16 < tddft.gjf > tddft.log
g16 < SP_stableopt_DKH4_x2c-TZVPall.gjf > SP_stableopt_DKH4_x2c-TZVPall.log
