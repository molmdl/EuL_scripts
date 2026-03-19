#!/bin/bash
#SBATCH -J EuL3_opt
##SBATCH -J tddft
#SBATCH -n 1
#SBATCH -c 48
#SBATCH -p cpu
##SBATCH -w node2
#SBATCH --mem=120GB
##SBATCH --mem-per-cpu=1250 #--mem=120GB

source /share/apps/intel/oneapi/setvars.sh
export g16root=/data/zhaojy/g16/
export GAUSS_SCRDIR=/data/nglokwan/.scr/
export GAUSS_EXEDIR=/data/zhaojy/g16/g16/bsd:/share/home/zhaojy/g16/g16
source /data/zhaojy/g16/g16/bsd/g16.profile

#g16 < optfreq_S0.gjf > optfreq_S0.log
#g16 < tddft.gjf > tddft.log
g16 < opt_s1.gjf > opt_s1.log
