#!/bin/bash
##SBATCH -J 120GB_optfreq
#SBATCH -J 50GB_sp
#SBATCH -n 1
#SBATCH -c 48
#SBATCH -p cpu
##SBATCH -w node1
##SBATCH --mem=120GB
#SBATCH --mem=50GB
##SBATCH --mem-per-cpu=1250 #--mem=120GB

source /share/apps/intel/oneapi/setvars.sh
export g16root=/data/zhaojy/g16/
export GAUSS_SCRDIR=/data/nglokwan/.scr/
export GAUSS_EXEDIR=/data/zhaojy/g16/g16/bsd:/share/home/zhaojy/g16/g16
source /data/zhaojy/g16/g16/bsd/g16.profile

#g16 < optfreq.gjf > optfreq.log
g16 < sp_solv.gjf > sp_solv.log
g16 < sp_gas.gjf > sp_gas.log

