#!/bin/bash
##SBATCH -J EuL3_opt
#SBATCH -J EuL3_opt
#SBATCH -n 1
#SBATCH -c 32
#SBATCH -p cpu
#SBATCH --mem=120GB
##SBATCH --mem-per-cpu=1250 #--mem=120GB

export g16root=/data/zhaojy/g16/
export GAUSS_SCRDIR=/data/nglokwan/.scr/
export GAUSS_EXEDIR=/data/zhaojy/g16/g16/bsd:/share/home/zhaojy/g16/g16
source /data/zhaojy/g16/g16/bsd/g16.profile

#g16 < optfreq_S0.gjf > optfreq_S0.log

g16 < ma_optfreq.gjf > ma_optfreq.log

