#!/bin/bash
#SBATCH -J EuL3_opt
##SBATCH -J tddft
#SBATCH -n 1
#SBATCH -c 36
#SBATCH -p cpu
#SBATCH -w node2
##SBATCH --exclusive
#SBATCH --mem-per-cpu=3000 #--mem=120GB
##SBATCH --mem-per-cpu=1250 #--mem=120GB

export g16root=/data/zhaojy/g16/
export GAUSS_SCRDIR=/data/nglokwan/.scr/
export GAUSS_EXEDIR=/data/zhaojy/g16/g16/bsd:/share/home/zhaojy/g16/g16
source /data/zhaojy/g16/g16/bsd/g16.profile

#g16 < optfreq_S0.gjf > optfreq_S0.log

#g16 < opt.gjf > opt.log
#g16 < freq.gjf > freq.log
g16 < td.gjf > td.log
