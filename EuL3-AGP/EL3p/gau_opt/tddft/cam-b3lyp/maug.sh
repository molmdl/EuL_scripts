#!/bin/bash
##SBATCH -J EuL3_opt
#SBATCH -J EuL3_opt
##SBATCH -J 1scf_cam-b3lyp_cc-pvtz
#SBATCH -n 1
#SBATCH -c 48
#SBATCH -p cpu
#SBATCH --mem=120GB
##SBATCH --mem-per-cpu=1250 #--mem=120GB

export g16root=/data/zhaojy/g16/
export GAUSS_SCRDIR=/data/nglokwan/.scr/
export GAUSS_EXEDIR=/data/zhaojy/g16/g16/bsd:/share/home/zhaojy/g16/g16
source /data/zhaojy/g16/g16/bsd/g16.profile

#g16 < optfreq_S0.gjf > optfreq_S0.log
#g16 < cam-b3lyp_opt.gjf > cam-b3lyp_opt.log
#g16 < cam-b3lyp_cc-pvtz_1scf.gjf > cam-b3lyp_cc-pvtz_1scf.log

#g16 < maug_SP_stableopt.gjf > maug_SP_stableopt.log
#
#g16 < maug_td.gjf > maug_td.log
g16 < maug_s1_opt.gjf > maug_s1_opt.log

