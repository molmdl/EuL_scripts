#!/bin/bash
#SBATCH -J EuL3_opt
##SBATCH -J tddft
##SBATCH -J 1scf_cam-b3lyp_cc-pvtz
#SBATCH -n 1
#SBATCH -c 48
#SBATCH -p cpu
#SBATCH -w node2
#SBATCH --mem-per-cpu=2500 #--mem=120GB
##SBATCH --mem-per-cpu=1250 #--mem=120GB

export g16root=/data/zhaojy/g16/
export GAUSS_SCRDIR=/data/nglokwan/.scr/
export GAUSS_EXEDIR=/data/zhaojy/g16/g16/bsd:/share/home/zhaojy/g16/g16
source /data/zhaojy/g16/g16/bsd/g16.profile

#g16 < optfreq_S0.gjf > optfreq_S0.log
#g16 < cam-b3lyp_opt.gjf > cam-b3lyp_opt.log
#g16 < cam-b3lyp_cc-pvtz_1scf.gjf > cam-b3lyp_cc-pvtz_1scf.log
#g16 < cam-b3lyp_optfreq.gjf > cam-b3lyp_optfreq.log

g16 < cam-b3lyp_td.gjf > cam-b3lyp_td.log

#ground state SP stable=opt

#g16 < opt_S1.gjf > opt_S1.log
