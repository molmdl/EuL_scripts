#!/bin/bash
#SBATCH --job-name=funnel #aagp1-E3P
#SBATCH --nodes=1            
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=22
#SBATCH --gres=gpu:1
#SBATCH -p rtx4090-short
##SBATCH -p main
##SBATCH -p workq
##SBATCH -p rtx4090

#source /share/apps/intel/oneapi/setvars.sh
#source /data/nglokwan/scripts/gmx2023.5_plumed2.9.3_gpu_shared_ompi.env
source /data/nglokwan/scripts/gmx2023.5_plumed2.9.3_gpu_shared.env

#source ~/scripts/gmx2023.5_plumed_gpu_py3.env

for i in $1 ; do
	gmx mdrun -deffnm fm_400ns_${i} -ntomp 22 -bonded gpu -nb gpu -update cpu -cpt 5 -plumed ../../aagp1_fm2_ext.dat -cpi fm_100ns_${i}.cpt -noappend 
done

echo

date
