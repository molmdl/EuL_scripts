#!/bin/bash
#SBATCH --job-name=fm7_E64Q_1.5us
#SBATCH --nodes=1            
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=22
##SBATCH --gres=gpu:1
##SBATCH -p rtx4090-short
##SBATCH -p main
##SBATCH -p workq
##SBATCH -p rtx4090
#SBATCH -p l40

#source /share/apps/intel/oneapi/setvars.sh
#source /data/nglokwan/scripts/gmx2023.5_plumed2.9.3_gpu_shared_ompi.env
source /data/nglokwan/scripts/gmx2023.5_plumed2.9.3_gpu_shared.env

#source ~/scripts/gmx2023.5_plumed_gpu_py3.env

export CUDA_VISIBLE_DEVICES=2

i=0
echo
gmx mdrun -deffnm fm_1500ns_${i} -ntomp 22 -bonded gpu -nb gpu -update cpu -cpt 5 -plumed ../../aagp2_E64Q_fm7_ext.dat  -cpi fm_1000ns_${i}.cpt -noappend

echo

date
