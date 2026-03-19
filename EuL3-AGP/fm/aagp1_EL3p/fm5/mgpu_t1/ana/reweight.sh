#!/bin/bash
#SBATCH --job-name=plumed_reweight
#SBATCH --nodes=1            
#SBATCH --ntasks=16
#SBATCH --cpus-per-task=1
#SBATCH -p cpu
##SBATCH -p workq
##SBATCH -p rtx4090

#source /share/apps/intel/oneapi/setvars.sh
#source /data/nglokwan/scripts/gmx2023.5_plumed2.9.3_gpu_shared_ompi.env
source /data/nglokwan/scripts/gmx2023.5_plumed2.9.3_gpu_shared.env

mpirun -np 16 plumed driver --mf_xtc ../fm_1000ns_?.xtc --plumed reweight.dat --mc ../../mcfile --kt 2.494339 --trajectory-stride 20 >> plumed_reweight.log 2>> plumed_reweight.log

