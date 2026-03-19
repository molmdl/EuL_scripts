#!/bin/bash
#SBATCH --job-name=fm7_AGP1_1us
#SBATCH --nodes=1            
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=22
#SBATCH --gres=gpu:1
##SBATCH -p rtx4090-short
##SBATCH -p main
#SBATCH -p workq
##SBATCH -p rtx4090
##SBATCH -p l40

#source /share/apps/intel/oneapi/setvars.sh
#source /data/nglokwan/scripts/gmx2023.5_plumed2.9.3_gpu_shared_ompi.env
source /data/nglokwan/scripts/gmx2023.5_plumed2.9.3_gpu_shared.env

#source ~/scripts/gmx2023.5_plumed_gpu_py3.env

i=$1

echo

gmx mdrun -deffnm fm_1000ns_${i} -ntomp 22 -bonded gpu -nb gpu -update cpu -cpt 1 -plumed ../../aagp1_fm7_ext.dat  -cpi fm_1000ns_${i}.cpt -nsteps 50000 # -noappend
gmx mdrun -deffnm fm_1000ns_${i} -ntomp 22 -bonded gpu -nb gpu -update cpu -cpt 1 -plumed ../../aagp1_fm7_ext.dat  -cpi fm_1000ns_${i}.cpt -nsteps 50000 # -noappend
gmx mdrun -deffnm fm_1000ns_${i} -ntomp 22 -bonded gpu -nb gpu -update cpu -cpt 1 -plumed ../../aagp1_fm7_ext.dat  -cpi fm_1000ns_${i}.cpt -nsteps 50000 # -noappend
gmx mdrun -deffnm fm_1000ns_${i} -ntomp 22 -bonded gpu -nb gpu -update cpu -cpt 1 -plumed ../../aagp1_fm7_ext.dat  -cpi fm_1000ns_${i}.cpt -nsteps 50000 # -noappend
gmx mdrun -deffnm fm_1000ns_${i} -ntomp 22 -bonded gpu -nb gpu -update cpu -cpt 1 -plumed ../../aagp1_fm7_ext.dat  -cpi fm_1000ns_${i}.cpt -nsteps 50000 # -noappend
gmx mdrun -deffnm fm_1000ns_${i} -ntomp 22 -bonded gpu -nb gpu -update cpu -cpt 1 -plumed ../../aagp1_fm7_ext.dat  -cpi fm_1000ns_${i}.cpt -nsteps 50000 # -noappend
rm fm_1000ns_${i}.gro \#f*gro*
gmx mdrun -deffnm fm_1000ns_${i} -ntomp 22 -bonded gpu -nb gpu -update cpu -cpt 1 -plumed ../../aagp1_fm7_ext.dat  -cpi fm_1000ns_${i}.cpt # -noappend

echo

date
