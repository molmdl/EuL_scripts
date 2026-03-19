#!/bin/bash
##SBATCH --job-name=fm6_AGP2_1us
#SBATCH --job-name=fm6_AGP2_1.5us
#SBATCH --nodes=1            
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=22
#SBATCH --gres=gpu:1
##SBATCH -p rtx4090-short
##SBATCH -p main
##SBATCH -p workq
##SBATCH -p rtx4090
#SBATCH -p l40

#source /share/apps/intel/oneapi/setvars.sh
#source /data/nglokwan/scripts/gmx2023.5_plumed2.9.3_gpu_shared_ompi.env
source /data/nglokwan/scripts/gmx2023.5_plumed2.9.3_gpu_shared.env

#source ~/scripts/gmx2023.5_plumed_gpu_py3.env

i=2

gmx mdrun -deffnm fm_1500ns_${i} -ntomp 22 -bonded gpu -nb gpu -update cpu -cpt 5 -plumed ../../aagp2_fm6_ext.dat -cpi fm_1000ns_${i}.cpt -noappend 
#gmx mdrun -deffnm fm_1500ns_${i} -ntomp 22 -bonded gpu -nb gpu -update cpu -cpt 2 -plumed ../../aagp2_fm6_ext.dat -cpi fm_1500ns_${i}.cpt -nsteps 20000 #-noappend 
#gmx mdrun -deffnm fm_1500ns_${i} -ntomp 22 -bonded gpu -nb gpu -update cpu -cpt 2 -plumed ../../aagp2_fm6_ext.dat -cpi fm_1500ns_${i}.cpt -nsteps 20000 #-noappend 
#gmx mdrun -deffnm fm_1500ns_${i} -ntomp 22 -bonded gpu -nb gpu -update cpu -cpt 2 -plumed ../../aagp2_fm6_ext.dat -cpi fm_1500ns_${i}.cpt -nsteps 20000 #-noappend 
#gmx mdrun -deffnm fm_1500ns_${i} -ntomp 22 -bonded gpu -nb gpu -update cpu -cpt 2 -plumed ../../aagp2_fm6_ext.dat -cpi fm_1500ns_${i}.cpt -nsteps 20000 #-noappend 
#gmx mdrun -deffnm fm_1500ns_${i} -ntomp 22 -bonded gpu -nb gpu -update cpu -cpt 2 -plumed ../../aagp2_fm6_ext.dat -cpi fm_1500ns_${i}.cpt -nsteps 20000 #-noappend 
#rm fm_1500ns_${i}.gro \#fm_1500ns_${i}.gro.*
#gmx mdrun -deffnm fm_1500ns_${i} -ntomp 22 -bonded gpu -nb gpu -update cpu -cpt 2 -plumed ../../aagp2_fm6_ext.dat -cpi fm_1500ns_${i}.cpt #-noappend 

echo

date
