#!/bin/bash
#SBATCH --job-name=fm7_E92V_1.5us
#SBATCH --nodes=1            
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=22
#SBATCH --gres=gpu:1
##SBATCH -p rtx4090-short
##SBATCH -p main
##SBATCH -p workq
#SBATCH -p rtx4090
##SBATCH -p l40

#source /share/apps/intel/oneapi/setvars.sh
#source /data/nglokwan/scripts/gmx2023.5_plumed2.9.3_gpu_shared_ompi.env
source /data/nglokwan/scripts/gmx2023.5_plumed2.9.3_gpu_shared.env

#source ~/scripts/gmx2023.5_plumed_gpu_py3.env

#i=$1
i=2

date
#echo sleep $((i*n))
#sleep $((i*n))
#echo

#gmx mdrun -deffnm fm_1500ns_${i} -ntomp 22 -bonded gpu -nb gpu -update cpu -cpt 2 -plumed ../../aagp2_E92V_fm7_ext.dat -cpi fm_1500ns_${i}.cpt -o fm_1500ns_2.part0006.trr -x fm_1500ns_2.part0006.xtc -e fm_1500ns_2.part0006.edr -g fm_1500ns_2.part0006.log
#gmx mdrun -deffnm fm_1500ns_${i} -ntomp 22 -bonded gpu -nb gpu -update cpu -cpt 1 -plumed ../../aagp2_E92V_fm7_ext.dat -cpi fm_1500ns_${i}.cpt -o fm_1500ns_2.part0007.trr -x fm_1500ns_2.part0007.xtc -e fm_1500ns_2.part0007.edr -g fm_1500ns_2.part0007.log
#gmx mdrun -deffnm fm_1500ns_${i} -ntomp 22 -bonded gpu -nb gpu -update cpu -cpt 1 -plumed ../../aagp2_E92V_fm7_ext.dat -cpi fm_1500ns_${i}.cpt -o fm_1500ns_2.part0008.trr -x fm_1500ns_2.part0008.xtc -e fm_1500ns_2.part0008.edr -g fm_1500ns_2.part0008.log
gmx mdrun -deffnm fm_1500ns_${i} -ntomp 22 -bonded gpu -nb gpu -update cpu -cpt 1 -plumed ../../aagp2_E92V_fm7_ext.dat -cpi fm_1500ns_${i}.cpt -o fm_1500ns_2.part0009.trr -x fm_1500ns_2.part0009.xtc -e fm_1500ns_2.part0009.edr -g fm_1500ns_2.part0009.log

echo

date
