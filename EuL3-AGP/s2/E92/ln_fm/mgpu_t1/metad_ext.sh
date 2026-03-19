#!/bin/bash
#SBATCH --job-name=Eu-E92_metad
#SBATCH --nodes=1            
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=22
#SBATCH --gres=gpu:1
#SBATCH -p rtx4090-short
##SBATCH -p main
##SBATCH -p workq
##SBATCH -p rtx4090
##SBATCH -p l40

#source /share/apps/intel/oneapi/setvars.sh
#source /data/nglokwan/scripts/gmx2023.5_plumed2.9.3_gpu_shared_ompi.env
source /data/nglokwan/scripts/gmx2023.5_plumed2.9.3_gpu_shared.env

#source ~/scripts/gmx2023.5_plumed_gpu_py3.env

i=$1
n=$((25+11+3))

#export CUDA_VISIBLE_DEVICES=$i

echo sleep $((i*n))
sleep $((i*n))
echo

gmx mdrun -deffnm metad_${i}_1us -ntomp 22 -bonded gpu -nb gpu -update cpu -cpt 5 -plumed ../metad_E92.dat  -cpi metad_${i}_1us.cpt --noappend

echo

date
