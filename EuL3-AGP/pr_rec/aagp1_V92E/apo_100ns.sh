#!/bin/bash
#SBATCH --job-name=V92E_apo
#SBATCH --nodes=1            
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=22
#SBATCH --gres=gpu:1
#SBATCH -p rtx4090-short
##SBATCH -p rtx4090
##SBATCH -p workq

#source /share/apps/intel/oneapi/setvars.sh
#source /share/apps/plumed-gromacs2023.5-gpu_mpi/bin/GMXRC
source ~/data/scripts/gmx2023.5_plumed2.9.3_gpu_shared.env

i=$1

gmx grompp -v -f apo_100ns.mdp -c rec_eq_${i}.gro -p topol.top -o apo_100ns_${i}.tpr -maxwarn 3
gmx mdrun -deffnm apo_100ns_${i} -ntmpi 1 -ntomp 22 -bonded gpu -nb gpu -update gpu -cpi 5

echo
echo date
