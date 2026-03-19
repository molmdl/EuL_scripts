#!/bin/bash
#SBATCH --job-name=apo
#SBATCH --nodes=1            
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=22
#SBATCH --gres=gpu:1
##SBATCH -p rtx4090-short
##SBATCH -p rtx4090
#SBATCH -p workq

source /share/apps/intel/oneapi/setvars.sh
source /share/apps/plumed-gromacs2023.5-gpu_mpi/bin/GMXRC

i=$1

gmx_mpi grompp -v -f apo_100ns.mdp -c rec_eq_${i}.gro -p topol.top -o apo_100ns_${i}.tpr -maxwarn 3
gmx_mpi mdrun -deffnm apo_100ns_${i} -ntomp 22 -bonded gpu -nb gpu -update gpu -cpi 5

echo
echo date
