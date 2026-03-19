#!/bin/bash
#SBATCH --job-name=E64V_pr
#SBATCH --nodes=1            
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=22
#SBATCH --gres=gpu:1
##SBATCH -p rtx4090-short
#SBATCH -p rtx4090
##SBATCH -p workq

#source /share/apps/intel/oneapi/setvars.sh
#source /share/apps/plumed-gromacs2023.5-gpu_mpi/bin/GMXRC

## Now try to minimize and do some pre-equilibration of the system
gmx grompp -v -f em.mdp -c ion.gro -p topol.top -o em.tpr 
gmx mdrun -deffnm em -ntomp 22 -ntmpi 1 
gmx grompp -v -f pr_pos.mdp -c em.gro -p topol.top -o pr_pos.tpr -r em.gro -maxwarn 3 
gmx mdrun -deffnm pr_pos -ntomp 22 -bonded gpu -nb gpu -update gpu -cpi 5 -ntmpi 1 
gmx grompp -v -f em.mdp -c pr_pos.gro -p topol.top -o pr_em.tpr -maxwarn 3 
gmx mdrun -deffnm pr_em -ntomp 22 -ntmpi 1 

for i in {0..3} ; do
	gmx grompp -v -f pr0.mdp -c pr_em.gro -p topol.top -o pr_${i}.tpr -maxwarn 3
	gmx mdrun -deffnm pr_${i} -ntomp 22 -bonded gpu -nb gpu -update gpu -cpi 5 -pme gpu -ntmpi 1 
	gmx grompp -v -f eq.mdp -c pr_${i}.gro -p topol.top -o rec_eq_${i}.tpr -maxwarn 3
	gmx mdrun -deffnm rec_eq_${i} -ntomp 22 -bonded gpu -nb gpu -update gpu -cpi 5 -pme gpu -ntmpi 1 
done

echo
echo date
