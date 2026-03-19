#!/bin/bash
#SBATCH --job-name=pr_metad
#SBATCH --nodes=1            
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=22
#SBATCH --gres=gpu:1
##SBATCH -p rtx4090-short
##SBATCH -p rtx4090
#SBATCH -p workq

source ~/scripts/gmx2023.5_plumed2.9.3_gpu_shared.env

## Now try to minimize and do some pre-equilibration of the system
gmx grompp -v -f ../em.mdp -c ion.gro -p sys.top -o em.tpr 
gmx mdrun -deffnm em -ntmpi 1 -ntomp 22
gmx grompp -v -f ../pr_pos.mdp -c em.gro -p sys.top -o pr_pos.tpr -r em.gro -maxwarn 3 
gmx mdrun -deffnm pr_pos -ntomp 22 -bonded gpu -nb gpu -update gpu -ntmpi 1 -cpi 5
gmx grompp -v -f ../em.mdp -c pr_pos.gro -p sys.top -o pr_em.tpr -maxwarn 3 
gmx mdrun -deffnm pr_em -ntmpi 1 -ntomp 22

for i in {0..3} ; do
	gmx grompp -v -f ../pr0.mdp -c pr_em.gro -p sys.top -o pr_${i}.tpr -maxwarn 3
	gmx mdrun -deffnm pr_${i} -ntmpi 1 -ntomp 22 -bonded gpu -nb gpu -update gpu -pme gpu -cpi 5
	gmx grompp -v -f ../eq.mdp -c pr_${i}.gro -p sys.top -o eq_${i}.tpr -maxwarn 3
	gmx mdrun -deffnm eq_${i} -ntmpi 1 -ntomp 22 -bonded gpu -nb gpu -update gpu -pme gpu -cpi 5
done

echo
echo date
