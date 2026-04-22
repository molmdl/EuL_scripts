#!/bin/bash
#SBATCH -J pr
#SBATCH -c 22
#SBATCH -n 1
#SBATCH -p workq
#SBATCH --gres=gpu:1

source ~/data/scripts/gmx2023.5_plumed2.9.3_gpu_shared.env
gmx mdrun -deffnm pr -ntmpi 1 -ntomp 22 -bonded gpu -nb gpu -update gpu -pme gpu -cpt 5 
#gmx grompp -f em.mdp -c pr.gro -p sys.top -o pr_em.tpr
#gmx mdrun -v -deffnm pr_em -ntmpi 1 -ntomp 8
for i in 0 1 2 3 ; do
	gmx grompp -f pr0.mdp -c pr.gro -p sys.top  -o pr_${i}.tpr -n index.ndx -maxwarn 1
	#gmx grompp -f pr0.mdp -c pr_em.gro -p sys.top  -o pr_${i}.tpr
	gmx mdrun -deffnm pr_$i -ntmpi 1 -ntomp 22 -bonded gpu -nb gpu -update gpu -pme gpu -cpt 5
done

