#!/bin/bash
#SBATCH -J prod_4x100ns
#SBATCH -c 8
#SBATCH -n 1
#SBATCH -p workq
#SBATCH --gres=gpu:1

for i in 0 1 2 3 ; do
	gmx grompp -f md.mdp -c pr_${i}.gro -p lig.top   -o prod_${i}.tpr
	gmx mdrun -deffnm prod_$i -ntmpi 1 -ntomp 8 -bonded gpu -nb gpu -update gpu -pme gpu -cpt 5
done

