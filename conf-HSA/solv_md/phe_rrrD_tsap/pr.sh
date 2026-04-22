#!/bin/bash
#SBATCH -J pr
#SBATCH -c 8
#SBATCH -n 1
#SBATCH -p workq
##SBATCH -p l40
#SBATCH --gres=gpu:1

gmx mdrun -deffnm pr -ntmpi 1 -ntomp 8 -bonded gpu -nb gpu -update gpu -pme gpu -cpt 5
#gmx mdrun -deffnm pr_0 -ntmpi 1 -ntomp 8 -bonded gpu -nb gpu -update gpu -pme gpu -cpt 5
for i in 0 1 2 3 ; do
	gmx grompp -f pr0.mdp -c pr.gro -p phe_rrrD_tsap.top  -o pr_${i}.tpr
	gmx mdrun -deffnm pr_$i -ntmpi 1 -ntomp 8 -bonded gpu -nb gpu -update gpu -pme gpu -cpt 5
done

