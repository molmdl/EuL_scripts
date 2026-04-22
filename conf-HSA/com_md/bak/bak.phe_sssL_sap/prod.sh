#!/bin/bash
#SBATCH -J prod_4x100ns
#SBATCH -c 22
#SBATCH -n 1
#SBATCH -p workq
##SBATCH -p rtx4090-short
##SBATCH -p l40
#SBATCH --gres=gpu:1

i=$1

gmx grompp -f md.mdp -c pr_${i}.gro -p sys.top   -o prod_${i}.tpr -n index.ndx -maxwarn 2
gmx mdrun -deffnm prod_$i -ntmpi 1 -ntomp 22 -bonded gpu -nb gpu -update gpu -pme gpu -cpt 5

