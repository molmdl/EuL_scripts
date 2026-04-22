#!/bin/bash
#SBATCH -J me_rrrL_sap_1
#SBATCH -c 22
#SBATCH -n 1
##SBATCH -p rtx4090-short
#SBATCH -p rtx4090
##SBATCH -p workq
##SBATCH -p l40
#SBATCH --gres=gpu:1

source ~/scripts/gmxMMPBSA.env

i=1

#gmx grompp -f em.mdp -c pr.gro -p sys.top -o pr_em.tpr
#gmx mdrun -v -deffnm pr_em -ntmpi 1 -ntomp 8
gmx grompp -f pr0.mdp -c pr.gro -p sys.top  -o pr_${i}.tpr -n index.ndx -maxwarn 1
#gmx grompp -f pr0.mdp -c pr_em.gro -p sys.top  -o pr_${i}.tpr
gmx mdrun -deffnm pr_$i -ntmpi 1 -ntomp 22 -bonded gpu -nb gpu -update gpu -pme gpu -cpt 5

gmx grompp -f md.mdp -c pr_${i}.gro -p sys.top   -o prod_${i}.tpr -n index.ndx -maxwarn 2
gmx mdrun -deffnm prod_$i -ntmpi 1 -ntomp 22 -bonded gpu -nb gpu -update gpu -pme gpu -cpt 5
