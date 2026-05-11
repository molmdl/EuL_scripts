#!/bin/bash
#SBATCH -J prep
#SBATCH -n 1
#SBATCH -c 22
##SBATCH -p workq
##SBATCH -p l40
#SBATCH -p rtx4090-short
#SBATCH --gres=gpu:1

gmx editconf -f com.gro -o box.gro -d 1 -bt dodecahedron -c
gmx solvate -p sys.top -cp box.gro -cs spc216 -o solv.gro
gmx grompp -f em.mdp -c solv.gro -p sys.top -o ion.tpr -maxwarn 2
echo SOL | gmx genion -s ion.tpr -p sys.top -neutral -nname CL -pname  NA -conc 0.15 -o ion.gro
gmx grompp -f em.mdp -c ion.gro -p sys.top -o em.tpr -maxwarn 2 
echo -e '"Protein" | "Other"  \n q' | gmx make_ndx -f em.tpr 
#gmx make_ndx -f lig.gro 
#gmx genrestr -f lig.gro 
gmx mdrun -deffnm em -ntmpi 1 -ntomp 22
gmx grompp -f pr.mdp -c em.gro -p sys.top  -r em.gro -o pr.tpr -maxwarn 2 
gmx mdrun -deffnm pr -ntmpi 1 -ntomp 22 -bonded gpu -nb gpu -update gpu -pme gpu -cpt 5 

