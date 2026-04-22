#!/bin/bash
gmx editconf -f com.gro -o box.gro -d 1 -bt dodecahedron -c
gmx solvate -p sys.top -cp box.gro -cs spc216 -o solv.gro
gmx grompp -f em.mdp -c solv.gro -p sys.top -o ion.tpr -maxwarn 2
echo SOL | gmx genion -s ion.tpr -p sys.top -neutral -nname CL -pname  NA -conc 0.15 -o ion.gro
gmx grompp -f em.mdp -c ion.gro -p sys.top -o em.tpr -maxwarn 2 
echo -e '"Protein" | "Other"  \n q' | gmx make_ndx -f em.tpr 
#gmx make_ndx -f lig.gro 
#gmx genrestr -f lig.gro 
gmx mdrun -deffnm em -ntmpi 1 -ntomp 8
gmx grompp -f pr.mdp -c em.gro -p sys.top  -r em.gro -o pr.tpr -maxwarn 2 
