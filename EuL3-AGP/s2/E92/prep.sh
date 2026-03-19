#!/bin/bash

#source ~/scripts/gmx2023.5_plumed2.9.3_gpu_shared.env
#gmx editconf -f ../ini/agp2_min1_dm1_f8_dry.gro -c -d 1 -o box.gro
gmx solvate -cp box.gro -cs spc216 -p sys.top -o solv.gro
gmx grompp -f ../em.mdp -c solv.gro -p sys.top -o ion.tpr -maxwarn 1
echo SOL | gmx genion -s ion.tpr  -p sys.top  -neutral -conc 0.15 -nname CL -pname NA -o ion.gro
