#!/bin/bash
RECNAME='2bxo'
#gmx_mpi solvate -cp test.gro -cs spc216.gro -o solv.gro -p sys.top 
#gmx_mpi grompp -f em.mdp -p sys.top -c solv.gro -o rec_addion.tpr -maxwarn 1
#echo SOL | gmx genion -s rec_addion -p sys.top -neutral -conc 0.15 -nname CL -pname NA -o ion.gro
gmx_mpi grompp -v -f em -c ion -p sys.top -o ${RECNAME}_lig_em #-n index.ndx
gmx_mpi mdrun -v -deffnm ${RECNAME}_lig_em -ntomp 1 -pin on -pinoffset 18
gmx_mpi grompp -v -f pr_pos -c ${RECNAME}_lig_em -p sys.top -o pr_pos_lig -r ${RECNAME}_lig_em #-n index.ndx
gmx_mpi mdrun -v -deffnm pr_pos_lig -ntomp 1 -pin on -pinoffset 18
gmx_mpi grompp -v -f em -c pr_pos_lig -p sys.top -o em_lig #-n index.ndx
gmx_mpi mdrun -v -deffnm em_lig -ntomp 1 -pin on -pinoffset 18
gmx_mpi grompp -f pr0.mdp -p sys.top -c em_lig.gro -o pr_0.tpr
gmx_mpi mdrun -v -deffnm pr_0 -ntomp 1 -pin on -pinoffset 18
