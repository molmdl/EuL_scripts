#!/bin/bash
#gmx grompp -v -f local_fm_test.mdp -p sys.top -c pr_0.gro -o local_fm_test.tpr
#gmx_mpi mdrun -v -deffnm local_fm_test -plumed plumed.dat -ntomp 16 -pinoffset 0 -pin on -nb gpu -bonded gpu
#gmx_mpi mdrun -deffnm fm_50ns -plumed plumed.dat -pinoffset 0 -pin on -nb gpu -bonded gpu -cpi fm_50ns.cpt -ntomp 12  #-v 
#gmx_mpi grompp -v -f fm_50ns.mdp -p sys.top -c pr_0.gro -o fm_50ns.tpr
#gmx_mpi mdrun -v -deffnm fm_50ns -plumed plumed.dat -ntomp 12 -pinoffset 0 -pin on -nb gpu -bonded gpu
gmx_mpi mdrun -deffnm fm_100ns_0 -plumed plumed.dat -pinoffset 0 -pin on -nb gpu -bonded gpu -ntomp 12  #-v 
