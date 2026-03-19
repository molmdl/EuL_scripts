#!/bin/bash
#SBATCH --job-name=pr
#SBATCH --nodes=1            
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=22
#SBATCH --gres=gpu:1
#SBATCH -p rtx4090-short
##SBATCH -p rtx4090
##SBATCH -p workq

source /data/nglokwan/scripts/gmx2023.5_plumed2.9.3_gpu_shared.env
#source /share/apps/intel/oneapi/setvars.sh
#source /share/apps/plumed-gromacs2023.5-gpu_mpi/bin/GMXRC

## Now try to minimize and do some pre-equilibration of the system
gmx grompp -v -f em.mdp -c ion.gro -p sys.top -o em.tpr 
gmx mdrun -deffnm em -ntomp 22 -ntmpi 1
gmx grompp -v -f pr_pos.mdp -c em.gro -p sys.top -o pr_pos.tpr -r em.gro -maxwarn 3  -n index.ndx
gmx mdrun -deffnm pr_pos -ntomp 22 -bonded gpu -nb gpu -update gpu -cpt 5  -ntmpi 1 
gmx grompp -v -f em.mdp -c pr_pos.gro -p sys.top -o pr_em.tpr -maxwarn 3 -n index.ndx 
gmx mdrun -deffnm pr_em -ntomp 22 -ntmpi 1 

##for i in {0..3} ; do
#for i in 0 ; do
#	gmx grompp -v -f pr0.mdp -c pr_em.gro -p sys.top -o pr_${i}.tpr -maxwarn 3 -n index.ndx 
#	mpirun -np 1 gmx mdrun -deffnm pr_${i} -ntomp 22 -bonded gpu -nb gpu -update cpu -cpt 5 -plumed ../aagp1_fm0.dat  
#	gmx grompp -v -f eq.mdp -c pr_${i}.gro -p sys.top -o rec_eq_${i}.tpr -maxwarn 3 -n index.ndx 
#	mpirun -np 1 gmx mdrun -deffnm rec_eq_${i} -ntomp 22 -bonded gpu -nb gpu -update cpu -cpt 5 -plumed ../aagp1_fm0.dat  
#done
#
##gmx mdrun -ntomp 22 -deffnm fm_300ns_${i} -plumed ../aagp1_fm0.dat -bonded gpu -nb gpu -cpt 5 -update cpu

echo

date
