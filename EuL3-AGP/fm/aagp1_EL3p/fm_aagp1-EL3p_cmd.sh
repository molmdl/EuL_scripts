#!/bin/bash
##SBATCH --job-name=openmpi_funnel #aagp1-E3P
#SBATCH --job-name=nompi_funnel #aagp1-E3P
#SBATCH --nodes=1            
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=22
#SBATCH --gres=gpu:1
#SBATCH -p rtx4090-short
##SBATCH -p main
##SBATCH -p workq
##SBATCH -p rtx4090

#source /share/apps/intel/oneapi/setvars.sh
#source /data/nglokwan/scripts/gmx2023.5_plumed2.9.3_gpu_shared_ompi.env
source /data/nglokwan/scripts/gmx2023.5_plumed2.9.3_gpu_shared.env

#source ~/scripts/gmx2023.5_plumed_gpu_py3.env

for i in $1 ; do
#        gmx_mpi grompp -v -f ../pr0.mdp -c ../pr_em.gro -p sys.top -o pr_${i}.tpr -maxwarn 3 -n ../index.ndx
#        gmx_mpi mdrun -deffnm pr_${i} -ntomp 22 -bonded gpu -nb gpu -update cpu -cpt 5 -plumed ../../aagp1_fm2.dat
#        gmx_mpi grompp -v -f ../eq.mdp -c pr_${i}.gro -p sys.top -o rec_eq_${i}.tpr -maxwarn 3 -n ../index.ndx
#        gmx_mpi mdrun -deffnm rec_eq_${i} -ntomp 22 -bonded gpu -nb gpu -update cpu -cpt 5 -plumed ../../aagp1_fm2.dat
#	gmx_mpi grompp -v -f ../prod_100ns.mdp -c rec_eq_${i}.gro -p sys.top -o fm_100ns_${i}.tpr -maxwarn 3 -n ../index.ndx 
#	gmx_mpi mdrun -deffnm fm_100ns_${i} -ntomp 22 -bonded gpu -nb gpu -update cpu -cpt 5 -plumed ../../aagp1_fm2.dat  
        gmx grompp -v -f ../pr0.mdp -c ../pr_em.gro -p sys.top -o pr_${i}.tpr -maxwarn 3 -n ../index.ndx
        gmx mdrun -deffnm pr_${i} -ntomp 22 -bonded gpu -nb gpu -update cpu -cpt 5 -plumed ../../aagp1_fm2.dat
        gmx grompp -v -f ../eq.mdp -c pr_${i}.gro -p sys.top -o rec_eq_${i}.tpr -maxwarn 3 -n ../index.ndx
        gmx mdrun -deffnm rec_eq_${i} -ntomp 22 -bonded gpu -nb gpu -update cpu -cpt 5 -plumed ../../aagp1_fm2.dat
	gmx grompp -v -f ../prod_100ns.mdp -c rec_eq_${i}.gro -p sys.top -o fm_100ns_${i}.tpr -maxwarn 3 -n ../index.ndx 
	gmx mdrun -deffnm fm_100ns_${i} -ntomp 22 -bonded gpu -nb gpu -update cpu -cpt 5 -plumed ../../aagp1_fm2.dat  
done

echo

date
