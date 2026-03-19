#!/bin/bash
#SBATCH --job-name=funnel_aagp1-E3P
#SBATCH --nodes=1            
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --gres=gpu:1
#SBATCH -p rtx4090
##SBATCH -p main
##SBATCH -p workq


source /share/apps/intel/oneapi/setvars.sh
source /share/home/nglokwan/install/sandbox/plumed/mpi.sh
module load /share/home/nglokwan/install/sandbox/plumed/lib/plumed/modulefile
source /share/apps/plumed-gromacs2023.5-gpu_mpi/bin/GMXRC
#source ~/scripts/gmx2022.6_plumed_gpu_py3.env

which mpirun

for i in 0 ; do
        gmx_mpi grompp -v -f pr0.mdp -c pr_em.gro -p sys.top -o pr_${i}.tpr -maxwarn 3 -n index.ndx
        mpirun -np 1 gmx_mpi mdrun -deffnm pr_${i} -ntomp 8 -bonded gpu -nb gpu -update cpu -cpt 5 -plumed ../aagp1_fm0.dat
        gmx_mpi grompp -v -f eq.mdp -c pr_${i}.gro -p sys.top -o rec_eq_${i}.tpr -maxwarn 3 -n index.ndx
        mpirun -np 1 gmx_mpi mdrun -deffnm rec_eq_${i} -ntomp 8 -bonded gpu -nb gpu -update cpu -cpt 5 -plumed ../aagp1_fm0.dat
	gmx_mpi grompp -v -f prod_100ns.mdp -c rec_eq_${i}.gro -p sys.top -o fm_100ns_${i}.tpr -maxwarn 3 -n index.ndx 
	mpirun -np 1 gmx_mpi mdrun -deffnm fm_100ns_${i} -ntomp 8 -bonded gpu -nb gpu -update cpu -cpt 5 -plumed ../aagp1_fm0.dat  
done

echo

date
