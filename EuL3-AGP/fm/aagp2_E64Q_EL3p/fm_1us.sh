#!/bin/bash
#SBATCH --job-name=fm7_E64Q_1us
#SBATCH --nodes=1            
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=22
#SBATCH --gres=gpu:1
##SBATCH -p rtx4090-short
##SBATCH -p main
##SBATCH -p workq
##SBATCH -p rtx4090
#SBATCH -p l40
##SBATCH -w gpu1

#source /share/apps/intel/oneapi/setvars.sh
#source /data/nglokwan/scripts/gmx2023.5_plumed2.9.3_gpu_shared_ompi.env
source /data/nglokwan/scripts/gmx2023.5_plumed2.9.3_gpu_shared.env

#source ~/scripts/gmx2023.5_plumed_gpu_py3.env

i=$1
n=$((25+11+3))

echo sleep $((i*n))
sleep $((i*n))
echo

#export CUDA_VISIBLE_DEVICES=$i

cp ../reference.pdb .
cp ../sys_fm.top sys.top

gmx grompp -v -f ../pr0.mdp -c ../pr_em.gro -p sys.top -o pr_${i}.tpr -maxwarn 3 -n ../index.ndx -r ../em.gro
gmx mdrun -deffnm pr_${i} -ntomp 22 -bonded gpu -nb gpu -update cpu -cpt 5 -plumed ../../aagp2_E64Q_fm7.dat
gmx grompp -v -f ../eq.mdp -c pr_${i}.gro -p sys.top -o rec_eq_${i}.tpr -maxwarn 3 -n ../index.ndx -r ../em.gro 
gmx mdrun -deffnm rec_eq_${i} -ntomp 22 -bonded gpu -nb gpu -update cpu -cpt 5 -plumed ../../aagp2_E64Q_fm7.dat
gmx grompp -v -f ../prod_1000ns.mdp -c rec_eq_${i}.gro -p sys.top -o fm_1000ns_${i}.tpr -maxwarn 3 -n ../index.ndx -r ../em.gro 
gmx mdrun -deffnm fm_1000ns_${i} -ntomp 22 -bonded gpu -nb gpu -update cpu -cpt 5 -plumed ../../aagp2_E64Q_fm7.dat  

echo

date
