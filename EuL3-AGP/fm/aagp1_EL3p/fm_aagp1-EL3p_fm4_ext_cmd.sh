#!/bin/bash
#SBATCH --job-name=FnlMtD_400ns
#SBATCH --nodes=1            
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=22
#SBATCH --gres=gpu:1
##SBATCH -p rtx4090-short
##SBATCH -p main
##SBATCH -p workq
#SBATCH -p rtx4090

#source /share/apps/intel/oneapi/setvars.sh
#source /data/nglokwan/scripts/gmx2023.5_plumed2.9.3_gpu_shared_ompi.env
source /data/nglokwan/scripts/gmx2023.5_plumed2.9.3_gpu_shared.env

#source ~/scripts/gmx2023.5_plumed_gpu_py3.env

i=$1

#gmx grompp -v -f ../pr0.mdp -c ../pr_em.gro -p sys.top -o pr_${i}.tpr -maxwarn 3 -n ../index.ndx
#gmx mdrun -deffnm pr_${i} -ntomp 22 -bonded gpu -nb gpu -update cpu -cpt 5 -plumed ../../aagp1_fm4.dat
#gmx grompp -v -f ../eq.mdp -c pr_${i}.gro -p sys.top -o rec_eq_${i}.tpr -maxwarn 3 -n ../index.ndx
#gmx mdrun -deffnm rec_eq_${i} -ntomp 22 -bonded gpu -nb gpu -update cpu -cpt 5 -plumed ../../aagp1_fm4.dat
#gmx grompp -v -f ../prod_400ns.mdp -c rec_eq_${i}.gro -p sys.top -o fm_400ns_${i}.tpr -maxwarn 3 -n ../index.ndx 
#gmx mdrun -deffnm fm_400ns_${i} -ntomp 22 -bonded gpu -nb gpu -update cpu -cpt 5 -plumed ../../aagp1_fm4.dat  

#gmx mdrun -deffnm fm_500ns_${i} -ntomp 22 -bonded gpu -nb gpu -update cpu -cpt 5 -plumed ../../aagp1_fm4_ext.dat -cpi fm_400ns_${i}.cpt -noappend 
#gmx mdrun -deffnm fm_700ns_${i} -ntomp 22 -bonded gpu -nb gpu -update cpu -cpt 5 -plumed ../../aagp1_fm4_ext.dat -cpi fm_500ns_${i}.cpt -noappend 
gmx mdrun -deffnm fm_1000ns_${i} -ntomp 22 -bonded gpu -nb gpu -update cpu -cpt 5 -plumed ../../aagp1_fm4_ext.dat -cpi fm_700ns_${i}.cpt -noappend 


echo

date
