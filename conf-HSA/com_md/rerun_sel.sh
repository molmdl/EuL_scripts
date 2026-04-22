#!/bin/bash

CWD=`pwd`

for line in `cat list4rerun.txt` ; do
	echo "$line"
	sys=`echo ${line} | awk -F',' '{print $1}'`
	trial=`echo ${line} | awk -F',' '{print $2}'`
	cd $sys
	sbatch <<- EOF
#!/bin/bash
#SBATCH -J ${sys}_${trial}
#SBATCH -c 22
#SBATCH -n 1
##SBATCH -p rtx4090-short
#SBATCH -p rtx4090
##SBATCH -p workq
##SBATCH -p l40
#SBATCH --gres=gpu:1

source ~/scripts/gmxMMPBSA.env

i=${trial}

#gmx grompp -f em.mdp -c pr.gro -p sys.top -o pr_em.tpr
#gmx mdrun -v -deffnm pr_em -ntmpi 1 -ntomp 8
gmx grompp -f pr0.mdp -c pr.gro -p sys.top  -o pr_\${i}.tpr -n index.ndx -maxwarn 1
#gmx grompp -f pr0.mdp -c pr_em.gro -p sys.top  -o pr_\${i}.tpr
gmx mdrun -deffnm pr_\$i -ntmpi 1 -ntomp 22 -bonded gpu -nb gpu -update gpu -pme gpu -cpt 5

gmx grompp -f md.mdp -c pr_\${i}.gro -p sys.top   -o prod_\${i}.tpr -n index.ndx -maxwarn 2
gmx mdrun -deffnm prod_\$i -ntmpi 1 -ntomp 22 -bonded gpu -nb gpu -update gpu -pme gpu -cpt 5

mkdir mmpbsa_\${i}
cd mmpbsa_\${i}
echo processing trj of trial \${i} ...
echo -e  'Protein_Other \n Protein_Other'| gmx trjconv -s ../prod_\${i}.tpr -f ../prod_\${i}.xtc -pbc mol -ur compact -o pbc.xtc -center -n ../index.ndx
echo -e 'Protein_Other \n Protein_Other' | gmx trjconv -s ../prod_\${i}.tpr -pbc cluster -f pbc.xtc -o com_traj.xtc -n ../index.ndx
rm pbc.xtc
echo copying files into mmpbsa_\${i} ...
cp ../prod_\${i}.tpr com.tpr
cp ../../../scripts/com/mmpbsa.?? .
sbatch mmpbsa.sh 
cd .. 

EOF
	cd $CWD
done
