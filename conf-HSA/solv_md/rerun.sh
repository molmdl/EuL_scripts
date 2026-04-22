#!/bin/bash

CWD=`pwd`

for line in `cat topol_list.txt` ; do
	echo "$line"
	sys=`echo ${line} | awk -F',' '{print $1}'`
	topol=`echo ${line} | awk -F',' '{print $2}'`
	cd $sys
	for trial in {0..3} ; do
		sbatch <<- EOF
#!/bin/bash
#SBATCH -J ${sys}_${trial}
#SBATCH -c 8
#SBATCH -n 1
##SBATCH -p rtx4090-short
##SBATCH -p rtx4090
#SBATCH -p workq
##SBATCH -p l40
##SBATCH --gres=gpu:1

source ~/scripts/gmxMMPBSA.env

i=${trial}

gmx grompp -f ${CWD}/rerun.mdp -c pr_\${i}.gro -p $topol -o rerun_\${i}.tpr  -maxwarn 2
gmx mdrun -deffnm rerun_\$i -ntmpi 1 -ntomp 8 -bonded cpu -nb cpu -update cpu -pme cpu -cpt 5 -rerun prod_\${i}.xtc

EOF
	done
	cd $CWD
done
