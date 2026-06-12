#!/bin/bash

for i in `cat $1 ` ; do
	cd $i
	sbatch <<- EOF
#!/bin/bash
#SBATCH -J xtb_omd
#SBATCH -n 1
#SBATCH -c 8
#SBATCH -p cpu

source /share/apps/intel/oneapi/setvars.sh
source ~/scripts/xtb.env
#export MKL_NUM_THREADS=16
#export OMP_NUM_THREADS=16

xtb input.xyz -omd --input ../xtb_md.inp --gfn2 --chrg 3 --gbsa h2o  > xtb_omd.log # --gbsa h2o  
EOF
	cd ..
done
