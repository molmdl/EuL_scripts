#!/bin/bash
#SBATCH -J sumhills
#SBATCH -n 1
#SBATCH -c 1
#SBATCH -p cpu

source ~/data/scripts/gmx2023.5_plumed2.9.3_gpu_shared.env

CWD=`pwd`

for s in E64 E92 ; do
	for i in {0..3} ; do
		if [[ -e "./${s}/mgpu_t${i}/HILLS" ]] ; then
			rm -r ${s}/mgpu_t${i}/sumhills 
			mkdir ${s}/mgpu_t${i}/sumhills 
			cd ${s}/mgpu_t${i}/sumhills 
			plumed sum_hills --hills ../HILLS --mintozero --stride 50000 --kt 2.494339 >> sumhills.log 2>> sumhills.log
			cd $CWD
			# no history
			rm -r ${s}/mgpu_t${i}/sumhills_nohist
			mkdir ${s}/mgpu_t${i}/sumhills_nohist
			cd ${s}/mgpu_t${i}/sumhills_nohist 
			plumed sum_hills --hills ../HILLS --mintozero --stride 100000 --kt 2.494339 --nohistory >> sumhills.log 2>> sumhills.log
		else
			echo HILLS not found for $s mgpu_t$i
		fi
	cd $CWD

	done
done
