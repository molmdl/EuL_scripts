#!/bin/bash

source ~/data/scripts/gmx2023.5_plumed2.9.3_gpu_shared.env

CWD=`pwd`

#for s in aagp1 aagp2 aagp1_V92E aagp2_E92V ; do 
#for s in aagp1_V92E aagp2_E92V ; do 
for s in aagp2_E92Q aagp2_E64V aagp2_E64Q ; do 
#for s in aagp2 ; do 
	#for i in {0..7} ; do
	for i in {0..3} 0_1 1_1 2_1 3_1 ; do
		if [[ -e "./${s}_EL3p/mgpu_t${i}/HILLS" ]] ; then
			rm -r ${s}_EL3p/mgpu_t${i}/sumhills 
			mkdir ${s}_EL3p/mgpu_t${i}/sumhills 
			cd ${s}_EL3p/mgpu_t${i}/sumhills 
			plumed sum_hills --hills ../HILLS --mintozero --stride 50000 --kt 2.494339 >> sumhills.log 2>> sumhills.log
		else
			echo HILLS not found for $s mgpu_t$i
		fi
		cd $CWD
	done
done
