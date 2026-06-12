#!/bin/bash
#for s in sssD sssL rrrL rrrD ; do for j in sap tsap ; do
#for s in sssD ; do for j in tsap ; do
#for s in sssD sssL rrrL ; do
#for s in phe_sssL_sap ; do
for s in ibp ibp_r ; do
	#cd phe_${s}_tsap
	#cd me_${s}_$j
	cd $s
	#python ../../scripts/com/bypass_angle_type3.py me_${s}_${j}.itp sys.top
	#python ../../scripts/com/bypass_angle_type3.py phe_${s}_tsap.itp sys.top
	#python ../../scripts/com/bypass_angle_type3.py lig_me_g.itp sys.top
	#python ../../scripts/com/bypass_angle_type3.py me_${s}.itp sys.top
	python ../../scripts/com/bypass_angle_type3.py ${s}.itp sys.top
	bash ../../scripts/com/trj4mmpbsa.sh 
	for i in {0..3} ; do 
		cd mmpbsa_$i 
		sbatch mmpbsa.sh 
		cd .. 
	done
	cd ..
done
