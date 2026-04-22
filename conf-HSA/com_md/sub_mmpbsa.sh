#!/bin/bash
#for s in sssD sssL rrrL rrrD ; do
for s in sssD sssL rrrL ; do
	#cd phe_${s}_tsap
	cd me_${s}_sap
	#python ../../scripts/com/bypass_angle_type3.py phe_${s}_tsap.itp sys.top
	python ../../scripts/com/bypass_angle_type3.py lig_me_g.itp sys.top
	#python ../../scripts/com/bypass_angle_type3.py me_${s}.itp sys.top
	bash ../../scripts/com/trj4mmpbsa.sh 
	for i in {0..3} ; do 
		cd mmpbsa_$i 
		sbatch mmpbsa.sh 
		cd .. 
	done
	cd ..
done
