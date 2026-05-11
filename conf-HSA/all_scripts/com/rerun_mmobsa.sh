#!/bin/bash
for line in `cat rerun_sel_top.txt` ; do
	echo "$line"
	sys=`echo ${line} | awk -F',' '{print $1}'`
	itp=`echo ${line} | awk -F',' '{print $2}'`
	#cd phe_${s}_tsap
	cd $sys
#	python ../../scripts/com/bypass_angle_type3.py $itp sys.top
#	bash ../../scripts/com/trj4mmpbsa.sh 
	for i in {0..3} ; do 
		cd mmpbsa_$i 
		cp ../../../scripts/com/mmpbsa.?? .
		sbatch mmpbsa.sh 
		cd .. 
	done
	cd ..
done
