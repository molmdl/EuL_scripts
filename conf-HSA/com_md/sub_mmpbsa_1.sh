#!/bin/bash
for s in `cat list.txt` ; do
	#cd phe_${s}_tsap
	cd $s
	for i in {0..3} ; do 
		cd mmpbsa_$i 
		cp ../../../scripts/com/mmpbsa.?? .
		for b in *SA* ; do
			mv $b bak.$b
		done
		sbatch mmpbsa.sh 
		cd .. 
	done
	cd ..
done
