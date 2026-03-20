#!/bin/bash
for l in lig2m lig2 ; do
	cd $l
	for i in `seq 0 20`; do 
		cd window$i
		sed s/WINID/${i}/g ../../prep_win_template.pbs > prep_win${i}.pbs
		qsub prep_win${i}.pbs
		cd ..
	done
	cd ..
done
#i=$1
##j=$2
#cd window$i
#sed s/WINID/${i}/g ../prep_config1_template.pbs > prep_win${i}.pbs #_gpu.pbs
##sed s/WINID/${i}/g ../pr_config1_template.pbs > prep_win${i}.pbs #_gpu.pbs
##sed -i s/GPUID/${j}/g prep_win${i}_gpu.pbs
#qsub prep_win${i}.pbs #_gpu.pbs
#cd ..
