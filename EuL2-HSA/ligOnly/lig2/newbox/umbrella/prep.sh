#!/bin/bash

for i in {0..20} ; do 
	cd window$i
	if [[ -e prep_win${i}.pbs ]] ; then
		mv prep_win${i}.pbs prep_${i}.pbs.bck
	fi
	sed s/WINID/${i}/g ../prep_config1_template.pbs > prep_win${i}.pbs
#	sed s/WINID/${i}/g ../prep_win_template.pbs > prep_win${i}.pbs	
	qsub prep_win${i}.pbs
#	qsub prep_win${i}_gpu.pbs
	cd ..
done
#i=$1
#j=$2
#cd window$i
#sed s/WINID/${i}/g ../prep_win_template.pbs > prep_win${i}_gpu.pbs
#sed -i s/GPUID/${j}/g prep_win${i}_gpu.pbs
##qsub prep_win${i}_gpu.pbs
#cd ..
