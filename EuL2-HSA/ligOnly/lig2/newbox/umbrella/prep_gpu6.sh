#!/bin/bash

#for i in `seq 0 26`; do 
#	cd window$i
#	sed s/WINID/${i}/g ../prep_win_template.pbs > prep_win${i}.pbs
#	qsub prep_win${i}.pbs
#	cd ..
#done
i=$1
j=$2
cd window$i
sed s/WINID/${i}/g ../prep_win_template_gpu6.pbs > prep_win${i}_gpu.pbs
sed -i s/GPUID/${j}/g prep_win${i}_gpu.pbs
#qsub prep_win${i}_gpu.pbs
cd ..
