#!/bin/bash
CWD=`pwd`

for line in `cat list4rerun.txt` ; do
        echo "${line}"
	sys=`echo ${line} | awk -F',' '{print $1}'`
	trial=`echo ${line} | awk -F',' '{print $2}'`
        cd $sys
	mkdir -p bak
	mv p*_${trial}* mmpbsa_$trial ./bak
	cd $CWD
done
