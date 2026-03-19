#!/bin/bash
CWD=`pwd`
for f in cv-time funnel2D-time ; do
	for i in `cat list.txt`; do 
		cd $i 
		inkscape.com ${f}.svg -o ${f}.png -b white
		cd $CWD 
	done
done
