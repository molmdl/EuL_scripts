#!/bin/bash
CWD=`pwd`
for i in `cat list_2cv.txt`; do 
	cd $i 
	inkscape.com cv-time.svg -o cv-time.png -b white
	inkscape.com funnel2D-time.svg -o funnel2D-time.png -b white
	cd $CWD 
done
