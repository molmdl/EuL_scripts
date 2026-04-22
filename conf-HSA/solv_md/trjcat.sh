#!/bin/bash
CWD=`pwd`

for i in `cat list.txt` ; do
	cd $i
	gmx trjcat -f prod_?.xtc -o solv_all.xtc -cat
	echo -e "MOL /n" | gmx trjconv -f solv_all.xtc -s prod_0.tpr -dump 0 -o v1.pdb
	cd $CWD
done
