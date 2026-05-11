#!/bin/bash

#python ../../scripts/com/bypass_angle_type3.py lig_g.itp sys.top
#python ../../scripts/com/bypass_angle_type3.py lig_me_g.itp sys.top
for i in {0..3} ; do
	mkdir mmpbsa_${i}
	cd mmpbsa_${i}
	echo processing trj of trial $i ...
	echo -e  'Protein_Other \n Protein_Other'| gmx trjconv -s ../prod_${i}.tpr -f ../prod_${i}.xtc -pbc mol -ur compact -o pbc.xtc -center -n ../index.ndx
	echo -e 'Protein_Other \n Protein_Other' | gmx trjconv -s ../prod_${i}.tpr -pbc cluster -f pbc.xtc -o com_traj.xtc -n ../index.ndx
	rm pbc.xtc
	echo copying files into mmpbsa_$i ...
	cp ../prod_${i}.tpr com.tpr
	cp ../../../scripts/com/mmpbsa.?? .
	cd ..
done

