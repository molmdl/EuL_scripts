#!/bin/bash
for i in `cat list_res.txt` ; do
	gmx_mpi -quiet energy -f rerun_test.edr -o Other_${i}.xvg <<-EOF
		Coul-SR:Other-Protein_${i}
		LJ-SR:Other-Protein_${i}
		0
		EOF
done
