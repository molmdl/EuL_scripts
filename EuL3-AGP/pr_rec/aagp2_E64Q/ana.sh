#!/bin/bash
source ~/data/scripts/gmx2023.5_plumed2.9.3_gpu_shared.env

for i in {0..3} ; do
	echo 1 1 | gmx trjconv -f ../apo_100ns_${i}.xtc -s ../apo_100ns_${i}.tpr -pbc mol -center -o pbcmol.xtc
	echo 4 1 | gmx trjconv -f pbcmol.xtc -s ../apo_100ns_${i}.tpr -fit rot+trans -o ../apo_100ns_${i}_pbc.xtc
	rm pbcmol.xtc 
	echo 1 | gmx trjconv -s ../apo_100ns_${i}.tpr -dump 0 -f ../apo_100ns_${i}_pbc.xtc -o ../apo_100ns_${i}_pbc_0.pdb
	# ana
	echo 4 4 | gmx rms -f ../apo_100ns_${i}_pbc.xtc -s ../apo_100ns_${i}.tpr -o apo_100ns_${i}_bb.rmsd
	echo 4 4 | gmx rmsf -s ../apo_100ns_${i}.tpr -f ../apo_100ns_${i}_pbc.xtc -oq apo100ns_${i}_rmsf_bfac.pdb -o apo100ns_${i}_bb_rmsf.xvg -res -dir apo100ns_${i}_bb_rmsf.log -oc apo100ns_${i}_bb_rmsf.xvg
	echo 4 2 | gmx rmsf -s ../apo_100ns_${i}.tpr -f ../apo_100ns_${i}_pbc.xtc -oq apo100ns_${i}_rmsf_noh_bfac.pdb -o apo100ns_${i}_noh_rmsf.xvg -res -dir apo100ns_${i}_noh_rmsf.log -oc apo100ns_${i}_noh_rmsf.xvg
done

ln -s ../apo_100ns_[0123]_pbc.xtc .
gmx trjcat -f ../apo_100ns_[0123]_pbc.xtc -o ../apo_100ns_merged.xtc -cat
echo 4 1 | gmx trjconv -f ../apo_100ns_merged.xtc -s ../apo_100ns_0.tpr -fit rot+trans -o ../apo_100ns_merged_fit.xtc 
echo 4 4 | gmx rms -f ../apo_100ns_merged.xtc -s ../apo_100ns_0.tpr -o apo_100ns_merged_bb.rmsd
echo 4 4 | gmx rmsf -s ../apo_100ns_0.tpr -f ../apo_100ns_merged.xtc -oq apo100ns_merged_rmsf_bfac.pdb -o apo100ns_merged_bb_rmsf.xvg -res -dir apo100ns_merged_bb_rmsf.log -oc apo100ns_merged_bb_rmsf.xvg
echo 4 2 | gmx rmsf -s ../apo_100ns_0.tpr -f ../apo_100ns_merged.xtc -oq apo100ns_merged_rmsf_noh_bfac.pdb -o apo100ns_mreged_noh_rmsf.xvg -res -dir apo100ns_mergd_noh_rmsf.log -oc apo100ns_merged_noh_rmsf.xvg
