#!/bin/bash
#source ~/scripts/gmx2023.5_plumed_gpu_py3.env 

for i in {0..3} ; do
	echo 1 1 | gmx trjconv -f pr_${i}.xtc -s pr_${i}.tpr -pbc mol -center -o pbcmol.xtc -b 10000 -ur compact
	echo 4 1 | gmx trjconv -f pbcmol.xtc -s pr_${i}.tpr -fit rot+trans -o apo_50ns_${i}_pbc.xtc
	rm pbcmol.xtc 
	echo 1 | gmx trjconv -s pr_${i}.tpr -dump 0 -f apo_50ns_${i}_pbc.xtc -o ../apo_50ns_${i}_pbc_0.pdb
	# ana
	echo 4 4 | gmx rms -f apo_50ns_${i}_pbc.xtc -s pr_${i}.tpr -o apo_50ns_${i}_bb.rmsd
	echo 4 4 | gmx rmsf -s pr_${i}.tpr -f apo_50ns_${i}_pbc.xtc -oq apo50ns_${i}_rmsf_bfac.pdb -o apo50ns_${i}_bb_rmsf.xvg -res -dir apo50ns_${i}_bb_rmsf.log -oc apo50ns_${i}_bb_rmsf.xvg
	echo 4 2 | gmx rmsf -s pr_${i}.tpr -f apo_50ns_${i}_pbc.xtc -oq apo50ns_${i}_rmsf_noh_bfac.pdb -o apo50ns_${i}_noh_rmsf.xvg -res -dir apo50ns_${i}_noh_rmsf.log -oc apo50ns_${i}_noh_rmsf.xvg
done

gmx trjcat -f ./apo_50ns_[0123]_pbc.xtc -o ./apo_50ns_merged.xtc -cat
echo 4 4 | gmx rms -f ./apo_50ns_merged.xtc -s ./pr_0.tpr -o apo_50ns_merged_bb.rmsd
echo 4 4 | gmx rmsf -s ./pr_0.tpr -f ./apo_50ns_merged.xtc -oq apo50ns_merged_rmsf_bfac.pdb -o apo50ns_merged_bb_rmsf.xvg -res -dir apo50ns_merged_bb_rmsf.log -oc apo50ns_merged_bb_rmsf.xvg
echo 4 2 | gmx rmsf -s ./pr_0.tpr -f ./apo_50ns_merged.xtc -oq apo50ns_merged_rmsf_noh_bfac.pdb -o apo50ns_mreged_noh_rmsf.xvg -res -dir apo50ns_mergd_noh_rmsf.log -oc apo50ns_merged_noh_rmsf.xvg

echo 4 1 | gmx trjconv  -f apo_50ns_merged.xtc -s pr_0.tpr -fit rot+trans -o apo_50ns_merged_fit.xtc
rm ./apo_50ns_merged.xtc
