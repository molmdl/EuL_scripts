#!/bin/bash
#SBATCH -J gnina
#SBATCH --gres=gpu:1
#SBATCH -n 1
#SBATCH -c 16
#SBATCH -p workq

#for l in phe me ; do
#for l in phe ; do
for l in  phe_rrrD me_rrrD ; do
	echo lig $l
	echo
	cd ${l}
	#cd ${l}_sssD
	for i in {0..9} ; do
		echo rec conf $i
		echo
		#gnina -r ../hsa${i}.pdb_ali.pdb -l ${l}_sssD.mol2 --autobox_ligand ../ref.pdb --autobox_add 12 --exhaustiveness 64 --num_modes 32 --addH off --stripH off --cpu 8 -o hsa${i}-${l}.sdf | tee -a hsa${i}-${l}.log 2>/dev/null
		gnina -r ../hsa${i}.pdb_ali.pdb -l ${l}.mol2 --autobox_ligand ../ref.pdb --autobox_add 12 --exhaustiveness 100 --num_modes 100 --addH off --stripH off --cpu 16 -o hsa${i}-${l}.sdf --min_rmsd_filter 3 --scoring ad4_scoring | tee -a hsa${i}-${l}.log 2>/dev/null
		echo
	done
	cd ..
	echo
done

#--out_flex rec-dzp_rec.pdb --full_flex_output --flexdist_ligand ref.pdb --flexdist 8 


