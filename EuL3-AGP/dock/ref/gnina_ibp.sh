#!/bin/bash
#SBATCH -J gnina
#SBATCH --gres=gpu:1
#SBATCH -n 1
#SBATCH -c 8
#SBATCH -p workq

echo redock
#gnina -r rec.pdb -l ibp.mol2 --autobox_ligand ref.pdb --autobox_add 12 --out_flex rec-ibp_rec.pdb --full_flex_output --flexdist_ligand ref.pdb --flexdist 8 --exhaustiveness 64 --num_modes 32 --addH off --stripH off --cpu 8 -o rec-ibp.pdb | tee -a rec-ibp.log 2> /dev/null
gnina -r rec.pdb -l ibp.mol2 --autobox_ligand ref.pdb --autobox_add 12 --exhaustiveness 64 --num_modes 32 --addH off --stripH off --cpu 8 -o rec-ibp.pdb | tee -a rec-ibp.log 2> /dev/null
gnina -r rec1.pdb -l ibp.mol2 --autobox_ligand ibp.pdb --autobox_add 12 --exhaustiveness 64 --num_modes 32 --stripH off --cpu 8 -o rec1-ibp.pdb | tee -a rec1-ibp.log 2> /dev/null

echo ensemble
for i in {0..9} ; do
	echo rec conf $i
	#gnina -r ../hsa${i}.pdb_ali.pdb -l ibp.mol2 --autobox_ligand ref.pdb --out_flex rec-ibp_rec.pdb --full_flex_output --flexdist_ligand ref.pdb --flexdist 8 --autobox_add 12 --exhaustiveness 64 --num_modes 32 --addH off --stripH off --cpu 8 -o hsa${i}-ibp.pdb | tee -a hsa${i}-ibp.log 2>/dev/null
	gnina -r ../hsa${i}.pdb_ali.pdb -l ibp.mol2 --autobox_ligand ref.pdb --autobox_add 12 --exhaustiveness 64 --num_modes 32 --addH off --stripH off --cpu 8 -o hsa${i}-ibp.pdb --pose_sort_order Energy | tee -a hsa${i}-ibp.log 2>/dev/null
done
