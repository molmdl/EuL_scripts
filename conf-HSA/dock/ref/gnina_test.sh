#!/bin/bash
#SBATCH -J gnina
#SBATCH --gres=gpu:1
#SBATCH -n 1
#SBATCH -c 8
#SBATCH -p workq

echo redock
#gnina -r rec.pdb -l dzp.mol2 --autobox_ligand ref.pdb --autobox_add 12 --out_flex rec-dzp_rec.pdb --full_flex_output --flexdist_ligand ref.pdb --flexdist 8 --exhaustiveness 64 --num_modes 32 --addH off --stripH off --cpu 8 -o rec-dzp.pdb | tee -a rec-dzp.log 2> /dev/null
gnina -r rec.pdb -l dzp.mol2 --autobox_ligand ref.pdb --autobox_add 12 --exhaustiveness 64 --num_modes 32 --addH off --stripH off --cpu 8 -o rec-dzp.pdb | tee -a rec-dzp.log 2> /dev/null

echo ensemble
for i in {0..9} ; do
	echo rec conf $i
	#gnina -r ../hsa${i}.pdb_ali.pdb -l dzp.mol2 --autobox_ligand ref.pdb --out_flex rec-dzp_rec.pdb --full_flex_output --flexdist_ligand ref.pdb --flexdist 8 --autobox_add 12 --exhaustiveness 64 --num_modes 32 --addH off --stripH off --cpu 8 -o hsa${i}-dzp.pdb | tee -a hsa${i}-dzp.log 2>/dev/null
	gnina -r ../hsa${i}.pdb_ali.pdb -l dzp.mol2 --autobox_ligand ref.pdb --autobox_add 12 --exhaustiveness 64 --num_modes 32 --addH off --stripH off --cpu 8 -o hsa${i}-dzp.pdb | tee -a hsa${i}-dzp.log 2>/dev/null
done
