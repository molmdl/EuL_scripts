#!/bin/bash
#SBATCH -J gnina_ref_ens
#SBATCH --gres=gpu:1
#SBATCH -n 1
#SBATCH -c 16
#SBATCH -p workq

cd dzp_rec
gnina -r ../rec.pdb -l dzp.mol2 --autobox_ligand ../ref.pdb --autobox_add 12 --exhaustiveness 100 --num_modes 100 --addH off --stripH off --cpu 16 -o rec-dzp.sdf --min_rmsd_filter 3 --scoring ad4_scoring --log rec-dzp.log 2>/dev/null
cd ..

cd ibp_rec1
gnina -r ../rec1.pdb -l ibp.mol2 --autobox_ligand ../ref.pdb --autobox_add 12 --exhaustiveness 100 --num_modes 100 --addH off --stripH off --cpu 16 -o rec1-dzp.sdf --min_rmsd_filter 3 --scoring ad4_scoring --log rec1-ibp.log 2>/dev/null
cd ..
