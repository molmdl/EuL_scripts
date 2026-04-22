#!/bin/bash

for l in phe_sssD_tsap phe_sssL_tsap phe_rrrD_tsap phe_rrrL_tsap ; do
	echo lig $l
	echo
	cd ${l}
	sbatch << EOF
#!/bin/bash
#SBATCH -J gnina_$l
#SBATCH --gres=gpu:1
#SBATCH -n 1
#SBATCH -c 16
#SBATCH -p workq


for i in {0..9} ; do
	echo rec conf \$i
	echo
	gnina -r ../hsa\${i}.pdb_ali.pdb -l ${l}.mol2 --autobox_ligand ../ref.pdb --autobox_add 12 --exhaustiveness 100 --num_modes 100 --addH off --stripH off --cpu 16 -o hsa\${i}-${l}.sdf --min_rmsd_filter 3 --scoring ad4_scoring | tee -a hsa\${i}-${l}.log 2>/dev/null
	echo
done

EOF
	cd ..
done

