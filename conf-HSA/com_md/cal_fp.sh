#!/bin/bash

source ~/scripts/gmxMMPBSA.env 

CWD=`pwd`

for i in `cat list.txt` ; do
	cd $i
	mkdir fp
	cd fp
	gmx trjcat -f ../mmpbsa_?/com_traj.xtc -o v1.xtc -cat
	cp ../mmpbsa_0/com.tpr .
	echo -e "Protein_Other \n" | gmx trjconv -f v1.xtc -s com.tpr -dump 0 -n ../index.ndx -o v1.pdb
	sbatch <<- EOF
#!/bin/bash
#SBATCH --job-name=fp_$i
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH -p cpu
#SBATCH --mem=32GB

source ~/scripts/miniconda_py3.env

export PYTHONPATH="../../../scripts/com:\$PYTHONPATH"

python ../../../scripts/com/fp.py \
  --topology v1.pdb \
  --trajectory v1.xtc \
  --receptor_sel "protein" \
  --ligand_sel "resname MOL" 
EOF
	cd $CWD
done
