echo -e  'Protein_Other \n Protein_Other'| gmx trjconv -s ../${1}.tpr -f ../${1}.xtc -pbc mol -ur compact -o pbc.xtc -center -n index.ndx
echo -e 'Protein_Other \n Protein_Other' | gmx trjconv -s ../${1}.tpr -pbc cluster -f pbc.xtc -o com_traj.xtc -n index.ndx
rm pbc.xtc

