gmx trjcat -f prod_?.xtc -cat -o trjcat.xtc
echo 0 | gmx trjconv -s prod_0.tpr -f trjcat.xtc -pbc mol -ur compact -o pbc.xtc
echo 0 | gmx trjconv -s prod_0.tpr -fit rot+trans -f pbc.xtc -o prod_4x100ns.xtc
rm trjcat.xtc pbc.xtc
