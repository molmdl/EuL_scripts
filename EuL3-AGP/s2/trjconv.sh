#source /opt/scripts/gromacs-2023.5_plumed2.9.3_py3.env 
source ~/data/scripts/gmx2023.5_plumed2.9.3_gpu_shared.env
gmx trjcat -f metad_?.xtc metad_?_1us.*xtc -o trjcat.xtc
echo 22 22 | gmx trjconv -s metad_?.tpr -f trjcat.xtc  -pbc cluster -n ../index.ndx -o pbc.xtc
echo 4 22 | gmx trjconv -s metad_?.tpr  -fit progressive -n ../index.ndx -f pbc.xtc -o v1.xtc
echo 22 | gmx trjconv -s metad_?.tpr  -dump 0  -n ../index.ndx -f v1.xtc -o v1.pdb
rm pbc.xtc trjcat.xtc
#gmx trjcat -f fm_1?00ns_?.*xtc -o fm_1500ns_${1}.xtc
#gmx trjcat -f fm_1?00ns_?.*xtc -o trjcat.xtc
#echo 22 22 | gmx trjconv -s fm_1500ns_?.tpr -f trjcat.xtc  -pbc cluster -n ../index.ndx -o pbc.xtc
#echo 4 22 | gmx trjconv -s fm_1500ns_?.tpr  -fit progressive -n ../index.ndx -f pbc.xtc -o v1.xtc
#echo 22 | gmx trjconv -s fm_1500ns_?.tpr  -dump 0  -n ../index.ndx -f v1.xtc -o v1.pdb
#rm pbc.xtc trjcat.xtc
