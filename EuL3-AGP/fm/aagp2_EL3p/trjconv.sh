#source /opt/scripts/gromacs-2023.5_plumed2.9.3_py3.env 
source ~/data/scripts/gmx2023.5_plumed2.9.3_gpu_shared.env
#echo 22 22 | gmx trjconv -s fm_1000ns_?.tpr -f fm_1000ns_?.xtc  -pbc cluster -n ../index.ndx -o pbc.xtc
#echo 4 22 | gmx trjconv -s fm_1000ns_?.tpr  -fit progressive -n ../index.ndx -f pbc.xtc -o v1.xtc
#echo 22 | gmx trjconv -s fm_1000ns_?.tpr  -dump 0  -n ../index.ndx -f v1.xtc -o v1.pdb
#rm pbc.xtc
gmx trjcat -f fm_1?00ns_?.*xtc -o fm_1500ns_${1}.xtc
echo 22 22 | gmx trjconv -s fm_1500ns_?.tpr -f fm_1500ns_?.xtc  -pbc cluster -n ../index.ndx -o pbc.xtc
echo 4 22 | gmx trjconv -s fm_1500ns_?.tpr  -fit progressive -n ../index.ndx -f pbc.xtc -o v1.xtc
echo 22 | gmx trjconv -s fm_1500ns_?.tpr  -dump 0  -n ../index.ndx -f v1.xtc -o v1.pdb
rm pbc.xtc
