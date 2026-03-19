echo 26 26 | gmx_mpi trjconv -s fm_500ns_?.tpr -f fm_500ns_?_merged.xtc -pbc mol -o pbc_mol.xtc -n index.ndx #-center
echo 26 26 | gmx_mpi trjconv -s fm_500ns_?.tpr -f pbc_mol.xtc  -pbc cluster -o pbc_clust.xtc -n index.ndx
echo 1 26 | gmx_mpi trjconv -s fm_300ns_?.tpr -f pbc_clust.xtc  -fit progressive -o v1.xtc -n index.ndx
echo 26 26 | gmx_mpi trjconv -s fm_500ns_?.tpr -f v1.xtc -pbc mol -o v1.pdb -n index.ndx -dump 0 #-center
rm pbc_*.xtc
