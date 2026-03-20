echo 1 0 | gmx_mpi trjconv -s prod_300ns.tpr -f prod_300ns.xtc -pbc mol -center -o pbcmol.xtc
echo 1 0 | gmx_mpi trjconv -s prod_300ns.tpr -pbc cluster -f pbcmol.xtc -o pbccluster.xtc 2> /dev/null
echo 1 0 | gmx_mpi trjconv -s prod_300ns.tpr -fit progressive -f pbccluster.xtc
rm pbccluster.xtc pbcmol.xtc
