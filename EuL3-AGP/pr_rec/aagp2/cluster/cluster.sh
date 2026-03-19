echo 1 1 | gmx_mpi cluster -method gromos -f ../apo_100ns_merged_fit.xtc -s ../apo_100ns_0.tpr -rmsmin 0.2 -cutoff 0.4 -sz -clid -cl clusters.xtc
echo 1 | gmx_mpi trjconv -s ../apo_100ns_0.tpr -f clusters.xtc -sep -pbc mol -o aagp2.gro
