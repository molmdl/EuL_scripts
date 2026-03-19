echo 1 1 | gmx cluster -method gromos -f ../prod_4x100ns.xtc -s ../pr_0.tpr -rmsmin 0.1 -cutoff 0.2 -sz -clid -cl clusters.xtc
echo 1 | gmx trjconv -s ../prod_0.tpr -f clusters.xtc -sep -pbc mol -o clust.gro
