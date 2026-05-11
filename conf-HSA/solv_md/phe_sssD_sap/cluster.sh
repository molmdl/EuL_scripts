echo 1 1 | gmx cluster -method gromos -f ./v1.xtc -s ./pr_0.tpr -rmsmin 0.1 -cutoff 0.1 -sz -clid -cl clusters.xtc
echo 1 | gmx trjconv -s ./pr_0.tpr -f clusters.xtc -sep -pbc mol -o phe_sssD_sap_0.pdb 
