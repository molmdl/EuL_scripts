#!/bin/bash

#
#cd me_rrrD_sap 
#mkdir mmpbsa_0
#cd mmpbsa_0
#echo processing trj of trial 0 ...
#echo -e  'Protein_Other \n Protein_Other'| gmx trjconv -s ../prod_0.tpr -f ../prod_0.xtc -pbc mol -ur compact -o pbc.xtc -center -n ../index.ndx
#echo -e 'Protein_Other \n Protein_Other' | gmx trjconv -s ../prod_0.tpr -pbc cluster -f pbc.xtc -o com_traj.xtc -n ../index.ndx
#rm pbc.xtc
#echo copying files into mmpbsa_0 ...
#cp ../prod_0.tpr com.tpr
#cp ../../../scripts/com/mmpbsa.?? .
#sbatch mmpbsa.sh 
#cd .. 
#
#cd ..
#
##
#cd me_rrrL_sap 
#mkdir mmpbsa_1
#cd mmpbsa_1
#echo processing trj of trial 1 ...
#echo -e  'Protein_Other \n Protein_Other'| gmx trjconv -s ../prod_1.tpr -f ../prod_1.xtc -pbc mol -ur compact -o pbc.xtc -center -n ../index.ndx
#echo -e 'Protein_Other \n Protein_Other' | gmx trjconv -s ../prod_1.tpr -pbc cluster -f pbc.xtc -o com_traj.xtc -n ../index.ndx
#rm pbc.xtc
#echo copying files into mmpbsa_1 ...
#cp ../prod_1.tpr com.tpr
#cp ../../../scripts/com/mmpbsa.?? .
#sbatch mmpbsa.sh 
#cd .. 
#
#cd ..
#
##
#cd phe_rrrL_tsap 
#mkdir mmpbsa_0
#cd mmpbsa_0
#echo processing trj of trial 0 ...
#echo -e  'Protein_Other \n Protein_Other'| gmx trjconv -s ../prod_0.tpr -f ../prod_0.xtc -pbc mol -ur compact -o pbc.xtc -center -n ../index.ndx
#echo -e 'Protein_Other \n Protein_Other' | gmx trjconv -s ../prod_0.tpr -pbc cluster -f pbc.xtc -o com_traj.xtc -n ../index.ndx
#rm pbc.xtc
#echo copying files into mmpbsa_0 ...
#cp ../prod_0.tpr com.tpr
#cp ../../../scripts/com/mmpbsa.?? .
#sbatch mmpbsa.sh 
#cd .. 
#
#cd ..
#
#
#cd me_sssD_sap 
#mkdir mmpbsa_3
#cd mmpbsa_3
#echo processing trj of trial 3 ...
#echo -e  'Protein_Other \n Protein_Other'| gmx trjconv -s ../prod_3.tpr -f ../prod_3.xtc -pbc mol -ur compact -o pbc.xtc -center -n ../index.ndx
#echo -e 'Protein_Other \n Protein_Other' | gmx trjconv -s ../prod_3.tpr -pbc cluster -f pbc.xtc -o com_traj.xtc -n ../index.ndx
#rm pbc.xtc
#echo copying files into mmpbsa_3 ...
#cp ../prod_3.tpr com.tpr
#cp ../../../scripts/com/mmpbsa.?? .
#sbatch mmpbsa.sh 
#cd .. 
#
#cd ..
#
##
#cd phe_rrrD_sap 
#mkdir mmpbsa_3
#cd mmpbsa_3
#echo processing trj of trial 3 ...
#echo -e  'Protein_Other \n Protein_Other'| gmx trjconv -s ../prod_3.tpr -f ../prod_3.xtc -pbc mol -ur compact -o pbc.xtc -center -n ../index.ndx
#echo -e 'Protein_Other \n Protein_Other' | gmx trjconv -s ../prod_3.tpr -pbc cluster -f pbc.xtc -o com_traj.xtc -n ../index.ndx
#rm pbc.xtc
#echo copying files into mmpbsa_3 ...
#cp ../prod_3.tpr com.tpr
#cp ../../../scripts/com/mmpbsa.?? .
#sbatch mmpbsa.sh 
#cd .. 
#
#cd ..
#
#

cd phe_sssD_tsap 
mkdir mmpbsa_0
cd mmpbsa_0
echo processing trj of trial 0 ...
echo -e  'Protein_Other \n Protein_Other'| gmx trjconv -s ../prod_0.tpr -f ../prod_0.xtc -pbc mol -ur compact -o pbc.xtc -center -n ../index.ndx
echo -e 'Protein_Other \n Protein_Other' | gmx trjconv -s ../prod_0.tpr -pbc cluster -f pbc.xtc -o com_traj.xtc -n ../index.ndx
rm pbc.xtc
echo copying files into mmpbsa_3 ...
cp ../prod_0.tpr com.tpr
cp ../../../scripts/com/mmpbsa.?? .
sbatch mmpbsa.sh 
cd .. 

mkdir mmpbsa_1
cd mmpbsa_1
echo processing trj of trial 0 ...
echo -e  'Protein_Other \n Protein_Other'| gmx trjconv -s ../prod_1.tpr -f ../prod_1.xtc -pbc mol -ur compact -o pbc.xtc -center -n ../index.ndx
echo -e 'Protein_Other \n Protein_Other' | gmx trjconv -s ../prod_1.tpr -pbc cluster -f pbc.xtc -o com_traj.xtc -n ../index.ndx
rm pbc.xtc
echo copying files into mmpbsa_3 ...
cp ../prod_1.tpr com.tpr
cp ../../../scripts/com/mmpbsa.?? .
sbatch mmpbsa.sh 
cd .. 

cd ..

#
