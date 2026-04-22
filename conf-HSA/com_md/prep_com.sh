#!/bin/bash

#for i in phe_sssD_sap ; do 
#for i in phe_sssD_sap phe_sssL_sap phe_rrrL_sap ; do 
#for i in me_sssD_sap me_sssL_sap me_rrrL_sap ; do 
#for i in me_sssD_sap me_sssL_sap me_rrrL_sap phe_rrrD_sap me_rrrD_sap ; do 
#for i in dzp ibp ; do
#for i in phe_sssD_tsap phe_sssL_tsap phe_rrrD_tsap phe_rrrL_tsap  ; do
#for i in phe_rrrL_tsap phe_rrrD_tsap ; do # phe_sssL_tsap ; do # phe_sssD_tsap ; do 
for i in me_sssD_sap me_sssL_sap me_rrrL_sap me_rrrD_sap ; do 
	cp ../scripts/com/*mdp $i 
	cd $i
	sbatch ../../scripts/com/prep.sh
	cd ..
done
