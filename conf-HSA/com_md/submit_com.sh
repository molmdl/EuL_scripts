#!/bin/bash

#for l in phe_sssD_sap phe_sssL_sap phe_rrrL_sap ; do 
#for l in me_sssD_sap me_sssL_sap me_rrrL_sap ; do 
#for l in me_sssD_sap me_sssL_sap me_rrrL_sap phe_rrrD_sap me_rrrD_sap ; do 
#for l in dzp ibp ; do
#for l in phe_sssD_tsap phe_sssL_tsap phe_rrrD_tsap phe_rrrL_tsap  ; do
#for l in phe_rrrL_tsap phe_rrrD_tsap ; do # phe_sssL_tsap ; do # phe_sssD_tsap ; do 
#for l in me_sssD_sap me_sssL_sap me_rrrL_sap me_rrrD_sap ; do 
for l in me_rrrD_sap ; do 
	cd $l
	#cp ../../scripts/com/*mdp .
	bash ../../scripts/com/pr_prod.sh $l
	cd ..
done
