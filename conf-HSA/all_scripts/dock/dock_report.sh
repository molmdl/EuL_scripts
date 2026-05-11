#!/bin/bash

i=$1
itp=$2

cd $i
ln -s ../hsa?.pdb_ali.gro .
#cp ../../solv_md/$/phe_${i}_tsap.itp .
python ../../scripts/dock/dock2com_2.2.1.py -i ${itp}.itp -s hsa?-${i}.sdf -t ${i}.mol2 --list-models --metric sasa > score_list.csv
rm hsa?.pdb_ali.gro
cd ..


