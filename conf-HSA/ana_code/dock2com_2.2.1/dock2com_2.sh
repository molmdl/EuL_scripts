#!/bin/bash

i=sssD
#cd phe_${i}_tsap
#ln -s ../hsa?.pdb_ali.gro .
#cp ../../solv_md/phe_${i}_tsap/phe_${i}_tsap.itp .
python3 dock2com_2.2.1.py -i phe_${i}_tsap.itp -s hsa5-phe_${i}_tsap.sdf -t phe_${i}_tsap.mol2 -r ../../rec/\#topol.top.1\# --ff-path ../../amber19SB_OL21_OL3_lipid17.ff/forcefield.itp --water-itp ../../amber19SB_OL21_OL3_lipid17.ff/opc3.itp --ions-itp ../../amber19SB_OL21_OL3_lipid17.ff/ions.itp --lig-gro best.gro --com-gro com.gro --rec-itp rec.itp --sys-top sys.top --model 23
#rm hsa?.pdb_ali.gro
#mkdir ../../com_md/phe_${i}_tsap/
#cp sys.top rec.itp com.gro best.gro phe_${i}_tsap.itp ../../rec/posre.itp  ../../com_md/phe_${i}_tsap/
#cd ..
