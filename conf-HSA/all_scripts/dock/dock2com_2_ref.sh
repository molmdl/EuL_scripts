#!/bin/bash

cd dzp
cp ../../ref/dzp/dzp.itp .
ln -s ../hsa?.pdb_ali.gro .
python ../../scripts/dock/dock2com_2.1.py -i dzp.itp -s hsa*-dzp.sdf -t dzp.mol2 -r ../../rec/\#topol.top.1\# --ff-path ../../amber19SB_OL21_OL3_lipid17.ff/forcefield.itp --water-itp ../../amber19SB_OL21_OL3_lipid17.ff/opc3.itp --ions-itp ../../amber19SB_OL21_OL3_lipid17.ff/ions.itp --lig-gro best.gro --com-gro com.gro --rec-itp rec.itp --sys-top sys.top --metric minimizedAffinity  
rm hsa?.pdb_ali.gro

mkdir ../../com_md/dzp/
cp sys.top rec.itp com.gro best.gro dzp.itp ../../rec/posre.itp ../../com_md/dzp/

cd ..

cd ibp
cp ../../ref/ibp/ibp.itp .
ln -s ../hsa?.pdb_ali.gro .
python ../../scripts/dock/dock2com_2.1.py -i ibp.itp -s hsa*-ibp.sdf -t ibp.mol2 -r ../../rec/\#topol.top.1\# --ff-path ../../amber19SB_OL21_OL3_lipid17.ff/forcefield.itp --water-itp ../../amber19SB_OL21_OL3_lipid17.ff/opc3.itp --ions-itp ../../amber19SB_OL21_OL3_lipid17.ff/ions.itp --lig-gro best.gro --com-gro com.gro --rec-itp rec.itp --sys-top sys.top --metric minimizedAffinity  
rm hsa?.pdb_ali.gro

mkdir ../../com_md/ibp/
cp sys.top rec.itp com.gro best.gro ibp.itp ../../rec/posre.itp  ../../com_md/ibp/

cd ..


