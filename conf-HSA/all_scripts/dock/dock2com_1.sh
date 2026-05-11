#!/bin/bash

cd phe_sssL
ln -s ../hsa?.pdb_ali.gro .
cp ../../solv_md/phe_sssL_sap/lig_g.itp .
python ../../scripts/dock/dock2com_1.py -i lig_g.itp -s hsa*-phe.sdf -t phe.mol2 -r ../../rec/\#topol.top.1\# --ff-path ../../amber19SB_OL21_OL3_lipid17.ff/forcefield.itp --water-itp ../../amber19SB_OL21_OL3_lipid17.ff/opc3.itp --ions-itp ../../amber19SB_OL21_OL3_lipid17.ff/ions.itp --lig-gro best.gro --com-gro com.gro --rec-itp rec.itp --sys-top sys.top 
rm hsa?.pdb_ali.gro

mkdir ../../com_md/phe_sssL_sap/
cp sys.top rec.itp com.gro best.gro lig_g.itp ../../com_md/phe_sssL_sap/

cd ..

cd phe_sssD
ln -s ../hsa?.pdb_ali.gro .
cp ../../solv_md/phe_sssD_sap_frmInvert/lig_g.itp .
python ../../scripts/dock/dock2com_1.py -i lig_g.itp -s hsa*-phe_sssD.sdf -t phe_sssD.mol2 -r ../../rec/\#topol.top.1\# --ff-path ../../amber19SB_OL21_OL3_lipid17.ff/forcefield.itp --water-itp ../../amber19SB_OL21_OL3_lipid17.ff/opc3.itp --ions-itp ../../amber19SB_OL21_OL3_lipid17.ff/ions.itp --lig-gro best.gro --com-gro com.gro --rec-itp rec.itp --sys-top sys.top 
rm hsa?.pdb_ali.gro

mkdir ../../com_md/phe_sssD_sap/
cp sys.top rec.itp com.gro best.gro lig_g.itp ../../com_md/phe_sssD_sap/

cd phe_rrrL
ln -s ../hsa?.pdb_ali.gro .
cp ../../solv_md/phe_rrrL_sap/lig_g.itp .
python ../../scripts/dock/dock2com_1.py -i lig_g.itp -s hsa*-phe.sdf -t phe.mol2 -r ../../rec/\#topol.top.1\# --ff-path ../../amber19SB_OL21_OL3_lipid17.ff/forcefield.itp --water-itp ../../amber19SB_OL21_OL3_lipid17.ff/opc3.itp --ions-itp ../../amber19SB_OL21_OL3_lipid17.ff/ions.itp --lig-gro best.gro --com-gro com.gro --rec-itp rec.itp --sys-top sys.top 
rm hsa?.pdb_ali.gro

mkdir ../../com_md/phe_rrrL_sap/
cp sys.top rec.itp com.gro best.gro lig_g.itp ../../com_md/phe_rrrL_sap/

cd ..

#cd ../../com_md/phe_sssD_sap/

