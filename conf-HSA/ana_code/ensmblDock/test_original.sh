#!/bin/bash
#echo original one-sdf cmd
#python mol2_reorder.py gro -i dzp.itp -s rec-dzp.sdf -o dzp_redock.gro -t dzp.mol2 --metric CNNscore
echo mol2_reorder.py
python mol2_reorder.py gro -i dzp.itp -s hsa*-dzp.sdf -o best.gro -t dzp.mol2 --metric CNNscore

echo dock2com
#cp \#topol.top.1\# topol_.top
python dock2com.py --rec-gro hsa7.pdb_ali.gro --lig-gro best.gro --rec-top \#topol.top.1\# --rec-itp rec.itp --lig-itp dzp.itp --sys-top sys.top --ff-path ../../amber19SB_OL21_OL3_lipid17.ff/forcefield.itp --water-itp ../../amber19SB_OL21_OL3_lipid17.ff/opc3.itp --ions-itp ../../amber19SB_OL21_OL3_lipid17.ff/ions.itp

