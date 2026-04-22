#!/bin/bash

python ../../ana_code/aln/mol2_reorder.py gro -i ibp.itp -s rec1-ibp.sdf -o ibp_redock.gro -t ibp.mol2 --metric CNNscore

python ../../ana_code/dock2com/dock2com.py --rec-gro rec1.gro --lig-gro ibp_redock.gro --rec-top topol.top --rec-itp rec.itp --lig-itp ibp.itp --sys-top sys.top --ff-path ../../amber19SB_OL21_OL3_lipid17.ff/forcefield.itp --water-itp ../../amber19SB_OL21_OL3_lipid17.ff/opc3.itp --ions-itp ../../amber19SB_OL21_OL3_lipid17.ff/ions.itp 
