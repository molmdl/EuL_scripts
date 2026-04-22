#!/bin/bash
# Test 1: With mol2 template (original command style)
# Test 2: Without mol2 template (ITP-only approach)
# Test 3: With --rec-pdb to convert PDB to GRO

echo "=== Test 1: With mol2 template ==="
python dock2com_1.py -i lig_g.itp -s hsa*-phe.sdf -t phe.mol2 -r ../../rec/\#topol.top.1\# --ff-path ../../amber19SB_OL21_OL3_lipid17.ff/forcefield.itp --water-itp ../../amber19SB_OL21_OL3_lipid17.ff/opc3.itp --ions-itp ../../amber19SB_OL21_OL3_lipid17.ff/ions.itp --lig-gro best.gro --com-gro com.gro --rec-itp rec.itp --sys-top sys.top --metric CNNscore

echo ""
echo "=== Test 2: Without mol2 template (ITP-only) ==="
python dock2com_1.py -i lig_g.itp -s hsa*-phe.sdf -r ../../rec/\#topol.top.1\# --ff-path ../../amber19SB_OL21_OL3_lipid17.ff/forcefield.itp --water-itp ../../amber19SB_OL21_OL3_lipid17.ff/opc3.itp --ions-itp ../../amber19SB_OL21_OL3_lipid17.ff/ions.itp --lig-gro best2.gro --com-gro com2.gro --rec-itp rec2.itp --sys-top sys2.top --metric CNNscore

echo ""
echo "=== Test 3: With --rec-pdb to convert PDB to GRO ==="
python dock2com_1.py -i lig_g.itp -s hsa*-phe.sdf -t phe.mol2 -r ../../rec/\#topol.top.1\# --ff-path ../../amber19SB_OL21_OL3_lipid17.ff/forcefield.itp --water-itp ../../amber19SB_OL21_OL3_lipid17.ff/opc3.itp --ions-itp ../../amber19SB_OL21_OL3_lipid17.ff/ions.itp --lig-gro best3.gro --com-gro com3.gro --rec-itp rec3.itp --sys-top sys3.top --metric CNNscore --rec-pdb hsa2.pdb_ali.pdb