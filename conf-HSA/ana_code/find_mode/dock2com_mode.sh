#!/bin/bash
# dock2com_mode.sh
# ================
# Re-run dock2com_2.2.1.py for each compound using the docking pose
# most similar to the MD target binding mode (selected by contact distance
# profile similarity, see find_similar_poses.py).
#
# Skipped:
#   phe_sssD       - original reference (starting structure for MD)
#   phe_rrrD_tsap  - already the current best by contact similarity (hsa6 m5)
#
# For each compound this script:
#   1. Moves the existing com_md/<dir> to ${COM_MD}/bak_new (backup)
#   2. Creates a fresh com_md/<dir>
#   3. Runs dock2com_2.2.1.py with the specific SDF and model
#   4. Copies outputs (best.gro, com.gro, sys.top, rec.itp, posre.itp, ligand.itp)
#      into the new com_md/<dir>

COM_MD=../../com_md
REC_TOP=../../rec/#topol.top.1#
FF_PATH=../../amber19SB_OL21_OL3_lipid17.ff/forcefield.itp
WATER_ITP=../../amber19SB_OL21_OL3_lipid17.ff/opc3.itp
IONS_ITP=../../amber19SB_OL21_OL3_lipid17.ff/ions.itp
POSRE=../../rec/posre.itp
DOCK_SCRIPT=../../scripts/dock/dock2com_2.2.1.py

mkdir -p ${COM_MD}/bak_new
# ---------------------------------------------------------------------------
# phe_sssL -> hsa6-phe.sdf model 6
# ---------------------------------------------------------------------------
cd phe_sssL
ln -s ../hsa?.pdb_ali.gro .
mv best.gro _best.gro
python ${DOCK_SCRIPT} \
    -i lig_g.itp \
    -s hsa6-phe.sdf \
    -t phe.mol2 \
    -r ${REC_TOP} \
    --ff-path ${FF_PATH} \
    --water-itp ${WATER_ITP} \
    --ions-itp ${IONS_ITP} \
    --lig-gro best.gro \
    --com-gro com.gro \
    --rec-itp rec.itp \
    --sys-top sys.top \
    --model 6
rm hsa?.pdb_ali.gro

mv ${COM_MD}/phe_sssL_sap ${COM_MD}/bak_new
mkdir ${COM_MD}/phe_sssL_sap
cp sys.top rec.itp com.gro best.gro lig_g.itp ${POSRE} ${COM_MD}/phe_sssL_sap/

cd ..

# ---------------------------------------------------------------------------
# phe_rrrD -> hsa6-phe_rrrD.sdf model 27
# ---------------------------------------------------------------------------
cd phe_rrrD
ln -s ../hsa?.pdb_ali.gro .
mv best.gro _best.gro
python ${DOCK_SCRIPT} \
    -i phe_rrrD.itp \
    -s hsa6-phe_rrrD.sdf \
    -t phe_rrrD.mol2 \
    -r ${REC_TOP} \
    --ff-path ${FF_PATH} \
    --water-itp ${WATER_ITP} \
    --ions-itp ${IONS_ITP} \
    --lig-gro best.gro \
    --com-gro com.gro \
    --rec-itp rec.itp \
    --sys-top sys.top \
    --model 27
rm hsa?.pdb_ali.gro

mv ${COM_MD}/phe_rrrD_sap ${COM_MD}/bak_new
mkdir ${COM_MD}/phe_rrrD_sap
cp sys.top rec.itp com.gro best.gro phe_rrrD.itp ${POSRE} ${COM_MD}/phe_rrrD_sap/

cd ..

# ---------------------------------------------------------------------------
# phe_rrrL -> hsa6-phe_rrrL.sdf model 6
# ---------------------------------------------------------------------------
cd phe_rrrL
ln -s ../hsa?.pdb_ali.gro .
mv best.gro _best.gro
python ${DOCK_SCRIPT} \
    -i lig_g.itp \
    -s hsa6-phe_rrrL.sdf \
    -t phe_rrrL.mol2 \
    -r ${REC_TOP} \
    --ff-path ${FF_PATH} \
    --water-itp ${WATER_ITP} \
    --ions-itp ${IONS_ITP} \
    --lig-gro best.gro \
    --com-gro com.gro \
    --rec-itp rec.itp \
    --sys-top sys.top \
    --model 6
rm hsa?.pdb_ali.gro

mv ${COM_MD}/phe_rrrL_sap ${COM_MD}/bak_new
mkdir ${COM_MD}/phe_rrrL_sap
cp sys.top rec.itp com.gro best.gro lig_g.itp ${POSRE} ${COM_MD}/phe_rrrL_sap/

cd ..

# ---------------------------------------------------------------------------
# phe_sssD_tsap -> hsa6-phe_sssD_tsap.sdf model 12
# ---------------------------------------------------------------------------
cd phe_sssD_tsap
ln -s ../hsa?.pdb_ali.gro .
mv best.gro _best.gro
python ${DOCK_SCRIPT} \
    -i phe_sssD_tsap.itp \
    -s hsa6-phe_sssD_tsap.sdf \
    -t phe_sssD_tsap.mol2 \
    -r ${REC_TOP} \
    --ff-path ${FF_PATH} \
    --water-itp ${WATER_ITP} \
    --ions-itp ${IONS_ITP} \
    --lig-gro best.gro \
    --com-gro com.gro \
    --rec-itp rec.itp \
    --sys-top sys.top \
    --model 12
rm hsa?.pdb_ali.gro

mv ${COM_MD}/phe_sssD_tsap ${COM_MD}/bak_new
mkdir ${COM_MD}/phe_sssD_tsap
cp sys.top rec.itp com.gro best.gro phe_sssD_tsap.itp ${POSRE} ${COM_MD}/phe_sssD_tsap/

cd ..

# ---------------------------------------------------------------------------
# phe_sssL_tsap -> hsa6-phe_sssL_tsap.sdf model 4
# ---------------------------------------------------------------------------
cd phe_sssL_tsap
ln -s ../hsa?.pdb_ali.gro .
mv best.gro _best.gro
python ${DOCK_SCRIPT} \
    -i phe_sssL_tsap.itp \
    -s hsa6-phe_sssL_tsap.sdf \
    -t phe_sssL_tsap.mol2 \
    -r ${REC_TOP} \
    --ff-path ${FF_PATH} \
    --water-itp ${WATER_ITP} \
    --ions-itp ${IONS_ITP} \
    --lig-gro best.gro \
    --com-gro com.gro \
    --rec-itp rec.itp \
    --sys-top sys.top \
    --model 4
rm hsa?.pdb_ali.gro

mv ${COM_MD}/phe_sssL_tsap ${COM_MD}/bak_new
mkdir ${COM_MD}/phe_sssL_tsap
cp sys.top rec.itp com.gro best.gro phe_sssL_tsap.itp ${POSRE} ${COM_MD}/phe_sssL_tsap/

cd ..

# ---------------------------------------------------------------------------
# phe_rrrL_tsap -> hsa6-phe_rrrL_tsap.sdf model 2
# ---------------------------------------------------------------------------
cd phe_rrrL_tsap
ln -s ../hsa?.pdb_ali.gro .
mv best.gro _best.gro
python ${DOCK_SCRIPT} \
    -i phe_rrrL_tsap.itp \
    -s hsa6-phe_rrrL_tsap.sdf \
    -t phe_rrrL_tsap.mol2 \
    -r ${REC_TOP} \
    --ff-path ${FF_PATH} \
    --water-itp ${WATER_ITP} \
    --ions-itp ${IONS_ITP} \
    --lig-gro best.gro \
    --com-gro com.gro \
    --rec-itp rec.itp \
    --sys-top sys.top \
    --model 2
rm hsa?.pdb_ali.gro

mv ${COM_MD}/phe_rrrL_tsap ${COM_MD}/bak_new
mkdir ${COM_MD}/phe_rrrL_tsap
cp sys.top rec.itp com.gro best.gro phe_rrrL_tsap.itp ${POSRE} ${COM_MD}/phe_rrrL_tsap/

cd ..

# ---------------------------------------------------------------------------
# me_sssD -> hsa5-me_sssD.sdf model 27
# ---------------------------------------------------------------------------
cd me_sssD
ln -s ../hsa?.pdb_ali.gro .
mv best.gro _best.gro
python ${DOCK_SCRIPT} \
    -i lig_me_g.itp \
    -s hsa5-me_sssD.sdf \
    -t me_sssD.mol2 \
    -r ${REC_TOP} \
    --ff-path ${FF_PATH} \
    --water-itp ${WATER_ITP} \
    --ions-itp ${IONS_ITP} \
    --lig-gro best.gro \
    --com-gro com.gro \
    --rec-itp rec.itp \
    --sys-top sys.top \
    --model 27
rm hsa?.pdb_ali.gro

mv ${COM_MD}/me_sssD_sap ${COM_MD}/bak_new
mkdir ${COM_MD}/me_sssD_sap
cp sys.top rec.itp com.gro best.gro lig_me_g.itp ${POSRE} ${COM_MD}/me_sssD_sap/

cd ..

# ---------------------------------------------------------------------------
# me_sssL -> hsa5-me.sdf model 19
# ---------------------------------------------------------------------------
cd me_sssL
ln -s ../hsa?.pdb_ali.gro .
mv best.gro _best.gro
python ${DOCK_SCRIPT} \
    -i lig_me_g.itp \
    -s hsa5-me.sdf \
    -t me.mol2 \
    -r ${REC_TOP} \
    --ff-path ${FF_PATH} \
    --water-itp ${WATER_ITP} \
    --ions-itp ${IONS_ITP} \
    --lig-gro best.gro \
    --com-gro com.gro \
    --rec-itp rec.itp \
    --sys-top sys.top \
    --model 19
rm hsa?.pdb_ali.gro

mv ${COM_MD}/me_sssL_sap ${COM_MD}/bak_new
mkdir ${COM_MD}/me_sssL_sap
cp sys.top rec.itp com.gro best.gro lig_me_g.itp ${POSRE} ${COM_MD}/me_sssL_sap/

cd ..

# ---------------------------------------------------------------------------
# me_rrrD -> hsa4-me_rrrD.sdf model 20
# ---------------------------------------------------------------------------
cd me_rrrD
ln -s ../hsa?.pdb_ali.gro .
mv best.gro _best.gro
python ${DOCK_SCRIPT} \
    -i me_rrrD.itp \
    -s hsa4-me_rrrD.sdf \
    -t me_rrrD.mol2 \
    -r ${REC_TOP} \
    --ff-path ${FF_PATH} \
    --water-itp ${WATER_ITP} \
    --ions-itp ${IONS_ITP} \
    --lig-gro best.gro \
    --com-gro com.gro \
    --rec-itp rec.itp \
    --sys-top sys.top \
    --model 20
rm hsa?.pdb_ali.gro

mv ${COM_MD}/me_rrrD_sap ${COM_MD}/bak_new
mkdir ${COM_MD}/me_rrrD_sap
cp sys.top rec.itp com.gro best.gro me_rrrD.itp ${POSRE} ${COM_MD}/me_rrrD_sap/

cd ..

# ---------------------------------------------------------------------------
# me_rrrL -> hsa3-me_rrrL.sdf model 30
# ---------------------------------------------------------------------------
cd me_rrrL
ln -s ../hsa?.pdb_ali.gro .
mv best.gro _best.gro
python ${DOCK_SCRIPT} \
    -i lig_me_g.itp \
    -s hsa3-me_rrrL.sdf \
    -t me_rrrL.mol2 \
    -r ${REC_TOP} \
    --ff-path ${FF_PATH} \
    --water-itp ${WATER_ITP} \
    --ions-itp ${IONS_ITP} \
    --lig-gro best.gro \
    --com-gro com.gro \
    --rec-itp rec.itp \
    --sys-top sys.top \
    --model 30
rm hsa?.pdb_ali.gro

mv ${COM_MD}/me_rrrL_sap ${COM_MD}/bak_new
mkdir ${COM_MD}/me_rrrL_sap
cp sys.top rec.itp com.gro best.gro lig_me_g.itp ${POSRE} ${COM_MD}/me_rrrL_sap/

cd ..
