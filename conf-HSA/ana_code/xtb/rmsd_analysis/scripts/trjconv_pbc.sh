#!/bin/bash
# trjconv_pbc.sh — PBC-correct solv_all.xtc and strip solvent+align for all 16 systems
#
# Following the trjconv.sh pipeline:
#   Step 1: gmx trjconv -pbc mol -ur compact  (make molecules whole)
#   Step 2: gmx trjconv -fit rot+trans          (remove overall rotation/translation)
#   Step 3: gmx trjconv -dump 0 (MOL only)       (extract MOL reference PDB)
#   Step 4: echo MOL | gmx trjconv               (strip solvent, keep MOL only)
#
# Output per system:
#   rmsd_analysis/solv_pbc/{sys}/solv_pbc.xtc   (PBC-fixed, fitted, MOL-only)
#   rmsd_analysis/solv_pbc/{sys}/ref.pdb        (frame-0 MOL reference)
#
# Usage:
#   bash rmsd_analysis/scripts/trjconv_pbc.sh
#   bash rmsd_analysis/scripts/trjconv_pbc.sh phe_sssD_sap

set -euo pipefail

SOLV_DIR="/share/home/nglokwan/dparker/dp_xinyi/solv_md"
OUT_DIR="/share/home/nglokwan/dparker/dp_xinyi/ana_code/xtb/rmsd_analysis/solv_pbc"

SYSTEMS=(
    me_rrrD_sap me_rrrD_tsap me_rrrL_sap me_rrrL_tsap
    me_sssD_sap me_sssD_tsap me_sssL_sap me_sssL_tsap
    phe_rrrD_sap phe_rrrD_tsap phe_rrrL_sap phe_rrrL_tsap
    phe_sssD_sap phe_sssD_tsap phe_sssL_sap phe_sssL_tsap
)

if [ $# -gt 0 ]; then
    SYSTEMS=("$@")
fi

for SYS in "${SYSTEMS[@]}"; do
    SYS_DIR="${SOLV_DIR}/${SYS}"
    OUT_SYS="${OUT_DIR}/${SYS}"
    mkdir -p "${OUT_SYS}"

    TPR="${SYS_DIR}/prod_0.tpr"
    XTC="${SYS_DIR}/solv_all.xtc"

    if [ ! -f "${TPR}" ] || [ ! -f "${XTC}" ]; then
        echo "SKIP ${SYS}: missing prod_0.tpr or solv_all.xtc"
        continue
    fi

    echo "=== Processing ${SYS} ==="

    # Step 1: Make molecules whole (-pbc mol) + compact (-ur compact)
    # Select group 0 (System) for output
    echo "  Step 1: PBC correction (-pbc mol -ur compact)..."
    echo 0 | gmx trjconv -s "${TPR}" -f "${XTC}" \
        -pbc mol -ur compact \
        -o "${OUT_SYS}/step1_pbc.xtc" 2>"${OUT_SYS}/step1.log"
    echo "  Done: $(wc -l < "${OUT_SYS}/step1.log") log lines"

    # Step 2: Remove overall rotation/translation (-fit rot+trans)
    # Fit on group 0 (System), output group 0 (System)
    echo "  Step 2: Fitting rotation+translation..."
    echo "0 0" | gmx trjconv -s "${TPR}" -f "${OUT_SYS}/step1_pbc.xtc" \
        -fit rot+trans \
        -o "${OUT_SYS}/step2_fit.xtc" 2>"${OUT_SYS}/step2.log"
    echo "  Done: $(wc -l < "${OUT_SYS}/step2.log") log lines"

    # Step 3: Extract frame-0 MOL-only reference PDB
    echo "  Step 3: Extract MOL reference PDB (frame 0)..."
    echo "MOL" | gmx trjconv -s "${TPR}" -f "${OUT_SYS}/step2_fit.xtc" \
        -dump 0 -o "${OUT_SYS}/ref.pdb" 2>"${OUT_SYS}/step3.log"
    echo "  Done: $(wc -l < "${OUT_SYS}/step3.log") log lines"

    # Step 4: Strip solvent — keep MOL only
    echo "  Step 4: Strip solvent (MOL only)..."
    echo "MOL" | gmx trjconv -s "${TPR}" -f "${OUT_SYS}/step2_fit.xtc" \
        -o "${OUT_SYS}/solv_pbc.xtc" 2>"${OUT_SYS}/step4.log"
    echo "  Done: $(wc -l < "${OUT_SYS}/step4.log") log lines"

    # Clean up intermediates
    rm -f "${OUT_SYS}/step1_pbc.xtc" "${OUT_SYS}/step2_fit.xtc"

    echo "  ${SYS} complete: ${OUT_SYS}/solv_pbc.xtc + ref.pdb"
    echo ""
done

echo "=== All done ==="
echo "Output structure per system:"
echo "  rmsd_analysis/solv_pbc/{sys}/solv_pbc.xtc  (PBC-fixed, fitted, MOL-only)"
echo "  rmsd_analysis/solv_pbc/{sys}/ref.pdb        (frame-0 MOL reference)"
