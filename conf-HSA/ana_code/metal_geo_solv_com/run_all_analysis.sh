#!/usr/bin/env bash
# run_all_analysis.sh
#
# Run metal coordination geometry analysis (SAP/TSAP) for all 12 systems,
# both solvated (solv_md) and receptor-complex (com_md) trajectories.
#
# Inputs per system
# -----------------
#   solv_md : <system>/prod_0.tpr  +  <system>/solv_all.xtc
#   com_md  : <system>/fp/com.tpr  +  <system>/fp/v1.xtc
#
# Outputs
# -------
#   <system>/solv_analysis/   – Part A + Part B outputs for solvated MD
#   <system>/com_analysis/    – Part A + Part B outputs for complex MD
#
# Usage
# -----
#   bash run_all_analysis.sh
#
# Run from the metal_geo_solv_com/ directory (where the script lives).

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR"

PYTHON="${PYTHON:-python}"
ANA="metal_geo_analysis_v2.py"

SYSTEMS=(
    me_rrrD_sap
    me_rrrL_sap
    me_sssD_sap
    me_sssL_sap
    phe_rrrD_sap
    phe_rrrD_tsap
    phe_rrrL_sap
    phe_rrrL_tsap
    phe_sssD_sap
    phe_sssD_tsap
    phe_sssL_sap
    phe_sssL_tsap
)

# ── Build argument lists ───────────────────────────────────────────────────────

SOLV_ARGS=()
for SYS in "${SYSTEMS[@]}"; do
    SOLV_ARGS+=(
        --system "$SYS"
        --tpr    "solv_md/${SYS}/prod_0.tpr"
        --xtc    "solv_md/${SYS}/solv_all.xtc"
        --outdir solv_analysis
    )
done

COM_ARGS=()
for SYS in "${SYSTEMS[@]}"; do
    COM_ARGS+=(
        --system "$SYS"
        --tpr    "com_md/${SYS}/fp/com.tpr"
        --xtc    "com_md/${SYS}/fp/v1.xtc"
        --outdir com_analysis
    )
done

# ── Run ───────────────────────────────────────────────────────────────────────

echo "========================================================"
echo "  Running SOLV analysis (12 systems)"
echo "========================================================"
"$PYTHON" "$ANA" "${SOLV_ARGS[@]}"

echo ""
echo "========================================================"
echo "  Running COM analysis (12 systems)"
echo "========================================================"
"$PYTHON" "$ANA" "${COM_ARGS[@]}"

echo ""
echo "All analyses complete."
echo "  Solvated results : <system>/solv_analysis/"
echo "  Complex  results : <system>/com_analysis/"
