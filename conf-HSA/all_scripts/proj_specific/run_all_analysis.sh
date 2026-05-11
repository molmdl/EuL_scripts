#!/bin/bash
# run_all_analysis.sh
#
# Submit metal coordination geometry analysis (v3) as individual SLURM jobs
# for each system and trajectory type.
#
# Each job analyzes one system (solv or com trajectory).
#
# v3 updates:
#   - Corrected TSAP twist angle (22.5° instead of 0°)
#   - Added twist angle classification
#   - Added chirality determination (Δ/Λ)
#   - Added SHAPE ChSM analysis (Part C)
#
# Usage
# -----
#   bash run_all_analysis.sh          # Submit all jobs
#   bash run_all_analysis.sh --dry    # Print commands without submitting
#

set -euo pipefail

DRY_RUN=false
if [ "${1:-}" = "--dry" ]; then
    DRY_RUN=true
    echo "DRY RUN - no jobs will be submitted"
    echo ""
fi

PYTHON="${PYTHON:-python}"
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ANA="${SCRIPT_DIR}/metal_geo_analysis_v3.py"

SYSTEMS=(
    me_rrrD_sap
    me_rrrL_sap
    me_sssD_sap
    me_sssL_sap
    me_rrrD_tsap
    me_rrrL_tsap
    me_sssD_tsap
    me_sssL_tsap
    phe_rrrD_sap
    phe_rrrD_tsap
    phe_rrrL_sap
    phe_rrrL_tsap
    phe_sssD_sap
    phe_sssD_tsap
    phe_sssL_sap
    phe_sssL_tsap
)

submit_job() {
    local sysname="$1"
    local trjtype="$2"
    local tpr="$3"
    local xtc="$4"
    local outdir="$5"
    
    local jobname="geo_${sysname}_${trjtype}"
    
    if [ "$DRY_RUN" = true ]; then
        echo "Would submit job: $jobname"
    else
        mkdir -p "$outdir"
        sbatch <<- EOF
#!/bin/bash
#SBATCH -J ${jobname}
#SBATCH -p workq
#SBATCH -c 1
#SBATCH -o ${outdir}/analysis_%j.out
#SBATCH -e ${outdir}/analysis_%j.err

set -euo pipefail

echo "========================================================"
echo "  System: ${sysname}"
echo "  Type:   ${trjtype}"
echo "  Job ID: \$SLURM_JOB_ID"
echo "  Date:   \$(date)"
echo "========================================================"

${PYTHON} ${ANA} \\
    --system "${sysname}" \\
    --tpr "${tpr}" \\
    --xtc "${xtc}" \\
    --skip-classification-chart \\
    --outdir "${outdir}"

echo ""
echo "Analysis complete: \$(date)"
echo "Output directory: ${outdir}"
EOF
        echo "Submitted: $jobname"
    fi
}

echo "========================================================"
echo "  Submitting SOLV analysis jobs (16 systems)"
echo "========================================================"

for SYS in "${SYSTEMS[@]}"; do
    TPR="solv_md/${SYS}/prod_0.tpr"
    XTC="solv_md/${SYS}/solv_all.xtc"
    OUTDIR="solv_analysis/${SYS}"
    
    if [ -f "$TPR" ] && [ -f "$XTC" ]; then
        submit_job "$SYS" "solv" "$TPR" "$XTC" "$OUTDIR"
    else
        echo "Skipping $SYS (solv): missing TPR or XTC"
    fi
done

echo ""
echo "========================================================"
echo "  Submitting COM analysis jobs (16 systems)"
echo "========================================================"

for SYS in "${SYSTEMS[@]}"; do
    TPR="com_md/${SYS}/fp/com.tpr"
    XTC="com_md/${SYS}/fp/v1.xtc"
    OUTDIR="com_analysis/${SYS}"
    
    if [ -f "$TPR" ] && [ -f "$XTC" ]; then
        submit_job "$SYS" "com" "$TPR" "$XTC" "$OUTDIR"
    else
        echo "Skipping $SYS (com): missing TPR or XTC"
    fi
done

echo ""
echo "========================================================"
if [ "$DRY_RUN" = true ]; then
    echo "  DRY RUN COMPLETE"
    echo "  Run without --dry to submit jobs"
else
    echo "  All jobs submitted"
    echo "  Check status with: squeue -u \$USER"
fi
echo "========================================================"
echo ""
echo "v3 New Features:"
echo "  - Part C: SHAPE ChSM analysis"
echo "  - Twist angle classification (primary geometry method)"
echo "  - Chirality determination (Δ/Λ)"
echo "  - Corrected TSAP definition (22.5° twist)"
