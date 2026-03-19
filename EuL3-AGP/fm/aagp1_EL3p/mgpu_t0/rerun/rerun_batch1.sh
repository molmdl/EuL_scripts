#!/bin/bash
#SBATCH --job-name=rerun_batch1
#SBATCH --partition=workq
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=4G
#SBATCH --time=24:00:00

TPR_FILE="../fm_1500ns_?.tpr"
TRAJ_FILE="../v1.xtc"
BATCH_NAME="batch1"
BATCH_NUM="0"
BATCH_SIZE="44"

# Grompp
echo "Running grompp for ${BATCH_NAME}..."
gmx grompp -f rerun_${BATCH_NAME}.mdp -c ${TPR_FILE} -n ../../index_res.ndx -o rerun_${BATCH_NAME}.tpr -maxwarn 2 -p ../sys.top

echo "Converting TPR..."
echo 22 | gmx convert-tpr -s rerun_${BATCH_NAME}.tpr -n ../../index_res.ndx -o rerun_${BATCH_NAME}.tpr

# Rerun
echo "Running mdrun -rerun..."
gmx mdrun -s rerun_${BATCH_NAME}.tpr -rerun ${TRAJ_FILE} -deffnm rerun_${BATCH_NAME} -ntmpi 1 -ntomp 16 -bonded cpu -nb cpu -update cpu

# Get all available energy terms (lists terms with numbers)
echo "Getting available energy terms..."
echo 0 | gmx energy -f rerun_${BATCH_NAME}.edr > energy_terms_${BATCH_NAME}.txt 2>&1

# Debug: show what we're searching
echo "Energy terms file content (first 50 lines):"
head -50 energy_terms_${BATCH_NAME}.txt

# Extract energies for each residue pair
echo "Extracting energy terms..."

# Get residues for this batch
START_LINE=$((BATCH_NUM * BATCH_SIZE + 1))
END_LINE=$(((BATCH_NUM + 1) * BATCH_SIZE))
mapfile -t BATCH_RESIDUES < <(sed -n "${START_LINE},${END_LINE}p" residues.list)

for RES in "${BATCH_RESIDUES[@]}"; do
    RES=$(echo "$RES" | tr -d ' ')  # Remove any spaces
    echo "Processing $RES..."
    
    # Get term numbers for this specific residue pair
    # Format in energy_terms file: "  23  LJ-14:Other-Other   24  Coul-SR:Other-Protein_MET_0"
    # Need to extract number RIGHT BEFORE the term
    COUL_NUM=$(grep -oE "[0-9]+[[:space:]]+Coul-SR:Other-${RES}" "energy_terms_${BATCH_NAME}.txt" | grep -oE "^[0-9]+")
    LJ_NUM=$(grep -oE "[0-9]+[[:space:]]+LJ-SR:Other-${RES}" "energy_terms_${BATCH_NAME}.txt" | grep -oE "^[0-9]+")
    
    echo "  Found term numbers - Coul: $COUL_NUM, LJ: $LJ_NUM"
    
    if [ -n "$COUL_NUM" ] && [ -n "$LJ_NUM" ]; then
        echo -e "${COUL_NUM}\n${LJ_NUM}\n0" | gmx energy -f rerun_${BATCH_NAME}.edr -o "rerun_${RES}.xvg"
    else
        echo "Warning: Could not find terms for $RES"
    fi
done

echo "Cleaning..."
rm -f rerun_${BATCH_NAME}.edr rerun_${BATCH_NAME}.tpr rerun_${BATCH_NAME}.mdp rerun_${BATCH_NAME}.xtc rerun_${BATCH_NAME}.log rerun_${BATCH_NAME}.gro \#*

echo "Done with ${BATCH_NAME}!"
