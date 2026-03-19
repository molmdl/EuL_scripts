#!/bin/bash
# Split into 4 runs for all energy groups

# Edit these placeholders
TPR_FILE="../fm_1500ns_?.tpr"
TRAJ_FILE="../v1.xtc"
PARTITION=workq
CPUS="16"
MEM="4G"

INDEX_FILE="../../index_res.ndx"

mkdir -p rerun
cd rerun

# Generate residues.list in working directory (trimmed)
grep 'Protein_' ${INDEX_FILE} | grep -v Other | awk '{print $2}' > residues.list

# Read residues into array
mapfile -t RESIDUES < residues.list
NUM_RES=${#RESIDUES[@]}
echo "Found $NUM_RES residues"

# Split into 4 batches
BATCH_SIZE=$(( (NUM_RES + 3) / 4 ))
echo "Splitting into batches of ~$BATCH_SIZE residues"

for BATCH_NUM in 0 1 2 3; do
    START_IDX=$((BATCH_NUM * BATCH_SIZE))
    END_IDX=$((START_IDX + BATCH_SIZE))
    if [ $END_IDX -gt $NUM_RES ]; then
        END_IDX=$NUM_RES
    fi
    
    BATCH_RESIDUES=("${RESIDUES[@]:$START_IDX:$BATCH_SIZE}")
    BATCH_NAME="batch$((BATCH_NUM + 1))"
    
    echo "Creating MDP for $BATCH_NAME (${#BATCH_RESIDUES[@]} residues)..."
    
    # Build energygrps line
    ENERGYGRPS="Other"
    for RES in "${BATCH_RESIDUES[@]}"; do
        ENERGYGRPS="${ENERGYGRPS} ${RES}"
    done
    
    # Create MDP for this batch
    cat > "rerun_${BATCH_NAME}.mdp" << MDP_EOF
integrator          =  md
dt                  =  0.002
nsteps              =  750000000
nstcomm             =  1
nstxout             =  0
nstvout             =  0
nstxtcout           =  0
energygrps = ENERGYGRPS_PLACEHOLDER

constraints = h-bonds
cutoff-scheme = Verlet
vdwtype = cutoff
vdw-modifier = force-switch
rlist = 1.2
rvdw = 1.2
rvdw-switch = 1.0
coulombtype = PME
rcoulomb = 1.2
DispCorr = no

Tcoupl              =  V-rescale
tc-grps             =  Protein non-Protein
tau_t               =  2.0   2.0
ref_t               =  300      300
Pcoupl              =  Parrinello-Rahman
tau_p               =  2.0
compressibility     =  4.5e-5
ref_p               =  1.0
gen_vel             =  no
gen_temp            =  300.0
gen_seed            =  -1

disre = simple
disre-fc = 500
disre-tau = 0
disre-mixed = no
nstdisreout = 500000

orire = yes
orire-fc = 1000
orire-tau = 0
orire-fitgrp = backbone
nstorireout = 500000
MDP_EOF
    
    sed -i "s|ENERGYGRPS_PLACEHOLDER|${ENERGYGRPS}|g" "rerun_${BATCH_NAME}.mdp"
done

echo "Submitting 4 SLURM jobs..."

# Submit jobs for each batch
for BATCH_NUM in 0 1 2 3; do
    BATCH_NAME="batch$((BATCH_NUM + 1))"
    
    # Create SLURM script for this batch
    cat > "rerun_${BATCH_NAME}.sh" << SLURM_EOF
#!/bin/bash
#SBATCH --job-name=rerun_${BATCH_NAME}
#SBATCH --partition=PARTITION_VAR
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=CPUS_VAR
#SBATCH --mem=MEM_VAR
#SBATCH --time=24:00:00

TPR_FILE="TPR_VAR"
TRAJ_FILE="TRAJ_VAR"
BATCH_NAME="BATCH_VAR"
BATCH_NUM="BATCH_NUM_VAR"
BATCH_SIZE="BATCH_SIZE_VAR"

# Grompp
echo "Running grompp for \${BATCH_NAME}..."
gmx grompp -f rerun_\${BATCH_NAME}.mdp -c \${TPR_FILE} -n ../../index_res.ndx -o rerun_\${BATCH_NAME}.tpr -maxwarn 2 -p ../sys.top

echo "Converting TPR..."
echo 22 | gmx convert-tpr -s rerun_\${BATCH_NAME}.tpr -n ../../index_res.ndx -o rerun_\${BATCH_NAME}.tpr

# Rerun
echo "Running mdrun -rerun..."
gmx mdrun -s rerun_\${BATCH_NAME}.tpr -rerun \${TRAJ_FILE} -deffnm rerun_\${BATCH_NAME} -ntmpi 1 -ntomp CPUS_VAR -bonded cpu -nb cpu -update cpu

# Get all available energy terms (lists terms with numbers)
echo "Getting available energy terms..."
echo 0 | gmx energy -f rerun_\${BATCH_NAME}.edr > energy_terms_\${BATCH_NAME}.txt 2>&1

# Debug: show what we're searching
echo "Energy terms file content (first 50 lines):"
head -50 energy_terms_\${BATCH_NAME}.txt

# Extract energies for each residue pair
echo "Extracting energy terms..."

# Get residues for this batch
START_LINE=\$((BATCH_NUM * BATCH_SIZE + 1))
END_LINE=\$(((BATCH_NUM + 1) * BATCH_SIZE))
mapfile -t BATCH_RESIDUES < <(sed -n "\${START_LINE},\${END_LINE}p" residues.list)

for RES in "\${BATCH_RESIDUES[@]}"; do
    RES=\$(echo "\$RES" | tr -d ' ')  # Remove any spaces
    echo "Processing \$RES..."
    
    # Get term numbers for this specific residue pair
    # Format in energy_terms file: "  23  LJ-14:Other-Other   24  Coul-SR:Other-Protein_MET_0"
    # Need to extract number RIGHT BEFORE the term
    COUL_NUM=\$(grep -oE "[0-9]+[[:space:]]+Coul-SR:Other-\${RES}" "energy_terms_\${BATCH_NAME}.txt" | grep -oE "^[0-9]+")
    LJ_NUM=\$(grep -oE "[0-9]+[[:space:]]+LJ-SR:Other-\${RES}" "energy_terms_\${BATCH_NAME}.txt" | grep -oE "^[0-9]+")
    
    echo "  Found term numbers - Coul: \$COUL_NUM, LJ: \$LJ_NUM"
    
    if [ -n "\$COUL_NUM" ] && [ -n "\$LJ_NUM" ]; then
        echo -e "\${COUL_NUM}\n\${LJ_NUM}\n0" | gmx energy -f rerun_\${BATCH_NAME}.edr -o "rerun_\${RES}.xvg"
    else
        echo "Warning: Could not find terms for \$RES"
    fi
done

echo "Cleaning..."
rm -f rerun_\${BATCH_NAME}.edr rerun_\${BATCH_NAME}.tpr rerun_\${BATCH_NAME}.mdp rerun_\${BATCH_NAME}.xtc rerun_\${BATCH_NAME}.log rerun_\${BATCH_NAME}.gro \#*

echo "Done with \${BATCH_NAME}!"
SLURM_EOF

    # Replace placeholders
    sed -i "s|PARTITION_VAR|${PARTITION}|g" "rerun_${BATCH_NAME}.sh"
    sed -i "s|CPUS_VAR|${CPUS}|g" "rerun_${BATCH_NAME}.sh"
    sed -i "s|MEM_VAR|${MEM}|g" "rerun_${BATCH_NAME}.sh"
    sed -i "s|TPR_VAR|${TPR_FILE}|g" "rerun_${BATCH_NAME}.sh"
    sed -i "s|TRAJ_VAR|${TRAJ_FILE}|g" "rerun_${BATCH_NAME}.sh"
    sed -i "s|BATCH_VAR|${BATCH_NAME}|g" "rerun_${BATCH_NAME}.sh"
    sed -i "s|BATCH_NUM_VAR|${BATCH_NUM}|g" "rerun_${BATCH_NAME}.sh"
    sed -i "s|BATCH_SIZE_VAR|${BATCH_SIZE}|g" "rerun_${BATCH_NAME}.sh"
    
    sbatch "rerun_${BATCH_NAME}.sh"
    
    echo "Submitted ${BATCH_NAME}"
done

echo "All 4 jobs submitted!"
