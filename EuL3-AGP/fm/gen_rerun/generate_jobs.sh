#!/bin/bash
# Generates individual SLURM job scripts for each residue

TPR_FILE="YOUR_TPR_FILE.tpr"
TRAJ_FILE="YOUR_TRAJ_FILE.xtc"
PARTITION="YOUR_PARTITION"
TIME="02:00:00"
CPUS="4"
MEM="4G"

mkdir -p job_scripts

while read -r residue; do
    jobname="rerun_${residue}"
    mdp_file="rerun_${residue}.mdp"
    tpr_out="rerun_${residue}.tpr"
    
    cat > "job_scripts/${jobname}.slurm" << JOBEOF
#!/bin/bash
#SBATCH --job-name=${jobname}
#SBATCH --output=${jobname}.out
#SBATCH --error=${jobname}.err
#SBATCH --partition=${PARTITION}
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=${CPUS}
#SBATCH --time=${TIME}
#SBATCH --mem=${MEM}

# Create MDP for this residue
sed "s/_RESIDUE_/${residue}/g" rerun.mdp > ${mdp_file}

# Grompp to create TPR
gmx grompp -f ${mdp_file} -c ${TPR_FILE} -n index_res.ndx -o ${tpr_out} -maxwarn 2

# Rerun calculation
gmx mdrun -s ${tpr_out} -rerun ${TRAJ_FILE} -deffnm ${jobname} -ntmpi 1 -ntomp ${CPUS}

# Cleanup temporary files
rm -f ${mdp_file} ${tpr_out} ${jobname}.trr ${jobname}.log ${jobname}.gro
JOBEOF
    
done < residues.list

echo "Generated $(ls job_scripts/*.slurm | wc -l) job scripts in job_scripts/"
echo "Edit TPR_FILE, TRAJ_FILE, PARTITION, TIME, CPUS, MEM variables in generate_jobs.sh if needed"
echo "Submit all jobs with: for job in job_scripts/*.slurm; do sbatch \$job; done"
