#!/bin/bash
# submit_rerun_jobs.sh

grep '^\\[ Protein_' index_res.ndx | sed 's/^\\[\\(.*\\)\\]/\\1/' > residues.list
N=\\$(wc -l < residues.list)
echo "Found \\${N} protein residues."
echo "Edit rerun_array.slurm placeholders (partition, cpus, time, mem, gromacs version, TPR, TRAJ), then:"
echo "sbatch --array=0-\\$((N-1)) rerun_array.slurm"

cat > rerun_array.slurm << 'ARRAY_EOF'
#!/bin/bash
#SBATCH --job-name=rerun_res
#SBATCH --output=rerun_%A_%a.out
#SBATCH --error=rerun_%A_%a.err

# PLACEHOLDERS - EDIT THESE:
#SBATCH --partition=cpu
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=4

# Load Gromacs

# Files - EDIT THESE:
TPR="md.tpr"  # Your .tpr file
TRAJ="traj.xtc"  # Your trajectory .xtc or .trr

res=\\$(sed -n "\\${SLURM_ARRAY_TASK_ID}p" residues.list)
mdp="rerun_\\${res}.mdp"
edr="rerun_\\${res}.edr"

sed "s/_RESIDUE_/\\${res}/g" rerun.mdp > "\\${mdp}"
gmx rerun -s "\\${TPR}" -f "\\${TRAJ}" -n index_res.ndx -mdp "\\${mdp}" -o "\\${edr}"
rm "\\${mdp}"
ARRAY_EOF
