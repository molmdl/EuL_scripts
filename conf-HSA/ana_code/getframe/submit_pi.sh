#!/bin/bash
#SBATCH --job-name=pi_interactions
#SBATCH --partition=cpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=24:00:00
#SBATCH --output=pi_interactions_%j.out
#SBATCH --error=pi_interactions_%j.err

cd /share/home/nglokwan/dparker/dp_xinyi/ana_code/getframe

python detect_pi_interactions.py

echo "Job completed at $(date)"
