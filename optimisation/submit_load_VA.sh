#!/bin/bash
#SBATCH --time=24:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=12
#SBATCH --mem=12G
#SBATCH -p long,vlong
#SBATCH --output=slurm-output-%j.log
#SBATCH --error=slurm-error-%j.log
 
# activate the veloxchem environment
source ~/.bashrc
conda activate aizynthmodels
 
# number of threads should match the SLURM specification
export OMP_NUM_THREADS=12
 
# start the calculation
python optimisation/load_VA_parts.py 
 
# end of script