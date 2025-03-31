#!/bin/bash
#SBATCH --job-name=fine-tune
#SBATCH --output="logs/fine_tune_chemformer-%j.log"
#SBATCH --nodes 1
##SBATCH --cpus-per-task=5
#SBATCH --ntasks=1
#SBATCH --mem=128gb
##days-hours:minutes:seconds
#SBATCH --time=2-00:00:00
#SBATCH -p medium-gpu
#SBATCH --gres=gpu:a100:1

source ~/.bashrc
conda activate aizynthmodels

config_template=$1
data=$2
model_path=$3
batch_size=$4
vocabulary_path=$5
output_directory=$6
config=$7

# for rate in `seq 0.001 0.001 0.01` ; do
# output_directory=results/hyperparam/LR_${rate}

cp $config_template $config

sed -i "s+DATA_PATH+${data}+g" ${config}
sed -i "s+MODEL_PATH+${model_path}+g" ${config}
sed -i "s+BATCH_SIZE+${batch_size}+g" ${config}
sed -i "s+VOCABULARY_PATH+${vocabulary_path}+g" ${config}
sed -i "s+OUTPUT_DIRECTORY+${output_directory}+g" ${config}

python -m aizynthmodels.chemformer.tools.fine_tune +config=${config}
