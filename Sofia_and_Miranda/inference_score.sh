#!/bin/bash
#SBATCH --array=0-19
#SBATCH --job-name=score-model
#SBATCH --output="logs/inference-score-%j.log"
#SBATCH --nodes 1
##SBATCH --cpus-per-task=4
#SBATCH --mem=64gb
##days-hours:minutes:seconds
#SBATCH --time=2-00:00:00
#SBATCH -p medium-gpu
#SBATCH --gres=gpu:a100:1
 
source ~/.bashrc
conda activate aizynthmodels
 
export HYDRA_FULL_ERROR=1
 
# Set i_chunk to 0 when running interactively
i_chunk=0 #$SLURM_ARRAY_TASK_ID
n_chunks=20
 
config_template=$1
data=$2
model_path=$3
vocabulary_path=$4
batch_size=$5
n_predictions=$6
task=$7
output_score_data=$8
output_predictions=$9
config=${10}${i_chunk}.yaml
 
echo config_template=${config_template}
echo config=${config}
 
cp $config_template ${config}
 
sed -i "s+DATA_PATH+${data}+g" ${config}
sed -i "s+MODEL_PATH+${model_path}+g" ${config}
sed -i "s+BATCH_SIZE+${batch_size}+g" ${config}
sed -i "s+VOCABULARY_PATH+${vocabulary_path}+g" ${config}
sed -i "s+NUM_PREDICTIONS+${n_predictions}+g" ${config}
sed -i "s+TASK+${task}+g" ${config}
sed -i "s+N_CHUNKS+${n_chunks}+g" ${config}
sed -i "s+I_CHUNK+${i_chunk}+g" ${config}
sed -i "s+OUTPUT_SCORE_DATA+${output_score_data}${i_chunk}.csv+g" ${config}
sed -i "s+OUTPUT_PREDICTIONS+${output_predictions}${i_chunk}.json+g" ${config}
 
python -m aizynthmodels.chemformer.tools.inference_score +config=${config}