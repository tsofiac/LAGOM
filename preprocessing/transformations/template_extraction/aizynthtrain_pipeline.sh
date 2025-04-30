#!/bin/bash
#SBATCH --time=2-00:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=50g
#SBATCH -p medium-gpu
#SBATCH --output="../logs/template-pipeline-%j.log"
#SBATCH --gres=gpu:1
#SBATCH --constraint='[ampere|volta]'

input_reactions=$1
output_reactions_mapped=$2
output_reaction_classes=$3
output_template_lib=$4

output_dir=$(dirname "$3")
echo ${PWD}
echo "Output dir: "${output_dir}

echo "input_reactions: "${input_reactions}
echo "output_reactions_mapped: "${output_reactions_mapped}
echo "output_reaction_classes: "${output_reaction_classes}
echo "output_template_lib: "${output_template_lib}

source ~/.bashrc 
conda activate /projects/mai/se_mai/users/klkf872_annie/miniconda3/envs/rxnmapper/

echo "Atom-mapping"
python scripts/mapping.py --input ${input_reactions} --output ${output_reactions_mapped}

conda activate aizynthmodels

export PATH=$PATH:/opt/scp/apps/gen-2019a/software/hazelnut/3.7.2-GCCcore-8.2.0-Python-3.7.2/bin/
export OE_LICENSE=/opt/scp/apps/system/software/oelicense/1.0/oe_license.seq1
export OE_DIR=/opt/scp/apps/gen-2019a/software/oetoolkits/2020.2.0-GCCcore-8.2.0-Python-3.7.2/openeye/toolkits
export LD_LIBRARY_PATH=/projects/mai/se_mai/casp/prod_models/repos/conda-env/lib/

echo "Reaction class prediction"
python -m rxnutils.pipeline.runner --data ${output_reactions_mapped} --pipeline configs/nm_pipeline.yml --output ${output_reaction_classes} --max-workers 1

echo "Transform format of csv-file"
python scripts/transform_data.py ${output_reaction_classes}

echo "Extract templates"
unset JUPYTER_PATH
python -m aizynthmodels.template_based.pipelines.template_pipeline run --config configs/template_pipeline_config.yaml --max-workers 32

mv chemformer_template_library.csv ${output_template_lib}