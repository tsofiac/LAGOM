# sbatch inference_score.sh \
./chemformer/inference_score.sh \
    inference_score.yaml \
    dataset/curated_data/gloryx_smiles_clean.csv \
    results/training/forward_prediction/version_0/checkpoints/epoch=10-step=88.ckpt \
    chemformer/bart_vocab_downstream.json \
    64 \
    1 \ 
    forward_prediction \
    results/evaluation/metrics_scores \
    results/evaluation/predictions \
    results/evaluation/inference_score


# sbatch inference_score.sh \
#     inference_score.yaml \
#     data.csv \
#     model_weights.ckpt \
#     chemformer/bart_vocab_downstream.json \
#     64 \
#     n_predictions \ 
#     forward_prediction \
#     output_score_data \
#     output_predictions 
