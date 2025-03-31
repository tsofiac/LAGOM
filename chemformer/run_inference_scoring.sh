# sbatch inference_score.sh \
./chemformer/inference_score.sh \
    chemformer/inference_score.yaml \
    dataset/finetune/combined_evaluation_finetune.csv \
    results/training/forward_prediction/version_61/checkpoints/epoch=4-step=135.ckpt \
    chemformer/bart_vocab_downstream.json \
    64 \
    20 \
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
