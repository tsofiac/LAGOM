# sbatch inference_score.sh \
./chemformer/inference_score.sh \
    chemformer/inference_score.yaml \
    dataset/finetune/combined_evaluation_finetune.csv \
    results/training/submitted/forward_prediction/version_5/checkpoints/epoch=67-step=816.ckpt \
    chemformer/bart_vocab_downstream.json \
    38 \
    20 \
    forward_prediction \
    results/evaluation/metrics_scores \
    results/evaluation/predictions \
    results/evaluation/inference_score

# 19 as batch size - ok!


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
