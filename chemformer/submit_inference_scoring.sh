sbatch chemformer/inference_score_submit.sh \
    chemformer/inference_score_submit.yaml \
    dataset/finetune/combined_evaluation_finetune.csv \
    results/training/forward_prediction/version_44/checkpoints/epoch=23-step=168.ckpt \
    chemformer/bart_vocab_downstream.json \
    38 \
    20 \
    forward_prediction \
    results/evaluation/base/metrics_scores \
    results/evaluation/base/predictions \
    results/evaluation/base/inference_score

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
