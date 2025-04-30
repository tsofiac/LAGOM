sbatch chemformer/inference_score_submit.sh \
    chemformer/inference_score_submit.yaml \
    dataset/alohomora/logp_evaluation.csv \
    results/training/submitted/forward_prediction/version_8/checkpoints/epoch=68-step=483.ckpt \
    chemformer/bart_vocab_downstream.json \
    38 \
    20 \
    forward_prediction \
    results/evaluation/final/metrics_scores \
    results/evaluation/final/predictions \
    results/evaluation/final/inference_score

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

# results/training/submitted/forward_prediction/version_4/checkpoints/epoch=110-step=1221.ckpt

# results/training/submitted/forward_prediction/version_2/checkpoints/epoch=81-step=574.ckpt