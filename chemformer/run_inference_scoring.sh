sbatch chemformer/inference_score.sh \
    chemformer/inference_score.yaml \
    dataset.csv \
    model_weights.ckpt \
    chemformer/bart_vocab_downstream.json \
    38 \
    20 \
    forward_prediction \
    results/evaluation/metrics_scores \
    results/evaluation/predictions \
    results/evaluation/inference_score

