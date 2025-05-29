sbatch chemformer/fine_tune.sh \
    chemformer/fine_tune_bart.yaml \
    dataset.csv \
    model_weights.ckpt \
    64 \
    chemformer/bart_vocab_downstream.json \
    results/training \
    results/training/config.yaml
