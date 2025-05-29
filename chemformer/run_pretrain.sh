sbatch chemformer/pretrain.sh \
    chemformer/pretrain_bart.yaml \
    dataset.csv \
    model_weights.ckpt \
    128 \
    0.5 \
    0.1 \
    chemformer/bart_vocab_downstream.json \
    results/pretraining \
    results/pretraining/config.yaml
