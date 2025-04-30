sbatch chemformer/pretrain.sh \
    chemformer/pretrain_bart.yaml \
    dataset/pretrain/mmp_new_split_finetune.csv \
    null \
    128 \
    0.5 \
    0.1 \
    chemformer/bart_vocab_downstream.json \
    results/pretraining/comb \
    results/pretraining/comb/config.yaml

# sbatch fine_tune.sh \
#     fine_tune_bart.yaml \
#     data.csv \
#     model_weights.ckpt \
#     64 \
#     bart_vocab.json \
#     results/training/ \
#     results/training/config.yaml
