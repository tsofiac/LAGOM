# sbatch fine_tune.sh \
./chemformer/fine_tune.sh \
    chemformer/fine_tune_bart.yaml \
    dataset/alohomora/csp3_finetune.csv \
    chemf_mmp_aug_mask_last \
    64 \
    chemformer/bart_vocab_downstream.json \
    results/training \
    results/training/config.yaml

# sbatch fine_tune.sh \
#     fine_tune_bart.yaml \
#     data.csv \
#     model_weights.ckpt \
#     64 \
#     bart_vocab.json \
#     results/training/ \
#     results/training/config.yaml
