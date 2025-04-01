    sbatch chemformer/fine_tune.sh \
    chemformer/pretrain_bart.yaml \
    dataset/finetune/mmp_all_finetune.csv\
    null \
    128 \
    chemformer/bart_vocab_downstream.json \
    results/pretraining \
    results/pretraining/config.yaml

# sbatch fine_tune.sh \
#     fine_tune_bart.yaml \
#     data.csv \
#     model_weights.ckpt \
#     64 \
#     bart_vocab.json \
#     results/training/ \
#     results/training/config.yaml
