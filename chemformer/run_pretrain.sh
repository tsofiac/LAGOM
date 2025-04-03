    sbatch chemformer/pretrain.sh \
    chemformer/pretrain_bart.yaml \
    dataset/finetune/mmp_100000rows_finetune.csv \
    step=1000000_mod_no_deepspeed.ckpt \
    500 \
    0 \
    0 \
    chemformer/bart_vocab_downstream.json \
    results/pretraining/base \
    results/pretraining/base/config.yaml

# sbatch fine_tune.sh \
#     fine_tune_bart.yaml \
#     data.csv \
#     model_weights.ckpt \
#     64 \
#     bart_vocab.json \
#     results/training/ \
#     results/training/config.yaml
