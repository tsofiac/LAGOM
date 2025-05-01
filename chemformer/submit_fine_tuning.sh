sbatch chemformer/fine_tune_submit.sh \
    chemformer/fine_tune_bart_submit.yaml \
    dataset/finetune/metatrans_finetune.csv \
    results/pretraining/ChemVA/forward_prediction/version_2/checkpoints/last.ckpt \ #CHECKA
    64 \
    chemformer/bart_vocab_downstream.json \
    results/training/meta/ \
    results/training/meta/config.yaml

# sbatch fine_tune.sh \
#     fine_tune_bart.yaml \
#     data.csv \
#     model_weights.ckpt \
#     64 \
#     bart_vocab.json \
#     results/training/ \
#     results/training/config.yaml
