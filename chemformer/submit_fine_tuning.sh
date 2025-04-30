sbatch chemformer/fine_tune_submit.sh \
    chemformer/fine_tune_bart_submit.yaml \
    dataset/alohomora/annotations_finetune.csv \
    pre-trained_models/chemVA_comb.ckpt \
    64 \
    chemformer/bart_vocab_downstream.json \
    results/training/submitted/ \
    results/training/submitted/config.yaml

# sbatch fine_tune.sh \
#     fine_tune_bart.yaml \
#     data.csv \
#     model_weights.ckpt \
#     64 \
#     bart_vocab.json \
#     results/training/ \
#     results/training/config.yaml
