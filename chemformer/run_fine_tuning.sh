# sbatch fine_tune.sh \
./chemformer/fine_tune.sh \
    chemformer/fine_tune_bart.yaml \
    dataset/dummy_data/unique_parents_metxbiodb_finetuning.csv \
    step=1000000_mod_no_deepspeed.ckpt \
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
