# sbatch fine_tune.sh \
./Sofia_and_Miranda/fine_tune.sh \
    Sofia_and_Miranda/fine_tune_bart.yaml \
    Sofia_and_Miranda/unique_parents_metxbiodb_finetuning.csv \
    step=1000000_mod_no_deepspeed.ckpt \
    64 \
    Sofia_and_Miranda/bart_vocab_downstream.json \
    Sofia_and_Miranda/results/training/ \
    Sofia_and_Miranda/results/training/config.yaml

# sbatch fine_tune.sh \
#     fine_tune_bart.yaml \
#     data.csv \
#     model_weights.ckpt \
#     64 \
#     bart_vocab.json \
#     results/training/ \
#     results/training/config.yaml
