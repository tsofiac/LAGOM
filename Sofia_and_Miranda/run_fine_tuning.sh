# sbatch fine_tune.sh \
./Sofia_and_Miranda/fine_tune.sh \
    Sofia_and_Miranda/fine_tune_bart.yaml \
    half_half_gloryx_reference_dataset.csv \
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
