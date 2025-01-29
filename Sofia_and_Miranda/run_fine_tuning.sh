# sbatch fine_tune.sh \
./Sofia_and_Miranda/fine_tune.sh \
    Sofia_and_Miranda/fine_tune_bart.yaml \
    Sofia_and_Miranda/testdata.csv \
    null \
    64 \
    tests/chemformer/data/simple_vocab.json \
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
