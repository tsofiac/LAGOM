# sbatch chemformer/fine_tune_submit.sh \
./chemformer/fine_tune_submit.sh \
    chemformer/fine_tune_bart_submit.yaml \
    dataset/alohomora/baseline_finetune.csv \
    pre-trained_models/step=1000000_mod_no_deepspeed.ckpt \
    64 \
    chemformer/bart_vocab_downstream.json \
    results/training/alohomora/ \
    results/training/alohomora/config.yaml

# sbatch fine_tune.sh \
#     fine_tune_bart.yaml \
#     data.csv \
#     model_weights.ckpt \
#     64 \
#     bart_vocab.json \
#     results/training/ \
#     results/training/config.yaml
