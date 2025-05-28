# sbatch chemformer/inference_score_submit.sh \
./chemformer/inference_score_submit.sh \
    chemformer/inference_score_submit.yaml \
    dataset/finetune/gloryx_finetune.csv \
    fine-tuned_models/metatrans_epoch=36-step=259.ckpt \
    chemformer/bart_vocab_downstream.json \
    38 \
    20 \
    forward_prediction \
    results/evaluation/alohomora/metrics_scores \
    results/evaluation/alohomora/predictions \
    results/evaluation/alohomora/inference_score

# 19 as batch size - ok!


# sbatch inference_score.sh \
#     inference_score.yaml \
#     data.csv \
#     model_weights.ckpt \
#     chemformer/bart_vocab_downstream.json \
#     64 \
#     n_predictions \ 
#     forward_prediction \
#     output_score_data \
#     output_predictions 

# results/training/submitted/forward_prediction/version_4/checkpoints/epoch=110-step=1221.ckpt

# results/training/submitted/forward_prediction/version_2/checkpoints/epoch=81-step=574.ckpt