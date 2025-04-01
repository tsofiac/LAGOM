import torch

# Specify the filename of the checkpoint
ckpt_filename = 'results/training/forward_prediction/version_0/checkpoints/epoch=9-step=70.ckpt'

# Load the checkpoint
try:
    checkpoint = torch.load(ckpt_filename, map_location=torch.device('cpu'), weights_only=True)
    print("Checkpoint loaded successfully.")
    
    print("Checkpoint keys:", checkpoint.keys())

    if 'hyper_parameters' in checkpoint:
        hyperparams = checkpoint['hyper_parameters']
        print("Hyper parameters:", hyperparams)


except FileNotFoundError:
    print(f"The file {ckpt_filename} was not found.")
except Exception as e:
    print(f"An error occurred while loading the checkpoint: {e}")