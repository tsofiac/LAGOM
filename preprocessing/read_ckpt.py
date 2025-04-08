import torch

# Specify the filename of the checkpoint
ckpt_filename = 'step=1000000_mod_no_deepspeed.ckpt'

# Load the checkpoint
try:
    # checkpoint = torch.load(ckpt_filename, map_location=torch.device('cpu'), weights_only=True) # for fine-tuned models
    checkpoint = torch.load(ckpt_filename, weights_only=False) # for chemf. pretrained model
    print("Checkpoint loaded successfully.")
    
    print("Checkpoint keys:", checkpoint.keys())

    if 'hyper_parameters' in checkpoint:
        hyperparams = checkpoint['hyper_parameters']
        print("Hyper parameters:", hyperparams)


except FileNotFoundError:
    print(f"The file {ckpt_filename} was not found.")
except Exception as e:
    print(f"An error occurred while loading the checkpoint: {e}")