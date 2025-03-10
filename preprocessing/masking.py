import pandas as pd
import logging
import random
from typing import Any, Dict, List, Optional, Tuple
import torch
from omegaconf import DictConfig


# ---------------SpanTokensMasker---------------

def apply_span_mask(smiles:str, mask_prob, span_lambda=2.0):
    tokens = list(smiles)
    curr_idx = 0
    masked = []
    token_mask = []

    mask_bools = [True, False]
    weights = [mask_prob, 1 - mask_prob]
    sampled_mask = random.choices(mask_bools, weights=weights, k=len(tokens))

    while curr_idx < len(tokens):
        #print(curr_idx)
        if sampled_mask[curr_idx]:
            # Sample length of mask from a Poisson distribution
            mask_len = torch.poisson(torch.tensor(span_lambda)).long().item()
            print(mask_len)
            masked.append("<MASK>")
            token_mask.append(True)
            if mask_len == 0:
                curr_idx += 1
            else:
                curr_idx += mask_len
        else:
            masked.append(tokens[curr_idx])
            token_mask.append(False)
            curr_idx += 1

        masked_smiles = detokenize(masked)

    return masked_smiles

def detokenize(tokens: List[str]) -> str:
        return ''.join(tokens)

# --------------ReplaceTokensMasker-----------------
# Different from the one in Chemformer

def apply_replace_mask(smiles:str, mask_prob):
    tokens = list(smiles)
    curr_idx = 0
    masked = []
    token_mask = []

    mask_bools = [True, False]
    weights = [mask_prob, 1 - mask_prob]
    sampled_mask = random.choices(mask_bools, weights=weights, k=len(tokens))

    while curr_idx < len(tokens):
        #print(curr_idx)
        if sampled_mask[curr_idx]:
            # Sample length of mask from a Poisson distribution
            masked.append("<MASK>")
            token_mask.append(True)
        else:
            masked.append(tokens[curr_idx])
            token_mask.append(False)
        curr_idx += 1    

        masked_smiles = detokenize(masked)

    return masked_smiles

# ----------------------------------------------
def mask_dataset(csv_file, mask_prob, masked_file, masker_type):

    df = pd.read_csv(csv_file)
    parent_smiles = df['parent_smiles']

    masked_data = []
    for index, row in df.iterrows():

        parent_smiles = row['parent_smiles']
        child_smiles = row['child_smiles']

        for _ in range(1):
            if masker_type == 'spanmask':
                smiles_masked = apply_span_mask(parent_smiles, mask_prob)
            elif masker_type == 'replacemask':
                smiles_masked = apply_replace_mask(parent_smiles, mask_prob)
            else:
                print("WARNING: ", masker_type, " is not a valid masker type")

            masked_row = {
                'parent_name': row['parent_name'],
                'child_name': row['child_name'],
                'parent_smiles': smiles_masked,
                'child_smiles': row['child_smiles'],
                'origin': row['origin'],
                'source': row['source'],
                'set': row['set']
            }
                
            masked_data.append(masked_row)

    masked_dataset = pd.DataFrame(masked_data)
    masked_dataset = pd.concat([df, masked_dataset])
    masked_dataset.to_csv(masked_file, index=False)


if __name__ == "__main__":
     
    mask_prob = 0.1
    masker_type = 'replacemask' #'replacemask'  'spanmask'

    csv_file = 'dataset/curated_data/combined_smiles.csv'
    masked_file = f'dataset/curated_data/{masker_type}_smiles_clean.csv'

    mask_dataset(csv_file, mask_prob, masked_file, masker_type)



   