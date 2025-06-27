import csv
import pandas as pd
from rdkit import Chem
import warnings

def is_valid_smiles(smiles):
    """Checks if a SMILES string is valid."""
    return Chem.MolFromSmiles(smiles) is not None

def load_metatrans(source_file, target_file, metatrans_csv, set_value):
    with open(source_file, 'r') as src, open(target_file, 'r') as tgt:
        source_smiles = src.readlines()
        target_smiles = tgt.readlines()

    if len(source_smiles) != len(target_smiles):
        raise ValueError("The source and target files have a different number of SMILES.")

    with open(metatrans_csv, 'w', newline='') as csvfile:
        writer = csv.writer(csvfile)
        
        writer.writerow(['parent_smiles', 'child_smiles', 'set'])
  
        # Write valid SMILES from both files into the CSV, removing spaces
        for parent, child in zip(source_smiles, target_smiles):
            parent_clean = parent.strip().replace(' ', '')
            child_clean = child.strip().replace(' ', '')
            
            # Validate SMILES before writing
            if is_valid_smiles(parent_clean) and is_valid_smiles(child_clean):
                writer.writerow([
                    parent_clean,
                    child_clean,
                    set_value  # Add the constant 'set' value
                ])
            else:
                warnings.warn(f"Invalid SMILES detected: Parent SMILES - {parent_clean}, Child SMILES - {child_clean}")

def combine_csv(train_csv, valid_csv, final_csv):
    df_train = pd.read_csv(train_csv)
    df_valid = pd.read_csv(valid_csv)
    
    # Concatenate the dataframes
    combined_df = pd.concat([df_train, df_valid], ignore_index=True)
    
    combined_df.to_csv(final_csv, index=False)

# Metatrans training data
source_file_train = 'dataset/raw_data/metatrans_source_train.txt'
target_file_train = 'dataset/raw_data/metatrans_target_train.txt'
metatrans_csv_train = 'dataset/extracted_data/metatrans_train.csv'

# Metatrans validation data
source_file_valid = 'dataset/raw_data/metatrans_source_valid.txt'
target_file_valid = 'dataset/raw_data/metatrans_target_valid.txt'
metatrans_csv_valid = 'dataset/extracted_data/metatrans_valid.csv'

metatrans_csv_final = 'dataset/extracted_data/metatrans_smiles.csv'

load_metatrans(source_file_train, target_file_train, metatrans_csv_train, set_value = 'train')
load_metatrans(source_file_valid, target_file_valid, metatrans_csv_valid, set_value = 'val')

combine_csv(metatrans_csv_train, metatrans_csv_valid, metatrans_csv_final)