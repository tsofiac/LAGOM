from rdkit import Chem, DataStructs
from rdkit.Chem import Descriptors, AllChem 
from rdkit.Chem.SaltRemover import SaltRemover
import torch
import pandas as pd
import os
from pathlib import Path
from collections import Counter
from standardize_smiles import standardize_smiles_collection
from sklearn.model_selection import GroupShuffleSplit


dataset_metx = 'dataset/curated_data/metxbiodb_smiles.csv'
dataset_drugbank = 'dataset/curated_data/drugbank_smiles.csv'

# df_metx = pd.read_csv(dataset_metx)
df_drugbank = pd.read_csv(dataset_drugbank)

# combined_df = pd.concat([df_metx, df_drugbank], ignore_index=True)

new_reactions = []
# Iterate through each row to create new reactions based on the 'origin'
for index, row in df_drugbank.iterrows():
    # Check if the reactant is different from the origin
    if row['parent_id'] != row['origin']:
        # Create a new reaction with the origin as the reactant and keep the product same
        new_reaction = {
            'parent_id': row['origin'],
            'child_id': row['child_id'],
            'origin': row['origin'],

        }
        # Append the new reaction to the augmented DataFrame
        new_reactions.append(new_reaction)

augmented_reactions = pd.DataFrame(new_reactions)

# Merge original data with the new augmented reactions
# final_df = pd.concat([combined_df, augmented_reactions], ignore_index=True)

# Print the final DataFrame to view the augmented data
augmented_reactions.to_csv('augmentation_drugbank.csv', index=False)

# parent_id,parent_name,child_id,child_name,enzymes,origin,parent_smiles,child_smiles