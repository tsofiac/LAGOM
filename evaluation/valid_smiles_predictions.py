import pandas as pd
import ast
from rdkit import Chem

def save_10_valid_smiles(input_file, output_file):

    df = pd.read_csv(input_file)
    
    for i in range(len(df)):
        sampled_molecules_i = ast.literal_eval(df.at[i, 'sampled_molecules'])

        # Filter for valid SMILES
        valid_smiles = [smile for smile in sampled_molecules_i if Chem.MolFromSmiles(smile) is not None]

        # Check for less than 10 valid SMILES
        if len(valid_smiles) < 10:
            print(f"Warning: Row {i} has less than 10 valid SMILES.")

        # Update the 'sampled_molecules' column with the valid SMILES
        df.at[i, 'sampled_molecules'] = valid_smiles[:10]  # Include at most 10 valid SMILES  

    df.to_csv(output_file, index=False)
   


input_file = 'evaluation/predictions/result_aug03.csv'
output_file = 'evaluation/predictions/10_result_aug03.csv'
save_10_valid_smiles(input_file, output_file)