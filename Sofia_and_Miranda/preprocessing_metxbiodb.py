from rdkit import Chem
import pandas as pd
from sklearn.model_selection import GroupShuffleSplit

def load_metbxbiodb():
    # Columns in raw data:
    # "biotid","substrate_name","substrate_cid","substrate_inchikey","substrate_inchi",
    # "enzyme","reaction_type","biotransformation_type","biosystem","prod_name","prod_cid",
    # "prod_inchikey","prod_inchi","reference"
    try:
        metxbiodb_df = pd.read_csv('raw_data/metxbiodb.csv')
        metxbiodb_df = metxbiodb_df[['substrate_name','substrate_inchi','prod_name','prod_inchi']]
        return metxbiodb_df
    
    except FileNotFoundError:
        print("Error: The file was not found. Please check the file path.")
    except Exception as e:
        print(f"An error occurred: {e}")

def metxbiodb_inchi_to_smiles(data):
    parent_child = []
    counter = 0
    for ind in data.index:
        try:
            parent_mol = Chem.inchi.MolFromInchi(data.loc[ind]["substrate_inchi"]) #metxbio_df[['substrate_name','substrate_inchi','prod_name','prod_inchi']]
            child_mol = Chem.inchi.MolFromInchi(data.loc[ind]["prod_inchi"])
            parent_smiles = Chem.rdmolfiles.MolToSmiles(parent_mol)
            child_smiles = Chem.rdmolfiles.MolToSmiles(child_mol)

            parent_name = data.loc[ind]["substrate_name"]
            child_name = data.loc[ind]["prod_name"]

            parent_child.append([parent_name, parent_smiles, child_name, child_smiles])
        except:
            counter += 1
            print("ERROR MISSING DATA, number: ", counter)
            print(data.loc[ind])
            
    smiles_metxbiodb_df = pd.DataFrame(parent_child)
    smiles_metxbiodb_df.columns = ['parent_name', 'parent_smiles', 'child_name', 'child_smiles']
  
    smiles_metxbiodb_df.to_csv('smiles_metxbiodb.csv', index=False)
    return smiles_metxbiodb_df

def modifiy_for_training(input_csv, output_csv = 'smiles_metxbiodb_finetuning.csv'):
    try:
        df = pd.read_csv(input_csv)
        
        # Create a new DataFrame with the desired structure
        new_df = pd.DataFrame({
            'products': df['child_smiles'],  # Map 'parent_smiles' to 'products'
            'reactants': df['parent_smiles'],  # Map 'child_smiles' to 'reactants'
            'set': ''  # Initialize 'set' column as empty
        })
        
        # Save the transformed DataFrame to the output CSV
        new_df.to_csv(output_csv, index=False)
        
        print(f"Transformation completed. The new file is saved as: {output_csv}")
    
    except FileNotFoundError:
        print(f"Error: The file {input_csv} was not found.")
    except KeyError as e:
        print(f"Error: Missing required column in the input data: {e}")
    except Exception as e:
        print(f"An unexpected error occurred: {e}")

def get_unique_parent_dataset(data): 
    unique_parent_dataset = []
    
    for parent in data["parent_smiles"].unique(): 
        children = data.loc[data['parent_smiles'] == parent]["child_smiles"]
        unique_parent_dataset.append((parent, list(children)[0]))

    unique_parent_df = pd.DataFrame(unique_parent_dataset)
    unique_parent_df.columns = ["parent_smiles", "child_smiles"]

    unique_parent_df.to_csv('unique_parents_metxbiodb.csv', index=False)
    return unique_parent_df

import pandas as pd
from sklearn.model_selection import train_test_split


def test_val_distribution(input_csv): #Input csv files has columns 'products, reactants, set', where the set column is empty
    
    df = pd.read_csv(input_csv)

    # Split the DataFrame into train (90%) and validation (10%) sets.
    # Random_state is for reproducibility
    train_df, val_df = train_test_split(df, test_size=0.1, random_state=42)

    # Assign "train" and "val" labels
    train_df['set'] = 'train'
    val_df['set'] = 'val'

    # Combine the train and validation DataFrames
    final_df = pd.concat([train_df, val_df]) #.reset_index(drop=True)

    # Save the DataFrame with the updated "set" column to a new CSV file. The old one will be overwritten
    final_df.to_csv(input_csv, sep='\t', index=False)

    print("The file ", input_csv, "is now divided into validation and training")




#metxbiodb_df = load_metbxbiodb()
#smiles_metxbiodb_df = metxbiodb_inchi_to_smiles(metxbiodb_df)
#modifiy_for_training('smiles_metxbiodb.csv')

#get_unique_parent_dataset(smiles_metxbiodb_df)
#modifiy_for_training('unique_parents_metxbiodb.csv', 'unique_parents_metxbiodb_finetuning.csv')
test_val_distribution('unique_parents_metxbiodb_finetuning.csv')