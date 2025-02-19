from rdkit import Chem
import pandas as pd
from sklearn.model_selection import GroupShuffleSplit
from sklearn.model_selection import train_test_split


def load_metxbiodb(file):
    # Columns in raw data:
    # "biotid","substrate_name","substrate_cid","substrate_inchikey","substrate_inchi",
    # "enzyme","reaction_type","biotransformation_type","biosystem","prod_name","prod_cid",
    # "prod_inchikey","prod_inchi","reference"
    try:
        metxbiodb_df = pd.read_csv(file)
        metxbiodb_df = metxbiodb_df[['substrate_name','substrate_inchi','prod_name','prod_inchi']]
        return metxbiodb_df
    
    except FileNotFoundError:
        print("Error: The file ", file,  " was  not found. Please check the file path.")
    except Exception as e:
        print(f"An error occurred: {e}")

def metxbiodb_inchi_to_smiles(data, output_file):
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
  
    smiles_metxbiodb_df.to_csv(output_file, index=False)


if __name__ == "__main__":

    dataset = 'dataset/raw_data/metxbiodb.csv'
    final = 'dataset/curated_data/metxbiodb_smiles.csv'

    metxbiodb = load_metxbiodb(dataset)
    metxbiodb_inchi_to_smiles(metxbiodb, final)

