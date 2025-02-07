from rdkit import Chem
import pandas as pd
from sklearn.model_selection import GroupShuffleSplit
from sklearn.model_selection import train_test_split

print("hello")

def load_metxbiodb():
    file = 'dataset/raw_data/metxbiodb.csv'
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

def modifiy_columns(df):
    try:
        
        # Create a new DataFrame with the desired structure
        new_df = pd.DataFrame({
            'products': df['child_smiles'],  # Map 'parent_smiles' to 'products'
            'reactants': df['parent_smiles'],  # Map 'child_smiles' to 'reactants'
            'set': ''  # Initialize 'set' column as empty
        })
        
        
        print(f"Transformation completed. The df now has columns: products, reactants, and set (empty)")
    
    except KeyError as e:
        print(f"Error: Missing required column in the input data: {e}")
    except Exception as e:
        print(f"An unexpected error occurred: {e}")

    return df

def get_unique_parent_df(df): #input df must have columns produtcts, reactants, set
    unique_parent_dataset = []
    
    for parent in df["reactants"].unique(): 
        children = df.loc[df['reactants'] == parent]["products"]
        unique_parent_dataset.append((parent, list(children)[0]))

    unique_parent_df = pd.DataFrame(unique_parent_dataset)
    unique_parent_df.columns = ["parent_smiles", "child_smiles"]

    unique_parent_df.to_csv('unique_parents_metxbiodb.csv', index=False)
    return unique_parent_df


def test_val_distribute(df): #Input df files has columns 'products, reactants, set', where the set column is empty

    # Split the DataFrame into train (90%) and validation (10%) sets.
    # Random_state is for reproducibility
    train_df, val_df = train_test_split(df, test_size=0.1, random_state=42)

    # Assign "train" and "val" labels
    train_df['set'] = 'train'
    val_df['set'] = 'val'

    # Combine the train and validation DataFrames
    final_df = pd.concat([train_df, val_df]).reset_index(drop=True)

    print("The df is now distributed into training and validation")

    return final_df

def df_to_csv(input_df, output_csv):
    input_df.to_csv(output_csv, sep='\t', index=False)



metxbiodb_df = load_metxbiodb() # Loads raw metxbiodata
smiles_metxbiodb_df = metxbiodb_inchi_to_smiles(metxbiodb_df) # Converts the data to smiles. Also creates csv file
smiles_metxbiodb_mod_columns_df = modifiy_columns(smiles_metxbiodb_df) # Outputs df ready for fine-tuning but the set-column is empty
metxbiodb_finetuning_ready_all_df = test_val_distribute(smiles_metxbiodb_mod_columns_df) # Distrubutes the datapoints into val and set

len_metxbiodb = len(metxbiodb_finetuning_ready_all_df)

df_to_csv(metxbiodb_finetuning_ready_all_df, 'dataset/dummy_data/smiles_metxbiodb_all')




### TO DO ###
#adapt unique parents and potentially cleaning









#metxbiodb_df = load_metbxbiodb()
#smiles_metxbiodb_df = metxbiodb_inchi_to_smiles(metxbiodb_df)
#modifiy_for_training('smiles_metxbiodb.csv')

#get_unique_parent_dataset(smiles_metxbiodb_df)
#modifiy_for_training('dataset/dummy_data/unique_parents_metxbiodb.csv', 'dataset/dummy_data/unique_parents_metxbiodb_finetuning.csv')
#test_val_distribution('unique_parents_metxbiodb_finetuning.csv')