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

def find_drug_origin_metxbiodb(data_file):
    df = pd.read_csv(data_file)

    # Initialize origins as parent_name or 'metxbiodb unknown' if NaN
    df['origin'] = df['parent_name'].fillna('metxbiodb unknown')

    max_iterations = len(df)  # Upper limit to avoid infinite loops
    for _ in range(max_iterations):
        update_made = False  # Flag to track if any updates are made
        for index, row in df.iterrows():
            current_origin = row['origin']

            # Check if the current origin is a child elsewhere, update if necessary
            for _, pot_parent_row in df.iterrows():
                if current_origin == pot_parent_row['child_name']:
                    new_origin = pot_parent_row['parent_name']
                    if current_origin != new_origin:
                        df.at[index, 'origin'] = new_origin
                        update_made = True

        # Exit the loop early if no updates are made
        if not update_made:
            break

    # Write the updated DataFrame back to a CSV
    df.to_csv(data_file, index=False)

def check_origin_vs_parent(data_file): #just for checking
    """
    Function to check if any row in the DataFrame has 'origin' not equal to 'parent_name'.
    
    Returns True if any such discrepancy exists and False otherwise.
    """
    # Read the DataFrame
    df = pd.read_csv(data_file)

    counter = 0

    # Check for discrepancies where 'origin' is not equal to 'parent_name'
    discrepancy_exists = False
    for index, row in df.iterrows():
        if row['origin'] != row['parent_name']:
            discrepancy_exists = True
            print(f"Discrepancy found at index {index}:")
            print(f"\tParent Name: {row['parent_name']}, Origin: {row['origin']}")
            counter += 1

    print("Nr of times where the parent and origin are not the same: ", counter)
    return discrepancy_exists

def check_origin_vs_parent_drugbank(data_file): #just for checking
    """
    Function to check if any row in the DataFrame has 'origin' not equal to 'parent_name'.
    
    Returns True if any such discrepancy exists and False otherwise.
    """
    # Read the DataFrame
    df = pd.read_csv(data_file)

    counter = 0

    # Check for discrepancies where 'origin' is not equal to 'parent_name'
    discrepancy_exists = False
    for index, row in df.iterrows():
        if row['origin'] != row['parent_id']:
            discrepancy_exists = True
            #print(f"Discrepancy found at index {index}:")
            #print(f"\tParent Name: {row['parent_name']}, Origin: {row['origin']}")
            counter += 1

    print("Nr of times where the parent and origin are not the same: ", counter)
    return discrepancy_exists

if __name__ == "__main__":

    dataset = 'dataset/raw_data/metxbiodb.csv'
    final = 'dataset/curated_data/metxbiodb_smiles.csv'

    metxbiodb = load_metxbiodb(dataset)
    metxbiodb_inchi_to_smiles(metxbiodb, final)
    find_drug_origin_metxbiodb(final)

    check_origin_vs_parent(final)

    # Comparing with drugbank
    # check_origin_vs_parent_drugbank('dataset/curated_data/drugbank_smiles.csv')

    '''
    MetXBioDB: 
    Nr of times where the parent and origin are not the same:  115

    DrugBank: 
    Nr of times where the parent and origin are not the same:  1448
    '''