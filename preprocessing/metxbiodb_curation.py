from rdkit import Chem
import pandas as pd
from sklearn.model_selection import GroupShuffleSplit
from sklearn.model_selection import train_test_split
from standardize_smiles_ours import standardize_molecule, standardize_smiles_collection, standardize_smiles_main

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

def modifiy_columns(df): # Returns a df ready for fine-tuning
    try:
        
        # Create a new DataFrame with the desired structure
        new_df = pd.DataFrame({
            'products': df['child_smiles'],  # Map 'parent_smiles' to 'products'
            'reactants': df['parent_smiles'],  # Map 'child_smiles' to 'reactants'
            'set': df['set']
        })
        
        
        print(f"Transformation completed. The df now has columns: products, reactants, and set (empty)")
    
    except KeyError as e:
        print(f"Error: Missing required column in the input data: {e}")
    except Exception as e:
        print(f"An unexpected error occurred: {e}")

    return new_df

def get_unique_parents(df):
    # Ensure df has 'parent_smiles' column, else raise an error
    if 'parent_smiles' not in df.columns:
        raise ValueError("The DataFrame must contain a 'parent_smiles' column.")

    unique_parent_dataset = []

    # Get a list of columns excluding 'parent_smiles'
    other_columns = [col for col in df.columns if col != 'parent_smiles']

    # Loop through unique parent SMILES
    for parent in df['parent_smiles'].unique():
        # Find the first row where 'parent_smiles' matches the current parent
        row = df[df['parent_smiles'] == parent].iloc[0]
        # Append the entire row
        unique_parent_dataset.append(row)

    # Create a new DataFrame from the unique dataset, preserving other columns
    unique_parent_df = pd.DataFrame(unique_parent_dataset).reset_index(drop=True)

    unique_parent_df = unique_parent_df[['parent_smiles'] + other_columns]

    return unique_parent_df


def test_val_distribute(df): #Input csv file. Function adds an extra column named 'set' and distributed into 'val' and 'set'

    # Add an empty column named 'set'
    df['set'] = None

    # Split the DataFrame into train (90%) and validation (10%) sets.
    # Random_state is for reproducibility
    train_df, val_df = train_test_split(df, test_size=0.1, random_state=26)

    # Assign "train" and "val" labels
    train_df['set'] = 'train'
    val_df['set'] = 'val'

    # Combine the train and validation DataFrames
    final_df = pd.concat([train_df, val_df]).reset_index(drop=True)

    print("The df is now distributed into training and validation")

    return final_df

def df_to_csv(input_df, output_csv): # outputs a csv with tabs
    input_df.to_csv(output_csv, sep='\t', index=False)


def curate_metxbiodb_main():
    metxbiodb_df = load_metxbiodb() # Loads raw metxbiodata
    smiles_metxbiodb_df = metxbiodb_inchi_to_smiles(metxbiodb_df) # Converts the data to smiles. Also creates csv file
    smiles_metxbiodb_set_df = test_val_distribute(smiles_metxbiodb_df) # adds a set column and distrubutes the datapoints into val and set
    smiles_metxbiodb_clean_df = standardize_smiles_main(smiles_metxbiodb_set_df) # Cleans data. Input df needs columns 'parent_smiles' and 'child_smiles'
    smiles_metxbiodb_mod_columns_df = modifiy_columns(smiles_metxbiodb_clean_df) # Outputs df ready for fine-tuning but the set-column is empty
    
    df_to_csv(smiles_metxbiodb_clean_df, 'dataset/dummy_data/preprocessed_metxbiodb/metxbiodb_clean.csv')
    df_to_csv(smiles_metxbiodb_mod_columns_df, 'dataset/dummy_data/preprocessed_metxbiodb/metxbiodb_clean_finetuning.csv')

    len_start = len(smiles_metxbiodb_df)
    len_clean = len(smiles_metxbiodb_clean_df)

    print("Curation with NO unique parents")
    print("Number of data points from start: ", len_start)
    print("Number of data points after cleaning: ", len_clean)

    val_count_start = (smiles_metxbiodb_set_df['set'] == 'val').sum()
    val_count_clean = (smiles_metxbiodb_clean_df['set'] == 'val').sum()

    print(f"The count of 'val' from start is: {val_count_start}")
    print(f"The count of 'val' after cleaning is: {val_count_clean}")
    print('The val distribution is now', val_count_clean/len_clean)


def curate_metxbiodb_unique_parents_main():
    metxbiodb_df = load_metxbiodb() # Loads raw metxbiodata
    smiles_metxbiodb_df = metxbiodb_inchi_to_smiles(metxbiodb_df) # Converts the data to smiles. Also creates csv file
    smiles_metxbiodb_set_df = test_val_distribute(smiles_metxbiodb_df) # adds a set column and distrubutes the datapoints into val and set
    smiles_metxbiodb_clean_df = standardize_smiles_main(smiles_metxbiodb_set_df) # Cleans data. Input df needs columns 'parent_smiles' and 'child_smiles'
    smiles_metxbiodb_clean_unique_df = get_unique_parents(smiles_metxbiodb_clean_df) # Gets a data set with only unique parents
    smiles_metxbiodb_mod_columns_df = modifiy_columns(smiles_metxbiodb_clean_unique_df) # Outputs df ready for fine-tuning but the set-column is empty
    
    df_to_csv(smiles_metxbiodb_mod_columns_df, 'dataset/dummy_data/preprocessed_metxbiodb/metxbiodb_clean_unique_parents_finetuning.csv')
    df_to_csv(smiles_metxbiodb_clean_unique_df, 'dataset/dummy_data/preprocessed_metxbiodb/metxbiodb_clean_unique_parents.csv')

    len_start = len(smiles_metxbiodb_df)
    len_clean = len(smiles_metxbiodb_clean_df)
    len_unique_parents = len(smiles_metxbiodb_clean_unique_df)

    print("Curation WITH unique parents")
    print("Number of data points from start: ", len_start)
    print("Number of data points after cleaning: ", len_clean)
    print("Number of data points with after cleaning and unique parents: ", len_unique_parents)
    
    val_count_start = (smiles_metxbiodb_set_df['set'] == 'val').sum()
    val_count_clean_unique = (smiles_metxbiodb_clean_unique_df['set'] == 'val').sum()

    print(f"The count of 'val' from start is: {val_count_start}")
    print(f"The count of 'val' after cleaning and unique parents is: {val_count_clean_unique}")

    print('The val distribution is now', val_count_clean_unique/len_unique_parents)

def curate_metxbiodb_unique_parents_no_clean_main():
    metxbiodb_df = load_metxbiodb() # Loads raw metxbiodata
    smiles_metxbiodb_df = metxbiodb_inchi_to_smiles(metxbiodb_df) # Converts the data to smiles. Also creates csv file
    smiles_metxbiodb_set_df = test_val_distribute(smiles_metxbiodb_df) # adds a set column and distrubutes the datapoints into val and set

    smiles_metxbiodb_unique_df = get_unique_parents(smiles_metxbiodb_set_df) # Gets a data set with only unique parents
    smiles_metxbiodb_mod_columns_df = modifiy_columns(smiles_metxbiodb_unique_df) # Outputs df ready for fine-tuning but the set-column is empty
    
    df_to_csv(smiles_metxbiodb_mod_columns_df, 'dataset/dummy_data/preprocessed_metxbiodb/metxbiodb_no_clean_unique_parents_finetuning.csv')
    df_to_csv(smiles_metxbiodb_unique_df, 'dataset/dummy_data/preprocessed_metxbiodb/metxbiodb_no_clean_unique_parents.csv')

    len_start = len(smiles_metxbiodb_df)
    len_unique_parents = len(smiles_metxbiodb_unique_df)

    print("Curation WITH unique parents")
    print("Number of data points from start: ", len_start)
    print("Number of data points unique parents: ", len_unique_parents)
    
    val_count_start = (smiles_metxbiodb_set_df['set'] == 'val').sum()
    val_count_unique = (smiles_metxbiodb_unique_df['set'] == 'val').sum()

    print(f"The count of 'val' from start is: {val_count_start}")
    print(f"The count of 'val' after cleaning and unique parents is: {val_count_unique}")

    print('The val distribution is now', val_count_unique/len_unique_parents)

#curate_metxbiodb_main()
#curate_metxbiodb_unique_parents_main()
curate_metxbiodb_unique_parents_no_clean_main()

#RESULTS:

'''
There are three data points with missing data. Thus the first step gets rid of these points. These are not counted

'''
#random state 26#
'''
Curation with cleaning and NO unique parents
Number of data points from start:  2130
Number of data points after cleaning:  2130
The count of 'val' from start is: 213
The count of 'val' after cleaning is: 213
The val distribution is now 0.1

Curation with cleaning AND unique parents
Number of data points from start:  2130
Number of data points after cleaning:  2130
Number of data points with after cleaning and unique parents:  1205
The count of 'val' from start is: 213
The count of 'val' after cleaning and unique parents is: 72
The val distribution is now 0.05975103734439834

Curation WITH unique parents and NO cleaning
Number of data points from start:  2130
Number of data points unique parents:  1251
The count of 'val' from start is: 213
The count of 'val' after cleaning and unique parents is: 76
The val distribution is now 0.060751398880895285

'''



#random state 13#
'''
Curation with NO unique parents
Number of data points from start:  2130
Number of data points after cleaning:  2130
The count of 'val' from start is: 213
The count of 'val' after cleaning is: 213

Curation WITH unique parents
Number of data points from start:  2130
Number of data points after cleaning:  2130
Number of data points with after cleaning and unique parents:  1205
The val distribution is now 0.06390041493775933
The count of 'val' from start is: 213
The count of 'val' after cleaning and unique parents is: 77
'''




#random state 42#
'''
Curation with NO unique parents
Number of data points from start:  2133
Number of data points after cleaning:  2130
The count of 'val' from start is: 213
The count of 'val' after cleaning is: 213
'''

#did unique parents first and then cleaning, so not how the code is now
'''
Curation WITH unique parents
Number of data points from start:  2133
Number of data points with only unique parents:  1251
Number of data points after cleaning:  1251
The count of 'val' from start is: 213
The count of 'val' after unique parents cleaning is: 63
'''
