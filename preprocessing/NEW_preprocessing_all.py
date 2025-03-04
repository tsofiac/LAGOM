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

# --------------------------Standardising SMILES and removing duplicates in each data set -----------------------------------

def standardize_smiles(input_file): #input df that has columns 'parent_smiles' and 'child_smiles'
    df = pd.read_csv(input_file)
    df['parent_smiles'] = standardize_smiles_collection(df['parent_smiles'], False) #'False' eliminates isomeres
    df['child_smiles'] = standardize_smiles_collection(df['child_smiles'], False)
    return df

def remove_duplicates(df, duplicates_data_file): #Removes duplicate reactions
 
    len_before = len(df)
    duplicates_df = df[df.duplicated(subset=['parent_smiles', 'child_smiles'], keep=False)]
    df = df.drop_duplicates(subset=['parent_smiles', 'child_smiles'], keep='first') 
    len_after = len(df)
    print("Total data points removed with remove_duplicates:", len_before - len_after)

    # Only save the duplicates DataFrame if it is not empty
    if not duplicates_df.empty:
        duplicates_df.to_csv(duplicates_data_file, index=False)

    return df

def non_equal_smiles(parent_smiles, child_smiles):
    # Returns True if the smiles are not equal, False otherwise.
    return parent_smiles != child_smiles

def remove_equal_parent_child(data, removed_data_file):

    total_removed = 0

    allowed_molecules = [non_equal_smiles(row['parent_smiles'], row['child_smiles']) for index, row in data.iterrows()]
    filtered_data = data[allowed_molecules]
    removed_data = data[[not allowed for allowed in allowed_molecules]]
    total_removed += allowed_molecules.count(False)

    if not removed_data.empty:
        removed_data.to_csv(removed_data_file, index=False)

    print(f"Total data points removed with remove_equal_parent_child: {total_removed}")

    return filtered_data

# --------------------------Combining datasets and removing duplicates with each other and with test data set -----------------------------------

def add_source_column(df, source):

    if source == 'drugbank':
        df['source'] = 'DrugBank'
    elif source == 'metxbiodb':
        df['source'] = 'MetXBioDB'
    elif source == 'gloryx':
        df['source'] = 'GLORYx' 
    else:
        print["Look over this."]
    
    return df
    
def combine_datasets(df1, df2, output_csv):
    selected_df1 = df1[['parent_name', 'child_name', 'parent_smiles', 'child_smiles', 'origin', 'source']].copy()
    selected_df2 = df2[['parent_name', 'child_name', 'parent_smiles', 'child_smiles', 'origin', 'source']].copy()

    combined_df = pd.concat([selected_df1, selected_df2], ignore_index=True)
    combined_df.to_csv(output_csv, index=False)

def remove_duplicates_combined(combined_csv, removed_duplicates_csv):  # Removes duplicate reactions and modifies source column
    df = pd.read_csv(combined_csv)
    len_before = len(df)
    
    # Create a mask that marks the duplicates
    df['is_duplicate'] = df.duplicated(subset=['parent_smiles', 'child_smiles'], keep=False)  # Keep all duplicates
    
    # Separate the duplicates (keep all versions)
    duplicates_df = df[df['is_duplicate']]
    
    # Modify the 'source' column for the rows that are kept
    df.loc[df['is_duplicate'], 'source'] = 'Both'  # Change 'source' value for rows that are kept

    # Remove the 'is_duplicate' column as it's no longer needed
    df = df.drop(columns='is_duplicate')
    duplicates_df = duplicates_df.drop(columns='is_duplicate')
    
    # Remove duplicates based on 'parent_smiles' and 'child_smiles', keeping the first occurrence
    df = df.drop_duplicates(subset=['parent_smiles', 'child_smiles'], keep='first')

    len_after = len(df)
    
    print("Total data points removed due to duplicates in both data sets:", len_before - len_after)

    df.to_csv(combined_csv, index=False)
    duplicates_df.to_csv(removed_duplicates_csv, index=False)


def compare_datasets(combined_csv, testdata_csv, removed_file):
    df1 = pd.read_csv(combined_csv)
    df2 = pd.read_csv(testdata_csv)
    # duplicate parent_smiles between df1 and df2
    duplicate_parent_smiles = set(df1['parent_smiles']).intersection(set(df2['parent_smiles']))
    # those not in the duplicates set
    non_duplicates_in_df1 = df1[~df1['parent_smiles'].isin(duplicate_parent_smiles)].copy()
    duplicates_df = df1[df1['parent_smiles'].isin(duplicate_parent_smiles)].copy()

    print('Number of duplicate parents reactions found: ', len(df1)-len(non_duplicates_in_df1))
    non_duplicates_in_df1.to_csv(combined_csv, index=False)
    duplicates_df.to_csv(removed_file, index=False)

# -------------------------Splitting the data-------------------------------
def test_val_distribute(data_file, val_size): #Input csv file. Function adds an extra column named 'set' and distributed into 'val' and 'set'
    df = pd.read_csv(data_file)
    # Add an empty column named 'set'
    df['set'] = None

    ## Using train_test_split
    # train_df, val_df = train_test_split(df, test_size=val_size, random_state=26)     # Random_state is for reproducibility
    
    # Assign "train" and "val" labels
    # train_df['set'] = 'train'
    # val_df['set'] = 'val'

    ## Using GroupShuffleSplit
    splitter = GroupShuffleSplit(test_size=val_size, n_splits=2, random_state=42)
    split = splitter.split(df, groups=df['origin']) 
    train_inds, val_inds = next(split)
    train_df = df.iloc[train_inds]
    val_df = df.iloc[val_inds]

    # Assign "train" and "val" labels
    train_df.loc[:, 'set'] = 'train'
    val_df.loc[:, 'set'] = 'val'

    # Combine the train and validation DataFrames
    final_df = pd.concat([train_df, val_df]).reset_index(drop=True)

    final_df.to_csv(data_file, index=False)

# --------------------------Filtering data--------------------------------------

def filter_data_on_both_sides(data_file, filter_method, removed_data_file): 
    data = pd.read_csv(data_file)

    # Ensure SMILES columns are strings
    data["child_smiles"] = data["child_smiles"].fillna('').astype(str)
    data["parent_smiles"] = data["parent_smiles"].fillna('').astype(str)

    total_removed = 0

    # Filtering based on child molecules
    allowed_children = [filter_method(metabolite) for metabolite in data["child_smiles"]]
    filtered_data = data[allowed_children]
    removed_data_children = data[[not allowed for allowed in allowed_children]]
    total_removed += allowed_children.count(False)

    # Filtering based on parent molecules
    allowed_parents = [filter_method(molecule) for molecule in filtered_data["parent_smiles"]]
    removed_data_parents = filtered_data[[not allowed for allowed in allowed_parents]]
    filtered_data = filtered_data[allowed_parents]
    total_removed += allowed_parents.count(False)

    removed_data = pd.concat([removed_data_children, removed_data_parents])

    filtered_data.to_csv(data_file, index=False)

    if not removed_data.empty:
        removed_data.to_csv(removed_data_file, index=False)

    print(f"Total data points removed with {filter_method.__name__}: {total_removed}")

def filter_data_on_one_side(data_file, filter_method, removed_data_file, if_parent = True): 
    name_property = "parent_smiles" if if_parent else "child_smiles"
    data = pd.read_csv(data_file)

    data[name_property] = data[name_property].fillna('').astype(str)

    total_removed = 0

    # Filtering based on parent molecules
    allowed_molecules = [filter_method(molecule) for molecule in data[name_property]]
    filtered_data = data[allowed_molecules]
    removed_data = data[[not allowed for allowed in allowed_molecules]]
    total_removed += allowed_molecules.count(False)

    filtered_data.to_csv(data_file, index=False)
    if not removed_data.empty:
        removed_data.to_csv(removed_data_file, index=False)

    print(f"Total data points removed with {filter_method.__name__}: {total_removed}")

# filtering methods
# ----------------------------------------------------------------
def valid_smiles(molecule): 
    try:
        return Chem.MolFromSmiles(molecule) is not None
    except:
        return False

def atoms_allowed_in_molecules(molecule): 
    try:
        atoms_to_include = ['C', 'N', 'S', 'O', 'H', 'F', 'I', 'P', 'Cl', 'Br']
        mol = Chem.MolFromSmiles(molecule)
        atoms = [atom.GetSymbol() for atom in mol.GetAtoms()]
        return set(atoms).issubset(set(atoms_to_include))
    except:
        return False

def molecule_allowed_based_on_weight(molecule, max_weight=750, min_weight=100): 
    try:
        mol_weight = Descriptors.ExactMolWt(Chem.MolFromSmiles(molecule))
        return min_weight <= mol_weight <= max_weight
    except:
        return False
# ----------------------------------------------------------------

def define_fingerprint_similarity(dataset):
    parent_smiles = dataset['parent_smiles'].tolist()
    child_smiles = dataset['child_smiles'].tolist()

    parent_mol = [Chem.MolFromSmiles(x) for x in parent_smiles]
    child_mol = [Chem.MolFromSmiles(x) for x in child_smiles]

    parent_fps = [AllChem.GetMorganFingerprintAsBitVect(x, radius=2, nBits=1024) for x in parent_mol]
    child_fps = [AllChem.GetMorganFingerprintAsBitVect(x, radius=2, nBits=1024) for x in child_mol]

    fingerprint_similarities = []
    for i in range(len(parent_smiles)):
        
        s = DataStructs.TanimotoSimilarity(parent_fps[i], child_fps[i]) 
        fingerprint_similarities.append(s)

    dataset['tanimoto'] = fingerprint_similarities

    return dataset

def fingerprints_allowed(similarity, min_similarity):
    if similarity < min_similarity: # or similarity == 1:
        return False
    return True

def filter_fingerprint_similarity(data_file, removed_data_file, min_similarity = 0.20):
    data = pd.read_csv(data_file)
    data = define_fingerprint_similarity(data)

    total_removed = 0

    allowed_molecules = [fingerprints_allowed(value, min_similarity) for value in data["tanimoto"]]
    filtered_data = data[allowed_molecules]
    removed_data = data[[not allowed for allowed in allowed_molecules]]
    total_removed += allowed_molecules.count(False)

    filtered_data.to_csv(data_file, index=False)
    if not removed_data.empty:
        removed_data.to_csv(removed_data_file, index=False)

    print(f"Total data points removed with fingerprint_similarity_filter: {total_removed}")

# ---------------------Reformat for Chemformer-------------------------------------

def get_unique_parents(input_file, output_file):
    df = pd.read_csv(input_file)
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

    val_count_end = (unique_parent_df['set'] == 'val').sum()

    print('Unique parent val distribution: ', val_count_end/len(unique_parent_df))

    unique_parent_df.to_csv(output_file, index=False)

def reformat_for_chemformer(input_file, output_file):
    df = pd.read_csv(input_file)

    df = df.rename(columns={
        "child_smiles": "products",
        "parent_smiles": "reactants",
    })

    df.to_csv(output_file, sep='\t', index=False)


    
if __name__ == "__main__":

    name = 'combined' # [ 'combined' 'drugbank' 'metxbiodb' 'augmented' ]
    preprocess_unique_parents = True

    val_size = 0.1
    min_similarity = 0.2
    
    clean_csv = f'dataset/curated_data/{name}_smiles_clean.csv'
    dataset_gloryx = 'dataset/curated_data/gloryx_smiles_clean.csv'

    removed_duplicates = f'dataset/removed_data/{name}_removed_duplicates.csv'
    removed_equal = f'dataset/removed_data/{name}_removed_equal.csv'
    removed_valid_smiles = f'dataset/removed_data/{name}_removed_valid_smiles.csv'
    removed_atoms_allowed = f'dataset/removed_data/{name}_removed_atoms_allowed.csv'
    removed_weights_allowed = f'dataset/removed_data/{name}_removed_weights_allowed.csv'
    removed_fingerprints = f'dataset/removed_data/{name}_removed_fingerprints.csv'

    finetune_csv = f'dataset/finetune/{name}_finetune.csv'
    unique = f'dataset/curated_data/{name}_unique_parents.csv'
    unique_finetune = f'dataset/finetune/{name}_unique_parents_finetune.csv'

    if name == 'combined':

        dataset_metx = 'dataset/curated_data/metxbiodb_smiles.csv'
        dataset_drugbank = 'dataset/curated_data/drugbank_smiles.csv'

        compare_removed_csv = 'dataset/removed_data/compare_removed_duplicates.csv'

        df_metx = standardize_smiles(dataset_metx)
        df_drugbank = standardize_smiles(dataset_drugbank)

        df_metx = remove_duplicates(df_metx, 'dataset/removed_data/metxbiodb_removed_duplicates.csv')
        df_drugbank = remove_duplicates(df_drugbank, 'dataset/removed_data/drugbank_removed_duplicates.csv') 

        df_metx = remove_equal_parent_child(df_metx, 'dataset/removed_data/metxbiodb_removed_equal.csv')
        df_drugbank = remove_equal_parent_child(df_drugbank, 'dataset/removed_data/drugbank_removed_equal.csv')

        df_metx = add_source_column(df_metx, 'metxbiodb')
        df_drugbank = add_source_column(df_drugbank, 'drugbank')

        combine_datasets(df_drugbank, df_metx, clean_csv)
        remove_duplicates_combined(clean_csv, removed_duplicates)
        # ----- to add augmented reactions to combined dataset --------
        df_clean = pd.read_csv(clean_csv)
        df_augmented = pd.read_csv('dataset/curated_data/augmented_data_clean.csv')
        # -------------------------------------------------------------
        combine_datasets(df_clean, df_augmented, clean_csv)
        compare_datasets(clean_csv, dataset_gloryx, compare_removed_csv)

        test_val_distribute(clean_csv, val_size)

        filter_data_on_both_sides(clean_csv, valid_smiles, removed_valid_smiles)
        filter_data_on_both_sides(clean_csv, atoms_allowed_in_molecules, removed_atoms_allowed)
        filter_data_on_one_side(clean_csv, molecule_allowed_based_on_weight, removed_weights_allowed, True)
        filter_fingerprint_similarity(clean_csv, removed_fingerprints, min_similarity)

        reformat_for_chemformer(clean_csv, finetune_csv)

        if preprocess_unique_parents:
            get_unique_parents(clean_csv, unique)
            reformat_for_chemformer(unique, unique_finetune)

    elif name == "augmented":
        
        augmentated_csv = 'dataset/curated_data/augmented_data_clean.csv'

        df_augmented = standardize_smiles('dataset/curated_data/augmented_data.csv')
        df_augmented = remove_duplicates(df_augmented, removed_duplicates)
        df_augmented = remove_equal_parent_child(df_augmented, removed_equal)

        df_augmented.to_csv(augmentated_csv, index=False)

        filter_data_on_both_sides(augmentated_csv, valid_smiles, removed_valid_smiles)
        filter_data_on_both_sides(augmentated_csv, atoms_allowed_in_molecules, removed_atoms_allowed)
        filter_data_on_one_side(augmentated_csv, molecule_allowed_based_on_weight, removed_weights_allowed, True)
        filter_fingerprint_similarity(augmentated_csv, removed_fingerprints, min_similarity)

    else:

        dataset = f'dataset/curated_data/{name}_smiles.csv'
        
        compare_removed_csv = f'dataset/removed_data/{name}_compare_removed_duplicates.csv'

        df = standardize_smiles(dataset)
        df = remove_duplicates(df, removed_duplicates)
        df = remove_equal_parent_child(df, removed_equal)

        df.to_csv(clean_csv, index=False)

        compare_datasets(clean_csv, dataset_gloryx, compare_removed_csv)
        test_val_distribute(clean_csv, val_size)

        filter_data_on_both_sides(clean_csv, valid_smiles, removed_valid_smiles)
        filter_data_on_both_sides(clean_csv, atoms_allowed_in_molecules, removed_atoms_allowed)
        filter_data_on_one_side(clean_csv, molecule_allowed_based_on_weight, removed_weights_allowed, True)
        filter_fingerprint_similarity(clean_csv, removed_fingerprints, min_similarity)
        
        reformat_for_chemformer(clean_csv, finetune_csv)

        if preprocess_unique_parents:
            get_unique_parents(clean_csv, unique)
            reformat_for_chemformer(unique, unique_finetune)




    #      metxbiodb_df = add_source_column('dataset/curated_data/metxbiodb_smiles_clean.csv', 'metxbiodb')
    # drugbank_df = add_source_column('dataset/curated_data/drugbank_smiles_clean.csv', 'drugbank')
    # gloryx_df = add_source_column('dataset/curated_data/gloryx_smiles_clean.csv', 'gloryx')

    # combine_datasets(drugbank_df, metxbiodb_df, combined_csv)
    # remove_duplicates(combined_csv, removed_csv)

    # compare_datasets(combined_csv, gloryx_df, compare_removed_csv)



    #      df.to_csv(output_file, index=False)

    # else: 
        