#### Imports ####
import pandas as pd
from rdkit import Chem
import json

from sklearn.model_selection import GroupShuffleSplit
from preprocessing.prev_master_thesis.utils import fingerprint_similarity_single, valid_smiles, molecule_allowed_based_on_weight, atoms_allowed_in_molecules, standardize_molecule

#### Code ####
def load_metxbiodb_raw():
    #### Columns in metxbiodb
    # "biotid","substrate_name","substrate_cid","substrate_inchikey",
    # "substrate_inchi","enzyme","reaction_type","biotransformation_type",
    # "biosystem","prod_name","prod_cid","prod_inchikey","prod_inchi","reference"
 
    metXBioDb_df = pd.read_csv('./dataset/raw_data/metxbiodb.csv')
    metXBioDb_df = metXBioDb_df[['substrate_name','substrate_inchi','prod_name','prod_inchi', 'enzyme','reaction_type']]

    return metXBioDb_df

def get_metabolites_datapoint(data): 
    new_datapoints = []
    metabolites = data["metabolites"]

    for metabolite in metabolites: 
        smiles = metabolite.get("smiles")
        gen = metabolite.get("generation")
        name = metabolite.get("metaboliteName")
        metabolites = metabolite.get("metabolites")
        parent_name = data["parent_name"]
        parent_smiles = data["parent_smiles"]

        new_datapoint = {"parent_name": parent_name, "parent_smiles": parent_smiles, "metabolite_smiles": smiles, "metabolite_generation": gen, "metabolite_name": name, "metabolites": metabolites}
        new_datapoints.append(new_datapoint)
    return new_datapoints

def get_next_generations_metabolites(df): 
    next_generation_raw = df.dropna(axis=0)[["metabolites", "metabolite_name", "metabolite_smiles"]]
    next_generation_raw.columns = ["metabolites", 'parent_name', 'parent_smiles']

    next_generation_list = []
    for index in next_generation_raw.index: 
        new_datapoints = get_metabolites_datapoint(next_generation_raw.loc[index])
        [next_generation_list.append(datapoint) for datapoint in new_datapoints]

    new_generation = pd.DataFrame(next_generation_list)
    return new_generation

def exist_new_generation(df): 
    new_generation = df.dropna(axis=0)[["metabolites"]]
    return len(new_generation) > 0

def load_gloryx_raw():
    with open('./dataset/raw_data/gloryx_test_dataset.json') as f:
        data = json.load(f)
        
        # get first generation
        gloryx_df = pd.json_normalize(data, "metabolites", ["drugName", "smiles"], record_prefix="metabolite.")
        gloryx_df.columns = ['metabolite_smiles', 'metabolite_generation', 'metabolite_name', 'metabolites', 'parent_name', 'parent_smiles']

        gloryx_combined_geneneration_df = gloryx_df

        # dataframe to get the rest of generations from 
        rest_of_generations_df = gloryx_df     
        gloryx_df = gloryx_df.drop(columns=["metabolites"])
        
        first_generation = gloryx_df
        # generate metabolites for each generation
        while exist_new_generation(rest_of_generations_df): 
            next_generation = get_next_generations_metabolites(rest_of_generations_df)
            rest_of_generations_df = next_generation

            next_generation = next_generation.drop(columns=["metabolites"])
            gloryx_df = pd.concat([gloryx_df, next_generation])
        
        return gloryx_df, gloryx_combined_geneneration_df, first_generation

def metxbiodb_inchi_to_smiles(data):
    parent_child = []

    for ind in data.index:
        try:
            parent_mol = Chem.inchi.MolFromInchi(data.loc[ind]["substrate_inchi"]) #metXBioDb_df[['substrate_name','substrate_inchi','prod_name','prod_inchi']]
            child_mol = Chem.inchi.MolFromInchi(data.loc[ind]["prod_inchi"])
            parent_smiles = Chem.rdmolfiles.MolToSmiles(parent_mol)
            child_smiles = Chem.rdmolfiles.MolToSmiles(child_mol)

            parent_name = data.loc[ind]["substrate_name"]
            child_name = data.loc[ind]["prod_name"]
            enzyme = data.loc[ind]["enzyme"]
            reaction_type = data.loc[ind]["reaction_type"]

            parent_child.append([parent_name, parent_smiles, child_name, child_smiles, enzyme, reaction_type])
        except:
            print(data.loc[ind])
            
    parent_child_df = pd.DataFrame(parent_child)
    parent_child_df.columns = ['parent_name', 'parent_smiles', 'metabolite_name', 'metabolite_smiles', 'enzymes', 'reaction_type']
    # from inchi to smiles 
    return parent_child_df

def filter_out_data(data, filter_method): 
    total_removed = 0
    allowed_metabolites = [filter_method(metabolite) for metabolite in data["metabolite_smiles"]]
    data = data[allowed_metabolites]
    total_removed += allowed_metabolites.count(False)

    allowed_molecules = [filter_method(molecule) for molecule in data["parent_smiles"]]
    data = data[allowed_molecules]
    total_removed += allowed_molecules.count(False)

    return data, total_removed

def clean_data(data, cleaning_method): 
    data["metabolite_smiles"] = data["metabolite_smiles"].apply(cleaning_method)
    data["parent_smiles"] = data["parent_smiles"].apply(cleaning_method)

    return data
    
# input should have metabolite_smiles and parent_smiles
def curate_data(data):
    data, invalid_smiles = filter_out_data(data, valid_smiles) 
    data, unlikely_removed = filter_out_data(data, atoms_allowed_in_molecules) 
    data, weight_removed = filter_out_data(data, molecule_allowed_based_on_weight)    
    data = clean_data(data, standardize_molecule)
    
    total_removed = unlikely_removed + weight_removed + invalid_smiles
    print("Removed total is: " + str(total_removed) + ". Invalid smiles: " + str(invalid_smiles) + ". Atoms based: " + str(unlikely_removed) + ". Weight based: " + str(weight_removed))
    return data

# Used to remove rows that does not have sufficient data. Such as smiles or inchi
def analyze_db_remove_bad_rows(df):
    initial_row_count = len(df)
    df.dropna(subset=["parent_smiles", "parent_inchi"], how='all', inplace=True)
    df.dropna(subset=["metabolite_smiles",'metabolite_inchi'],how='all' ,inplace=True)
    final_row_count = len(df)
    dropped_count = initial_row_count - final_row_count
    print("Dropped count: ",dropped_count)
    #df.to_csv(input_file, index=False)
    return df

# Compares two dataframes and removes duplicates. Returns the source with removed duplicates 
def compare_dataset_remove_duplicates(source_df, target_df):
    equal_count = 0
    for index_src, row_src in source_df.iterrows():
        for _, row_tgt in target_df.iterrows():
            if (row_src['parent_smiles'] == row_tgt['parent_smiles'] and row_src['metabolite_smiles'] == row_tgt['metabolite_smiles']):
                equal_count+=1
                source_df = source_df.drop([index_src])
                break
      
    print("number of equal elements: ",equal_count)
    return source_df

# Parses H
def get_unique_parent_dataset(data): 
    unique_parent_dataset = []
    
    for parent in data["parent_smiles"].unique(): 
        children = data.loc[data['parent_smiles'] == parent]["metabolite_smiles"]
        unique_parent_dataset.append((parent, list(children)[0]))

    unique_parent_df = pd.DataFrame(unique_parent_dataset)
    unique_parent_df.columns = ["parent_smiles", "metabolite_smiles"]
    return unique_parent_df

def save_dataset(dataset, filename):
    df = pd.DataFrame(dataset)
    df.to_csv(filename,index=False)

# Combines two datasets 
def combine_datasets(dataset1, dataset2, remove_duplicates = True):
    # If one want to remove the duplicates in the two ds before combining them.
    if(remove_duplicates):
        # Compares the two datasets. First argument will be the one which rows will be removed. 
        dataset2 = compare_dataset_remove_duplicates(dataset2, dataset1)     
    dataset_combined = pd.concat([dataset1, dataset2])
    return dataset_combined


def prepare_metxbiodb(fps_cutoff_upper=1): 
    # metxbiodb inchi to smiles
    data = load_metxbiodb_raw()
    smiles_data = metxbiodb_inchi_to_smiles(data)
    filtered_data = curate_data(smiles_data)
      
    # Add fingerprints 
    fingerprints = []
    for _, row in filtered_data.iterrows():
            fps = round(fingerprint_similarity_single([row['parent_smiles']], [row['metabolite_smiles']])[0],2)
            fingerprints.append(fps)
    filtered_data = filtered_data.assign(fingerprint_similarity=fingerprints)
    before_count = len(filtered_data)
    filtered_data = filtered_data[filtered_data['fingerprint_similarity'] < fps_cutoff_upper]
    print(f"Based on {fps_cutoff_upper}, removed these many rows: ", before_count - len(filtered_data))
    
    # Smaller dataset
    smaller_dataset = filtered_data.head(50)
    save_dataset(smaller_dataset, 'dataset/processed_data/metxbiodb_smiles_small.csv')

    # Unique parents
    unique_parent_data = get_unique_parent_dataset(smaller_dataset)
    save_dataset(unique_parent_data, 'dataset/processed_data/metxbiodb_unique_parents.csv')

    save_dataset(filtered_data, 'dataset/processed_data/metxbiodb_smiles.csv')
    return filtered_data

def prepare_drugbank_duplicates(): 
    data_smiles = pd.read_csv(f'dataset/processed_data/drugbank_smiles_duplicates_saved.csv')
    filtered_data = curate_data(data_smiles)
    save_dataset(filtered_data, 'dataset/processed_data/drugbank_smiles_duplicates_saved.csv')

def prepare_gloryx(): 
    gloryx_df, gloryx_combined_geneneration_df, first_generation_gloryx_df = load_gloryx_raw()
    gloryx_df = curate_data(gloryx_df)
    gloryx_combined_geneneration_df = curate_data(gloryx_combined_geneneration_df)
    first_generation_gloryx_df = clean_data(first_generation_gloryx_df, standardize_molecule)
    save_dataset(gloryx_df, 'dataset/processed_data/gloryx_smiles.csv')
    save_dataset(gloryx_combined_geneneration_df, 'dataset/raw_data/gloryx_combined_generations_map.csv')
    save_dataset(first_generation_gloryx_df, 'dataset/test/gloryx_smiles_first_generation.csv') # NOTE: not curated since same as biotransformer etc

    return gloryx_df, gloryx_combined_geneneration_df, first_generation_gloryx_df
   
# prepares drugbank smiles file. Source comes from raw data and a finalized curated version is saved in processed data
def prepare_drugbank(fps_cutoff_lower = 0, fps_cutoff_upper=1):

    parent_child_smiles_finalized = "dataset/raw_data/drugbank_pairs_smiles_finalized.csv"
    drugbank_finalized = "dataset/processed_data/drugbank_smiles.csv"

    drugbank_df_raw = pd.read_csv(parent_child_smiles_finalized)
    # Checks to filter out the raw dataset    
    ## Remove bad rows
    drugbank_df_raw = analyze_db_remove_bad_rows(drugbank_df_raw)
    ## cannonicalize and remove bad data
    drugbank_df_raw = curate_data(drugbank_df_raw)
    ## Remove duplicates within the file
    drugbank_df_final = drugbank_df_raw.drop_duplicates(subset=['parent_smiles','metabolite_smiles'], keep='first')
    # Drop glucoronic acid
    drugbank_df_final = drugbank_df_final.drop(drugbank_df_final[drugbank_df_final['metabolite_name'] == "Glucuronic acid"].index)
    # Add fingerprint similarity
    fingerprints = []
    for _, row in drugbank_df_final.iterrows():
            fps = round(fingerprint_similarity_single([row['parent_smiles']], [row['metabolite_smiles']])[0],2)
            fingerprints.append(fps)
    drugbank_df_final = drugbank_df_final.assign(fingerprint_similarity=fingerprints)
    before_count = len(drugbank_df_final)
    drugbank_df_final = drugbank_df_final[drugbank_df_final['fingerprint_similarity'] < fps_cutoff_upper]
    
    if(fps_cutoff_lower > 0):
        drugbank_with_cutoff_df = drugbank_df_final[drugbank_df_final['fingerprint_similarity'] >= fps_cutoff_lower]
        print(f"Based on FPS: {fps_cutoff_lower} and {fps_cutoff_upper}, removed these many rows: ", before_count - len(drugbank_with_cutoff_df))

        drugbank_with_cutoff_df.to_csv(drugbank_finalized, index=False)
        # Save and return a cuttofed version of the dataset
        save_dataset(drugbank_with_cutoff_df, drugbank_finalized)
        return drugbank_with_cutoff_df
    # Save and return finalized dataframe
    save_dataset(drugbank_df_final, drugbank_finalized)
    return drugbank_df_final

# Augments data by creating new pairs with second or latter generational transformations
# Example of its logic
# A -> B -> C => 
# A -> B, gen 0
# A -> C, gen 1
# B -> C, gen 0
# It adds fingerprint similarity as a column.
# It appends this data to the metabolic dataset as a new one. 
def augment_data_add_to_new_dataset(input_path):
    df = pd.read_csv(input_path)
    new_pairs_matrix = []
   
   # Go over all unique parent smiles. 
    unique_parents = df['parent_smiles'].unique().tolist()
    for _,parent in enumerate(unique_parents):
        
        # Add all of the parents metabolites to a list we want to go over
        child_smiles_list = df[df['parent_smiles'] == parent]['metabolite_smiles'].tolist()
        child_smiles_list = list(zip(child_smiles_list,[0]*len(child_smiles_list)))
        child_iter = 0

        # This loop traverse every child of the parent. Each childs child will be mapped as a new reaction with the parent. I.e. parent -> grandchild 
        # The grandchilds are also added to the list while the checked child is removed
        while(len(child_smiles_list) > 0):
            active_child, generation = child_smiles_list.pop()
            active_childs_children =  df[df['parent_smiles'] == active_child]['metabolite_smiles'].tolist()
            # If a found grandchild already is in the list we check then we exlude it. If the grandchild is identical to child or parent. This is to not create an infinite loop 
            active_childs_children = [child for child in active_childs_children if child != active_child and child != parent]
            # Each of the found childs are given an generation index
            active_childs_children = list(zip(active_childs_children,[int(generation+1)]*len(active_childs_children))) 
            
            if(child_iter > 30):
                break
            
            # We add fingerprint similarity for the parent and the new child reaction
            for acc,gen in active_childs_children:
                fps = round(fingerprint_similarity_single([parent], [acc])[0],2)
                new_pairs_matrix.append([parent, acc, gen, fps ])

            child_smiles_list += active_childs_children
            child_iter+=1
    # parent_name, parent_smiles, parent_inchi, metabolite_name, metabolite_smiles, metabolite_inchi, enzymes, reaction_type
    augmented_df = pd.DataFrame(new_pairs_matrix,columns=['parent_smiles', 'metabolite_smiles','generation' ,'fingerprint_similarity'],dtype='object')
    augmented_df = augmented_df.drop_duplicates(subset=['parent_smiles','metabolite_smiles'], keep="first")
    complete_data = []
    
    # We map each new reaction with the original dataset to restore valuable information. I.e name, inchi
    for i,row in augmented_df.iterrows():
        parent_values = df[df['parent_smiles'] == row['parent_smiles']].values[0]
        metabolite_values = df[df['metabolite_smiles'] == row['metabolite_smiles']].values[0]
        parent_name = parent_values[0]
        parent_inchi = parent_values[2]
        metabolite_name = metabolite_values[3]
        metabolite_inchi = metabolite_values[5]
        enzymes = ""
        reaction_type = ""
        complete_data.append(
            [parent_name, 
             row['parent_smiles'], 
             parent_inchi, 
             metabolite_name, 
             row['metabolite_smiles'], 
             metabolite_inchi, 
             enzymes,
             reaction_type,
             row['generation'], 
             row['fingerprint_similarity']])
    complete_data_df = pd.DataFrame(complete_data, 
                                    columns=['parent_name', 
                                             'parent_smiles', 
                                             'parent_inchi', 
                                             'metabolite_name', 
                                             'metabolite_smiles',
                                             'metabolite_inchi', 
                                             'enzymes', 
                                             'reaction_type', 
                                             'generation', 
                                             'fingerprint_similarity'])
    # Add FPS to the metabolic dataset and generation value 0 
    df = df.assign(generation=[0]*len(df))
    metabolic_fps = []
    for i, row in df.iterrows():
        fps = round(fingerprint_similarity_single([row['parent_smiles']], [row['metabolite_smiles']])[0],2)
        metabolic_fps.append(fps)

    df = df.assign(fingerprint_similarity=metabolic_fps)
    # Combine the two datasets
    df = df.drop_duplicates(subset=['parent_smiles','metabolite_smiles'], keep='first')

    combined_df = combine_datasets(df,complete_data_df)
    save_dataset(complete_data_df,'dataset/train/augmented_smiles.csv' )
    save_dataset(combined_df,'dataset/train/metabolic_and_augmented_smiles.csv' )

############### PRETRAINING DATA ##################################
def prepare_matched_molecular_pair(): 
    matched_molecular_pair_df = pd.read_csv("dataset/raw_data/paired_mmp.csv")
    filtered_data = curate_data(matched_molecular_pair_df)
    save_dataset(filtered_data, 'dataset/train/matched_molecular_pairs.csv')
    return filtered_data

def get_smaller_mmp(): 
    df = pd.read_csv('dataset/processed_data/matched_molecular_pairs.csv')
    include = [sim > 0.35 for sim in fingerprint_similarity_single(df["parent_smiles"], df["metabolite_smiles"])]
    temp = include.count(True)
    print(f"Fingerprint over 0.35 {temp}")
    save_dataset(df[include], 'dataset/processed_data/matched_molecular_pairs_smaller.csv')

def get_smaller_mmp_random(size): 
    df = pd.read_csv('dataset/processed_data/matched_molecular_pairs.csv')
    save_dataset(df.sample(size), 'dataset/processed_data/matched_molecular_pairs_smaller.csv')

################ ENSEMBLE DATA ####################################
# Compares two dataframes and get duplicates.
def get_duplicates(source_df, target_df):
    duplicates = []
    for _, row_src in source_df.iterrows():
        for _, row_tgt in target_df.iterrows():
            if (row_src['parent_smiles'] == row_tgt['parent_smiles'] and row_src['metabolite_smiles'] == row_tgt['metabolite_smiles']):
                duplicates.append(list(row_src))
                break 
    df = pd.DataFrame(duplicates)
    df.columns = ["parent_name", "parent_smiles", "metabolite_name", "metabolite_smiles", "enzymes", "reaction_types", "fingerprint_similarity"]
    return df

def divide_duplicates(): 
    metxbiodb = pd.read_csv('dataset/train/metxbiodb_smiles.csv')
    drugbank = pd.read_csv('dataset/train/drugbank_smiles.csv')
    duplicates = get_duplicates(metxbiodb, drugbank)
    
    splitter = GroupShuffleSplit(test_size=0.5, random_state=42)
    split = splitter.split(duplicates, groups=duplicates['parent_smiles'])
    metxbiodb_inds, drugbank_inds = next(split)

    metxbiodb_dups = duplicates.iloc[metxbiodb_inds] 
    
    drugbank_dups = duplicates.iloc[drugbank_inds] 

    metxbiodb = compare_dataset_remove_duplicates(metxbiodb, duplicates)
    drugbank = compare_dataset_remove_duplicates(drugbank, duplicates)

    metxbiodb = combine_datasets(metxbiodb, metxbiodb_dups)
    drugbank = combine_datasets(drugbank, drugbank_dups)

    return metxbiodb, drugbank

def main():
  
    #folders
    train_folder = "dataset/train"
    test_folder = "dataset/test"
    processed_folder ="dataset/processed_data"
    
    # run configuration 
    # NOTE: Configure these flags when running this script. 
    # It is recommended to run the preprocessing on their own so you can monitor the results sooner.
    prepare_train_test_data = True
    prepare_mmp = False
    make_smaller_mmp = False
    prepare_augmented_train_data = False
    prepare_ensemble_data = False
    prepare_unique_data = False
            
    if(prepare_train_test_data):
        # prepare GloryX, MetxBio, Drugbank

        _ ,gloryx_combined_geneneration_df, gloryx_first_generation = prepare_gloryx()
        drugbank_df = prepare_drugbank(fps_cutoff_upper = 1,fps_cutoff_lower = 0.15)
        metxbiodb_df = prepare_metxbiodb(fps_cutoff_upper = 1)    
        # Remove data from metx and drugbank that is already in gloryx
        metxbiodb_df = compare_dataset_remove_duplicates(metxbiodb_df, gloryx_first_generation)
        drugbank_df =  compare_dataset_remove_duplicates(drugbank_df, gloryx_first_generation)

        # Combine the datasets into one metabolic set
        metxbiodb_drugbank_df = combine_datasets(drugbank_df, metxbiodb_df)
        
        # We split this set into train and test data. This based on parent name.  
        ## NOTE: DEFAULT RANDOM STATE IS 42 TO MAKE SURE WE USE SAME SPLIT
        splitter = GroupShuffleSplit(test_size=0.1, n_splits=2, random_state=42)
        split = splitter.split(metxbiodb_drugbank_df, groups=metxbiodb_drugbank_df['parent_smiles'])
        train_inds, test_inds = next(split)
        metabolic_train_df = metxbiodb_drugbank_df.iloc[train_inds]
        metabolic_test_df = metxbiodb_drugbank_df.iloc[test_inds]
        
        save_dataset(metabolic_train_df,f'{train_folder}/metabolic_smiles.csv')
        save_dataset(metabolic_test_df,f'{test_folder}/metabolic_smiles.csv')
        
        # remove data shared with test and metxbiodb and drugbank
        metxbiodb_df = compare_dataset_remove_duplicates(metxbiodb_df, metabolic_test_df)
        drugbank_df = compare_dataset_remove_duplicates(drugbank_df, metabolic_test_df)

        save_dataset(metxbiodb_df,f'{train_folder}/metxbiodb_smiles.csv')
        save_dataset(drugbank_df,f'{train_folder}/drugbank_smiles.csv')

        if(prepare_unique_data):
            metxbiodb_drugbank_unique_df = metxbiodb_drugbank_df.drop_duplicates(subset=['parent_smiles'], keep='first')
            # We split this set into train and test data. This based on parent name.  
            ## NOTE: DEFAULT RANDOM STATE IS 42 TO MAKE SURE WE USE SAME SPLIT
            splitter = GroupShuffleSplit(test_size=0.1, n_splits=2, random_state=42)
            split = splitter.split(metxbiodb_drugbank_unique_df, groups=metxbiodb_drugbank_unique_df['parent_smiles'])
            train_inds, test_inds = next(split)
            metabolic_train_unique_df = metxbiodb_drugbank_unique_df.iloc[train_inds]
            metabolic_test_unique_df = metxbiodb_drugbank_unique_df.iloc[test_inds]
            
            save_dataset(metabolic_train_unique_df,f'{processed_folder}/metabolic_unique_train_smiles.csv')
            save_dataset(metabolic_test_unique_df,f'{processed_folder}/metabolic_unique_test_smiles.csv')

    if(prepare_mmp):
        prepare_matched_molecular_pair()

    if(make_smaller_mmp): 
        get_smaller_mmp_random(500000)

    if(prepare_augmented_train_data):
        augment_data_add_to_new_dataset(f'{train_folder}/metabolic_smiles.csv')

    if(prepare_ensemble_data):
        new_metx, new_drugb = divide_duplicates()
        save_dataset(new_metx,f'{train_folder}/metxbiodb_smiles_dups.csv')
        save_dataset(new_drugb,f'{train_folder}/drugbank_smiles_dups.csv')

    
if __name__ == "__main__":    
   main()


