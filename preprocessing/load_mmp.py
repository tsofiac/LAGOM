import pandas as pd
import numpy as np


##########################LOAD###########################

def get_chembl(filepath): 
    file = open(filepath, "r")
    content = file.readlines()
    file.close()

    result = [line.strip("\n").split("\t") for line in content]
    columns = result[0]
    result = result[1:]
    df = pd.DataFrame(result)
    df.columns = columns

    return df

def create_datapoint(chembldb, row): 
    v_analog_id = row[0] 
    v_analog = row[1]
    chembl_id_metabolite = row[2]

    metabolite = chembldb[np.where(chembldb[:,0] == chembl_id_metabolite)]
    if len(metabolite) == 0 : 
        return None
    metabolite_smiles = metabolite[0,1]
    return [v_analog_id, v_analog, chembl_id_metabolite, metabolite_smiles]

def create_datapoints(mmps, chembldb): 
    # mmps: VIRTUAL_ANALOG_IDX,VIRTUAL_ANALOG,CHEMBL_COMPOUND_ID,CHEMBL_TARGET_IDs
    # chembld: chembl_id	canonical_smiles	standard_inchi	standard_inchi_key 
    result = [create_datapoint(chembldb, row) for row in mmps]
    result = list(filter(lambda item: item is not None, result))
    return result

def create_mmp(mmp_df, chembldb, output_file): 
    total = []
    
    total = create_datapoints(mmp_df.to_numpy(), chembldb.to_numpy())
    mmp_smiles = pd.DataFrame(total)
    mmp_smiles.columns = ["virtual_analog_id", "parent_smiles", "chembl_id_metabolite", "metabolite_smiles"]
    mmp_smiles.to_csv(output_file, index=False)
    return mmp_smiles

'''
mmp_split = np.array_split(mmp_df.to_numpy(), 20)
for current in list(range(12,20)): 
    print(f"Create_dataset of {current}")
    current_mmp = mmp_split[current]
    result = create_datapoints(current_mmp, chembldb.to_numpy())
    total = total + result

    mmp_smiles = pd.DataFrame(total)
    mmp_smiles.columns = ["virtual_analog_id", "parent_smiles", "chembl_id_metabolite", "metabolite_smiles"]
    output_file = "dataset/raw_data/all_paired_mmp" + str(current) + ".csv"
    mmp_smiles.to_csv(output_file, index=False)
'''


def main():    
    molecular_match_pairs_filepath = "dataset/raw_data/1297204_Virtual_Analogs.dat"
    chembldb_filepath = "dataset/raw_data/chembl_35_chemreps.txt" # CHEMBL35
    output_file = "dataset/curated_data/paired_mmp.csv"

    mmp_df = pd.read_csv(molecular_match_pairs_filepath)
    chembldb = get_chembl(chembldb_filepath)
    create_mmp(mmp_df, chembldb, output_file)
    
if __name__ == "__main__":
    main()





############### PRETRAINING DATA ##################################
# def prepare_matched_molecular_pair(): 
#     matched_molecular_pair_df = pd.read_csv("dataset/raw_data/paired_mmp.csv")
#     filtered_data = curate_data(matched_molecular_pair_df)
#     save_dataset(filtered_data, 'dataset/train/matched_molecular_pairs.csv')
#     return filtered_data

# def get_smaller_mmp(): 
#     df = pd.read_csv('dataset/processed_data/matched_molecular_pairs.csv')
#     include = [sim > 0.35 for sim in fingerprint_similarity_single(df["parent_smiles"], df["metabolite_smiles"])]
#     temp = include.count(True)
#     print(f"Fingerprint over 0.35 {temp}")
#     save_dataset(df[include], 'dataset/processed_data/matched_molecular_pairs_smaller.csv')

# def get_smaller_mmp_random(size): 
#     df = pd.read_csv('dataset/processed_data/matched_molecular_pairs.csv')
#     save_dataset(df.sample(size), 'dataset/processed_data/matched_molecular_pairs_smaller.csv')