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

    print("You now have a df")

    return df

def create_datapoint(chembldb, row): 
    v_analog_id = row[0] 
    v_analog = row[1]
    chembl_id_metabolite = row[2]

    # if v_analog_id % 10000 == 0:
    #     print(f'v_analog_id {v_analog_id}')

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

# def create_datapoints(mmps, chembldb):
#     result = []
#     for row in mmps:
#         datapoint = create_datapoint(chembldb, row)
#         # if datapoint == "STOP":
#         #     break  # Exit loop when the signal is received
#         if datapoint is not None:
#             result.append(datapoint)
#     return result

def create_mmp(mmp_df, chembldb, output_file): 
    total = []
    
    total = create_datapoints(mmp_df.to_numpy(), chembldb.to_numpy())
    mmp_smiles = pd.DataFrame(total)
    mmp_smiles.columns = ["virtual_analog_id", "parent_smiles", "chembl_id_metabolite", "child_smiles"]
    mmp_smiles.to_csv(output_file, index=False)
    return mmp_smiles



def main():    

    # # Division into 10 parts
    # section_size = 1101304

    # start_row = 10
    # end_row = 15

    # start_row = 0
    # end_row = 1101304

    # start_row = 1101305
    # end_row = 2202607

    # start_row = 2202608
    # end_row = 3303912

    # start_row = 3303913
    # end_row = 4405216 

    # start_row = 4405217
    # end_row = 5506519

    start_row = 5506520
    end_row = 6607824

    # start_row = 6607825
    # end_row = 7709127

    # start_row = 7709128
    # end_row = 8810431

    # start_row = 8810432
    # end_row = 9911735

    # start_row = 9911736
    # end_row = 11013037

    print("rows:", start_row, "to", end_row)

    molecular_match_pairs_filepath = "dataset/raw_data/1297204_Virtual_Analogs.dat"
    chembldb_filepath = "dataset/raw_data/chembl_35_chemreps.txt" # CHEMBL35
    output_file = f"dataset/curated_data/paired_mmp_rows_{start_row}_to_{end_row}.csv"

    mmp_df = pd.read_csv(molecular_match_pairs_filepath)
    chembldb = get_chembl(chembldb_filepath)

    print("length", len(mmp_df))

    mmp_df_subset = mmp_df.iloc[start_row:end_row]

    create_mmp(mmp_df_subset, chembldb, output_file)

    # mmp_split = np.array_split(mmp_df.to_numpy(), 40)
    # for current in list(range(12,20)): 
    #     print(f"Create_dataset of {current}")
    #     current_mmp = mmp_split[current]
    #     result = create_datapoints(current_mmp, chembldb.to_numpy())
    #     total = total + result

    #     mmp_smiles = pd.DataFrame(total)
    #     mmp_smiles.columns = ["virtual_analog_id", "parent_smiles", "chembl_id_metabolite", "child_smiles"]
    #     output_file = "dataset/curated_data/paired_mmp" + str(current) + ".csv"
    #     mmp_smiles.to_csv(output_file, index=False)
        
if __name__ == "__main__":
    main()

    # # Define start and end rows
    # start_row = 4382477
    # end_row = 5408311

    # df = pd.read_csv('dataset/curated_data/paired_mmp_all.csv')

    # df = df.iloc[start_row:end_row]
    # # Save the selected rows to a new CSV file
    # df.to_csv(f'dataset/curated_data/paired_mmp_rows_{start_row}_to_{end_row}.csv', index=False)
    # # data_df.to_csv('dataset/curated_data/paired_mmp_rows_5506520_to_6607824.csv', index=False)




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