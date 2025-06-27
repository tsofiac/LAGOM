import pandas as pd
import numpy as np


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

    metabolite = chembldb[np.where(chembldb[:, 0] == chembl_id_metabolite)]
    if len(metabolite) == 0:
        return None
    metabolite_smiles = metabolite[0, 1]
    return [v_analog_id, v_analog, chembl_id_metabolite, metabolite_smiles]


def create_datapoints(VAs, chembldb):
    # VAs: VIRTUAL_ANALOG_IDX,VIRTUAL_ANALOG,CHEMBL_COMPOUND_ID,CHEMBL_TARGET_IDs
    # chembld: chembl_id	canonical_smiles	standard_inchi	standard_inchi_key
    result = [create_datapoint(chembldb, row) for row in VAs]
    result = list(filter(lambda item: item is not None, result))

    return result


def create_VA(VA_df, chembldb, output_file):
    total = []

    total = create_datapoints(VA_df.to_numpy(), chembldb.to_numpy())
    VA_smiles = pd.DataFrame(total)
    VA_smiles.columns = [
        "virtual_analog_id",
        "parent_smiles",
        "chembl_id_metabolite",
        "child_smiles",
    ]
    VA_smiles.to_csv(output_file, index=False)
    return VA_smiles


def main():
    # # Division into 10 parts

    # start_row = 0
    # end_row = 1101304

    # start_row = 1101304
    # end_row = 2202607

    # start_row = 2202607
    # end_row = 3303912

    # start_row = 3303912
    # end_row = 4405216

    # start_row = 4405216
    # end_row = 5506519

    # start_row = 5506519
    # end_row = 6607824

    # start_row = 6607824
    # end_row = 7709127

    # start_row = 7709127
    # end_row = 8810431

    # start_row = 8810431
    # end_row = 9911735

    start_row = 9911735
    end_row = None  # 11013037

    print("rows:", start_row, "to", end_row)

    molecular_match_pairs_filepath = "dataset/raw_data/1297204_Virtual_Analogs.dat"
    chembldb_filepath = "dataset/raw_data/chembl_35_chemreps.txt"  # CHEMBL35
    output_file = f"dataset/curated_data/VA_rows_{start_row}_to_{end_row}.csv"

    VA_df = pd.read_csv(molecular_match_pairs_filepath)
    chembldb = get_chembl(chembldb_filepath)

    VA_df_subset = VA_df.iloc[start_row:end_row]

    create_VA(VA_df_subset, chembldb, output_file)


if __name__ == "__main__":
    main()
