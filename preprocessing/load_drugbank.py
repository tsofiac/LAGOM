import xml.etree.ElementTree as ET
import pandas as pd
from rdkit import Chem


def sdf_to_csv(sdf_file, csv_file, is_drug_information=False):
    name_property = "GENERIC_NAME" if is_drug_information else "NAME"
    data = []
    suppl = Chem.SDMolSupplier(sdf_file)
    for mol in suppl:
        if mol is None:
            continue

        drugbank_id = (
            mol.GetProp("DRUGBANK_ID") if mol.HasProp("DRUGBANK_ID") else "NaN"
        )
        inchi = Chem.inchi.MolToInchi(mol) if mol else "NaN"
        smiles = Chem.MolToSmiles(mol) if mol else "NaN"
        name = mol.GetProp(name_property) if mol.HasProp(name_property) else "NaN"

        data.append(
            {"dbid": drugbank_id, "name": name, "smiles": smiles, "inchi": inchi}
        )

    df = pd.DataFrame(data, columns=["dbid", "name", "smiles", "inchi"])
    df.to_csv(csv_file, index=False)


def reformat_external_csv(input_file, output_file):
    db_external_df = pd.read_csv(input_file)
    db_external_df = db_external_df.rename(
        columns={
            "DrugBank ID": "dbid",
            "Name": "name",
            "SMILES": "smiles",
            "InChI": "inchi",
        }
    )
    db_external_df = db_external_df[["dbid", "name", "smiles", "inchi"]]
    db_external_df.to_csv(output_file, index=False)


def remove_bad_drug_metabolite_rows(input_file):
    df = pd.read_csv(input_file)
    initial_row_count = len(df)
    df.dropna(subset=["dbid", "name"], how="all", inplace=True)
    df.dropna(subset=["smiles", "inchi"], how="all", inplace=True)
    final_row_count = len(df)
    dropped_count = initial_row_count - final_row_count
    print("Missing data: ", dropped_count)
    df.to_csv(input_file, index=False)


def count_nan_smiles(csv_file):
    df = pd.read_csv(csv_file)
    nan_count = df["smiles"].isna().sum()
    print("Missing SMILES: ", nan_count)


def inchi_to_smiles(csv_file):
    df = pd.read_csv(csv_file)
    for index, row in df[df["smiles"].isna()].iterrows():
        try:
            mol = Chem.inchi.MolFromInchi(row["inchi"])
            if mol is not None:
                smiles = Chem.MolToSmiles(mol)
                df.at[index, "smiles"] = smiles
            else:
                print(
                    f"Conversion failed for entry {row['dbid']} with InChI: {row['inchi']}"
                )

        except Exception as e:
            print(f"ERROR processing row {index}, DBID: {row['dbid']}, error: {e}")
    df.to_csv(csv_file, index=False)


def combine_datasets(first_file, second_file, output_file):
    first_df = pd.read_csv(first_file, on_bad_lines="skip")
    second_df = pd.read_csv(second_file, on_bad_lines="skip")

    first_df.drop(columns=["inchi"], inplace=True, errors="ignore")
    second_df.drop(columns=["inchi"], inplace=True, errors="ignore")

    output_df = pd.concat([first_df, second_df], ignore_index=True)

    output_df.to_csv(output_file, index=False)


def remove_duplicates(csv_file):
    df = pd.read_csv(csv_file)
    initial_row_count = len(df)

    df = df[
        df["dbid"].isnull()
        | ~df[df["dbid"].notnull()].duplicated(subset="dbid", keep="first")
    ]
    df = df[
        df["name"].isnull()
        | ~df[df["name"].notnull()].duplicated(subset="name", keep="first")
    ]
    final_row_count = len(df)
    dropped_count = initial_row_count - final_row_count
    print("Number of duplicates removed: ", dropped_count)

    df.to_csv(csv_file, index=False)


class Reaction:
    def __init__(self):
        self.left_side_name = ""
        self.left_side_id = ""
        self.right_side_name = ""
        self.right_side_id = ""
        self.enzymes = []

    def has_data(self):
        return any([getattr(self, attr) for attr in vars(self)])


def generate_reaction_pairs(input_file, output_file):
    context = ET.iterparse(input_file, events=("start", "end"))
    parent_child_ids = []
    count = 0

    local_reaction = Reaction()
    is_in_local_reaction = False
    is_in_enzymes = False
    right_or_left_reaction = "left"
    for i, (event, element) in enumerate(context):
        if i % 1000000 == 0:
            print(f"iteration {i}")

        # Split the tag to remove garbage
        if element.tag == "root":
            tag = "root"
        else:
            tag = element.tag.split("}")[1]

        if event == "start":
            # We have entered a local reaction.
            if tag == "reaction":
                is_in_local_reaction = True
                # If local reaction has data then it is complete or it needs to be investigated
                if local_reaction.has_data():
                    # We need at least an id or a name on both sides of the reaction
                    # To map data
                    if (
                        local_reaction.left_side_id == ""
                        and local_reaction.left_side_name == ""
                    ) or (
                        local_reaction.right_side_id == ""
                        and local_reaction.right_side_name == ""
                    ):
                        count += 1
                    # If all data exists then we add it to the list
                    else:
                        # We join the list of enzymes to single string
                        enzyme_string = ", ".join(local_reaction.enzymes)

                        parent_child_ids.append(
                            [
                                local_reaction.left_side_id,
                                local_reaction.left_side_name,
                                local_reaction.right_side_id,
                                local_reaction.right_side_name,
                                enzyme_string,
                            ]
                        )
                    # We reset here since we want to look for the new reaction
                    local_reaction = Reaction()

            # To keep track on what data we are looking at
            if tag == "left-element" and is_in_local_reaction:
                right_or_left_reaction = "left"
            if tag == "right-element" and is_in_local_reaction:
                right_or_left_reaction = "right"
            if tag == "drugbank-id" and is_in_local_reaction:
                if right_or_left_reaction == "left":
                    local_reaction.left_side_id = element.text
                elif right_or_left_reaction == "right":
                    local_reaction.right_side_id = element.text
            if tag == "name" and is_in_local_reaction:
                if right_or_left_reaction == "left":
                    local_reaction.left_side_name = element.text
                elif right_or_left_reaction == "right":
                    local_reaction.right_side_name = element.text
            # Get the names of the enzymes
            if tag == "enzyme" and is_in_local_reaction:
                right_or_left_reaction = "none"
                is_in_enzymes = True
            if is_in_enzymes:
                if tag == "name":
                    if element.text is not None:
                        local_reaction.enzymes.append(element.text)
                    is_in_enzymes = False
            # We know its the end of the reaction if this tag shows
            if tag == "snp-effects":
                is_in_enzymes = False
                is_in_local_reaction = False

        elif event == "end":
            # Clean up the element when the end tag is encountered
            element.clear()

    print("Number of rows with missing data:", count)
    df = pd.DataFrame(
        parent_child_ids,
        columns=["parent_id", "parent_name", "child_id", "child_name", "enzymes"],
    )
    df.to_csv(output_file, index=False)


def add_missing_ID(data_file):
    df = pd.read_csv(data_file)

    missing_parent_id_rows = df[df["parent_id"].isnull()]
    missing_child_id_rows = df[df["child_id"].isnull()]

    df_no_missing_id = df[df["parent_id"].notnull() & df["child_id"].notnull()]

    for index, missing_row in missing_parent_id_rows.iterrows():
        parent_name = missing_row["parent_name"]

        matching_rows = df_no_missing_id[df_no_missing_id["parent_name"] == parent_name]

        if not matching_rows.empty:
            df.at[index, "parent_id"] = matching_rows.iloc[0]["parent_id"]
        else:
            matching_rows = df_no_missing_id[
                df_no_missing_id["child_name"] == parent_name
            ]

            if not matching_rows.empty:
                df.at[index, "parent_id"] = matching_rows.iloc[0]["child_id"]

    for index, missing_row in missing_child_id_rows.iterrows():
        child_name = missing_row["child_name"]

        matching_rows = df_no_missing_id[df_no_missing_id["child_name"] == child_name]
        if not matching_rows.empty:
            df.at[index, "child_id"] = matching_rows.iloc[0]["child_id"]
        else:
            matching_rows = df_no_missing_id[
                df_no_missing_id["parent_name"] == child_name
            ]
            if not matching_rows.empty:
                df.at[index, "child_id"] = matching_rows.iloc[0]["parent_id"]

    df.to_csv(data_file, index=False)


def find_drug_origin_drugbank(data_file):
    original_df = pd.read_csv(data_file)
    original_df["origin"] = "unknown"
    only_db_parents_df = original_df[
        ~original_df["parent_id"].str.contains("DBMET").fillna(False)
    ]  # All parents that not contain DBMET
    only_dbmet_parents_df = original_df[
        original_df["parent_id"].str.contains("DBMET").fillna(False)
    ]  # All parents that contain DBMET

    only_db_parents_df.loc[:, "origin"] = only_db_parents_df["parent_id"]

    both_db_rows = only_db_parents_df[
        ~only_db_parents_df["child_id"].str.contains("DBMET").fillna(False)
    ]
    db_and_dbmet_rows = only_db_parents_df[
        only_db_parents_df["child_id"].str.contains("DBMET").fillna(False)
    ]

    # DB-DB fix ----------------------------------------
    db_child_ids = both_db_rows["child_id"].tolist()

    both_db_child_also_parent = both_db_rows[
        both_db_rows["parent_id"].isin(db_child_ids)
    ]

    for index, row in both_db_child_also_parent.iterrows():
        origin_parent_id = both_db_rows[both_db_rows["child_id"] == row["parent_id"]][
            "origin"
        ]
        if not origin_parent_id.empty:
            both_db_child_also_parent.at[index, "origin"] = origin_parent_id.values[0]

    both_db_rows.update(both_db_child_also_parent)

    # DB-DBMET fix ----------------------------------------
    db_child_ids = both_db_rows["child_id"].tolist()
    db_parent_with_other_origin = db_and_dbmet_rows[
        db_and_dbmet_rows["parent_id"].isin(db_child_ids)
    ]

    for index, row in db_parent_with_other_origin.iterrows():
        origin_parent_id = both_db_rows[both_db_rows["child_id"] == row["parent_id"]][
            "origin"
        ]
        if not origin_parent_id.empty:
            db_parent_with_other_origin.at[index, "origin"] = origin_parent_id.values[0]

    db_and_dbmet_rows = db_and_dbmet_rows[
        ~db_and_dbmet_rows["parent_id"].isin(db_child_ids)
    ]

    only_db_parents_df = pd.concat([both_db_rows, db_parent_with_other_origin])
    only_db_parents_df = pd.concat([only_db_parents_df, db_and_dbmet_rows])

    # DBMET-DBMET fix ----------------------------------------
    num_new_rows = 1
    while num_new_rows != 0:
        child_ids = only_db_parents_df["child_id"].tolist()
        # Store rows where parents are child in list
        metabolite_parent_with_drug_origin_df = only_dbmet_parents_df[
            only_dbmet_parents_df["parent_id"].isin(child_ids)
        ]

        # Assign 'origin' for metabolite_parent_with_drug_origin_df
        for index, row in metabolite_parent_with_drug_origin_df.iterrows():
            origin_parent_id = only_db_parents_df[
                only_db_parents_df["child_id"] == row["parent_id"]
            ]["origin"]
            if not origin_parent_id.empty:
                metabolite_parent_with_drug_origin_df.at[index, "origin"] = (
                    origin_parent_id.values[0]
                )

        # Drop all children that are not children of drugs or drug metabolites
        only_dbmet_parents_df = only_dbmet_parents_df[
            ~only_dbmet_parents_df["parent_id"].isin(child_ids)
        ]

        only_db_parents_df = pd.concat(
            [only_db_parents_df, metabolite_parent_with_drug_origin_df]
        )
        num_new_rows = len(metabolite_parent_with_drug_origin_df)

    combined_data = pd.concat([only_db_parents_df, only_dbmet_parents_df])

    # Make sure there are no NaN values in 'origin'
    combined_data["origin"] = combined_data["origin"].fillna("unknown")

    # Give all 'unknown' origins a unique id
    unknown_counter = 0
    for index, row in combined_data.iterrows():
        if row["origin"] == "unknown":
            unknown_counter += 1
            combined_data.at[index, "origin"] = f"unknown{unknown_counter}"

    print("Number of reactions with unknown origin: ", unknown_counter)
    combined_data.to_csv(data_file, index=False)


def map_smiles_to_reaction_pairs(
    reaction_pairs_file, structures_file, output_file, is_drug=False
):
    smiles_property = "parent_smiles" if is_drug else "child_smiles"
    id_property = "parent_id" if is_drug else "child_id"
    name_property = "parent_name" if is_drug else "child_name"

    reaction_pairs_df = pd.read_csv(reaction_pairs_file)
    structures_df = pd.read_csv(structures_file)

    reaction_pairs_df[smiles_property] = None

    # Fill the smiles property based on id_property and name_property
    for index, row in reaction_pairs_df.iterrows():
        # Check for id_property match
        if pd.notna(row[id_property]):
            matched_row = structures_df[structures_df["dbid"] == row[id_property]]
            if not matched_row.empty:
                reaction_pairs_df.at[index, smiles_property] = matched_row[
                    "smiles"
                ].values[0]

        # Check for name_property match if smiles not found yet
        elif pd.notna(row[name_property]):
            matched_row = structures_df[structures_df["name"] == row[name_property]]
            if not matched_row.empty:
                reaction_pairs_df.at[index, smiles_property] = matched_row[
                    "smiles"
                ].values[0]

    reaction_pairs_df.to_csv(output_file, index=False)


def clean_smiles(input_file):
    df = pd.read_csv(input_file)
    initial_row_count = len(df)
    # Remove rows where either 'parent_smiles' or 'child_smiles' is missing
    df.dropna(subset=["parent_smiles", "child_smiles"], how="any", inplace=True)
    final_row_count = len(df)
    dropped_count = initial_row_count - final_row_count
    print("Number of dropped rows: ", dropped_count)

    df["source"] = "DrugBank"
    df.to_csv(input_file, index=False)


if __name__ == "__main__":
    get_drug_structures = True
    get_metabolite_structures = True
    get_external_structures = True
    get_reaction_pairs = True
    combine_all_structures = True
    extend_dataset = True

    ## Drug structure - drugbank_drug_structures.sdf
    drugbank_drug_structures = "dataset/raw_data/drugbank_drug_structures.sdf"
    parsed_drug_structures = "dataset/extracted_data/drugbank_drug_structures.csv"
    if get_drug_structures:
        sdf_to_csv(
            drugbank_drug_structures, parsed_drug_structures, is_drug_information=True
        )
        remove_bad_drug_metabolite_rows(parsed_drug_structures)
        count_nan_smiles(parsed_drug_structures)
        # Missing data:: 0
        # Missing SMILES: 0

    ## Metabolite structure - drugbank_metabolite_structures.sdf
    drugbank_metabolite_structures = (
        "dataset/raw_data/drugbank_metabolite_structures.sdf"
    )
    parsed_metabolite_structures = (
        "dataset/extracted_data/drugbank_metabolite_structures.csv"
    )
    if get_metabolite_structures:
        sdf_to_csv(
            drugbank_metabolite_structures,
            parsed_metabolite_structures,
            is_drug_information=False,
        )
        remove_bad_drug_metabolite_rows(parsed_metabolite_structures)
        count_nan_smiles(parsed_metabolite_structures)
        # Missing data:: 0
        # Missing SMILES: 0

    ## External drugs - drugbank_external_structures.csv
    drugbank_external_structures = "dataset/raw_data/drugbank_external_structures.csv"
    parsed_external_stuctures = (
        "dataset/extracted_data/drugbank_external_structures_cleaned.csv"
    )
    if get_external_structures:
        reformat_external_csv(drugbank_external_structures, parsed_external_stuctures)
        remove_bad_drug_metabolite_rows(parsed_external_stuctures)
        count_nan_smiles(parsed_external_stuctures)
        inchi_to_smiles(parsed_external_stuctures)
        count_nan_smiles(parsed_external_stuctures)
        # Missing data:: 850
        # Missing SMILES: 3
        # Missing SMILES: 0

    # Get the reaction pairs - drugbank_full_database.xml
    drugbank_full_database = "dataset/raw_data/drugbank_full_database.xml"
    drugbank_reaction_pairs = "dataset/extracted_data/drugbank_reaction_pairs.csv"
    if get_reaction_pairs:
        generate_reaction_pairs(drugbank_full_database, drugbank_reaction_pairs)
        add_missing_ID(drugbank_reaction_pairs)
        find_drug_origin_drugbank(drugbank_reaction_pairs)
        # Missing data: 45
        # Number of reactions with unknown origin:  0

    # Combine all strucutres into one file
    drugbank_full_structures = "dataset/extracted_data/drugbank_full_structures.csv"
    if combine_all_structures:
        combine_datasets(
            parsed_drug_structures,
            parsed_metabolite_structures,
            drugbank_full_structures,
        )
        combine_datasets(
            drugbank_full_structures,
            parsed_external_stuctures,
            drugbank_full_structures,
        )
        remove_duplicates(drugbank_full_structures)
        # Number of duplicates removed:  12392

    # Extend reaction pairs with SMILES
    mapped_smiles_to_reactions = "dataset/extracted_data/drugbank_smiles.csv"
    if extend_dataset:
        map_smiles_to_reaction_pairs(
            drugbank_reaction_pairs,
            drugbank_full_structures,
            mapped_smiles_to_reactions,
            True,
        )
        map_smiles_to_reaction_pairs(
            mapped_smiles_to_reactions,
            drugbank_full_structures,
            mapped_smiles_to_reactions,
            False,
        )
        clean_smiles(mapped_smiles_to_reactions)
        # Missing data:  556
