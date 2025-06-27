import json
import pandas as pd
from standardize_smiles import standardize_smiles_collection


def load_gloryx(json_file, comma_file, comma_unique_file, tab_file):
    with open(json_file, "r", encoding="utf8") as file:
        reader = file.read()
        # replace \ to \\ for correct reading of file
        reader = reader.replace("\\", "\\\\")
        data = json.loads(reader)

    child_smiles = []
    parent_smiles = []
    child_name = []
    parent_name = []
    generation_column = []
    set_column = []

    def process_metabolites(metabolites, drug_smiles, drug_name):
        for metabolite in metabolites:
            metabolite_smiles = metabolite["smiles"]
            metabolite_name = metabolite.get("metaboliteName", None)
            generation = metabolite.get("generation", None)

            child_smiles.append(metabolite_smiles)
            parent_smiles.append(drug_smiles)
            child_name.append(metabolite_name)
            parent_name.append(drug_name)
            generation_column.append(generation)
            set_column.append("test")

            if "metabolites" in metabolite:
                process_metabolites(
                    metabolite["metabolites"], metabolite_smiles, metabolite_name
                )

    for entry in data:
        drug_name = entry["drugName"]
        drug_smiles = entry["smiles"]

        if "metabolites" in entry:
            process_metabolites(entry["metabolites"], drug_smiles, drug_name)

    df = pd.DataFrame(
        {
            "parent_name": parent_name,
            "child_name": child_name,
            "parent_smiles": parent_smiles,
            "child_smiles": child_smiles,
            "generation": generation_column,
            "set": set_column,
        }
    )

    df["source"] = "GLORYx"

    df["parent_smiles"] = standardize_smiles_collection(df["parent_smiles"], False)
    df["child_smiles"] = standardize_smiles_collection(df["child_smiles"], False)

    df = df[df["generation"] == 1]  # Only first generation reactions
    df.to_csv(comma_file, index=False)
    df = df.drop_duplicates(subset="parent_smiles", keep="first")  # Only unique parent
    df.to_csv(comma_unique_file, index=False)

    df = df.rename(
        columns={
            "child_smiles": "products",
            "parent_smiles": "reactants",
        }
    )

    df.to_csv(tab_file, sep="\t", index=False)


if __name__ == "__main__":
    gloryx = "dataset/raw_data/gloryx_test_dataset.json"
    gloryx_comma = "dataset/curated_data/gloryx_smiles_clean.csv"
    gloryx_unique = "dataset/curated_data/gloryx_unique_smiles_clean.csv"
    gloryx_tab = "dataset/finetune/gloryx_finetune.csv"

    load_gloryx(gloryx, gloryx_comma, gloryx_unique, gloryx_tab)
