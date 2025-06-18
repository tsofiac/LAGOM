#!/usr/bin/env python3
import pandas as pd


def get_unique_parents(input_file):
    # For curated_data file
    df = pd.read_csv(input_file)

    # For finetune file
    # df = pd.read_csv(input_file, sep='\t')

    if "parent_smiles" in df.columns:
        unique_parents = df["parent_smiles"].unique()

    elif "reactants" in df.columns:
        unique_parents = df["reactants"].unique()

    return len(unique_parents)


if __name__ == "__main__":
    # input_file = 'dataset/raw_data/mmp_new_split_finetune.csv'
    input_file = "dataset/curated_data/combined_smiles_clean.csv"

    unique = get_unique_parents(input_file)

    print(unique)
