#!/usr/bin/env python3
import pandas as pd

def get_unique_parents(input_file):

    # For curated_data file
    # df = pd.read_csv(input_file)

    # For finetune file
    df = pd.read_csv(input_file, sep='\t')

    # df = df[df['source'].isin(['DrugBank'])]

    # unique_parent_dataset = []

    # print(df.columns)

    if 'parent_smiles' in df.columns:

        unique_parents = df['parent_smiles'].unique()

        # for parent in df['parent_smiles'].unique():
        #     row = df[df['parent_smiles'] == parent].iloc[0]
        #     unique_parent_dataset.append(row)
    elif 'reactants' in df.columns:

        unique_parents = df['reactants'].unique()
        # for parent in df['reactants'].unique():
        #     row = df[df['reactants'] == parent].iloc[0]
        #     unique_parent_dataset.append(row)

    # unique_parent_df = pd.DataFrame(unique_parent_dataset).reset_index(drop=True)
    # unique_parent_df = unique_parent_df[df.columns]
    # print('here 2')

    return len(unique_parents)


if __name__ == "__main__":
    input_file = 'dataset/raw_data/mmp_new_split_finetune.csv'
    # input_file = 'dataset/curated_data/gloryx_smiles_clean.csv'

    unique = get_unique_parents(input_file)

    print(unique)