import pandas as pd
import numpy as np
from collections import Counter

import sys
import os
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))
from preprocessing.NEW_preprocessing_all import add_possible_products

def split_metabolites_4(input_csv, split1_csv, split2_csv, split3_csv, split4_csv):

    np.random.seed(6)

    df_split1 = []
    df_split2 = []
    df_split3 = []
    df_split4 = []

    # Read the CSV file into a DataFrame
    df = pd.read_csv(input_csv)

    # Filter the DataFrame to include only rows where the 'set'-column is 'train'
    df_train = df[df['set'] == 'train']
    df_val = df[df['set'] == 'val']
    df_test = df[df['set'] == 'test']

    # Group by 'parent_smiles' and get all rows for each unique parent
    for parent, parent_rows in df_train.groupby('parent_smiles'):
        num_rows = len(parent_rows)

        if num_rows == 1:
            df_split1.append(parent_rows.iloc[0])
            df_split2.append(parent_rows.iloc[0])
            df_split3.append(parent_rows.iloc[0])
            df_split4.append(parent_rows.iloc[0])

        elif num_rows == 2:
            df_split1.append(parent_rows.iloc[0])
            df_split2.append(parent_rows.iloc[0])
            df_split3.append(parent_rows.iloc[1])
            df_split4.append(parent_rows.iloc[1])

        elif num_rows == 3:
            df_split1.append(parent_rows.iloc[0])
            df_split2.append(parent_rows.iloc[1])
            df_split3.append(parent_rows.iloc[2])
            random_index = np.random.randint(0, num_rows)
            df_split4.append(parent_rows.iloc[random_index])

        elif num_rows > 3:
            equally_divisble_rows = num_rows - (num_rows % 4)
            parent_rows_list = list(parent_rows.iterrows())
            for i in range(equally_divisble_rows):
                _, row = parent_rows_list[i]
                if i % 4 == 0:
                    df_split1.append(row)
                elif i % 4 == 1:
                    df_split2.append(row)
                elif i % 4 == 2:
                    df_split3.append(row)
                elif i % 4 == 3:
                    df_split4.append(row)
            if num_rows>equally_divisble_rows:
                split_list = [1,2,3,4]
                for i in range(equally_divisble_rows, num_rows):
                    _, row = parent_rows_list[i]
                    if not split_list:
                        split_list = [1,2,3,4]
                    random_split_index = np.random.choice(split_list)
                    split_list.remove(random_split_index)
                    if random_split_index == 1:
                        df_split1.append(row)
                    elif random_split_index == 2:
                        df_split2.append(row)
                    elif random_split_index == 3:
                        df_split3.append(row)
                    elif random_split_index == 4:
                        df_split4.append(row)
            
    # Convert lists to DataFrames
    df_split1 = pd.DataFrame(df_split1)
    df_split2 = pd.DataFrame(df_split2)
    df_split3 = pd.DataFrame(df_split3)
    df_split4 = pd.DataFrame(df_split4)

    # Concatenate df_val and df_train to all splits
    df_split1 = pd.concat([df_split1, df_val], ignore_index=True)
    df_split2 = pd.concat([df_split2, df_val], ignore_index=True)
    df_split3 = pd.concat([df_split3, df_val], ignore_index=True)
    df_split4 = pd.concat([df_split4, df_val], ignore_index=True)

    df_split1 = pd.concat([df_split1, df_test], ignore_index=True)
    df_split2 = pd.concat([df_split2, df_test], ignore_index=True)
    df_split3 = pd.concat([df_split3, df_test], ignore_index=True)
    df_split4 = pd.concat([df_split4, df_test], ignore_index=True)

    # Print lengths of each DataFrame
    print(f'The wholed dataframe contains {len(df)} rows')
    print(f'Split 1 contains {len(df_split1)} rows')
    print(f'Split 2 contains {len(df_split2)} rows')
    print(f'Split 3 contains {len(df_split3)} rows')
    print(f'Split 4 contains {len(df_split4)} rows')

    # Convert to csv
    df_split1.to_csv(split1_csv, index=False)
    df_split2.to_csv(split2_csv, index=False)
    df_split3.to_csv(split3_csv, index=False)
    df_split4.to_csv(split4_csv, index=False)

    # return df_split1, df_split2, df_split3, df_split4

def reformat_for_chemformer(input_file, output_file):
    df = pd.read_csv(input_file)

    df = df.rename(columns={
        "child_smiles": "products",
        "parent_smiles": "reactants",
    })

    df.to_csv(output_file, sep='\t', index=False)

def split_on_tb_clusters(input_file, n_splits, name):
    df = pd.read_csv(input_file)

    df_train = df[df['set'] == 'train']
    df_val = df[df['set'] == 'val']
    df = pd.concat([df_train, df_val])

    clusters = df['TB_Cluster']

    # Count the number of rows in each cluster
    cluster_counts = Counter(clusters)
     # Sort clusters by size (largest first)
    sorted_clusters = sorted(cluster_counts.items(), key=lambda x: x[1], reverse=True)

    # Initialise the splits
    splits = [pd.DataFrame() for _ in range(n_splits)]

    # Distribute clusters to the split with the least number of rows
    for cluster_id, _ in sorted_clusters:
        cluster_rows = df[df['TB_Cluster'] == cluster_id]
        # Find the split with the least number of rows
        split_lengths = [len(split) for split in splits]
        split_index = split_lengths.index(min(split_lengths))
        # Add the cluster rows to that split
        splits[split_index] = pd.concat([splits[split_index], cluster_rows])

    for i, split in enumerate(splits):
        csv = f'dataset/curated_data/tb_split_{i + 1}_of_{n_splits}_{name}.csv'
        finetune_csv = f'dataset/finetune/tb_split_{i + 1}_of_{n_splits}_{name}_finetune.csv'
        split.to_csv(csv, index=False)
        reformat_for_chemformer(csv, finetune_csv)
        add_possible_products(finetune_csv)
        print(f"Split {i + 1} contains {len(split)} rows.")

    
    return splits

if __name__ == "__main__":

    # Define split type. 4 splits will be generated
    split_type = 'children' #'random' 'parents' 'children'

    if split_type == 'parents':
        tb_file = 'dataset/curated_data/tb_output_parent.csv'
        split_on_tb_clusters(tb_file, 4, 'parent')
        # Reformating for Chemformer fine-tuning is done in the function

    elif split_type == 'children':
        tb_file = 'dataset/curated_data/tb_output_child.csv'
        split_on_tb_clusters(tb_file, 4, 'child')
        # Reformating for Chemformer fine-tuning is done in the function

    elif split_type == 'random':

        input_csv = 'dataset/curated_data/combined_smiles_clean.csv'

        split1_csv = 'dataset/curated_data/4_split1_combined_smiles_clean.csv'
        split2_csv = 'dataset/curated_data/4_split2_combined_smiles_clean.csv'
        split3_csv = 'dataset/curated_data/4_split3_combined_smiles_clean.csv'
        split4_csv = 'dataset/curated_data/4_split4_combined_smiles_clean.csv'

        fine_tune_split1 = 'dataset/finetune/4_split1_random_finetune.csv'
        fine_tune_split2 = 'dataset/finetune/4_split2_random_finetune.csv'
        fine_tune_split3 = 'dataset/finetune/4_split3_random_finetune.csv'
        fine_tune_split4 = 'dataset/finetune/4_split4_random_finetune.csv'

        split_metabolites_4(input_csv, split1_csv, split2_csv, split3_csv, split4_csv)

        reformat_for_chemformer(split1_csv, fine_tune_split1)
        reformat_for_chemformer(split2_csv, fine_tune_split2)
        reformat_for_chemformer(split3_csv, fine_tune_split3)
        reformat_for_chemformer(split4_csv, fine_tune_split4)

        add_possible_products(fine_tune_split1)
        add_possible_products(fine_tune_split2)
        add_possible_products(fine_tune_split3)
        add_possible_products(fine_tune_split4)