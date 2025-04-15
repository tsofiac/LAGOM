import pandas as pd
import numpy as np

def count_metabolites(input_csv): # Counts metabolites in training data
    # Read the CSV file into a DataFrame
    df = pd.read_csv(input_csv)

    # Filter the DataFrame to include only rows where 'set' is 'train'
    df_train = df[df['set'] == 'train']

    # Group by 'parent_smiles' and count the number of metabolites for each parent
    metabolite_counts = df_train.groupby('parent_smiles').size()

    # Count how many parents have a specific number of metabolites
    count_of_metabolite_counts = metabolite_counts.value_counts().sort_index()

    # Print the results
    for metabolites, count in count_of_metabolite_counts.items():
        print(f"Nr of drugs with {metabolites} metabolite(s): {count}")

    # # print rows
    # # Filter to get only parents with exactly 5 metabolites
    # parents_with_five_metabolites = metabolite_counts[metabolite_counts == 1]

    # # Print parent identifiers with exactly 5 metabolites
    # print("Parents with 1 metabolites:")
    # for parent in parents_with_five_metabolites.index:
    #     print(parent)

def split_metabolites(input_csv, split1_csv, split2_csv, split3_csv, split4_csv):

    np.random.seed(6)

    df_split1 = []
    df_split2 = []
    df_split3 = []
    df_split4 = []

    # Read the CSV file into a DataFrame
    df = pd.read_csv(input_csv)

    # Filter the DataFrame to include only rows where 'set' is 'train'
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
    print(len(df))
    print(len(df_split1))
    print(len(df_split2))
    print(len(df_split3))
    print(len(df_split4))

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

input_csv = 'dataset/curated_data/combined_smiles_clean.csv'

split1_csv = 'dataset/curated_data/split1_combined_smiles_clean.csv'
split2_csv = 'dataset/curated_data/split2_combined_smiles_clean.csv'
split3_csv = 'dataset/curated_data/split3_combined_smiles_clean.csv'
split4_csv = 'dataset/curated_data/split4_combined_smiles_clean.csv'

fine_tune_split1 = 'dataset/finetune/split1_finetune.csv'
fine_tune_split2 = 'dataset/finetune/split2_finetune.csv'
fine_tune_split3 = 'dataset/finetune/split3_finetune.csv'
fine_tune_split4 = 'dataset/finetune/split4_finetune.csv'

count_metabolites(input_csv)
split_metabolites(input_csv, split1_csv, split2_csv, split3_csv, split4_csv)

reformat_for_chemformer(split1_csv, fine_tune_split1)
reformat_for_chemformer(split2_csv, fine_tune_split2)
reformat_for_chemformer(split3_csv, fine_tune_split3)
reformat_for_chemformer(split4_csv, fine_tune_split4)