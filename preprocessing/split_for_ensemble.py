import pandas as pd
import numpy as np
from rdkit import Chem, DataStructs
from scipy.cluster.hierarchy import linkage, fcluster
from scipy.spatial.distance import squareform
#from add_possible_products import add_possible_products

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

def split_metabolites_4(input_csv, split1_csv, split2_csv, split3_csv, split4_csv):

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

def split_metabolites_6(input_csv, split1_csv, split2_csv, split3_csv, split4_csv, split5_csv, split6_csv):

    np.random.seed(4)

    # Initialize lists for each split
    df_split1_list = []
    df_split2_list = []
    df_split3_list = []
    df_split4_list = []
    df_split5_list = []
    df_split6_list = []

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
            # Duplicate the single row across all splits
            df_split1_list.append(parent_rows.iloc[0])
            df_split2_list.append(parent_rows.iloc[0])
            df_split3_list.append(parent_rows.iloc[0])
            df_split4_list.append(parent_rows.iloc[0])
            df_split5_list.append(parent_rows.iloc[0])
            df_split6_list.append(parent_rows.iloc[0])

        elif num_rows == 2:
            # Distribute two rows across six splits, with duplication
            df_split1_list.append(parent_rows.iloc[0])
            df_split2_list.append(parent_rows.iloc[0])
            df_split3_list.append(parent_rows.iloc[0])
            df_split4_list.append(parent_rows.iloc[1])
            df_split5_list.append(parent_rows.iloc[1])
            df_split6_list.append(parent_rows.iloc[1])

        elif num_rows == 3:
            # Distribute three rows across six splits, with duplication
            df_split1_list.append(parent_rows.iloc[0])
            df_split2_list.append(parent_rows.iloc[1])
            df_split3_list.append(parent_rows.iloc[2])
            df_split4_list.append(parent_rows.iloc[0])
            df_split5_list.append(parent_rows.iloc[1])
            df_split6_list.append(parent_rows.iloc[2])

        elif num_rows == 4:
            df_split1_list.append(parent_rows.iloc[0])
            df_split2_list.append(parent_rows.iloc[1])
            df_split3_list.append(parent_rows.iloc[2])
            df_split4_list.append(parent_rows.iloc[3])
            random_index1 = np.random.randint(0, num_rows)
            random_index2 = np.random.randint(0, num_rows)
            df_split5_list.append(parent_rows.iloc[random_index1])
            df_split6_list.append(parent_rows.iloc[random_index2])

        elif num_rows == 5:
            df_split1_list.append(parent_rows.iloc[0])
            df_split2_list.append(parent_rows.iloc[1])
            df_split3_list.append(parent_rows.iloc[2])
            df_split4_list.append(parent_rows.iloc[3])
            df_split5_list.append(parent_rows.iloc[4])
            random_index = np.random.randint(0, num_rows)
            df_split6_list.append(parent_rows.iloc[random_index])

        elif num_rows > 5:
            equally_divisible_rows = num_rows - (num_rows % 6)
            parent_rows_list = list(parent_rows.iterrows())
            for i in range(equally_divisible_rows):
                _, row = parent_rows_list[i]
                if i % 6 == 0:
                    df_split1_list.append(row)
                elif i % 6 == 1:
                    df_split2_list.append(row)
                elif i % 6 == 2:
                    df_split3_list.append(row)
                elif i % 6 == 3:
                    df_split4_list.append(row)
                elif i % 6 == 4:
                    df_split5_list.append(row)
                elif i % 6 == 5:
                    df_split6_list.append(row)

            # Add remaining rows into a random split without duplicates
            split_list = [1, 2, 3, 4, 5, 6]
            if num_rows > equally_divisible_rows:
                for i in range(equally_divisible_rows, num_rows):
                    _, row = parent_rows_list[i]

                    if not split_list:
                        split_list = [1, 2, 3, 4, 5, 6]

                    random_split_index = np.random.choice(split_list)
                    split_list.remove(random_split_index)

                    if random_split_index == 1:
                        df_split1_list.append(row)
                    elif random_split_index == 2:
                        df_split2_list.append(row)
                    elif random_split_index == 3:
                        df_split3_list.append(row)
                    elif random_split_index == 4:
                        df_split4_list.append(row)
                    elif random_split_index == 5:
                        df_split5_list.append(row)
                    elif random_split_index == 6:
                        df_split6_list.append(row)

    # Convert lists to DataFrames
    df_split1 = pd.DataFrame(df_split1_list)
    df_split2 = pd.DataFrame(df_split2_list)
    df_split3 = pd.DataFrame(df_split3_list)
    df_split4 = pd.DataFrame(df_split4_list)
    df_split5 = pd.DataFrame(df_split5_list)
    df_split6 = pd.DataFrame(df_split6_list)

    df_split1 = pd.concat([df_split1, df_val], ignore_index=True)
    df_split2 = pd.concat([df_split2, df_val], ignore_index=True)
    df_split3 = pd.concat([df_split3, df_val], ignore_index=True)
    df_split4 = pd.concat([df_split4, df_val], ignore_index=True)
    df_split5 = pd.concat([df_split5, df_val], ignore_index=True)
    df_split6 = pd.concat([df_split6, df_val], ignore_index=True)

    df_split1 = pd.concat([df_split1, df_test], ignore_index=True)
    df_split2 = pd.concat([df_split2, df_test], ignore_index=True)
    df_split3 = pd.concat([df_split3, df_test], ignore_index=True)
    df_split4 = pd.concat([df_split4, df_test], ignore_index=True)
    df_split5 = pd.concat([df_split5, df_test], ignore_index=True)
    df_split6 = pd.concat([df_split6, df_test], ignore_index=True)

    # Print lengths of each DataFrame
    print(len(df_split1))
    print(len(df_split2))
    print(len(df_split3))
    print(len(df_split4))
    print(len(df_split5))
    print(len(df_split6))

    # Convert to CSV
    df_split1.to_csv(split1_csv, index=False)
    df_split2.to_csv(split2_csv, index=False)
    df_split3.to_csv(split3_csv, index=False)
    df_split4.to_csv(split4_csv, index=False)
    df_split5.to_csv(split5_csv, index=False)
    df_split6.to_csv(split6_csv, index=False)

    # return df_split1, df_split2, df_split3, df_split4, df_split5, df_split6

def reformat_for_chemformer(input_file, output_file):
    df = pd.read_csv(input_file)

    df = df.rename(columns={
        "child_smiles": "products",
        "parent_smiles": "reactants",
    })

    df.to_csv(output_file, sep='\t', index=False)

   #add_possible_products(output_file)

def split_tanimoto_score(child_smiles, n_splits): #input is a list of smiles


    child_mol = [Chem.MolFromSmiles(x) for x in child_smiles]
    child_fps = [AllChem.GetMorganFingerprintAsBitVect(x, radius=2, nBits=1024) for x in child_mol]

    num_molecules = len(child_smiles)
    similarity_matrix = np.zeros((num_molecules, num_molecules)) 

    for i in range(num_molecules):
        for j in range(i+1, num_molecules):
            sim = DataStructs.TanimotoSimilarity(child_fps[i], child_fps[j]) 
            similarity_matrix[i,j] = sim
            similarity_matrix[j,i] = sim

    # Convert similarity matrix to distance matrix for clustering
    distance_matrix = 1 - similarity_matrix
    condensed_distance_matrix = squareform(distance_matrix, checks=False)

    # Perform hierarchical clustering
    linkage_matrix = linkage(condensed_distance_matrix, method='average')
    
    # Generate cluster labels
    cluster_labels = fcluster(linkage_matrix, n_splits, criterion='maxclust')

    # Organize molecules into splits based on cluster labels
    splits = [[] for _ in range(n_splits)]
    for molecule_idx, cluster_label in enumerate(cluster_labels):
        splits[cluster_label - 1].append(child_smiles[molecule_idx])

    return splits

splits = split_tanimoto_score(['CCO', 'CCN', 'CCC', 'COC'], 2)
for i, split in enumerate(splits):
    print(f"Split {i + 1}: {split}")


# if __name__ == "__main__":

#     splits = '6' #'4' #'6'

#     input_csv = 'dataset/curated_data/combined_smiles_clean.csv'

#     split1_csv = f'dataset/curated_data/{splits}_split1_combined_smiles_clean.csv'
#     split2_csv = f'dataset/curated_data/{splits}_split2_combined_smiles_clean.csv'
#     split3_csv = f'dataset/curated_data/{splits}_split3_combined_smiles_clean.csv'
#     split4_csv = f'dataset/curated_data/{splits}_split4_combined_smiles_clean.csv'
#     split5_csv = f'dataset/curated_data/{splits}_split5_combined_smiles_clean.csv'
#     split6_csv = f'dataset/curated_data/{splits}_split6_combined_smiles_clean.csv'

#     fine_tune_split1 = f'dataset/finetune/{splits}_split1_finetune.csv'
#     fine_tune_split2 = f'dataset/finetune/{splits}_split2_finetune.csv'
#     fine_tune_split3 = f'dataset/finetune/{splits}_split3_finetune.csv'
#     fine_tune_split4 = f'dataset/finetune/{splits}_split4_finetune.csv'
#     fine_tune_split5 = f'dataset/finetune/{splits}_split5_finetune.csv'
#     fine_tune_split6 = f'dataset/finetune/{splits}_split6_finetune.csv'

#     count_metabolites(input_csv)

#     if splits == '4':

#         split_metabolites_4(input_csv, split1_csv, split2_csv, split3_csv, split4_csv)

#         reformat_for_chemformer(split1_csv, fine_tune_split1)
#         reformat_for_chemformer(split2_csv, fine_tune_split2)
#         reformat_for_chemformer(split3_csv, fine_tune_split3)
#         reformat_for_chemformer(split4_csv, fine_tune_split4)

#     if splits == '6':
#         split_metabolites_6(input_csv, split1_csv, split2_csv, split3_csv, split4_csv, split5_csv, split6_csv)

#         reformat_for_chemformer(split1_csv, fine_tune_split1)
#         reformat_for_chemformer(split2_csv, fine_tune_split2)
#         reformat_for_chemformer(split3_csv, fine_tune_split3)
#         reformat_for_chemformer(split4_csv, fine_tune_split4)
#         reformat_for_chemformer(split5_csv, fine_tune_split5)
#         reformat_for_chemformer(split6_csv, fine_tune_split6)