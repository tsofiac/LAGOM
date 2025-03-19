import pandas as pd

#combined_mmp_file = 'dataset/curated_data/paired_mmp_all.csv' 

# df1 = pd.read_csv('dataset/curated_data/paired_mmp_rows_0_to_1101304.csv')
# df2 = pd.read_csv('dataset/curated_data/paired_mmp_rows_1101305_to_2202607.csv')
# df3 = pd.read_csv('dataset/curated_data/paired_mmp_rows_2202608_to_3303912.csv')
# df4 = pd.read_csv('dataset/curated_data/paired_mmp_rows_3303912_to_4405216.csv')
# df5 = pd.read_csv('dataset/curated_data/paired_mmp_rows_4405217_to_5506519.csv')
# df6 = pd.read_csv('dataset/curated_data/paired_mmp_rows_5506520_to_6607824.csv')
# df7 = pd.read_csv('dataset/curated_data/paired_mmp_rows_6607825_to_7709127.csv')
# df8 = pd.read_csv('dataset/curated_data/paired_mmp_rows_7709128_to_8810431.csv')
# df9 = pd.read_csv('dataset/curated_data/paired_mmp_rows_8810432_to_9911735.csv')
# df10 = pd.read_csv('dataset/curated_data/paired_mmp_rows_9911736_to_11013037.csv')

combined_mmp_file = 'dataset/curated_data/mmp_atoms_allowed_all_smiles_clean.csv'

df1 = pd.read_csv('dataset/curated_data/mmp_atoms_allowed_1_smiles_clean.csv')
df2 = pd.read_csv('dataset/curated_data/mmp_atoms_allowed_2_smiles_clean.csv')
df3 = pd.read_csv('dataset/curated_data/mmp_atoms_allowed_3_smiles_clean.csv')
df4 = pd.read_csv('dataset/curated_data/mmp_atoms_allowed_4_smiles_clean.csv')
df5 = pd.read_csv('dataset/curated_data/mmp_atoms_allowed_5_smiles_clean.csv')
df6 = pd.read_csv('dataset/curated_data/mmp_atoms_allowed_6_smiles_clean.csv')
df7 = pd.read_csv('dataset/curated_data/mmp_atoms_allowed_7_smiles_clean.csv')
df8 = pd.read_csv('dataset/curated_data/mmp_atoms_allowed_8_smiles_clean.csv')
df9 = pd.read_csv('dataset/curated_data/mmp_atoms_allowed_9_smiles_clean.csv')
df10 = pd.read_csv('dataset/curated_data/mmp_atoms_allowed_10_smiles_clean.csv')


import pandas as pd

def randomly_select_rows(csv_file, small_file, n):
    """
    Randomly select n rows from a CSV file.

    Parameters:
    - csv_file (str): The path to the CSV file.
    - n (int): The number of rows to randomly select.

    Returns:
    - pd.DataFrame: A DataFrame containing the randomly selected rows.
    """
    # Load the CSV file into a DataFrame
    df = pd.read_csv(csv_file, low_memory = False)
    
    # Randomly sample `n` rows from the DataFrame without replacement
    random_sample = df.sample(n=n, random_state=1)

    random_sample.to_csv(small_file, index=False)

def reformat_for_chemformer(input_file, output_file):
    df = pd.read_csv(input_file)

    df = df.rename(columns={
        "child_smiles": "products",
        "parent_smiles": "reactants",
    })

    df.to_csv(output_file, sep='\t', index=False)


df_combined = pd.concat([df1, df2, df3, df4, df5, df6, df7, df8, df9, df10])
df_combined.to_csv(combined_mmp_file, index=False)

n = 100000
small_file = f'dataset/curated_data/mmp_atoms_allowed_{n}rows_smiles_clean.csv'
randomly_select_rows(combined_mmp_file, small_file, n)

finetune_file = f'dataset/finetune/mmp_{n}rows_finetune.csv'
reformat_for_chemformer(small_file, finetune_file)

