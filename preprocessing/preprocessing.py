from rdkit import Chem
import pandas as pd
from sklearn.model_selection import GroupShuffleSplit
from sklearn.model_selection import train_test_split
from standardize_smiles import standardize_smiles_collection
# from standardize_smiles_ours import standardize_molecule, standardize_smiles_collection, standardize_smiles_main


def test_val_distribute(input_file, output_file, val_size): #Input csv file. Function adds an extra column named 'set' and distributed into 'val' and 'set'
    df = pd.read_csv(input_file)
    # Add an empty column named 'set'
    df['set'] = None

    # Random_state is for reproducibility
    train_df, val_df = train_test_split(df, test_size=val_size, random_state=26)

    # Assign "train" and "val" labels
    train_df['set'] = 'train'
    val_df['set'] = 'val'

    # Combine the train and validation DataFrames
    final_df = pd.concat([train_df, val_df]).reset_index(drop=True)

    final_df.to_csv(output_file, index=False)


def get_unique_parents(input_file, additional_file):
    df = pd.read_csv(input_file)
    # Ensure df has 'parent_smiles' column, else raise an error
    if 'parent_smiles' not in df.columns:
        raise ValueError("The DataFrame must contain a 'parent_smiles' column.")

    unique_parent_dataset = []

    # Get a list of columns excluding 'parent_smiles'
    other_columns = [col for col in df.columns if col != 'parent_smiles']

    # Loop through unique parent SMILES
    for parent in df['parent_smiles'].unique():
        # Find the first row where 'parent_smiles' matches the current parent
        row = df[df['parent_smiles'] == parent].iloc[0]
        # Append the entire row
        unique_parent_dataset.append(row)

    # Create a new DataFrame from the unique dataset, preserving other columns
    unique_parent_df = pd.DataFrame(unique_parent_dataset).reset_index(drop=True)
    unique_parent_df = unique_parent_df[['parent_smiles'] + other_columns]

    val_count_end = (unique_parent_df['set'] == 'val').sum()

    print('Unique parent val distribution: ', val_count_end/len(unique_parent_df))

    unique_parent_df.to_csv(input_file, index=False)
    unique_parent_df.to_csv(additional_file, columns=['parent_smiles','parent_name','child_name','child_smiles'],index=False)


def reformat_for_chemformer(input_file):
    df = pd.read_csv(input_file)
    try:
        new_df = pd.DataFrame({
            'products': df['child_smiles'],
            'reactants': df['parent_smiles'],
            'set': df['set']
        })
        
    except KeyError as e:
        print(f"Error: Missing required column in the input data: {e}")
    except Exception as e:
        print(f"An unexpected error occurred: {e}")

    new_df.to_csv(input_file, sep='\t', index=False)



def small_dataset(input_file, output_file, size):
    # Read the CSV file into a DataFrame
    df = pd.read_csv(input_file)
    
    # Check if size is less than or equal to the number of rows in the DataFrame
    if size > len(df):
        raise ValueError("Size is greater than the number of rows in the dataset.")
    
    # Select a random sample of rows
    random_sample = df.sample(n=size)
    
    # Write the random sample to an output CSV file
    random_sample.to_csv(output_file, index=False)




if __name__ == "__main__":

    val_size = 0.1
    preprocess = True
    preprocess_unique_parents = True

    size = 733
    get_small_dataset = False

    name = 'drugbank' # [ 'drugbank' 'metxbiodb' ]
    
    # SPECIFY WHICH DATASET TO USE!
    dataset = f'dataset/curated_data/{name}_smiles_clean.csv'

    finetune = f'dataset/finetune/{name}_finetune.csv'
    finetune_small = f'dataset/finetune/{name}_finetune_small.csv'
    if preprocess:
        test_val_distribute(dataset, finetune, val_size)
        reformat_for_chemformer(finetune)
        if get_small_dataset:
            small_dataset(finetune, finetune_small, size)

    unique = f'dataset/curated_data/{name}_unique_parents.csv'
    unique_finetune = f'dataset/finetune/{name}_unique_parents_finetune.csv'
    unique_finetune_small = f'dataset/finetune/{name}_unique_parents_finetune_small.csv'
    if preprocess_unique_parents:
        test_val_distribute(dataset, unique_finetune, val_size)
        get_unique_parents(unique_finetune, unique)
        reformat_for_chemformer(unique_finetune)
        if get_small_dataset:
            small_dataset(unique_finetune, unique_finetune_small, size)

