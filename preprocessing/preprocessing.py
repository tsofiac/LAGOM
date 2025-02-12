from rdkit import Chem
import pandas as pd
from sklearn.model_selection import GroupShuffleSplit
from sklearn.model_selection import train_test_split
from standardize_smiles import standardize_smiles_collection
# from standardize_smiles_ours import standardize_molecule, standardize_smiles_collection, standardize_smiles_main


def test_val_distribute(input_file, output_file): #Input csv file. Function adds an extra column named 'set' and distributed into 'val' and 'set'
    df = pd.read_csv(input_file)
    # Add an empty column named 'set'
    df['set'] = None

    # Split the DataFrame into train (90%) and validation (10%) sets.
    # Random_state is for reproducibility
    train_df, val_df = train_test_split(df, test_size=0.1, random_state=26)

    # Assign "train" and "val" labels
    train_df['set'] = 'train'
    val_df['set'] = 'val'

    # Combine the train and validation DataFrames
    final_df = pd.concat([train_df, val_df]).reset_index(drop=True)

    print("The df is now distributed into training and validation")

    final_df.to_csv(output_file, index=False)


# From standardize_smiles_ours
def standardize_smiles_main(input_file): #input df that has columns 'parent_smiles' and 'child_smiles'
    df = pd.read_csv(input_file)
    df['parent_smiles'] = standardize_smiles_collection(df['parent_smiles'], False) #'False' eliminates isomeres
    df['child_smiles'] = standardize_smiles_collection(df['child_smiles'], False)
    df.to_csv(input_file, index=False)


def get_unique_parents(input_file):
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
    unique_parent_df.to_csv(input_file, index=False)


def modifiy_columns(input_file, output_file):
    df = pd.read_csv(input_file)
    try:
        new_df = pd.DataFrame({
            'products': df['child_smiles'],
            'reactants': df['parent_smiles'],
            'set': df['set']
        })
        
        print(f"Transformation completed. The df now has columns: products, reactants, and set (empty)")
    
    except KeyError as e:
        print(f"Error: Missing required column in the input data: {e}")
    except Exception as e:
        print(f"An unexpected error occurred: {e}")

    new_df.to_csv(output_file, sep='\t', index=False)


if __name__ == "__main__":

    preprocess_clean = True
    preprocess_clean_unique_parents = True
    preprocess_unique_parents_no_clean = True

    name = 'metxbiodb' # [ 'drugbank' 'metxbiodb' ]
    dataset = f'dataset/curated_data/{name}_smiles.csv'
    

    clean = f'dataset/curated_data/{name}_clean.csv'
    clean_finetune = f'dataset/finetune/{name}_clean_finetune.csv'
    if preprocess_clean:
        test_val_distribute(dataset, clean)
        standardize_smiles_main(clean)
        modifiy_columns(clean, clean_finetune)


    clean_unique = f'dataset/curated_data/{name}_clean_unique_parents.csv'
    clean_unique_finetune = f'dataset/finetune/{name}_clean_unique_parents_finetune.csv'
    if preprocess_clean_unique_parents:
        test_val_distribute(dataset, clean_unique)
        standardize_smiles_main(clean_unique)
        get_unique_parents(clean_unique)
        modifiy_columns(clean_unique, clean_unique_finetune)


    unique = f'dataset/curated_data/{name}_unique_parents.csv'
    unique_finetune = f'dataset/finetune/{name}_unique_parents_finetune.csv' 
    if preprocess_unique_parents_no_clean:
        test_val_distribute(dataset, unique)
        get_unique_parents(unique)
        modifiy_columns(unique, unique_finetune)



# def df_to_csv_tab(input_df, output_csv): # outputs a csv with tabs
#     input_df.to_csv(output_csv, sep='\t', index=False)

# def df_to_csv(input_df, output_csv):
#     input_df.to_csv(output_csv, index=False)


# # --- RUN FUNCTIONS --- CHANGE FILE NAMES
# def preprocess_clean_main(df): #Input df needs columns 'parent_smiles' and 'child_smiles'
    
#     df_distributed = test_val_distribute(df) # adds a set column and distrubutes the datapoints into val and set
#     df_clean = standardize_smiles_main(df_distributed) # Cleans data. Input df needs columns 'parent_smiles' and 'child_smiles'
#     df_finetune = modifiy_columns(df_clean) # Outputs df ready for fine-tuning but the set-column is empty
    
#     df_to_csv(df_clean, 'dataset/curated_data/drugbank_clean.csv')
#     df_to_csv_tab(df_finetune, 'dataset/curated_data/drugbank_clean_finetune.csv')

#     len_start = len(df)
#     len_clean = len(df_clean)

#     print("Just cleaning conducted")
#     print("Number of data points from start: ", len_start)
#     print("Number of data points after cleaning: ", len_clean)

#     val_count_start = (df_distributed['set'] == 'val').sum()
#     val_count_clean = (df_clean['set'] == 'val').sum()

#     print(f"The count of 'val' from start is: {val_count_start}")
#     print(f"The count of 'val' after cleaning is: {val_count_clean}")
#     print('The val distribution is now', val_count_clean/len_clean)

# def preprocess_clean_unique_parents_main(df): #Input df needs columns 'parent_smiles' and 'child_smiles'
    
#     df_distributed = test_val_distribute(df) # adds a set column and distrubutes the datapoints into val and set
#     df_clean = standardize_smiles_main(df_distributed) # Cleans data. Input df needs columns 'parent_smiles' and 'child_smiles'
#     df_clean_unique = get_unique_parents(df_clean)
#     df_finetune = modifiy_columns(df_clean_unique) # Outputs df ready for fine-tuning but the set-column is empty
    
#     df_to_csv(df_clean_unique, 'dataset/curated_data/drugbank_clean_unique_parents.csv')
#     df_to_csv_tab(df_finetune, 'dataset/curated_data/drugbank_clean_unique_parents_finetune.csv')

#     len_start = len(df)
#     len_end = len(df_clean_unique)

#     print("Cleaning conducted. Unique parents")
#     print("Number of data points from start: ", len_start)
#     print("Number of data points after cleaning and deleting any parent duplicates: ", len_end)

#     val_count_start = (df_distributed['set'] == 'val').sum()
#     val_count_end = (df_clean_unique['set'] == 'val').sum()

#     print(f"The count of 'val' from start is: {val_count_start}")
#     print(f"The count of 'val' after cleaning is: {val_count_end}")
#     print('The val distribution is now', val_count_end/len_end)

# def preprocess_unique_parents_no_clean_main(df): #Input df needs columns 'parent_smiles' and 'child_smiles'
    
#     df_distributed = test_val_distribute(df) # adds a set column and distrubutes the datapoints into val and set
#     df_unique = get_unique_parents(df_distributed)
#     df_finetune = modifiy_columns(df_unique) # Outputs df ready for fine-tuning but the set-column is empty
    
#     df_to_csv(df_unique, 'dataset/curated_data/drugbank_unique_parents.csv')
#     df_to_csv_tab(df_finetune, 'dataset/curated_data/drugbank_unique_parents_finetune.csv')

#     len_start = len(df)
#     len_end = len(df_unique)

#     print("Unique parents - No cleaning conducted")
#     print("Number of data points from start: ", len_start)
#     print("Number of data points after not deleting any parent duplicates: ", len_end)

#     val_count_start = (df_distributed['set'] == 'val').sum()
#     val_count_end = (df_unique['set'] == 'val').sum()

#     print(f"The count of 'val' from start is: {val_count_start}")
#     print(f"The count of 'val' after cleaning is: {val_count_end}")
#     print('The val distribution is now', val_count_end/len_end)


# drugbank_df = pd.read_csv('dataset/curated_data/drugbank_smiles.csv')

# preprocess_clean_main(drugbank_df)
# preprocess_clean_unique_parents_main(drugbank_df)
# preprocess_unique_parents_no_clean_main(drugbank_df)