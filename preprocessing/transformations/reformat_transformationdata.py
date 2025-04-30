import pandas as pd

def change_set_column(transformation_csv, original_csv, new_transformation_csv):
    # Read the CSV files into DataFrames
    transformations_df = pd.read_csv(transformation_csv)
    original_df = pd.read_csv(original_csv)
    
    # Identify columns that are common to both DataFrames, excluding 'transformation' column
    common_columns = transformations_df.columns.intersection(original_df.columns).tolist()
    common_columns.remove('set')

    # Merge on all common columns
    merged_df = original_df.merge(
        transformations_df[common_columns + ['transformation']],
        on=common_columns,
        how='left'
    )
    
    # Save the updated DataFrame to a new CSV file
    merged_df.to_csv(new_transformation_csv, index=False)

def remove_multiple_fragments(transformation_csv, filtered_transformation_csv, removed_data_csv):

    transformations_df = pd.read_csv(transformation_csv)

    #if there is a full stop (.) in the column 'transformation', remove row from df 
    filtered_df = transformations_df[~transformations_df['transformation'].str.contains('\.', na=False)]
    removed_df = transformations_df[transformations_df['transformation'].str.contains('\.', na=False)]

    filtered_df.to_csv(filtered_transformation_csv, index=False)
    removed_df.to_csv(removed_data_csv, index = False)

    print('filtering complete')

def reformat_for_chemformer(input_file, output_file):
    df = pd.read_csv(input_file)

    df = df.rename(columns={
        "child_smiles": "products",
        "parent_smiles": "reactants",
    })

    df.to_csv(output_file, sep='\t', index=False)

def create_transform_evaluation(filtered_transformation_csv, evaluation_transform_csv):
    df = pd.read_csv(filtered_transformation_csv)
    test_df = df[df['set'] == 'test']
    test_df.to_csv(evaluation_transform_csv, index=False)

def get_unique_parents(input_file, output_file):
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
    unique_parent_df = unique_parent_df[df.columns]

    unique_parent_df.to_csv(output_file, index=False)


transformation_csv = 'dataset/curated_data/combined_smiles_clean_transform_oldset.csv'
new_transformation_csv ='dataset/curated_data/combined_smiles_clean_transform_newset.csv'
original_csv = 'dataset/curated_data/combined_smiles_clean.csv'

filtered_transformation_csv = 'dataset/curated_data/combined_smiles_clean_transform_filtered.csv'
removed_data_csv ='dataset/curated_data/filtered/filtered_transformations.csv'

finetune_csv = 'dataset/finetune/combined_transformations_finetune.csv'

evaluation_all_csv = 'dataset/curated_data/combined_evaluation_transformations.csv'
evaluation_unique_csv = 'dataset/curated_data/combined_evaluation_transformations_unique.csv'
evaluation_finetune_csv = 'dataset/finetune/combined_evaluation_transformations_finetune.csv'


change_set_column(transformation_csv, original_csv, new_transformation_csv)
remove_multiple_fragments(new_transformation_csv, filtered_transformation_csv, removed_data_csv)
reformat_for_chemformer(filtered_transformation_csv, finetune_csv)

# evaluation file
create_transform_evaluation(filtered_transformation_csv, evaluation_all_csv)
get_unique_parents(evaluation_all_csv, evaluation_unique_csv)
reformat_for_chemformer(evaluation_unique_csv, evaluation_finetune_csv)