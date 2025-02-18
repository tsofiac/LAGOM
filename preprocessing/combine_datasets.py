import pandas as pd

def add_source_column(csv, source): #source is either 'drugbank' or 'metxbiodb'
    
    df = pd.read_csv(csv)

    if source == 'drugbank':
        df['source'] = 'DrugBank'
    elif source == 'metxbiodb':
        df['source'] = 'MetXBioDB'
    else:
        print["Input source is not set as either 'drugbank' or 'metxbiodb'. Look over this."]
    
    return df
    
def combine_datasets(metxbiodb_df, drugbank_df, output_csv):
    # Select only the relevant columns
    drugbank_selected_df = drugbank_df[['parent_name', 'parent_smiles', 'child_name', 'child_smiles', 'source']].copy()
    metxbiodb_selected_df = metxbiodb_df[['parent_name', 'parent_smiles', 'child_name', 'child_smiles', 'source']].copy()
    # Concatenate the two DataFrames
    combined_df = pd.concat([metxbiodb_selected_df, drugbank_selected_df], ignore_index=True)

    # Save the combined DataFrame to a new CSV file
    combined_df.to_csv(output_csv, index=False)

    print("Combined CSV file created successfully!")

def remove_duplicates(combined_csv, removed_duplicates_csv):  # Removes duplicate reactions and modifies source column
    df = pd.read_csv(combined_csv)
    len_before = len(df)
    
    # Create a mask that marks the duplicates
    df['is_duplicate'] = df.duplicated(subset=['parent_smiles', 'child_smiles'], keep=False)  # Keep all duplicates
    
    # Separate the duplicates (keep all versions)
    duplicates_df = df[df['is_duplicate']]
    
    # Modify the 'source' column for the rows that are kept
    df.loc[~df['is_duplicate'], 'source'] = 'Both'  # Change 'source' value for rows that are kept

    # Remove the 'is_duplicate' column as it's no longer needed
    df = df.drop(columns='is_duplicate')
    
    # Remove duplicates based on 'parent_smiles' and 'child_smiles', keeping the first occurrence
    df = df.drop_duplicates(subset=['parent_smiles', 'child_smiles'], keep='first')

    len_after = len(df)
    
    print("Total data points removed due to duplicates:", len_before - len_after)

    df.to_csv(combined_csv, index=False)
    duplicates_df.to_csv(removed_duplicates_csv, index=False)

   
# def get_duplicates(source_df, target_df):
#     source_df = source_df[['parent_name', 'parent_smiles', 'child_name', 'child_smiles', 'source']].copy()
#     target_df = target_df[['parent_name', 'parent_smiles', 'child_name', 'child_smiles', 'source']].copy()

#     duplicates = []
#     for _, row_src in source_df.iterrows():
#         for _, row_tgt in target_df.iterrows():
#             if (row_src['parent_smiles'] == row_tgt['parent_smiles'] and row_src['child_smiles'] == row_tgt['child_smiles']):
#                 duplicates.append(list(row_src))
#                 duplicates.append(list(row_tgt))
#                 break 
#     df_duplicates = pd.DataFrame(duplicates)

#     df_duplicates.to_csv('data_analysis/duplicates.csv', index=False)

#     return df_duplicates

if __name__ == "__main__":

    output_csv = 'dataset/curated_data/combined_smiles_clean.csv'

    metxbiodb_df = add_source_column('dataset/curated_data/metxbiodb_smiles_clean.csv', 'metxbiodb')
    drugbank_df = add_source_column('dataset/curated_data/drugbank_smiles_clean.csv', 'drugbank')

    combine_datasets(metxbiodb_df, drugbank_df, output_csv)

    remove_duplicates(output_csv, 'dataset/removed_data/combined_duplicates.csv')