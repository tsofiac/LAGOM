import pandas as pd

def add_source_column(csv, source):
    
    df = pd.read_csv(csv)

    if source == 'drugbank':
        df['source'] = 'DrugBank'
    elif source == 'metxbiodb':
        df['source'] = 'MetXBioDB'
    elif source == 'gloryx':
        df['source'] = 'GLORYx' 
    else:
        print["Look over this."]
    
    return df
    
def combine_datasets(df1, df2, output_csv=None):
    selected_df1 = df1[['parent_name', 'child_name', 'parent_smiles', 'child_smiles', 'source']].copy()
    selected_df2 = df2[['parent_name', 'child_name', 'parent_smiles', 'child_smiles', 'source']].copy()

    combined_df = pd.concat([selected_df1, selected_df2], ignore_index=True)

    # Save the combined DataFrame to a new CSV file
    if output_csv is None:
        return combined_df
    else:
        combined_df.to_csv(output_csv, index=False)
        return None


def remove_duplicates(combined_csv, removed_duplicates_csv):  # Removes duplicate reactions and modifies source column
    df = pd.read_csv(combined_csv)
    len_before = len(df)
    
    # Create a mask that marks the duplicates
    df['is_duplicate'] = df.duplicated(subset=['parent_smiles', 'child_smiles'], keep=False)  # Keep all duplicates
    
    # Separate the duplicates (keep all versions)
    duplicates_df = df[df['is_duplicate']]
    
    # Modify the 'source' column for the rows that are kept
    df.loc[df['is_duplicate'], 'source'] = 'Both'  # Change 'source' value for rows that are kept

    # Remove the 'is_duplicate' column as it's no longer needed
    df = df.drop(columns='is_duplicate')
    duplicates_df = duplicates_df.drop(columns='is_duplicate')
    
    # Remove duplicates based on 'parent_smiles' and 'child_smiles', keeping the first occurrence
    df = df.drop_duplicates(subset=['parent_smiles', 'child_smiles'], keep='first')

    len_after = len(df)
    
    print("Total data points removed due to duplicates:", len_before - len_after)

    df.to_csv(combined_csv, index=False)
    duplicates_df.to_csv(removed_duplicates_csv, index=False)


def compare_datasets(combined_csv, df2, removed_file):

    df1 = pd.read_csv(combined_csv)

    combined_df = combine_datasets(df1, df2)

    combined_df['is_duplicate'] = combined_df.duplicated(subset=['parent_smiles', 'child_smiles'], keep=False) 

    duplicates_df = combined_df[combined_df['is_duplicate']]

    combined_df = combined_df.drop(columns='is_duplicate')
    duplicates_df = duplicates_df.drop(columns='is_duplicate')

    combined_df = combined_df.drop_duplicates(subset=['parent_smiles', 'child_smiles'], keep='last')

    combined_df = combined_df[combined_df['source'] != 'GLORYx']

    print('Number of duplicates found: ', int(len(duplicates_df)/2))

    combined_df.to_csv(combined_csv, index=False)
    duplicates_df.to_csv(removed_file, index=False)



if __name__ == "__main__":

    combined_csv = 'dataset/curated_data/combined_smiles_clean.csv'
    removed_csv = 'dataset/removed_data/combined_removed_duplicates.csv'
    compare_removed_csv = 'dataset/removed_data/compare_removed_duplicates.csv'

    metxbiodb_df = add_source_column('dataset/curated_data/metxbiodb_smiles_clean.csv', 'metxbiodb')
    drugbank_df = add_source_column('dataset/curated_data/drugbank_smiles_clean.csv', 'drugbank')
    gloryx_df = add_source_column('dataset/curated_data/gloryx_smiles_clean.csv', 'gloryx')

    combine_datasets(metxbiodb_df, drugbank_df, combined_csv)
    remove_duplicates(combined_csv, removed_csv)

    compare_datasets(combined_csv, gloryx_df, compare_removed_csv)

    