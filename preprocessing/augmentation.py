import pandas as pd


def augment_drugbank(input_file):
    df_drugbank = pd.read_csv(input_file)

    new_reactions = []
    for index, row in df_drugbank.iterrows():
        # Check if the reactant is different from the origin
        if row['parent_id'] != row['origin']:
            # Create a new reaction with the origin as the reactant and keep the product same
            new_reaction = {
                'parent_id': row['origin'],
                'child_id': row['child_id'],
                'child_name': row['child_name'],
                'child_smiles': row['child_smiles'],
                'origin': row['origin'],
            }
            new_reactions.append(new_reaction)
    augmented_reactions = pd.DataFrame(new_reactions)

    unique_df_drugbank = df_drugbank[['parent_id', 'parent_name', 'parent_smiles']].drop_duplicates(subset='parent_id')
    augmented_reactions = pd.merge(
        augmented_reactions,
        unique_df_drugbank[['parent_id', 'parent_name', 'parent_smiles']],
        on='parent_id',
        how='left'
    )
    augmented_reactions['source'] = 'DrugBank Augmented'
    desired_column_order = [
        'parent_id', 'parent_name', 'child_id', 
        'child_name', 'parent_smiles', 'child_smiles',
        'origin', 'source'
    ]
    augmented_reactions = augmented_reactions[desired_column_order]
    return augmented_reactions

def augment_metxbiodb(input_file):
    df_metxbiodb = pd.read_csv(input_file)

    new_reactions = []
    for index, row in df_metxbiodb.iterrows():
        # Check if the reactant is different from the origin
        if row['parent_name'] != row['origin']:
            # Create a new reaction with the origin as the reactant and keep the product same
            new_reaction = {
                'parent_name': row['origin'],
                'child_name': row['child_name'],
                'child_smiles': row['child_smiles'],
                'origin': row['origin'],
            }
            new_reactions.append(new_reaction)
    augmented_reactions = pd.DataFrame(new_reactions)

    unique_df_metxbiodb = df_metxbiodb[['parent_name', 'parent_smiles']].drop_duplicates(subset='parent_name')
    augmented_reactions = pd.merge(
        augmented_reactions,
        unique_df_metxbiodb[['parent_name', 'parent_smiles']],
        on='parent_name',
        how='left'
    )
    augmented_reactions['source'] = 'MetXBioDB Augmented'
    desired_column_order = [
        'parent_name', 'child_name', 
        'parent_smiles', 'child_smiles',
        'origin', 'source'
    ]
    augmented_reactions = augmented_reactions[desired_column_order]
    return augmented_reactions

def join(df1, df2, output_file):
    combined = pd.concat([df1,df2], ignore_index=True)
    combined.to_csv(output_file, index=False)



if __name__ == "__main__":

    metxbiodb = 'dataset/curated_data/metxbiodb_smiles.csv'
    drugbank = 'dataset/curated_data/drugbank_smiles.csv'

    augmented_data = 'dataset/curated_data/augmented_data.csv'

    augmented_drugbank = augment_drugbank(drugbank)
    augmented_metxbiodb = augment_metxbiodb(metxbiodb)

    join(augmented_drugbank, augmented_metxbiodb, augmented_data)

    