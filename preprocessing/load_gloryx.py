import json
import pandas as pd


def load_gloryx(input_file, output_file):
    # Load JSON data from the input file
    with open(input_file, 'r') as f:
        data = json.load(f)

    # Prepare lists to hold the extracted data
    child_smiles = []
    parent_smiles = []
    child_name = []
    parent_name = []
    set_column = []

    def process_metabolites(metabolites, drug_smiles, drug_name):
        for metabolite in metabolites:
            metabolite_smiles = metabolite['smiles']
            metabolite_name = metabolite.get('metaboliteName', None)
            
            # Append the data to the lists
            child_smiles.append(metabolite_smiles)
            parent_smiles.append(drug_smiles)
            child_name.append(metabolite_name)
            parent_name.append(drug_name)
            set_column.append('test')
            
            # Process any nested metabolites
            if 'metabolites' in metabolite:
                process_metabolites(metabolite['metabolites'], metabolite_smiles, drug_name)


    # Iterate over each entry in the JSON data
    for entry in data:
        drug_name = entry['drugName']
        drug_smiles = entry['smiles']
        
        # Process each metabolite
        if 'metabolites' in entry:
            process_metabolites(entry['metabolites'], drug_smiles, drug_name)

    # Create a DataFrame with the collected data
    df = pd.DataFrame({
        'parent_name': parent_name,
        'child_name': child_name,
        'parent_smiles': parent_smiles,
        'child_smiles': child_smiles,
        'set': set_column
    })

    # Calculate number of unique parents by parent_name
    unique_parents_by_name = df['parent_name'].nunique()
    print(f"Number of unique parent names: {unique_parents_by_name}")

    # Alternatively, calculate number of unique parents by parent_smiles
    unique_parents_by_smiles = df['parent_smiles'].nunique()
    print(f"Number of unique parent SMILES: {unique_parents_by_smiles}")


    # Write the DataFrame to a TSV file (tab-separated values)
    df.to_csv(output_file, index=False)


gloryx = 'dataset/raw_data/gloryx_test_dataset.json'
gloryx_loaded = 'dataset/curated_data/gloryx_test_dataset.csv'

load_gloryx(gloryx, gloryx_loaded)
