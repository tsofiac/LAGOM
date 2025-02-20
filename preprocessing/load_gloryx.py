import json
import pandas as pd


def load_gloryx(input_file, output_file):
    with open(input_file, 'r', encoding="utf8") as file:
        reader = file.read()
        # replace \ to \\ for correct reading of file
        reader = reader.replace('\\', '\\\\')
        data = json.loads(reader)

    child_smiles = []
    parent_smiles = []
    child_name = []
    parent_name = []
    generation_column = []
    set_column = []

    def process_metabolites(metabolites, drug_smiles, drug_name):
        for metabolite in metabolites:
            metabolite_smiles = metabolite['smiles']
            metabolite_name = metabolite.get('metaboliteName', None)
            generation = metabolite.get('generation', None)
            
            child_smiles.append(metabolite_smiles)
            parent_smiles.append(drug_smiles)
            child_name.append(metabolite_name)
            parent_name.append(drug_name)
            generation_column.append(generation)
            set_column.append('test')
            
            # Process any nested metabolites
            if 'metabolites' in metabolite:
                process_metabolites(metabolite['metabolites'], metabolite_smiles, metabolite_name)

    for entry in data:
        drug_name = entry['drugName']
        drug_smiles = entry['smiles']
        
        if 'metabolites' in entry:
            process_metabolites(entry['metabolites'], drug_smiles, drug_name)

    df = pd.DataFrame({
        'parent_name': parent_name,
        'child_name': child_name,
        'parent_smiles': parent_smiles,
        'child_smiles': child_smiles,
        'generation': generation_column,
        'set': set_column
    })

    generation_1_df = df[df['generation'] == 1]
    generation_1_df.to_csv(output_file, index=False)


gloryx = 'dataset/raw_data/gloryx_test_dataset.json'
gloryx_loaded = 'dataset/curated_data/gloryx_smiles.csv'

load_gloryx(gloryx, gloryx_loaded)
