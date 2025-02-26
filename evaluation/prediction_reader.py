import pandas as pd
import json

def json_to_csv(json_file, csv_file):

    with open(json_file, 'r') as file:
        data = json.load(file)

    data_entries = data.get('data', [])

    indices = []
    log_lhs_list = []
    sampled_molecules_list = []
    target_smiles_list = []

    for entry in data_entries:
        index = entry.get('index')
        log_lhs = entry.get('log_lhs', [])
        sampled_molecules = entry.get('sampled_molecules', [])
        target_smiles = entry.get('target_smiles', '')

        repeated_index = [index] * len(target_smiles) #* nr_of_samples
        # ----------------------------------------------------------------------------------
        # nr_of_samples = len(sampled_molecules[0])
        # repeated_log_lhs = [sublist for sublist in log_lhs for _ in range(nr_of_samples)]
        # ----------------------------------------------------------------------------------
        
        indices.extend(repeated_index)
        log_lhs_list.extend(log_lhs)
        target_smiles_list.extend(target_smiles)
        sampled_molecules_list.extend(sampled_molecules)

    df = pd.DataFrame({
        'index': indices,
        'log_lhs': log_lhs_list,
        'target_smiles': target_smiles_list,
        'sampled_molecules': sampled_molecules_list
    })

    df.to_csv(csv_file, index=False)



json_file = 'results/evaluation/predictions0.json'
csv_file = f"evaluation/predictions/prediction.csv"

json_to_csv(json_file, csv_file)


