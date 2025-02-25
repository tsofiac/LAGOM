import pandas as pd
import json
import ast

def json_to_csv(json_file, csv_file):

    # Load the data as a Python dictionary
    with open(json_file, 'r') as file:
        data = json.load(file)

    # Access the data entries
    data_entries = data.get('data', [])

    # Prepare lists to store each column's cumulative data
    indices = []
    log_lhs_list = []
    sampled_molecules_list = []
    target_smiles_list = []

    # Iterate through each entry and collect data for the DataFrame
    for entry in data_entries:
        index = entry.get('index')
        log_lhs = entry.get('log_lhs', [])
        sampled_molecules = entry.get('sampled_molecules', [])
        target_smiles = entry.get('target_smiles', '')

        # Append data to cumulative lists
        indices.append(index)
        log_lhs_list.append(log_lhs)
        
        # Flatten the nested list of sampled_molecules
        sampled_molecules_flat = [item[0] for item in sampled_molecules]
        sampled_molecules_list.append(';'.join(sampled_molecules_flat))

        target_smiles_list.append(target_smiles)

    # Create DataFrame from the cumulative lists
    df = pd.DataFrame({
        'index': indices,
        'log_lhs': log_lhs_list,
        'sampled_molecules': sampled_molecules_list,
        'target_smiles': target_smiles_list
    })

    # Write the entire DataFrame to a single CSV file
    df.to_csv(csv_file, index=False)

def reformat_csv(input, i):

    df = pd.read_csv(input)

    print('Number of indices: ', len(df))

    sampled_row = df['sampled_molecules'][i]
    sampled_row_list = sampled_row.split(';')

    target_row = df['target_smiles'][i]
    target_row_list = ast.literal_eval(target_row)

    log_row = df['log_lhs'][i]
    log_row_list = ast.literal_eval(log_row)

    df_new = pd.DataFrame({
            'log_lhs': log_row_list,
            'target_smiles': target_row_list,
            'sampled_molecules': sampled_row_list
        })

    df_new.to_csv(input, index=False)



# Choose prediction file
json_file = 'results/evaluation/predictions0.json'
# Choose index between 0-2
i = 0

csv_file = f"evaluation/predictions/prediction{i}.csv"
json_to_csv(json_file, csv_file)
reformat_csv(csv_file, i)