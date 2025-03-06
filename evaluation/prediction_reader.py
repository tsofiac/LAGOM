import pandas as pd
import json
import ast
import sys
import os
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
from preprocessing.standardize_smiles import standardize_smiles_collection




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


def present_result(input1,input2,output):

    test_df = pd.read_csv(input1)
    prediction_df = pd.read_csv(input2)

    # Define a custom function to aggregate while preserving order
    def agg_order_preserved(series):
        return list(series)

    # Group by 'parent_name' and 'parent_smiles', keeping the order preserved
    grouped = test_df.groupby(['parent_name', 'parent_smiles'], sort=False).agg({
        'child_name': agg_order_preserved,  # Custom aggregate function
        'child_smiles': agg_order_preserved # Custom aggregate function
    }).reset_index()


    if len(grouped) == len(prediction_df):
        grouped['sampled_molecules'] = prediction_df['sampled_molecules']
    else:
        raise ValueError("The number of rows in 'grouped' does not match the number of rows in 'prediction_df'.")

    # Write the result to a new CSV file
    grouped.to_csv(output, index=False)

def score_result(input_file):

    df = pd.read_csv(input_file)

    parent_name = df['parent_name']
    sampled_molecules = df['sampled_molecules']
    child_name = df['child_name']
    child_smiles = df['child_smiles']
    
    score1 = 0
    score3 = 0
    score5 = 0
    score10 = 0
    scatter_x = []
    scatter_y = []
    for i in range(len(parent_name)):
        sampled_molecules_i = ast.literal_eval(sampled_molecules[i])
        child_smiles_i = ast.literal_eval(child_smiles[i])
        child_name_i = ast.literal_eval(child_name[i])

        sampled_molecules_i = standardize_smiles_collection(sampled_molecules_i, False)
        
        # print(f'For {parent_name[i]}: ')
        count_top1 = 0
        count_top3 = 0
        count_top5 = 0
        count_top10 = 0
        for j in range(len(child_smiles_i)):
            for k in range(len(sampled_molecules_i)):
                if child_smiles_i[j] == sampled_molecules_i[k]:
                    # print(f'\t"{child_name_i[j]}" matches with sampled molecule {k+1}')
                    if k < 10:
                        count_top10 += 1
                    if k < 5:
                        count_top5 += 1
                    if k < 3:
                        count_top3 += 1
                    if k < 1:
                        count_top1 += 1
                    break

        # print(f'{parent_name[i]}: {count_top10} / {len(child_smiles_i)}')

        scatter_x.append(len(child_smiles_i))
        scatter_y.append(count_top10)

        if count_top1 > 0:
            score1 += 1 
        if count_top3 > 0:
            score3 += 1 
        if count_top5 > 0:
            score5 += 1 
        if count_top10 > 0:
            score10 += 1 

        # total = len(child_smiles_i) if len(child_smiles_i) < 1 else 1
        # score1 = score1 + (count_top1 / total)
        # total = len(child_smiles_i) if len(child_smiles_i) < 3 else 3
        # score3 = score3 + (count_top3 / total)
        # total = len(child_smiles_i) if len(child_smiles_i) < 5 else 5
        # score5 = score5 + (count_top5 / total)
        # total = len(child_smiles_i) if len(child_smiles_i) < 10 else 10
        # score10 = score10 + (count_top10 / total)

    score1 = score1 / len(parent_name)
    score3 = score3 / len(parent_name)
    score5 = score5 / len(parent_name)
    score10 = score10 / len(parent_name)

    print('Score top1: ', score1)
    print('Score top3: ', score3)
    print('Score top5: ', score5)
    print('Score top10: ', score10)

    print(scatter_x)
    print(scatter_y)



gloryx = 'dataset/curated_data/gloryx_smiles_clean.csv'
json_file = 'results/evaluation/predictions0.json'

csv_file = "evaluation/predictions/prediction_version4.csv"
result_file = "evaluation/predictions/result_version4.csv"

json_to_csv(json_file, csv_file)
present_result(gloryx, csv_file, result_file)

score_result(result_file)

