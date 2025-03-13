import pandas as pd
import json
import ast
import sys
import os
from rdkit import Chem
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

        repeated_index = [index] * len(target_smiles) 

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


def present_result(input1,input2):

    test_df = pd.read_csv(input1)
    prediction_df = pd.read_csv(input2)

    def agg_order_preserved(series):
        return list(series)

    grouped = test_df.groupby(['parent_name', 'parent_smiles'], sort=False).agg({
        'child_name': agg_order_preserved, 
        'child_smiles': agg_order_preserved 
    }).reset_index()

    if len(grouped) == len(prediction_df):
        grouped['sampled_molecules'] = prediction_df['sampled_molecules']
    else:
        raise ValueError("The number of rows in 'grouped' does not match the number of rows in 'prediction_df'.")

    grouped.to_csv(input2, index=False)


def save_n_valid_smiles(input_file, max_metabolites=12):

    df = pd.read_csv(input_file)

    parent_smiles = df['parent_smiles']
    
    count = 0
    for i in range(len(parent_smiles)):
        sampled_molecules_i = ast.literal_eval(df.at[i, 'sampled_molecules'])
        parent_smiles_i = parent_smiles[i]

        valid_smiles = []
        count_dup = 0
        count_parent = 0
        for smiles in sampled_molecules_i:
            if Chem.MolFromSmiles(smiles) is not None:
                if smiles not in valid_smiles:
                    if smiles != parent_smiles_i:
                        valid_smiles.append(smiles)
                    else:
                        count_parent += 1
                else:
                    count_dup += 1

        # valid_smiles = [smile for smile in sampled_molecules_i if Chem.MolFromSmiles(smile) is not None]

        print("nr of valid smiles: ",len(valid_smiles))
        print('nr of dup: ', count_dup)
        print('nr of parent dup: ', count_parent)

        if len(valid_smiles) < max_metabolites:
            count += 1

        df.at[i, 'sampled_molecules'] = valid_smiles[:max_metabolites]

    print('Nr of drugs with less than 12 valid SMILES: ', count)
    df.to_csv(input_file, index=False)

def top_k_accuracy(parent_name, sampled_molecules, child_smiles, child_name):

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

    return score1, score3, score5, score10, scatter_x, scatter_y

def precision_score(parent_name, sampled_molecules, child_smiles):

    TP = 0
    FP = 0
    for i in range(len(parent_name)):
        sampled_molecules_i = ast.literal_eval(sampled_molecules[i])
        child_smiles_i = ast.literal_eval(child_smiles[i])

        sampled_molecules_i = standardize_smiles_collection(sampled_molecules_i, False)

        TP_i = 0
        for j in range(len(child_smiles_i)):
            for k in range(len(sampled_molecules_i)):
                if k <= j and child_smiles_i[j] == sampled_molecules_i[k]:
                    TP_i += 1
            
        TP = TP + TP_i
        FP = FP + (len(child_smiles_i) - TP_i)
    
    return TP / (TP + FP)


def score_result(input_file):

    df = pd.read_csv(input_file)

    parent_name = df['parent_name']
    sampled_molecules = df['sampled_molecules']
    child_name = df['child_name']
    child_smiles = df['child_smiles']

    score1, score3, score5, score10, scatter_x, scatter_y = top_k_accuracy(parent_name, sampled_molecules, child_smiles, child_name)

    print(f'Score top1: {score1:.3f}')
    print(f'Score top3: {score3:.3f}')
    print(f'Score top5: {score5:.3f}')
    print(f'Score top10: {score10:.3f}')

    precision = precision_score(parent_name, sampled_molecules, child_smiles)

    print(f'Precision: {precision:.3f}')

    print(scatter_x)
    print(scatter_y)
    




gloryx = 'dataset/curated_data/gloryx_smiles_clean.csv'
json_file = 'results/evaluation/predictions0.json'

name = 'version23'

csv_file = f"evaluation/predictions/result_{name}.csv"

json_to_csv(json_file, csv_file)
present_result(gloryx, csv_file)
save_n_valid_smiles(csv_file, 12)
score_result(csv_file)


'''
With parent dup, dup and smiles
Score top1: 0.405
Score top3: 0.568
Score top5: 0.568
Score top10: 0.595
Precision: 0.184
'''

'''
With dup and valid smiles
Score top1: 0.135
Score top3: 0.541
Score top5: 0.568
Score top10: 0.595
Precision: 0.132
'''

'''
With only valid smiles
Score top1: 0.135
Score top3: 0.541
Score top5: 0.568
Score top10: 0.595
Precision: 0.132
'''