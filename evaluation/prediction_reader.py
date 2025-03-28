import pandas as pd
import json
import ast
import sys
import os
from rdkit import Chem
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
from preprocessing.standardize_smiles import standardize_smiles_collection, standardize_molecule
import math


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

    grouped = test_df.groupby(['parent_smiles'], sort=False).agg({
        'parent_name': 'first',  # Get the first entry for 'parent_name' in each group
        'child_name': agg_order_preserved, 
        'child_smiles': agg_order_preserved 
    }).reset_index()

    if len(grouped) == len(prediction_df):
        grouped['sampled_molecules'] = prediction_df['sampled_molecules']
    else:
        raise ValueError("The number of rows in 'grouped' does not match the number of rows in 'prediction_df'.")

    grouped.to_csv(input2, index=False)


def save_n_valid_smiles(input_file, max_metabolites=10):

    df = pd.read_csv(input_file)

    parent_smiles = df['parent_smiles']
    
    total_valid_smiles = 0
    total_predictions = 0
    for i in range(len(parent_smiles)):
        sampled_molecules_i = ast.literal_eval(df.at[i, 'sampled_molecules'])
        total_predictions += len(sampled_molecules_i)
        parent_smiles_i = parent_smiles[i]

        valid_smiles = []
        count_dup = 0
        count_parent = 0
        for smiles in sampled_molecules_i:
            if Chem.MolFromSmiles(smiles) is not None:
                total_valid_smiles += 1
                smiles = standardize_molecule(smiles)
                if smiles not in valid_smiles:
                    if smiles != parent_smiles_i:
                        valid_smiles.append(smiles)
                    else:
                        count_parent += 1
                else:
                    count_dup += 1

        print('nr of dup: ', count_dup)
        print('nr of parent dup: ', count_parent)
        print("nr of valid smiles left: ",len(valid_smiles))

        df.at[i, 'sampled_molecules'] = valid_smiles[:max_metabolites]

    validity = total_valid_smiles / total_predictions
    print('\nValidity: ', validity)

    df.to_csv(input_file, index=False)

def specify_df(df, specification = None):

    parent_name = df['parent_name']
    child_smiles = df['child_smiles']

    if specification == 1:

        mask = []
        for i in range(len(parent_name)):
            child_smiles_i = ast.literal_eval(child_smiles[i])

            if len(child_smiles_i) == 1:
                mask.append(True)
            else:
                mask.append(False)

        filtered_df = df[mask]
        filtered_df = filtered_df.reset_index(drop=True)

        return filtered_df
    
    elif specification == 2:

        mask = []
        for i in range(len(parent_name)):
            child_smiles_i = ast.literal_eval(child_smiles[i])

            if len(child_smiles_i) > 1:
                mask.append(True)
            else:
                mask.append(False)

        filtered_df = df[mask]
        filtered_df = filtered_df.reset_index(drop=True)

        return filtered_df
    
    elif specification == 3:

        mask = []
        for i in range(len(parent_name)):
            child_smiles_i = ast.literal_eval(child_smiles[i])

            if len(child_smiles_i) > 2:
                mask.append(True)
            else:
                mask.append(False)

        filtered_df = df[mask]
        filtered_df = filtered_df.reset_index(drop=True)

        return filtered_df

    else:

        return df


# def at_least_one_metabolite(parent_name, sampled_molecules, child_smiles, child_name):

#     score1 = 0
#     score3 = 0
#     score5 = 0
#     score10 = 0
#     scatter_x = []
#     scatter_y = []
#     for i in range(len(parent_name)):
#         sampled_molecules_i = ast.literal_eval(sampled_molecules[i])
#         child_smiles_i = ast.literal_eval(child_smiles[i])
#         # child_name_i = ast.literal_eval(child_name[i])
        
#         # print(f'For {parent_name[i]}: ')
#         count_top1 = 0
#         count_top3 = 0
#         count_top5 = 0
#         count_top10 = 0
#         for j in range(len(child_smiles_i)):
#             for k in range(len(sampled_molecules_i)):
#                 if child_smiles_i[j] == sampled_molecules_i[k]:
#                     # print(f'\t"{child_name_i[j]}" matches with sampled molecule {k+1}')
#                     if k < 10:
#                         count_top10 += 1
#                     if k < 5:
#                         count_top5 += 1
#                     if k < 3:
#                         count_top3 += 1
#                     if k < 1:
#                         count_top1 += 1
#                     break

#         # print(f'{parent_name[i]}: {count_top10} / {len(child_smiles_i)}')
#         scatter_x.append(len(child_smiles_i))
#         scatter_y.append(count_top10)

#         if count_top1 > 0:
#             score1 += 1 
#         if count_top3 > 0:
#             score3 += 1 
#         if count_top5 > 0:
#             score5 += 1 
#         if count_top10 > 0:
#             score10 += 1 

#          # total = len(child_smiles_i) if len(child_smiles_i) < 1 else 1
#          # score1 = score1 + (count_top1 / total)
#          # total = len(child_smiles_i) if len(child_smiles_i) < 3 else 3
#          # score3 = score3 + (count_top3 / total)
#          # total = len(child_smiles_i) if len(child_smiles_i) < 5 else 5
#          # score5 = score5 + (count_top5 / total)
#          # total = len(child_smiles_i) if len(child_smiles_i) < 10 else 10
#          # score10 = score10 + (count_top10 / total)

#     score1 = score1 / len(parent_name)
#     score3 = score3 / len(parent_name)
#     score5 = score5 / len(parent_name)
#     score10 = score10 / len(parent_name)

#     return score1, score3, score5, score10, scatter_x, scatter_y

# def top_k_accuracy(parent_name, sampled_molecules, child_smiles):

#     score1 = 0
#     score3 = 0
#     score5 = 0
#     score10 = 0
#     scatter_x = []
#     scatter_y = []
#     nr_of_metabolites = 0
#     for i in range(len(parent_name)):
#         sampled_molecules_i = ast.literal_eval(sampled_molecules[i])
#         child_smiles_i = ast.literal_eval(child_smiles[i])

#         nr_of_metabolites += len(child_smiles_i)
#         for j in range(len(child_smiles_i)):
#             for k in range(len(sampled_molecules_i)):
#                 if child_smiles_i[j] == sampled_molecules_i[k]:
#                     if k < 10:
#                         score10 += 1
#                     if k < 5:
#                         score5 += 1
#                     if k < 3:
#                         score3 += 1
#                     if k < 1:
#                         score1 += 1
#                     break

#     print(score1)
#     print(score3)
#     print(score5)
#     print(score10)
#     print(nr_of_metabolites)
#     score1 = score1 / nr_of_metabolites
#     score3 = score3 / nr_of_metabolites
#     score5 = score5 / nr_of_metabolites
#     score10 = score10 / nr_of_metabolites

#     return score1, score3, score5, score10

def count_valid_smiles(parent_name, sampled_molecules, max_metabolites):

    count = 0
    zero_count = 0
    for i in range(len(parent_name)):
        sampled_molecules_i = ast.literal_eval(sampled_molecules[i])

        if len(sampled_molecules_i) < max_metabolites:
            count += 1
        if len(sampled_molecules_i) < 1:
            zero_count += 1

    print(f'\nNr of drugs with less than {max_metabolites} valid SMILES: ', count)
    print('Nr of drugs with zero valid SMILES: ', zero_count)

def count_correct_metabolites(parent_name, sampled_molecules, child_smiles):

    top1 = []
    top3 = []
    top5 = []
    top10 = []
    reference = []
    predictions = []
    for i in range(len(parent_name)):
        sampled_molecules_i = ast.literal_eval(sampled_molecules[i])
        child_smiles_i = ast.literal_eval(child_smiles[i])

        count_top1 = 0
        count_top3 = 0
        count_top5 = 0
        count_top10 = 0
        for j in range(len(child_smiles_i)):
            for k in range(len(sampled_molecules_i)):
                if child_smiles_i[j] == sampled_molecules_i[k]:
                    if k < 10:
                        count_top10 += 1
                    if k < 5:
                        count_top5 += 1
                    if k < 3:
                        count_top3 += 1
                    if k < 1:
                        count_top1 += 1
                    break

        top1.append(count_top1)
        top3.append(count_top3)
        top5.append(count_top5)
        top10.append(count_top10)
        reference.append(len(child_smiles_i))
        predictions.append(len(sampled_molecules_i))

    return top1, top3, top5, top10, reference, predictions

def at_least_one_metabolite(top1, top3, top5, top10, reference):

    score1 = 0
    for i in range(len(top1)):
        if top1[i] >= 1:
            score1 += 1
    score1 = score1 / len(top1)

    score3 = 0
    for i in range(len(top3)):
        if top3[i] >= 1:
            score3 += 1
    score3 = score3 / len(top3)

    score5 = 0
    for i in range(len(top5)):
        if top5[i] >= 1:
            score5 += 1
    score5 = score5 / len(top5)

    score10 = 0
    for i in range(len(top10)):
        if top10[i] >= 1:
            score10 += 1
    score10 = score10 / len(top10)

    return score1, score3, score5, score10

def at_least_half_metabolites(top1, top3, top5, top10, reference):

    score1 = 0
    for i in range(len(top1)):
        if top1[i] >= math.ceil(0.5*reference[i]):# or top1[i] == 1:
            score1 += 1
    score1 = score1 / len(top1)

    score3 = 0
    for i in range(len(top3)):
        if top3[i] >= math.ceil(0.5*reference[i]):# or top3[i] == 3:
            score3 += 1
    score3 = score3 / len(top3)

    score5 = 0
    for i in range(len(top5)):
        if top5[i] >= math.ceil(0.5*reference[i]):# or top5[i] == 5:
            score5 += 1
    score5 = score5 / len(top5)

    score10 = 0
    for i in range(len(top10)):
        if top10[i] >= math.ceil(0.5*reference[i]):# or top10[i] == 10:
            score10 += 1
    score10 = score10 / len(top10)

    return score1, score3, score5, score10

def all_metabolites(top1, top3, top5, top10, reference):

    score1 = 0
    for i in range(len(top1)):
        if top1[i] == reference[i]:# or top1[i] == 1:
            score1 += 1
    score1 = score1 / len(top1)

    score3 = 0
    for i in range(len(top3)):
        if top3[i] == reference[i]:# or top3[i] == 3:
            score3 += 1
    score3 = score3 / len(top3)

    score5 = 0
    for i in range(len(top5)):
        if top5[i] == reference[i]:# or top5[i] == 5:
            score5 += 1
    score5 = score5 / len(top5)

    score10 = 0
    for i in range(len(top10)):
        if top10[i] == reference[i]:# or top10[i] == 10:
            score10 += 1
    score10 = score10 / len(top10)

    return score1, score3, score5, score10

def precision_score(parent_name, sampled_molecules, child_smiles):

    TP = 0
    FP = 0
    for i in range(len(parent_name)):
        sampled_molecules_i = ast.literal_eval(sampled_molecules[i])
        child_smiles_i = ast.literal_eval(child_smiles[i])

        TP_i = 0
        for j in range(len(child_smiles_i)):
            for k in range(len(sampled_molecules_i)):
                if k <= j and child_smiles_i[j] == sampled_molecules_i[k]:
                    TP_i += 1
            
        TP = TP + TP_i
        FP = FP + (len(child_smiles_i) - TP_i)
    
    return TP / (TP + FP)

def recall_score(parent_name, sampled_molecules, child_smiles):

    TP = 0
    FN = 0
    for i in range(len(parent_name)):
        sampled_molecules_i = ast.literal_eval(sampled_molecules[i])
        child_smiles_i = ast.literal_eval(child_smiles[i])

        TP_i = 0
        for j in range(len(child_smiles_i)):
            for k in range(len(sampled_molecules_i)):
                if child_smiles_i[j] == sampled_molecules_i[k]:
                    TP_i += 1
                    break

        TP = TP + TP_i
        FN = FN + (len(child_smiles_i) - TP_i)
    
    return TP / (TP + FN)

def recall_and_precision(top10, reference, predictions):

    TP = sum(top10)
    FN = sum(reference) - TP
    FP = sum(predictions) - TP

    recall = TP / (TP + FN)
    precision = TP / (TP + FP)

    return recall, precision


def score_result(input_file, max_metabolites, specification):

    df = pd.read_csv(input_file)

    df = specify_df(df, specification)

    parent_name = df['parent_name']
    sampled_molecules = df['sampled_molecules']
    child_smiles = df['child_smiles']

    count_valid_smiles(parent_name, sampled_molecules, max_metabolites)

    top1, top3, top5, top10, reference, predictions = count_correct_metabolites(parent_name, sampled_molecules, child_smiles)

    print('\t')
    print(f'Total identified metabolites: {sum(top10)} / {sum(reference)}')
    print('Total number of predictions: ', sum(predictions))

    score1, score3, score5, score10 = at_least_one_metabolite(top1, top3, top5, top10, reference)

    print('\nAt least one metabolite: ')
    print(f'Score top1: {score1:.3f}')
    print(f'Score top3: {score3:.3f}')
    print(f'Score top5: {score5:.3f}')
    print(f'Score top10: {score10:.3f}')

    score1, score3, score5, score10 = at_least_half_metabolites(top1, top3, top5, top10, reference)

    print('\nAt least half metabolites: ')
    print(f'Score top1: {score1:.3f}')
    print(f'Score top3: {score3:.3f}')
    print(f'Score top5: {score5:.3f}')
    print(f'Score top10: {score10:.3f}')

    score1, score3, score5, score10 = all_metabolites(top1, top3, top5, top10, reference)

    print('\nAll metabolites: ')
    print(f'Score top1: {score1:.3f}')
    print(f'Score top3: {score3:.3f}')
    print(f'Score top5: {score5:.3f}')
    print(f'Score top10: {score10:.3f}')

    precision = precision_score(parent_name, sampled_molecules, child_smiles)

    print('\t')
    print(f'Our precision: {precision:.3f}')

    recall, precision = recall_and_precision(top10, reference, predictions)

    print(f'Precision: {precision:.3f}')
    print(f'Recall: {recall:.3f}')

    print('\t')
    print(reference)
    print(top10)
    


testset = 'dataset/curated_data/combined_evaluation.csv' # max: 10
# testset = 'dataset/curated_data/gloryx_smiles_clean.csv' # gloryx -- max: 12
json_file = 'results/evaluation/predictions0.json'

name = 'version57'
specification = None # 1 (only_child) 2 (more than 1) 3 (more than 2) None (all)
max_metabolites = 10

csv_file = f"evaluation/predictions/result_{name}.csv"

json_to_csv(json_file, csv_file)
present_result(testset, csv_file)
save_n_valid_smiles(csv_file, max_metabolites)
score_result(csv_file, max_metabolites, specification)

