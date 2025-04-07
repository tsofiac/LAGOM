import pandas as pd
import json
import ast
import sys
import os
from rdkit import Chem, DataStructs
from rdkit.Chem import AllChem 
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

def present_result(input1,input2, output):

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
        grouped['log_lhs'] = prediction_df['log_lhs']
    else:
        raise ValueError("The number of rows in 'grouped' does not match the number of rows in 'prediction_df'.")

    grouped.to_csv(output, index=False)

# def calculate_fingerprint_similarity(parent, sampled):

#     parent_mol = Chem.MolFromSmiles(parent)
#     sampled_mol = Chem.MolFromSmiles(sampled)

#     parent_fps = AllChem.GetMorganFingerprintAsBitVect(parent_mol, radius=2, nBits=1024)
#     sampled_fps = AllChem.GetMorganFingerprintAsBitVect(sampled_mol, radius=2, nBits=1024)

#     fingerprint_similarity = DataStructs.TanimotoSimilarity(parent_fps, sampled_fps) 
   
#     return fingerprint_similarity

def save_valid_smiles(input_file):

    df = pd.read_csv(input_file)

    parent_smiles = df['parent_smiles']
    
    total_valid_smiles = 0
    total_predictions = 0
    mean_predictions = 0
    for i in range(len(parent_smiles)):
        sampled_molecules_i = ast.literal_eval(df.at[i, 'sampled_molecules'])
        log_lhs_i = ast.literal_eval(df.at[i, 'log_lhs'])
        total_predictions += len(sampled_molecules_i)
        parent_smiles_i = parent_smiles[i]

        valid_smiles = []
        valid_log_lhs = []
        count_dup = 0
        count_parent = 0
        for j, smiles in enumerate(sampled_molecules_i):
            if Chem.MolFromSmiles(smiles) is not None:
                total_valid_smiles += 1
                smiles = standardize_molecule(smiles)
                if smiles not in valid_smiles:
                    if smiles != parent_smiles_i:
                        valid_smiles.append(smiles)
                        valid_log_lhs.append(log_lhs_i[j])
                    else:
                        count_parent += 1
                else:
                    count_dup += 1

        print('nr of dup: ', count_dup)
        print('nr of parent dup: ', count_parent)
        print("nr of valid smiles left: ",len(valid_smiles))

        mean_predictions += len(valid_smiles)

        df.at[i, 'sampled_molecules'] = valid_smiles
        df.at[i, 'log_lhs'] = valid_log_lhs

    validity = total_valid_smiles / total_predictions
    mean_validity = mean_predictions / len(parent_smiles)
    print(f'\nValidity: {validity:.3f}')
    print(f'Mean number of SMILES per drug: {mean_validity:.3f} / {len(sampled_molecules_i)}')

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
    
def concat_multiple_predictions(input1, input2, output):

    df1 = pd.read_csv(input1)
    df2 = pd.read_csv(input2)

    combined_sampled_molecules = []

    parent_name = df1['parent_name']
    sampled_molecules_df1 = df1['sampled_molecules']
    sampled_molecules_df2 = df2['sampled_molecules']

    for i in range(len(parent_name)):
        sampled_molecules_df1_i = ast.literal_eval(sampled_molecules_df1[i])
        sampled_molecules_df2_i = ast.literal_eval(sampled_molecules_df2[i])

        combined_sampled_molecules_i = list(set(sampled_molecules_df1_i).union(set(sampled_molecules_df2_i)))
        combined_sampled_molecules.append(combined_sampled_molecules_i)

    df1['sampled_molecules'] = combined_sampled_molecules

    df1.to_csv(output, index=False)

def count_correct_metabolites(input_file, max_metabolites, specification):

    df = pd.read_csv(input_file)

    df = specify_df(df, specification)

    parent_name = df['parent_name']
    child_smiles = df['child_smiles']
    sampled_molecules = df['sampled_molecules']
    df['sampled_boolean'] = None

    count = 0
    zero_count = 0
    top1 = []
    top3 = []
    top5 = []
    top10 = []
    all = []
    reference = []
    predictions = []
    for i in range(len(parent_name)):
        sampled_molecules_i = ast.literal_eval(sampled_molecules[i])
        child_smiles_i = ast.literal_eval(child_smiles[i])

        if len(sampled_molecules_i) < max_metabolites:
            count += 1
        if len(sampled_molecules_i) < 1:
            zero_count += 1

        count_top1 = 0
        count_top3 = 0
        count_top5 = 0
        count_top10 = 0
        count_all = 0
        sampled_boolean = [False] * len(sampled_molecules_i)
        for j in range(len(child_smiles_i)):
            for k in range(len(sampled_molecules_i)):
                if child_smiles_i[j] == sampled_molecules_i[k]:
                    count_all += 1
                    sampled_boolean[k] = True
                    if k < 10:
                        count_top10 += 1
                    if k < 5:
                        count_top5 += 1
                    if k < 3:
                        count_top3 += 1
                    if k < 1:
                        count_top1 += 1

        top1.append(count_top1)
        top3.append(count_top3)
        top5.append(count_top5)
        top10.append(count_top10)
        all.append(count_all)
        reference.append(len(child_smiles_i))
        predictions.append(len(sampled_molecules_i))

        df.at[i, 'sampled_boolean'] = sampled_boolean

    df.to_csv(input_file, index=False)
    print(f'\nNr of drugs with less than {max_metabolites} valid SMILES: ', count)
    print('Nr of drugs with zero valid SMILES: ', zero_count)
    return top1, top3, top5, top10, all, reference, predictions

def at_least_one_metabolite(top1, top3, top5, top10, all, reference):

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

    score_all = 0
    for i in range(len(all)):
        if all[i] >= 1:
            score_all += 1
    score_all = score_all / len(all)

    return score1, score3, score5, score10, score_all

def at_least_half_metabolites(top1, top3, top5, top10, all, reference):

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

    score_all = 0
    for i in range(len(all)):
        if all[i] >= math.ceil(0.5*reference[i]):# or top10[i] == 10:
            score_all += 1
    score_all = score_all / len(all)

    return score1, score3, score5, score10, score_all

def all_metabolites(top1, top3, top5, top10, all, reference):

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

    score_all = 0
    for i in range(len(all)):
        if all[i] == reference[i]:# or top10[i] == 10:
            score_all += 1
    score_all = score_all / len(all)

    return score1, score3, score5, score10, score_all

def precision_score(input_file):

    df = pd.read_csv(input_file)

    parent_name = df['parent_name']
    child_smiles = df['child_smiles']
    sampled_molecules = df['sampled_molecules']

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


def recall_and_precision(top10, reference, predictions):

    TP = sum(top10)
    FN = sum(reference) - TP
    FP = sum(predictions) - TP

    recall = TP / (TP + FN)
    precision = TP / (TP + FP)

    return recall, precision


def score_result(input_file, max_metabolites, specification):

    top1, top3, top5, top10, all, reference, predictions = count_correct_metabolites(input_file, max_metabolites, specification)

    print('\t')
    print(f'Total identified metabolites: {sum(all)} / {sum(reference)}')
    print(f'Total number of predictions: {sum(predictions)}')

    score1, score3, score5, score10, score_all = at_least_one_metabolite(top1, top3, top5, top10, all, reference)

    print('\nAt least one metabolite: ')
    print(f'Score top1: {score1:.3f}')
    print(f'Score top3: {score3:.3f}')
    print(f'Score top5: {score5:.3f}')
    print(f'Score top10: {score10:.3f}')
    print(f'Score all: {score_all:.3f}')

    score1, score3, score5, score10, score_all = at_least_half_metabolites(top1, top3, top5, top10, all, reference)

    print('\nAt least half metabolites: ')
    print(f'Score top1: {score1:.3f}')
    print(f'Score top3: {score3:.3f}')
    print(f'Score top5: {score5:.3f}')
    print(f'Score top10: {score10:.3f}')
    print(f'Score all: {score_all:.3f}')

    score1, score3, score5, score10, score_all = all_metabolites(top1, top3, top5, top10, all, reference)

    print('\nAll metabolites: ')
    print(f'Score top1: {score1:.3f}')
    print(f'Score top3: {score3:.3f}')
    print(f'Score top5: {score5:.3f}')
    print(f'Score top10: {score10:.3f}')
    print(f'Score all: {score_all:.3f}')

    precision = precision_score(input_file)

    print('\t')
    print(f'Our precision: {precision:.3f}')

    recall, precision = recall_and_precision(all, reference, predictions)

    print(f'Precision: {precision:.3f}')
    print(f'Recall: {recall:.3f}')

    print('\t')
    print(reference)
    print(all)
    


testset = 'dataset/curated_data/combined_evaluation.csv' # max: 10
# testset = 'dataset/curated_data/gloryx_smiles_clean.csv' # gloryx -- max: 12
json_predictions = 'results/evaluation/predictions0.json'

status = 'score' # 'score' 'combine' 'new'
name = 'version7_BS32'
specification = 0 # 0 (all) 1 (only_child) 2 (more than 1) 3 (more than 2) 
max_metabolites = 10

csv_predictions = f"evaluation/predictions/predictions_{name}.csv"
csv_result = f"evaluation/result/result_{name}.csv"

if status == 'new':
    json_to_csv(json_predictions, csv_predictions)
    present_result(testset, csv_predictions, csv_result)
    save_valid_smiles(csv_result)
    score_result(csv_result, max_metabolites, specification)

elif status == 'score':
    present_result(testset, csv_predictions, csv_result)
    save_valid_smiles(csv_result)
    score_result(csv_result, max_metabolites, specification)

elif status == 'combine':
    csv_comb = f"evaluation/result/result_comb_{name}.csv"
    # concat_multiple_predictions("evaluation/result/result_test1.csv", "evaluation/result/result_test2.csv", csv_comb)
    # score_result(csv_comb, max_metabolites, specification)
    concat_multiple_predictions(csv_comb, "evaluation/result/result_test0.csv", csv_comb)
    score_result(csv_comb, max_metabolites, specification)
else:
    print('Wrong status')

