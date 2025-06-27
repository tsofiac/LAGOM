import pandas as pd
import json
import ast
import sys
import os
from rdkit import Chem, DataStructs

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))
from preprocessing.standardize_smiles import standardize_molecule
import math
import numpy as np
from rdkit.Chem import AllChem


def json_to_csv(json_file, csv_file):
    with open(json_file, "r") as file:
        data = json.load(file)

    data_entries = data.get("data", [])

    indices = []
    log_lhs_list = []
    sampled_molecules_list = []
    target_smiles_list = []

    for entry in data_entries:
        index = entry.get("index")
        log_lhs = entry.get("log_lhs", [])
        sampled_molecules = entry.get("sampled_molecules", [])
        target_smiles = entry.get("target_smiles", "")

        repeated_index = [index] * len(target_smiles)

        indices.extend(repeated_index)
        log_lhs_list.extend(log_lhs)
        target_smiles_list.extend(target_smiles)
        sampled_molecules_list.extend(sampled_molecules)

    df = pd.DataFrame(
        {
            "index": indices,
            "log_lhs": log_lhs_list,
            "target_smiles": target_smiles_list,
            "sampled_molecules": sampled_molecules_list,
        }
    )

    df.to_csv(csv_file, index=False)


def present_result(input1, input2, output):
    test_df = pd.read_csv(input1)
    prediction_df = pd.read_csv(input2)

    def agg_order_preserved(series):
        return list(series)

    grouped = (
        test_df.groupby(["parent_smiles"], sort=False)
        .agg(
            {
                "parent_name": "first",  # Get the first entry for 'parent_name' in each group
                "child_name": agg_order_preserved,
                "child_smiles": agg_order_preserved,
            }
        )
        .reset_index()
    )

    if len(grouped) == len(prediction_df):
        grouped["sampled_molecules"] = prediction_df["sampled_molecules"]
        grouped["log_lhs"] = prediction_df["log_lhs"]
    else:
        raise ValueError(
            "The number of rows in 'grouped' does not match the number of rows in 'prediction_df'."
        )

    grouped.insert(0, "index", prediction_df["index"])

    grouped.to_csv(output, index=False)


def calculate_fingerprint_similarity(parent, sampled):
    parent_mol = Chem.MolFromSmiles(parent)
    sampled_mol = Chem.MolFromSmiles(sampled)

    parent_fps = AllChem.GetMorganFingerprintAsBitVect(parent_mol, radius=2, nBits=1024)
    sampled_fps = AllChem.GetMorganFingerprintAsBitVect(
        sampled_mol, radius=2, nBits=1024
    )

    fingerprint_similarity = DataStructs.TanimotoSimilarity(parent_fps, sampled_fps)

    return fingerprint_similarity


def save_valid_smiles(input_file, batches):
    df = pd.read_csv(input_file)

    main_df = pd.DataFrame()

    validity = []
    mean_validity = []
    for i in range(batches):
        batch_df = df.loc[df["index"] == i]
        batch_df.reset_index(drop=True, inplace=True)

        parent_smiles = batch_df["parent_smiles"]

        total_valid_smiles = 0
        total_predictions = 0
        mean_predictions = 0
        for i in range(len(parent_smiles)):
            sampled_molecules_i = ast.literal_eval(batch_df.at[i, "sampled_molecules"])
            log_lhs_i = ast.literal_eval(batch_df.at[i, "log_lhs"])
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

            print("nr of dup: ", count_dup)
            print("nr of parent dup: ", count_parent)
            print("nr of valid smiles left: ", len(valid_smiles))

            mean_predictions += len(valid_smiles)

            batch_df.at[i, "sampled_molecules"] = valid_smiles
            batch_df.at[i, "log_lhs"] = valid_log_lhs

        validity.append(total_valid_smiles / total_predictions)
        mean_validity.append(mean_predictions / len(parent_smiles))

        main_df = pd.concat([main_df, batch_df], ignore_index=True)

    print(validity)
    validity_mean, validity_var = mean_and_variance(validity)
    mean_validity_mean, mean_validity_var = mean_and_variance(mean_validity)
    print(f"\nValidity: {validity_mean:.3f} $\pm$ {validity_var:.3f}")
    print(
        f"Mean number of SMILES per drug: {mean_validity_mean:.3f} $\pm$ {mean_validity_var:.3f} / {len(sampled_molecules_i)}"
    )

    main_df.to_csv(input_file, index=False)


def specify_df(df, specification=None):
    parent_name = df["parent_name"]
    child_smiles = df["child_smiles"]

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


def concat_multiple_predictions(input_files, output, n=None):
    dataframes = [pd.read_csv(input_file) for input_file in input_files]

    combined_sampled_molecules = []

    parent_name = dataframes[0]["parent_name"]

    for i in range(len(parent_name)):
        combined_sampled_molecules_i = set()

        for df in dataframes:
            sampled_molecules_i = ast.literal_eval(df["sampled_molecules"][i])

            if n is not None:
                sampled_molecules_i = sampled_molecules_i[:n]

            combined_sampled_molecules_i.update(sampled_molecules_i)

        combined_sampled_molecules.append(list(combined_sampled_molecules_i))

    result_df = dataframes[0].copy()
    result_df["sampled_molecules"] = combined_sampled_molecules

    result_df.to_csv(output, index=False)


def count_metabolites(input_file):
    df = pd.read_csv(input_file)

    def count_molecules(row):
        try:
            molecules_list = ast.literal_eval(row)
            return len(molecules_list)
        except Exception as e:
            print(f"Error converting row to list: {e}")
            return 0

    df["molecule_count"] = df["sampled_molecules"].apply(count_molecules)

    average_molecules = round(df["molecule_count"].mean(), 2)

    print(f"The average number of sampled molecules per drug is: {average_molecules}")


def count_correct_metabolites(input_file, batch, specification, fingerprint=None):
    df = pd.read_csv(input_file)

    df = specify_df(df, specification)

    if batch is not None:
        df = df.loc[df["index"] == batch]
        df.reset_index(drop=True, inplace=True)

    parent_name = df["parent_name"]
    child_smiles = df["child_smiles"]
    sampled_molecules = df["sampled_molecules"]

    top1 = []
    top1_pred = []
    top3 = []
    top3_pred = []
    top5 = []
    top5_pred = []
    top10 = []
    top10_pred = []
    all = []
    all_pred = []
    reference = []

    for i in range(len(parent_name)):
        sampled_molecules_i = ast.literal_eval(sampled_molecules[i])
        top1_pred.append(len(sampled_molecules_i[0:1]))
        top3_pred.append(len(sampled_molecules_i[0:3]))
        top5_pred.append(len(sampled_molecules_i[0:5]))
        top10_pred.append(len(sampled_molecules_i[0:10]))
        child_smiles_i = ast.literal_eval(child_smiles[i])

        count_top1 = 0
        count_top3 = 0
        count_top5 = 0
        count_top10 = 0
        count_all = 0
        for j in range(len(child_smiles_i)):
            for k in range(len(sampled_molecules_i)):
                if fingerprint is not None:
                    if (
                        calculate_fingerprint_similarity(
                            child_smiles_i[j], sampled_molecules_i[k]
                        )
                        >= fingerprint
                    ):
                        count_all += 1
                        if k < 10:
                            count_top10 += 1
                        if k < 5:
                            count_top5 += 1
                        if k < 3:
                            count_top3 += 1
                        if k < 1:
                            count_top1 += 1
                else:
                    if child_smiles_i[j] == sampled_molecules_i[k]:
                        count_all += 1
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
        all_pred.append(len(sampled_molecules_i))
        reference.append(len(child_smiles_i))

    return (
        top1,
        top1_pred,
        top3,
        top3_pred,
        top5,
        top5_pred,
        top10,
        top10_pred,
        all,
        all_pred,
        reference,
    )


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
        if top1[i] >= math.ceil(0.5 * reference[i]):
            score1 += 1
    score1 = score1 / len(top1)

    score3 = 0
    for i in range(len(top3)):
        if top3[i] >= math.ceil(0.5 * reference[i]):
            score3 += 1
    score3 = score3 / len(top3)

    score5 = 0
    for i in range(len(top5)):
        if top5[i] >= math.ceil(0.5 * reference[i]):
            score5 += 1
    score5 = score5 / len(top5)

    score10 = 0
    for i in range(len(top10)):
        if top10[i] >= math.ceil(0.5 * reference[i]):
            score10 += 1
    score10 = score10 / len(top10)

    score_all = 0
    for i in range(len(all)):
        if all[i] >= math.ceil(0.5 * reference[i]):
            score_all += 1
    score_all = score_all / len(all)

    return score1, score3, score5, score10, score_all


def all_metabolites(top1, top3, top5, top10, all, reference):
    score1 = 0
    for i in range(len(top1)):
        if top1[i] == reference[i]:
            score1 += 1
    score1 = score1 / len(top1)

    score3 = 0
    for i in range(len(top3)):
        if top3[i] == reference[i]:
            score3 += 1
    score3 = score3 / len(top3)

    score5 = 0
    for i in range(len(top5)):
        if top5[i] == reference[i]:
            score5 += 1
    score5 = score5 / len(top5)

    score10 = 0
    for i in range(len(top10)):
        if top10[i] == reference[i]:
            score10 += 1
    score10 = score10 / len(top10)

    score_all = 0
    for i in range(len(all)):
        if all[i] == reference[i]:
            score_all += 1
    score_all = score_all / len(all)

    return score1, score3, score5, score10, score_all


def precision_score(input_file):
    df = pd.read_csv(input_file)

    parent_name = df["parent_name"]
    child_smiles = df["child_smiles"]
    sampled_molecules = df["sampled_molecules"]

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


def recall_and_precision(correct, predictions, reference):
    TP = sum(correct)
    FP = sum(predictions) - TP
    FN = sum(reference) - TP

    recall = TP / (TP + FN)
    precision = TP / (TP + FP)
    F1 = 2 * (TP) / (2 * TP + FP + FN)

    return recall, precision, F1


def mean_and_variance(list):
    np_list = np.array(list)
    mean = np.mean(np_list)
    variance = np.var(np_list, ddof=1)  # ddof=1 for sample variance

    return mean, variance


def score_result(
    input_file, batches, print_result=False, fingerprint=None, specification=0
):
    (
        top1,
        top1_pred,
        top3,
        top3_pred,
        top5,
        top5_pred,
        top10,
        top10_pred,
        all,
        all_pred,
        reference,
    ) = count_correct_metabolites(input_file, None, specification, fingerprint)

    if print_result:
        print("\t")
        print(f"Total identified metabolites: {sum(all)} / {sum(reference)}")
        print(f"Total number of predictions: {sum(all_pred)}")
        print(f"Mean number of predictions per drug: {np.mean(all_pred)}")

        print("\t")
        print(
            f"Total identified metabolites in top-10: {sum(top10)} / {sum(reference)}"
        )
        print(f"Total number of predictions in top-10: {sum(top10_pred)}")
        print(f"Mean number of predictions per drug in top-10: {np.mean(top10_pred)}")

    recall1_list = []
    precision1_list = []
    recall3_list = []
    precision3_list = []
    recall5_list = []
    precision5_list = []
    recall10_list = []
    precision10_list = []
    f1_10_list = []
    recall_all_list = []
    precision_all_list = []
    f1_all_list = []
    score1_one_list = []
    score3_one_list = []
    score5_one_list = []
    score10_one_list = []
    score_all_one_list = []
    score1_all_list = []
    score3_all_list = []
    score5_all_list = []
    score10_all_list = []
    score_all_all_list = []
    for i in range(batches):
        (
            top1,
            top1_pred,
            top3,
            top3_pred,
            top5,
            top5_pred,
            top10,
            top10_pred,
            all,
            all_pred,
            reference,
        ) = count_correct_metabolites(input_file, i, specification, fingerprint)
        score1_one, score3_one, score5_one, score10_one, score_all_one = (
            at_least_one_metabolite(top1, top3, top5, top10, all, reference)
        )
        score1_all, score3_all, score5_all, score10_all, score_all_all = (
            all_metabolites(top1, top3, top5, top10, all, reference)
        )

        score1_one_list.append(score1_one)
        score3_one_list.append(score3_one)
        score5_one_list.append(score5_one)
        score10_one_list.append(score10_one)
        score_all_one_list.append(score_all_one)

        score1_all_list.append(score1_all)
        score3_all_list.append(score3_all)
        score5_all_list.append(score5_all)
        score10_all_list.append(score10_all)
        score_all_all_list.append(score_all_all)

        recall1, precision1, _ = recall_and_precision(top1, top1_pred, reference)
        recall1_list.append(recall1)
        precision1_list.append(precision1)

        recall3, precision3, _ = recall_and_precision(top3, top3_pred, reference)
        recall3_list.append(recall3)
        precision3_list.append(precision3)

        recall5, precision5, _ = recall_and_precision(top5, top5_pred, reference)
        recall5_list.append(recall5)
        precision5_list.append(precision5)

        recall10, precision10, f1_10 = recall_and_precision(
            top10, top10_pred, reference
        )
        recall10_list.append(recall10)
        precision10_list.append(precision10)
        f1_10_list.append(f1_10)

        recall_all, precision_all, f1_all = recall_and_precision(
            all, all_pred, reference
        )
        recall_all_list.append(recall_all)
        precision_all_list.append(precision_all)
        f1_all_list.append(f1_all)

    if print_result:
        score1_one_mean, score1_one_var = mean_and_variance(score1_one_list)
        score10_one_mean, score10_one_var = mean_and_variance(score10_one_list)

        score10_all_mean, score10_all_var = mean_and_variance(score10_all_list)

        score_all_one_mean, score_all_one_var = mean_and_variance(score_all_one_list)
        score_all_all_mean, score_all_all_var = mean_and_variance(score_all_all_list)

        recall10_mean, recall10_var = mean_and_variance(recall10_list)
        precision10_mean, precision10_var = mean_and_variance(precision10_list)
        f1_10_mean, f1_10_var = mean_and_variance(f1_10_list)

        recall_all_mean, recall_all_var = mean_and_variance(recall_all_list)
        precision_all_mean, precision_all_var = mean_and_variance(precision_all_list)
        f1_all_mean, f1_all_var = mean_and_variance(f1_all_list)

        print("\nAt least one metabolite: ")
        print(f"Score1: {score1_one_mean:.2f} $\pm$ {score1_one_var:.2f}")
        print(f"Score10: {score10_one_mean:.2f} $\pm$ {score10_one_var:.2f}")
        print(f"Score all: {score_all_one_mean:.2f} $\pm$ {score_all_one_var:.2f}")
        print("\nAll metabolites: ")
        print(f"Score10: {score10_all_mean:.2f} $\pm$ {score10_all_var:.2f}")
        print(f"Score all: {score_all_all_mean:.2f} $\pm$ {score_all_all_var:.2f}")
        print("\t")
        print(f"Precision @ 10: {precision10_mean:.2f} $\pm$ {precision10_var:.2f}")
        print(f"Recall @ 10: {recall10_mean:.2f} $\pm$ {recall10_var:.2f}")
        print(f"F1 @ 10: {f1_10_mean:.2f} $\pm$ {f1_10_var:.2f}")
        print(
            f"Precision @ all: {precision_all_mean:.2f} $\pm$ {precision_all_var:.2f}"
        )
        print(f"Recall @ all: {recall_all_mean:.2f} $\pm$ {recall_all_var:.2f}")
        print(f"F1 @ all: {f1_all_mean:.2f} $\pm$ {f1_all_var:.2f}")

    return (
        [recall1_list, recall3_list, recall5_list, recall10_list, recall_all_list],
        [
            precision1_list,
            precision3_list,
            precision5_list,
            precision10_list,
            precision_all_list,
        ],
        [
            score1_one_list,
            score3_one_list,
            score5_one_list,
            score10_one_list,
            score_all_one_list,
        ],
        [
            score1_all_list,
            score3_all_list,
            score5_all_list,
            score10_all_list,
            score_all_all_list,
        ],
    )


if __name__ == "__main__":
    benchmark = False  # True if GLORYx, False if Evaluation
    status = "score"  # 'score' 'combine' 'new'
    name = "chemVA-Met_base_05"

    specification = 0  # 0 (all) 1 (only_child) 2 (more than 1) 3 (more than 2)
    fingerprint = (
        None  # 1 (similarity = 1) 0.8 (similarity >= 0.8) None (exact SMILES string)
    )

    # If combine:
    ensemble_list = [
        "evaluation/scores/ensemble/result_ChemVA-Met_Ssplit1.csv",
        "evaluation/scores/ensemble/result_ChemVA-Met_Ssplit2.csv",
        "evaluation/scores/ensemble/result_ChemVA-Met_Ssplit3.csv",
        "evaluation/scores/ensemble/result_ChemVA-Met_Ssplit4.csv",
    ]
    samples_per_model = 5

    if benchmark:
        bs = 1
        testset = "dataset/curated_data/gloryx_smiles_clean.csv"
    else:
        bs = 4
        testset = "dataset/curated_data/LAGOM_evaluation.csv"

    json_predictions = "results/evaluation/scores/predictions0.json"
    csv_predictions = f"evaluation/scores/predictions_{name}.csv"
    csv_result = f"evaluation/scores/result_{name}.csv"

    if status == "new":
        json_to_csv(json_predictions, csv_predictions)
        present_result(testset, csv_predictions, csv_result)
        save_valid_smiles(csv_result, bs)
        score_result(csv_result, bs, True, fingerprint, specification)

    elif status == "score":
        present_result(testset, csv_predictions, csv_result)
        save_valid_smiles(csv_result, bs)
        score_result(csv_result, bs, True, fingerprint, specification)

    elif status == "combine":
        csv_comb = f"evaluation/scores/result_{name}_{samples_per_model}_per_model.csv"
        concat_multiple_predictions(ensemble_list, csv_comb, samples_per_model)
        score_result(csv_comb, bs, True, fingerprint, specification)
        count_metabolites(csv_comb)
    else:
        print("Wrong status")
