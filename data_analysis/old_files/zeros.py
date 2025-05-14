import pandas as pd
import csv

def find_datapoints(zero_indices, predictions_csv, zeros_csv):
    # Load the predictions CSV into a DataFrame
    predictions_df = pd.read_csv(predictions_csv)

    # Filter the DataFrame to get rows corresponding to the zero_indices
    filtered_df = predictions_df.iloc[zero_indices]
    
    # filtered_df = filtered_df[['parent_smiles', 'child_smiles']]

    # Save the new DataFrame to the zeros CSV file
    filtered_df.to_csv(zeros_csv, index=False)
     

predictions_csv = 'dataset/curated_data/combined_evaluation_unique.csv'
zeros_csv = 'dataset/curated_data/zero_predict.csv'

def present_result(input1,input2, output):

    test_df = pd.read_csv(input1)
    prediction_df = pd.read_csv(input2)

    def agg_order_preserved(series):
        return list(series)

    grouped = test_df.groupby(['parent_smiles'], sort=False).agg({
        'parent_name': 'first',  # Get the first entry for 'parent_name' in each group
        'child_name': agg_order_preserved, 
        'child_smiles': agg_order_preserved,
        'enzymes': agg_order_preserved 
    }).reset_index()

    if len(grouped) == len(prediction_df):
        grouped['sampled_molecules'] = prediction_df['sampled_molecules']
    else:
        raise ValueError("The number of rows in 'grouped' does not match the number of rows in 'prediction_df'.")

    grouped.to_csv(output, index=False)



#zero_indices = allzero(v1, v2, v3, v4)
#find_datapoints(zero_indices, predictions_csv, zeros_csv)

#testset = 'dataset/curated_data/combined_evaluation.csv'
#present_result(testset, 'dataset/curated_data/zero_predict.csv','dataset/curated_data/zero_predict_full.csv') #funkar ej

zeros_csv = 'dataset/curated_data/zero_predict.csv'
evaluation_csv = 'dataset/curated_data/combined_evaluation.csv'
output_csv = 'dataset/curated_data/zero_predict_all_children.csv'


def find_all_children(zeros_csv, evaluation_csv, ouput_csv):
    zeros_df = pd.read_csv(zeros_csv)
    evaluation_df = pd.read_csv(evaluation_csv)
    output_df = pd.DataFrame()

    for index, row in zeros_df.iterrows():
        parent_smiles_zero = row['parent_smiles']

        # Filter evaluation_df where parent_smiles equals parent_smiles_zero
        matched_rows = evaluation_df[evaluation_df['parent_smiles'] == parent_smiles_zero]

        # Append matched rows to output_df
        output_df = pd.concat([output_df, matched_rows], ignore_index=True)

    output_df.to_csv(output_csv, index = False)

# find_all_children(zeros_csv, evaluation_csv, output_csv)

def number_of_children(csv):
    # Read the CSV into a DataFrame
    df = pd.read_csv(csv)

    # Count how many times each 'parent_smiles' occurs
    child_counts = df['parent_smiles'].value_counts()

    # Create a summary of occurrence frequencies
    occurrence_summary = child_counts.value_counts().sort_index()

    # Print the results
    print("Number of children summary:")
    for number_of_children, frequency in occurrence_summary.items():
        print(f"There are {frequency} parents with {number_of_children} child(ren).")

# csv_to_count = 'dataset/curated_data/zero_predict_all_children.csv'
# csv_to_count = 'dataset/curated_data/combined_evaluation.csv'
csv_to_count = 'dataset/curated_data/combined_smiles_clean.csv'

number_of_children(csv_to_count)