import pandas as pd

def create_csv(input1,input2,output):

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

    

    # Write the result to a new CSV file
    grouped.to_csv(output, index=False)


gloryx = 'dataset/curated_data/gloryx_smiles_clean.csv'
pred = 'evaluation/predictions/prediction_03-03_unique.csv'
out = 'test.csv'

create_csv(gloryx, pred, out)