import pandas as pd


def add_possible_products(input_file):
    df = pd.read_csv(input_file, sep='\t') 

    def agg_order_preserved(series):
        return list(series)

    grouped = df.groupby(['reactants'], sort=False).agg({
        'parent_name': 'first',  # Get the first entry for 'parent_name' in each group
        'child_name': agg_order_preserved, 
        'products': agg_order_preserved 
    }).reset_index()

    grouped.rename(columns={'products': 'possible_products'}, inplace=True)

    df = df.merge(grouped[['reactants', 'possible_products']], on='reactants', how='left')

    df['possible_products'] = df['possible_products'].apply(lambda x: '.'.join(x))

    df.to_csv(input_file, sep='\t', index=False)


file = 'dataset/finetune/combined_finetune.csv'
add_possible_products(file)
