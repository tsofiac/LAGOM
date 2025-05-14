import pandas as pd

### See if parents and metabolites are same ###

def count_matching_smiles(data_csv):
    # Load the CSV file into a DataFrame
    df = pd.read_csv(data_csv)

    num_matching_smiles = 0
    matching_smiles = []

    for index, row in df.iterrows():
        if row['parent_smiles'] == row['child_smiles']:
            num_matching_smiles += 1
            matching_smiles.append(row)

    print('Number of equal parent and child:', num_matching_smiles)
    print(matching_smiles)
    
print("DRUGBANK")
count_matching_smiles('dataset/curated_data/drugbank_clean_unique_parents.csv')

print()
print("METXBIODB")
count_matching_smiles('dataset/preprocessed_metxbiodb/metxbiodb_clean_unique_parents.csv')
