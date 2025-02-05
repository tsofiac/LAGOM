import json
import pandas as pd

# Load the JSON data from the file
with open('gloryx_reference_dataset.json', 'r') as json_file:
    data = json.load(json_file)

# Prepare lists to hold the extracted data
products = []
reactants = []
set_column = []

# Iterate over each entry in the JSON data
for entry in data:
    parent_smiles = entry['Parent molecule']['SMILES']
  
    # Iterate over metabolites for each parent molecule
    for metabolite in entry['Metabolites']:
        metabolite_smiles = metabolite['SMILES']
  
        # Append the SMILES and 'train' set to the lists
        products.append(metabolite_smiles)
        reactants.append(parent_smiles)
        set_column.append('train')

# Calculate the cutoff index for the last 10% to be labeled as 'val'
cutoff_index = int(len(products) * 0.9)

# Update the last 10% of the set_column to 'val'
set_column[cutoff_index:] = ['val'] * (len(products) - cutoff_index)

# Create a DataFrame with the collected data
df1 = pd.DataFrame({
    'products': products,
    'reactants': reactants,
    'set': set_column
})

# Write the DataFrame to a TSV file (tab-separated values)
df1.to_csv('10percent_gloryx_reference_dataset.csv', sep='\t', index=False)

# Create df2 for the half-and-half dataset
half_index = len(products) // 2

# Select half of the products and reactants for 'train'
half_train_products = products[:half_index]
half_train_reactants = reactants[:half_index]
half_train_set = ['train'] * half_index

# Duplicate the same for 'val'
half_val_products = products[:half_index]
half_val_reactants = reactants[:half_index]
half_val_set = ['val'] * half_index

# Combine the halves to form the full DataFrame
df2 = pd.DataFrame({
    'products': half_train_products + half_val_products,
    'reactants': half_train_reactants + half_val_reactants,
    'set': half_train_set + half_val_set
})

# Write the second DataFrame to a TSV file (tab-separated)
df2.to_csv('half_half_gloryx_reference_dataset.csv', sep='\t', index=False)
