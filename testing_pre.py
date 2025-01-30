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
df = pd.DataFrame({
    'products': products,
    'reactants': reactants,
    'set': set_column
})

# Write the DataFrame to a TSV file (tab-separated values)
df.to_csv('gloryx_reference_dataset.csv', sep='\t', index=False)