import pandas as pd

# Load the MetXBio CSV
metxbiodb_csv = 'dataset/curated_data/metxbiodb_smiles.csv'
metxbiodb_df = pd.read_csv(metxbiodb_csv)

# Add a source column with the value 'MetXBioDB'
metxbiodb_df['source'] = 'MetXBioDB'

# Load the DrugBank CSV
drugbank_csv = 'dataset/curated_data/drugbank_smiles.csv'
drugbank_df = pd.read_csv(drugbank_csv)

# Select only the relevant columns from the DrugBank DataFrame
drugbank_selected_df = drugbank_df[['parent_name', 'parent_smiles', 'child_name', 'child_smiles']].copy()
drugbank_selected_df['source'] = 'DrugBank'

# Concatenate the two DataFrames
combined_df = pd.concat([metxbiodb_df, drugbank_selected_df], ignore_index=True)

# Save the combined DataFrame to a new CSV file
combined_df.to_csv('dataset/curated_data/combined_smiles.csv', index=False)

print("Combined CSV file created successfully!")