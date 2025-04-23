import pandas as pd

csv = 'dataset/curated_data/combined_smiles_clean.csv'
df = pd.read_csv(csv)

# Create the new 'ReactionSmilesClean' column
df['ReactionSmilesClean'] = df.apply(lambda row: f"{row['parent_smiles']}>>{row['child_smiles']}", axis=1)

output_csv = 'preprocessing/transformations/reactions_smiles_clean.csv'  # Replace with your desired output file path

# Save the DataFrame to a new tab-separated CSV
df.to_csv(output_csv, sep='\t', index=False, columns=['ReactionSmilesClean'])

print(f"New file saved to {output_csv}")
