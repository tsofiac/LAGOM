
import pandas as pd

csv = 'dataset/finetune/combined_finetune.csv'

df = pd.read_csv(csv, sep='\t')

# Check which rows in the 'set' column are empty
empty_rows = df[df['set'].isnull()]

# Display the empty rows
print(empty_rows)


#513 rows