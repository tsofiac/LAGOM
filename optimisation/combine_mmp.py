import pandas as pd

def combine_filtered_mmp_files(combined_mmp_file):
    # Base directory for storing or referencing the CSV files
    base_dir = '/projects/cc/se_users/carlsson_ksmq649/MasterThesis/dataset/mmp'
    
    # Define row ranges that are to be used for naming
    row_ranges = [
        (0, 1101304),
        (1101304, 2202607),
        (2202607, 3303912),
        (3303912, 4405216),
        (4405216, 5506519),
        (5506519, 6607824),
        (6607824, 7709127),
        (7709127, 8810431),
        (8810431, 9911735),
        (9911735, 'end')
    ]
    
    # List to hold the DataFrames
    dataframes = []

    for start_row, end_row in row_ranges:
        # Construct file name
        name = f'mmp_{start_row}_{end_row}'
        file_path = f'{base_dir}/{name}_smiles_clean.csv'
        
        # Read each CSV as a DataFrame and append to the list
        try:
            df = pd.read_csv(file_path)
            dataframes.append(df)
        except FileNotFoundError:
            print(f"File not found: {file_path}")
            continue
    
    # Concatenate all DataFrames in the list
    df_combined = pd.concat(dataframes, ignore_index=True)
    
    # Save the combined DataFrame to a new CSV file
    df_combined.to_csv(combined_mmp_file, index=False)

combined_mmp_file = 'dataset/mmp/mmp_all_smiles_clean.csv' 
combine_filtered_mmp_files(combined_mmp_file)