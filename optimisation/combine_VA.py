import pandas as pd


def combine_filtered_VA_files(combined_VA_file):
    # Base directory for storing or referencing the CSV files
    base_dir = "../dataset/VA"

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
        (9911735, "end"),
    ]

    dataframes = []

    for start_row, end_row in row_ranges:
        name = f"VA_{start_row}_{end_row}"
        file_path = f"{base_dir}/{name}_smiles_clean.csv"

        try:
            df = pd.read_csv(file_path)
            dataframes.append(df)
        except FileNotFoundError:
            print(f"File not found: {file_path}")
            continue

    # Concatenate all DataFrames in the list
    df_combined = pd.concat(dataframes, ignore_index=True)

    df_combined.to_csv(combined_VA_file, index=False)


combined_VA_file = "dataset/VA/VA_smiles_clean.csv"
combine_filtered_VA_files(combined_VA_file)
