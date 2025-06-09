
import pandas as pd


def compare_datasets(gloryx_csv, metatrans_csv, removed_file, gloryx_final_csv):
    df1 = pd.read_csv(gloryx_csv)
    df2 = pd.read_csv(metatrans_csv)
    # duplicate parent_smiles between df1 and df2
    duplicate_parent_smiles = set(df1['parent_smiles']).intersection(set(df2['parent_smiles']))
    # those not in the duplicates set
    non_duplicates_in_df1 = df1[~df1['parent_smiles'].isin(duplicate_parent_smiles)].copy()
    duplicates_df = df1[df1['parent_smiles'].isin(duplicate_parent_smiles)].copy()

    print('Number of overlapping reactions removed: ', len(df1)-len(non_duplicates_in_df1))
    non_duplicates_in_df1.to_csv(gloryx_final_csv, index=False)
    duplicates_df.to_csv(removed_file, index=False)

def get_unique_parents(file, tab_file):
    df = pd.read_csv(file)
    df = df.drop_duplicates(subset="parent_smiles", keep="first")  # Only unique parent

    df = df.rename(
        columns={
            "child_smiles": "products",
            "parent_smiles": "reactants",
        }
    )

    df.to_csv(tab_file, sep="\t", index=False)

gloryx_csv = 'dataset/curated_data/gloryx_smiles_clean.csv'
metatrans_csv = 'dataset/extracted_data/metatrans_smiles.csv'
removed_file = 'for_article/gloryx_metatrans_overlap.csv' 
gloryx_final_csv = 'for_article/gloryx_removed_metatrans.csv' 
gloryx_finetune_csv = 'for_article/gloryx_article_finetune.csv'

print("Comparing datasets to remove overlapping reactions from GLORYx")
compare_datasets(gloryx_csv, metatrans_csv, removed_file, gloryx_final_csv)
print("Creating Chemformer-friendly file")
get_unique_parents(gloryx_final_csv, gloryx_finetune_csv)