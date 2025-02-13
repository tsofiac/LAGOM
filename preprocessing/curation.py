from rdkit import Chem
import pandas as pd
from sklearn.model_selection import GroupShuffleSplit
from sklearn.model_selection import train_test_split
from standardize_smiles import standardize_smiles_collection


def standardize_smiles_main(input_file,output_file): #input df that has columns 'parent_smiles' and 'child_smiles'
    df = pd.read_csv(input_file)
    df['parent_smiles'] = standardize_smiles_collection(df['parent_smiles'], False) #'False' eliminates isomeres
    df['child_smiles'] = standardize_smiles_collection(df['child_smiles'], False)
    df.to_csv(output_file, index=False)


if __name__ == "__main__":


    name = 'drugbank' # [ 'drugbank' 'metxbiodb' ]

    dataset = f'dataset/curated_data/{name}_smiles.csv'

    standard = f'dataset/curated_data/{name}_smiles_standard.csv'

    standardize_smiles_main(dataset,standard)
