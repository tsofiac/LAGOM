import numpy as np
from rdkit import Chem
from rdkit.Chem.SaltRemover import SaltRemover
from chembl_structure_pipeline.standardizer import standardize_mol
import pandas as pd


def standardize_molecule(smiles, isomericSmiles=True):
    
    '''
    Standardization of SMILES by standardiser;
    '''
    
    try:
        rm = SaltRemover()
        mol = Chem.MolFromSmiles(smiles)
        mol = rm.StripMol(mol, dontRemoveEverything=True)
        if len(Chem.MolToSmiles(mol).split('.')) > 1:
            salt = Chem.MolToSmiles(mol)
            frag_dict = {len(k): k for k in salt.split('.')}
            max_frag = frag_dict[max(list(frag_dict.keys()))]
            mol = Chem.MolFromSmiles(max_frag)
        mol = standardize_mol(mol)
        clean_smiles = Chem.MolToSmiles(mol, isomericSmiles=isomericSmiles, kekuleSmiles=True) #removing (False) or keeping (True) stereochemistry, kekuleSmiles added
        if Chem.MolFromSmiles(clean_smiles) is None:
            clean_smiles = None  
    except:
        clean_smiles = np.nan
    return clean_smiles

def standardize_smiles_collection(smiles_list, isomericSmiles=True):
    
    '''
    Standardization of SMILES collection by standardiser;
    '''
    
    lookup = {}
    return_smiles = []
    for smiles in smiles_list:
        if smiles in lookup:
            standardized_smiles = lookup[smiles]
        else:
            standardized_smiles = standardize_molecule(smiles, isomericSmiles=isomericSmiles)
            lookup[smiles] = standardized_smiles
        return_smiles.append(standardized_smiles)
        
    return return_smiles


def standardize_smiles(input_file): #input df that has columns 'parent_smiles' and 'child_smiles'
    df = pd.read_csv(input_file)
    if 'parent_smiles' in df.columns:
        df['parent_smiles'] = standardize_smiles_collection(df['parent_smiles'], False) #'False' eliminates isomeres
        df['child_smiles'] = standardize_smiles_collection(df['child_smiles'], False)
    else:
        df['reactants'] = standardize_smiles_collection(df['reactants'], False) #'False' eliminates isomeres
        df['products'] = standardize_smiles_collection(df['products'], False)

    df.to_csv(input_file, index=False)



def standardize_smiles_tab(input_file): #input df that has columns 'parent_smiles' and 'child_smiles'
    df = pd.read_csv(input_file, sep='\t')
    if 'parent_smiles' in df.columns:
        df['parent_smiles'] = standardize_smiles_collection(df['parent_smiles'], False) #'False' eliminates isomeres
        df['child_smiles'] = standardize_smiles_collection(df['child_smiles'], False)
    else:
        df['reactants'] = standardize_smiles_collection(df['reactants'], False) #'False' eliminates isomeres
        df['products'] = standardize_smiles_collection(df['products'], False)

    df.to_csv(input_file, sep='\t', index=False)


combined = 'dataset/curated_data/combined_smiles_clean.csv'
combined_eval = 'dataset/curated_data/combined_evaluation.csv'
combined_eval_unique = 'dataset/curated_data/combined_evaluation_unique.csv'

combine_finetune = 'dataset/finetune/combined_finetune.csv'
combined_eval_finetune = 'dataset/finetune/combined_evaluation_finetune.csv'

standardize_smiles(combined)
# standardize_smiles(combined_eval)
# standardize_smiles(combined_eval_unique)

# standardize_smiles_tab(combine_finetune)
# standardize_smiles_tab(combined_eval_finetune)

