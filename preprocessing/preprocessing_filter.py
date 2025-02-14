from rdkit import Chem, DataStructs
from rdkit.Chem import Descriptors, AllChem 
from rdkit.Chem.SaltRemover import SaltRemover
import torch
import pandas as pd
import os
from pathlib import Path
from collections import Counter
from standardize_smiles import standardize_smiles_collection


def standardize_smiles(input_file,output_file): #input df that has columns 'parent_smiles' and 'child_smiles'
    df = pd.read_csv(input_file)
    df['parent_smiles'] = standardize_smiles_collection(df['parent_smiles'], False) #'False' eliminates isomeres
    df['child_smiles'] = standardize_smiles_collection(df['child_smiles'], False)
    df.to_csv(output_file, index=False)

# saving removed data
def filter_data_on_both_sides(data_file, filter_method, removed_data_file): 
    data = pd.read_csv(data_file)
    total_removed = 0

    # Filtering based on child molecules
    allowed_children = [filter_method(metabolite) for metabolite in data["child_smiles"]]
    filtered_data = data[allowed_children]
    removed_data_children = data[[not allowed for allowed in allowed_children]]
    total_removed += allowed_children.count(False)

    # Filtering based on parent molecules
    allowed_parents = [filter_method(molecule) for molecule in filtered_data["parent_smiles"]]
    removed_data_parents = filtered_data[[not allowed for allowed in allowed_parents]]
    filtered_data = filtered_data[allowed_parents]
    total_removed += allowed_parents.count(False)

    removed_data = pd.concat([removed_data_children, removed_data_parents])

    filtered_data.to_csv(data_file, index=False)
    removed_data.to_csv(removed_data_file, index=False)

    print(f"Total data points removed: {total_removed}")


def filter_data_on_one_side(data_file, filter_method, removed_data_file, if_parent = True): 
    name_property = "parent_smiles" if if_parent else "child_smiles"
    data = pd.read_csv(data_file)

    total_removed = 0

    # Filtering based on parent molecules
    allowed_molecules = [filter_method(molecule) for molecule in data[name_property]]
    filtered_data = data[allowed_molecules]
    removed_data = data[[not allowed for allowed in allowed_molecules]]
    total_removed += allowed_molecules.count(False)

    filtered_data.to_csv(data_file, index=False)
    removed_data.to_csv(removed_data_file, index=False)

    print(f"Total data points removed: {total_removed}")

# filtering methods
# ----------------------------------------------------------------
def valid_smiles(molecule): 
    return Chem.MolFromSmiles(molecule) is not None

def atoms_allowed_in_molecules(molecule): 
    atoms_to_include = ['C', 'N', 'S', 'O', 'H', 'F', 'I', 'P', 'Cl', 'Br'] # B and Si removed in contrast to last year
    mol = Chem.MolFromSmiles(molecule)
    atoms = [atom.GetSymbol() for atom in mol.GetAtoms()]
    return set(atoms).issubset(set(atoms_to_include)) #Returns true if all atoms in the molecule are included in the specified set

def molecule_allowed_based_on_weight(molecule, max_weight=700, min_weight=100): 
    mol_weight = Descriptors.ExactMolWt(Chem.MolFromSmiles(molecule))
    if mol_weight <= max_weight and mol_weight >= min_weight: 
        return True
    return False 
# ----------------------------------------------------------------

def define_finger_print_similarity(dataset):
    parent_smiles = dataset['parent_smiles'].tolist()
    child_smiles = dataset['child_smiles'].tolist()

    parent_mol = [Chem.MolFromSmiles(x) for x in parent_smiles]
    child_mol = [Chem.MolFromSmiles(x) for x in child_smiles]

    parent_fps = [AllChem.GetMorganFingerprintAsBitVect(x, radius=2, nBits=1024) for x in parent_mol]
    child_fps = [AllChem.GetMorganFingerprintAsBitVect(x, radius=2, nBits=1024) for x in child_mol]

    fingerprint_similarities = []
    for i in range(len(parent_smiles)):
        
        s = DataStructs.TanimotoSimilarity(parent_fps[i], child_fps[i]) 
        fingerprint_similarities.append(s)

    dataset['tanimoto'] = fingerprint_similarities

    return dataset

def fingerprints_allowed(similarity, min_similarity):
    if similarity < min_similarity or similarity == 1:
        return False
    return True

def finger_print_similarity_filter(data_file, removed_data_file, min_similarity = 0.15):
    data = pd.read_csv(data_file)
    data = define_finger_print_similarity(data)

    total_removed = 0

    allowed_molecules = [fingerprints_allowed(value, min_similarity) for value in data["tanimoto"]]
    filtered_data = data[allowed_molecules]
    removed_data = data[[not allowed for allowed in allowed_molecules]]
    total_removed += allowed_molecules.count(False)

    filtered_data.to_csv(data_file, index=False)
    removed_data.to_csv(removed_data_file, index=False)

    print(f"Total data points removed in fingerprint: {total_removed}")




# ---------------------------------

# name = 'metxbiodb'
# df = pd.read_csv(f'dataset/curated_data/{name}_clean_unique_parents.csv')

# filtered_data, removed_data, total_removed = finger_print_similarity_filter(df)
# filtered_data.to_csv(f'dataset/curated_data/filtered/filtered_fingerprint_{name}_clean_unique_parents.csv', index=False)

# print("Total data points removed: ", total_removed)
# print(removed_data)


# ---------------------------------

'''
name = 'metxbiodb' # [ 'drugbank' 'metxbiodb' ]'
filter_method = valid_smiles #remember to change file name

df = pd.read_csv(f'dataset/curated_data/{name}_clean_unique_parents.csv')

filtered_data, removed_data, total_removed=filter_out_data_on_parent_and_child(df, filter_method)

filtered_data.to_csv(f'dataset/curated_data/filtered/filtered_valid_smiles_{name}_clean_unique_parents.csv', index=False)

print("Total data points removed: ", total_removed)
print(removed_data)
'''

if __name__ == "__main__":


    name = 'metxbiodb' # [ 'drugbank' 'metxbiodb' ]

    dataset = f'dataset/curated_data/{name}_smiles.csv'
    clean = f'dataset/curated_data/{name}_smiles_clean.csv'

    removed_valid_smiles = f'dataset/removed_data/{name}_removed_valid_smiles.csv'
    removed_atoms_allowed = f'dataset/removed_data/{name}_removed_atoms_allowed.csv'
    removed_weights_allowed = f'dataset/removed_data/{name}_removed_weights_allowed.csv'
    removed_fingerprints = f'dataset/removed_data/{name}_removed_fingerprints.csv'

    standardize_smiles(dataset, clean)
    filter_data_on_both_sides(clean, valid_smiles, removed_valid_smiles)
    filter_data_on_both_sides(clean, atoms_allowed_in_molecules, removed_atoms_allowed)
    filter_data_on_both_sides(clean, molecule_allowed_based_on_weight, removed_weights_allowed)
    finger_print_similarity_filter(clean, removed_fingerprints)