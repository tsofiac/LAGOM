from rdkit import Chem, DataStructs
from rdkit.Chem import Descriptors, AllChem 
from rdkit.Chem.SaltRemover import SaltRemover
import torch
import pandas as pd
import os
from pathlib import Path
from collections import Counter

# input should have child_smiles and parent_smiles
# def filter_out_data_on_parent_and_child(data, filter_method): 
#     total_removed = 0

#     # Filtering based on metabolites
#     allowed_metabolites = [filter_method(metabolite) for metabolite in data["child_smiles"]]
#     data = data[allowed_metabolites]
#     total_removed += allowed_metabolites.count(False)

#     # Filtering based on parent molecules
#     allowed_molecules = [filter_method(molecule) for molecule in data["parent_smiles"]]
#     data = data[allowed_molecules]
#     total_removed += allowed_molecules.count(False)

#     return data, total_removed

# saving removed data
def filter_out_data_on_parent_and_child(data, filter_method): 
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

    return filtered_data, removed_data, total_removed

def filter_out_data_on_parent(data, filter_method): 
    total_removed = 0

    # Filtering based on parent molecules
    allowed_molecules = [filter_method(molecule) for molecule in data["parent_smiles"]]
    data = data[allowed_molecules]
    total_removed += allowed_molecules.count(False)

    return data, total_removed

def molecule_allowed_based_on_weight(molecule, max_weight=700, min_weight=100): 
    mol_weight = Descriptors.ExactMolWt(Chem.MolFromSmiles(molecule))
    if mol_weight <= max_weight and mol_weight >= min_weight: 
        return True
    return False 

def valid_smiles(molecule): 
    return Chem.MolFromSmiles(molecule) is not None

def atoms_allowed_in_molecules(molecule): 
    atoms_to_include = ['C', 'N', 'S', 'O', 'H', 'F', 'I', 'P', 'Cl', 'Br'] # B and Si removed in contrast to last year
    mol = Chem.MolFromSmiles(molecule)
    atoms = [atom.GetSymbol() for atom in mol.GetAtoms()]
    return set(atoms).issubset(set(atoms_to_include)) #Returns true if all atoms in the molecule are included in the specified set


name = 'metxbiodb' # [ 'drugbank' 'metxbiodb' ]'
filter_method = atoms_allowed_in_molecules #remember to change file name

df = pd.read_csv(f'dataset/curated_data/{name}_clean_unique_parents.csv')

filtered_data, removed_data, total_removed=filter_out_data_on_parent_and_child(df, filter_method)

filtered_data.to_csv(f'dataset/curated_data/filtered/filtered_atoms_{name}_clean_unique_parents.csv', index=False)

print("Total data points removed: ", total_removed)
print(removed_data)