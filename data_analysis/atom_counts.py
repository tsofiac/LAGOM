from rdkit import Chem, DataStructs
from rdkit.Chem import Descriptors, AllChem 
from rdkit.Chem.SaltRemover import SaltRemover
import torch
import pandas as pd
import os
from pathlib import Path
from collections import Counter

def define_atoms(molecule):
    mol = Chem.MolFromSmiles(molecule)
    atoms = [atom.GetSymbol() for atom in mol.GetAtoms()]
    atom_counts = Counter(atoms)
    return atom_counts

def analyse_atoms(data, column):
    total_atoms_count = Counter()
    
    for molecule in data[column]:
        atom_counts = define_atoms(molecule)
        total_atoms_count.update(atom_counts)
    
    return total_atoms_count

def print_atom_counts(atom_counts):
    for atom, count in atom_counts.items():
        print(f"{atom}: {count}")

name = 'drugbank' # [ 'drugbank' 'metxbiodb' ]
dataset = f'dataset/curated_data/{name}_smiles.csv'
data = pd.read_csv(dataset)

parents_count = analyse_atoms(data, 'parent_smiles')
children_count = analyse_atoms(data, 'child_smiles')

print("Parent SMILES atom counts:")
print_atom_counts(parents_count)

print("\nChild SMILES atom counts::")
print_atom_counts(children_count)