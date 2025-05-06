from rdkit import Chem
import numpy as np
from chembl_structure_pipeline.standardizer import standardize_mol

# Canonicalise
smiles = "CC(=O)OC1=CC=CC=C1C(O)=O"
mol = Chem.MolFromSmiles(smiles)
mol = standardize_mol(mol)
smiles = Chem.MolToSmiles(mol, canonical=True, isomericSmiles=False)
print(smiles)

# Randomise to non-canonicalise
mol = Chem.MolFromSmiles(smiles)
atom_order = list(range(mol.GetNumAtoms()))
np.random.shuffle(atom_order)
mol_rand = Chem.RenumberAtoms(mol, atom_order)
new_smiles = Chem.MolToSmiles(mol_rand, canonical=False)
print(new_smiles)

# Verify same smiles
new_mol = Chem.MolFromSmiles(new_smiles)
new_smiles = Chem.MolToSmiles(new_mol, canonical=True)
print(new_smiles)