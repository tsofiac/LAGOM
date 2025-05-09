from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdFMCS
import random

def find_mcs_and_generate_smiles(smiles1, smiles2, shuffle_start=False):
    # Create RDKit molecule objects
    mol1 = Chem.MolFromSmiles(smiles1)
    mol2 = Chem.MolFromSmiles(smiles2)

    # Find the Maximum Common Substructure (MCS)
    mcs = rdFMCS.FindMCS([mol1, mol2])
    mcs_mol = Chem.MolFromSmarts(mcs.smartsString)

    # Get the atom indices of the MCS in both molecules
    match1 = mol1.GetSubstructMatch(mcs_mol)
    match2 = mol2.GetSubstructMatch(mcs_mol)
    
    min_match = min([len(match1), len(match2)])

    # Choose a random starting atom from the MCS
    if shuffle_start:
        start_atom_index = random.randint(0, min_match - 1)
        start_atom1 = match1[start_atom_index]
        start_atom2 = match2[start_atom_index]
    else:
        start_atom1 = match1[0]
        start_atom2 = match2[0]

    # Generate new SMILES with the same starting point
    new_smiles1 = Chem.MolToSmiles(mol1, rootedAtAtom=start_atom1, canonical=False, isomericSmiles=False)
    new_smiles2 = Chem.MolToSmiles(mol2, rootedAtAtom=start_atom2, canonical=False, isomericSmiles=False)

    return new_smiles1, new_smiles2, mcs.smartsString

# Example usage
smiles1 = "CC(=O)OC1=CC=CC=C1C(=O)O"  # Aspirin
smiles2 = "CC(=O)OC1=CC=C(C=C1)C(=O)O"  # p-Acetylsalicylic acid

# Function to demonstrate multiple runs
def demonstrate_multiple_runs(n=5):
    print(f"Demonstrating {n} runs with different starting points:")
    for i in range(n):
        ordered_smiles1, ordered_smiles2, _ = find_mcs_and_generate_smiles(smiles1, smiles2, shuffle_start=True)
        print(f"Run {i+1}:")
        print(f"  Molecule 1: {ordered_smiles1}")
        print(f"  Molecule 2: {ordered_smiles2}")
        print()

# Demonstrate multiple runs
demonstrate_multiple_runs()

# Without shuffling (for comparison)
print("Without shuffling:")
normal_smiles1, normal_smiles2, mcs_smarts = find_mcs_and_generate_smiles(smiles1, smiles2, shuffle_start=False)
print(f"Molecule 1: {normal_smiles1}")
print(f"Molecule 2: {normal_smiles2}")
print(f"\nMCS SMARTS: {mcs_smarts}")