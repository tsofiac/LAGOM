from rdkit import Chem
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

    min_match = min(len(match1), len(match2))

    # Choose a starting atom from the MCS
    if shuffle_start:
        start_atom_index = random.randint(0, min_match - 1)
    else:
        start_atom_index = 0

    # Generate new SMILES with the MCS first
    new_smiles1 = generate_mcs_first_smiles(mol1, match1, start_atom_index)
    new_smiles2 = generate_mcs_first_smiles(mol2, match2, start_atom_index)

    return new_smiles1, new_smiles2, mcs.smartsString

def generate_mcs_first_smiles(mol, mcs_match, start_atom_index):
    # Create a new molecule with only the MCS atoms
    mcs_mol = Chem.RWMol()
    atom_map = {}

    for i, idx in enumerate(mcs_match):
        atom = mol.GetAtomWithIdx(idx)
        new_atom = Chem.Atom(atom.GetAtomicNum())
        new_idx = mcs_mol.AddAtom(new_atom)
        atom_map[idx] = new_idx

    # Add bonds within the MCS
    for bond in mol.GetBonds():
        begin_idx = bond.GetBeginAtomIdx()
        end_idx = bond.GetEndAtomIdx()
        if begin_idx in atom_map and end_idx in atom_map:
            mcs_mol.AddBond(atom_map[begin_idx], atom_map[end_idx], bond.GetBondType())

    # Generate SMILES for the MCS part
    mcs_smiles = Chem.MolToSmiles(mcs_mol, rootedAtAtom=atom_map[mcs_match[start_atom_index]], canonical=False, isomericSmiles=False)

    # Add non-MCS atoms and bonds
    for atom in mol.GetAtoms():
        if atom.GetIdx() not in mcs_match:
            new_atom = Chem.Atom(atom.GetAtomicNum())
            new_idx = mcs_mol.AddAtom(new_atom)
            atom_map[atom.GetIdx()] = new_idx

    for bond in mol.GetBonds():
        begin_idx = bond.GetBeginAtomIdx()
        end_idx = bond.GetEndAtomIdx()
        if not (begin_idx in mcs_match and end_idx in mcs_match):
            mcs_mol.AddBond(atom_map[begin_idx], atom_map[end_idx], bond.GetBondType())

    # Generate final SMILES
    final_smiles = Chem.MolToSmiles(mcs_mol, rootedAtAtom=atom_map[mcs_match[start_atom_index]], canonical=False, isomericSmiles=False)

    return final_smiles

# Example usage
smiles1 = 'CN1CCC23C4C1Cc1c2c(c(O)cc1)OC3C(O)C=C4' # parent
smiles2 = 'CN1CCC23C4C1Cc1c2c(c(O)cc1)OC3C(OC1OC(C(=O)O)C(O)C(O)C1O)C=C4' # child

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