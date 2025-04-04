from rdkit import Chem
from rdkit.Chem import rdFMCS
import numpy as np
import random
from random import shuffle
from rxnutils.chem.utils import split_smiles_from_reaction
import pandas as pd
from standardize_smiles import standardize_smiles_collection


# ------------------- NOT USING -----------------------------------------
def augment_smiles(smiles, augment_prob): 

    smiles_aug = []
    for smi in smiles:
        if augment_prob < np.random.rand():
            smiles_aug.append(smi)
            continue

        mols = list(map(Chem.MolFromSmiles, split_smiles_from_reaction(smi)))
        for _ in range(3):
            try:
                smi_new = [Chem.MolToSmiles(mol, doRandom=True) for mol in mols]
                shuffle(smi_new)
                smiles_aug.append(".".join(smi_new))
                break
            except Exception as e:
                print(f"Augmentation failed for {smi} with error: {e}")
        else:
            smiles_aug.append(smi)
            print(f"Augmentation failed three times for {smi}, returning unaugmented original")

    return smiles_aug
# -------------------------------------------------------------------------

def randomize_mol_restricted(mol):
    # Standard shuffle surprisingly leads to 35% slower code.
    atom_order = list(range(mol.GetNumAtoms()))
    np.random.shuffle(atom_order)
    return Chem.RenumberAtoms(mol, atom_order)

def augment_smiles_restricted(smiles):

    mol1 = Chem.MolFromSmiles(smiles[0])
    mol2 = Chem.MolFromSmiles(smiles[1])
    for _ in range(5):

        mol_rand1 = randomize_mol_restricted(mol1)
        mol_rand2 = randomize_mol_restricted(mol2)
        new_smiles1 = Chem.MolToSmiles(mol_rand1, canonical=False)
        new_smiles2 = Chem.MolToSmiles(mol_rand2, canonical=False)
        new_smiles = [new_smiles1, new_smiles2]
        if new_smiles[0] != smiles[0] and new_smiles[1] != smiles[1]:
            aug_smiles = new_smiles
            break
    else:
        aug_smiles = [None, None]
        print(f"Augmentation failed five times for {smiles}")

    return aug_smiles

def find_mcs_and_generate_smiles(smiles, shuffle_start=False):

    mol1 = Chem.MolFromSmiles(smiles[0])
    mol2 = Chem.MolFromSmiles(smiles[1])

    for _ in range(5):

        mcs = rdFMCS.FindMCS([mol1, mol2])
        mcs_mol = Chem.MolFromSmarts(mcs.smartsString)

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
        new_smiles = [new_smiles1, new_smiles2]
        if new_smiles[0] != smiles[0] and new_smiles[1] != smiles[1]:
            aug_smiles = new_smiles
            break

    else:
        aug_smiles = [None, None]
        print(f"Augmentation failed five times for {smiles}")

    
    return aug_smiles

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

def find_mcs_and_generate_smiles_strict(smiles, shuffle_start=False):

    mol1 = Chem.MolFromSmiles(smiles[0])
    mol2 = Chem.MolFromSmiles(smiles[1])

    for _ in range(5):

        mcs = rdFMCS.FindMCS([mol1, mol2])
        mcs_mol = Chem.MolFromSmarts(mcs.smartsString)

        match1 = mol1.GetSubstructMatch(mcs_mol)
        match2 = mol2.GetSubstructMatch(mcs_mol)
        
        min_match = min([len(match1), len(match2)])

        # Choose a starting atom from the MCS
        if shuffle_start:
            start_atom_index = random.randint(0, min_match - 1)
        else:
            start_atom_index = 0

        # Generate new SMILES with the MCS first
        new_smiles1 = generate_mcs_first_smiles(mol1, match1, start_atom_index)
        new_smiles2 = generate_mcs_first_smiles(mol2, match2, start_atom_index)
        new_smiles = [new_smiles1, new_smiles2]
        if new_smiles[0] != smiles[0] and new_smiles[1] != smiles[1]:
            aug_smiles = new_smiles
            break

    else:
        aug_smiles = [None, None]
        print(f"Augmentation failed five times for {smiles}")

    return aug_smiles


def augment_dataset(csv_file, augmented_file, nr_of_reactions):

    df = pd.read_csv(csv_file)
    df = df[df['set'] == 'train']

    parent_smiles = df['parent_smiles']

    augmented_data = []
    for index, row in df.iterrows():
        print(index)
        parent_smiles = row['parent_smiles']
        child_smiles = row['child_smiles']

        for _ in range(nr_of_reactions):

            # smiles_aug = augment_smiles_restricted([parent_smiles, child_smiles])
            # smiles_aug = find_mcs_and_generate_smiles([parent_smiles, child_smiles], shuffle_start=True)
            smiles_aug = find_mcs_and_generate_smiles_strict([parent_smiles, child_smiles], shuffle_start=True)
            if None not in smiles_aug:
                    
                augmented_row = {
                    'parent_name': row['parent_name'],
                    'child_name': row['child_name'],
                    'parent_smiles': smiles_aug[0],
                    'child_smiles': smiles_aug[1],
                    'origin': row['origin'],
                    'source': row['source'],
                    'set': row['set']
                }
                    
                augmented_data.append(augmented_row)

    augmented_dataset = pd.DataFrame(augmented_data)

    augmented_dataset.to_csv(augmented_file, index=False)


def reformat_for_chemformer(input1_file, input2_file, output_file):
    df1 = pd.read_csv(input1_file)
    df2 = pd.read_csv(input2_file)

    df2['source'] = 'Randomised'

    df = pd.concat([df1, df2], ignore_index=True)

    df = df.rename(columns={
        "child_smiles": "products",
        "parent_smiles": "reactants",
    })

    df.to_csv(output_file, sep='\t', index=False)


combined_csv = 'dataset/curated_data/combined_smiles_clean.csv'
augmented_csv = 'dataset/curated_data/randomised.csv'
combined_finetune = 'dataset/finetune/combined_finetune.csv'


# augment_dataset(combined_csv, augmented_csv, 1)

specific_csv = 'dataset/curated_data/randomised_mcs.csv'
reformat_for_chemformer(combined_csv, specific_csv, combined_finetune)

# combined_df = pd.read_csv(combined_csv)
# parent_smiles = combined_df['parent_smiles']
# print(parent_smiles[0])

# aug_smiles = augment_smiles(parent_smiles, 1)
# print(aug_smiles[0])
# aug_smiles_sd = standardize_smiles_collection(aug_smiles, False)
# print(aug_smiles_sd[0])