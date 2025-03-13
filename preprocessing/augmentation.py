from rdkit import Chem
import numpy as np
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

def augment_smiles_restricted(smiles, augment_prob):

    smiles_aug = []

    augment_prob_og = augment_prob
    augment_prob = 1.0  # To avoid double application

    for smi in smiles:
        if augment_prob < np.random.rand():
            smiles_aug.append(smi)
            continue

        mol = Chem.MolFromSmiles(smi)
        for _ in range(3):
            try:
                mol_rand = randomize_mol_restricted(mol)
                smiles_aug.append(Chem.MolToSmiles(mol_rand, canonical=False))
                break
            except Exception as e:
                print(f"Augmentation failed for {smi} with error: {e}")
        else:
            smiles_aug.append(smi)
            print(f"Augmentation failed three times for {smi}, returning unaugmented original")

    augment_prob = augment_prob_og
    return smiles_aug

def augment_dataset(csv_file, augment_prob, augmented_file, nr_of_reactions):

    df = pd.read_csv(csv_file)
    parent_smiles = df['parent_smiles']

    augmented_data = []
    for index, row in df.iterrows():

        parent_smiles = row['parent_smiles']
        child_smiles = row['child_smiles']

        for _ in range(nr_of_reactions):

            if 1-augment_prob < np.random.rand():
                smiles_aug = augment_smiles_restricted([parent_smiles, child_smiles], augment_prob)
            else:
                smiles_aug = [parent_smiles, child_smiles]

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


def reformat_for_chemformer(input_file, output_file):
    df = pd.read_csv(input_file)

    df = df.rename(columns={
        "child_smiles": "products",
        "parent_smiles": "reactants",
    })

    df.to_csv(output_file, sep='\t', index=False)


combined_csv = 'dataset/curated_data/combined_smiles_clean.csv'
augmented_csv = 'dataset/curated_data/randomised.csv'

finetune = 'dataset/finetune/randomised_finetune.csv'

augment_dataset(combined_csv, 1, augmented_csv, 2)
reformat_for_chemformer(augmented_csv, finetune)

# combined_df = pd.read_csv(combined_csv)
# parent_smiles = combined_df['parent_smiles']
# print(parent_smiles[0])

# aug_smiles = augment_smiles(parent_smiles, 1)
# print(aug_smiles[0])
# aug_smiles_sd = standardize_smiles_collection(aug_smiles, False)
# print(aug_smiles_sd[0])