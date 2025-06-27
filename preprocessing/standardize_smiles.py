import numpy as np
from rdkit import Chem
from rdkit.Chem.SaltRemover import SaltRemover
from chembl_structure_pipeline.standardizer import standardize_mol


def standardize_molecule(smiles, isomericSmiles=True):
    """
    Standardization of SMILES by standardiser;
    """

    try:
        rm = SaltRemover()
        mol = Chem.MolFromSmiles(smiles)
        mol = rm.StripMol(mol, dontRemoveEverything=True)
        if len(Chem.MolToSmiles(mol).split(".")) > 1:
            salt = Chem.MolToSmiles(mol)
            frag_dict = {len(k): k for k in salt.split(".")}
            max_frag = frag_dict[max(list(frag_dict.keys()))]
            mol = Chem.MolFromSmiles(max_frag)
        mol = standardize_mol(mol)
        clean_smiles = Chem.MolToSmiles(
            mol, isomericSmiles=isomericSmiles
        )  # removing (False) or keeping (True) stereochemistry
        if Chem.MolFromSmiles(clean_smiles) is None:
            clean_smiles = None
    except Exception:
        clean_smiles = np.nan
    return clean_smiles


def standardize_smiles_collection(smiles_list, isomericSmiles=True):
    """
    Standardization of SMILES collection by standardiser;
    """

    lookup = {}
    return_smiles = []
    for smiles in smiles_list:
        if smiles in lookup:
            standardized_smiles = lookup[smiles]
        else:
            standardized_smiles = standardize_molecule(
                smiles, isomericSmiles=isomericSmiles
            )
            lookup[smiles] = standardized_smiles
        return_smiles.append(standardized_smiles)

    return return_smiles
