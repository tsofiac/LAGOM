from rdkit import Chem, DataStructs
from rdkit.Chem import Descriptors, AllChem 
from rdkit.Chem.SaltRemover import SaltRemover
import torch
import pandas as pd
import os
from pathlib import Path

def set_up_device(): 
    # The device which will run the model. CUDA OR CPU
    device = "cuda" if torch.cuda.is_available() else "mps" if torch.has_mps or torch.backends.mps.is_available() else "cpu"
    print("Using device: ", device)
    
    if (device == 'cuda'):
        print(f"Device name: {torch.cuda.get_device_name(device.index)}")
        print(f"Device memory: {torch.cuda.get_device_properties(device.index).total_memory / 1024 ** 3} GB")
    elif (device == 'mps'):
        print(f"Device name: <mps>")
    device = torch.device(device)
    return device

####
def standardize_molecule(smiles): 
    rm = SaltRemover()
    mol = Chem.MolFromSmiles(smiles)
    # if all is salt, leave the last one. If two last is the same, leaves both
    mol = rm.StripMol(mol, dontRemoveEverything=True) 

    # if aromatic bonds is left, remove the smallest of them
    if len(Chem.MolToSmiles(mol).split('.')) > 1: 
        salt = Chem.MolToSmiles(mol)
        frag_dict = {len(k): k for k in salt.split('.')}
        max_frag = frag_dict[max(list(frag_dict.keys()))]
        mol = Chem.MolFromSmiles(max_frag)    
    smiles = Chem.MolToSmiles(mol, isomericSmiles=False, kekuleSmiles=True) 
    return smiles

def canonicalise_smiles(smiles): 
    mol = Chem.MolFromSmiles(smiles)
    canonical = Chem.MolToSmiles(mol, isomericSmiles=False)
    return canonical

def atoms_allowed_in_molecules(molecule): 
    atoms_to_include = ['C', 'N', 'S', 'O', 'H', 'F', 'I', 'P', 'B', 'Cl', 'Br', 'Si']
    mol = Chem.MolFromSmiles(molecule)
    atoms = [atom.GetSymbol() for atom in mol.GetAtoms()]
    return set(atoms).issubset(set(atoms_to_include))

def molecule_allowed_based_on_weight(molecule, max_weight=1000, min_weight=100): 
    mol_weight = Descriptors.ExactMolWt(Chem.MolFromSmiles(molecule))
    if mol_weight <= max_weight and mol_weight >= min_weight: 
        return True
    return False 

def valid_smiles(molecule): 
    return Chem.MolFromSmiles(molecule) is not None

def validate_molecule_match(expected, predicted): 
    if not valid_smiles(predicted): 
        return False
    if not valid_smiles(expected): 
        return False
    expected_smiles = canonicalise_smiles(expected)
    predicted_smiles = canonicalise_smiles(predicted)

    return expected_smiles == predicted_smiles
   # Chem.RDKFingerprint(mol1)

def fingerprint_similarity_single(correct, predicted): 
    parents = [Chem.MolFromSmiles(smiles) for smiles in correct if valid_smiles(smiles)]
    metabolites = [Chem.MolFromSmiles(smiles) for smiles in predicted if valid_smiles(smiles)]

    fpgen = AllChem.GetRDKitFPGenerator()
    fps_parents = [fpgen.GetFingerprint(x) for x in parents]
    fps_metabolites = [fpgen.GetFingerprint(x) for x in metabolites]

    similarities = []
    for (fpp, fpm) in zip(fps_parents, fps_metabolites): 
        sim = DataStructs.TanimotoSimilarity(fpp, fpm)
        
        similarities.append(sim)
    return similarities

def fingerprint_similarities_multiple(correct, predicted_lists): 
    cor = [Chem.MolFromSmiles(smiles) for smiles in correct if valid_smiles(smiles)]
    precs_lists = [[Chem.MolFromSmiles(smiles) for smiles in predicted if valid_smiles(smiles)] for predicted in predicted_lists]

    fpgen = AllChem.GetRDKitFPGenerator()
    fps_cor = [fpgen.GetFingerprint(x) for x in cor]
    fps_precs = [[fpgen.GetFingerprint(x) for x in prec] for prec in precs_lists]
    print(type(fps_precs))

    similarities = []
    for fps in fps_precs: 
        sims = []
        for (fpc, fpp) in zip(fps_cor, fps): 
            sim = DataStructs.TanimotoSimilarity(fpc, fpp)
            
            sims.append(sim)
        similarities.append(sims)
    return similarities

# Saves a log file from a training run into a designated folder src/data/
def save_output(data, col_name, model_name, epoch, filestep = None):
    df = pd.DataFrame({col_name:data, "epoch": epoch})
    Path(f"data/{col_name}").mkdir(parents=True, exist_ok=True)
    file_path = f'data/{col_name}/{model_name}.csv'
    file_exists = os.path.isfile(file_path)
    
    if filestep is not None: 
        Path(f"data/{col_name}/{filestep}").mkdir(parents=True, exist_ok=True)
        file_path = f'data/{col_name}/{filestep}/{model_name}.csv'
        file_exists = os.path.isfile(file_path)

    if(file_exists):
        existing_df = pd.read_csv(file_path)
        combined_df = pd.concat([existing_df, df], ignore_index=True)
        combined_df.to_csv(file_path, index=False)
    else:
        df.to_csv(file_path, index=False)