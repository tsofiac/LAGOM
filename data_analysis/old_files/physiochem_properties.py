import argparse
from rdkit import Chem
# from rdkit.Chem import PandasTools
from rdkit.Chem import Descriptors
from rdkit.Chem import rdMolDescriptors
from rdkit.Chem import Lipinski
from rdkit.Chem import Crippen
import pandas as pd

'''
Script to calculate physicochemical properties of molecules:
get_physicochemical_properties to calculate 15 descriptors:
number of nitrogen atoms, number of oxygen atoms, number of chiral centers, molecular weight, number of heavy atoms,
number of hydrogen bond acceptors, number of hydrogen bond donors, logP, topological polar surface area, number of
aromatic atoms, number of rings, fraction of Csp3 atoms, 
number of sulfur atoms, number of halogen atoms, molar refractivity 
'''
def manual_add_molecule_column(molDF, smiles_column, molecule_column='Molecule'):
    """
    Manually converts SMILES to RDKit molecule objects and adds them as a new column.
    
    :param molDF: DataFrame containing SMILES strings.
    :param smiles_column: Name of the column containing SMILES.
    :param molecule_column: Name of the column to be created/used for RDKit molecule objects.
    """
    molDF[molecule_column] = molDF[smiles_column].apply(Chem.MolFromSmiles)

def get_physicochemical_properties(molDF,smiles_column):
    '''
    Applies all property calculations to the ring systems of the dataframe and stores each property in a new column

    :param molDF: dataframe with ring systems as SMILES in the column 'ringSmiles'
    :return: a dataframe with ring system molecules and their properties
    '''
    # PandasTools.AddMoleculeColumnToFrame(molDF, smiles_column, 'Molecule')
    manual_add_molecule_column(molDF, smiles_column)
    print('Start calculcating properties.')
    molDF['N'] = molDF['Molecule'].apply(get_molecule_composition, args=(7,))
    molDF['O'] = molDF['Molecule'].apply(get_molecule_composition, args=(8,))
    molDF['chiral'] = molDF['Molecule'].apply(get_nof_chiral_centers)
    molDF['MW'] = molDF['Molecule'].apply(get_MW)
    molDF['heavy_atoms'] = molDF['Molecule'].apply(num_heavy_atoms)
    molDF['h_acc'] = molDF['Molecule'].apply(num_of_h_acceptors_and_donors, args=(True,))
    molDF['h_don'] = molDF['Molecule'].apply(num_of_h_acceptors_and_donors, args=(False,))
    molDF['logP'] = molDF['Molecule'].apply(get_logp)
    molDF['TPSA'] = molDF['Molecule'].apply(get_TPSA)
    molDF['numAro'] = molDF['Molecule'].apply(num_aromatic_atoms)
    molDF['numRings'] = molDF['Molecule'].apply(num_rings)
    molDF['frac_csp3'] = molDF['Molecule'].apply(fraction_csp3)    
    molDF['S'] = molDF['Molecule'].apply(get_molecule_composition, args=(16,))
    molDF['nHalogens'] = molDF['Molecule'].apply(num_halogens)
    molDF['MR'] = molDF['Molecule'].apply(get_mr)


def get_molecule_composition(mol, requestedAtomicNum):
    '''
    Counts the number of atoms of a given element in the ring system

    :param mol: the ring system molecule
    :param requestedAtomicNum: atomic number of the element for which the occurrence should be counted
    :return: the number of atoms of an element
    '''
    counter = 0
    for atom in mol.GetAtoms():
        atomicNum = atom.GetAtomicNum()
        if atomicNum == requestedAtomicNum:
            counter += 1
    return counter


def get_nof_chiral_centers(mol):
    return len(Chem.FindMolChiralCenters(mol, includeUnassigned=True))


def get_MW(mol):
    return round(Descriptors.MolWt(mol), 3)


def num_heavy_atoms(mol):
    return Lipinski.HeavyAtomCount(mol)


def num_of_h_acceptors_and_donors(mol, acc=True):
    if acc:
        return Lipinski.NumHAcceptors(mol)
    else:
        return Lipinski.NumHDonors(mol)


def get_logp(mol):
    return round(Crippen.MolLogP(mol), 3)


def get_TPSA(mol):
    return round(Descriptors.TPSA(mol), 3)


def num_aromatic_atoms(mol):
    numAromaticAtoms = 0
    for atom in mol.GetAtoms():
        if atom.GetIsAromatic():
            numAromaticAtoms += 1
    return numAromaticAtoms


def sum_formal_charge(mol):
    formalCharge = 0
    for atom in mol.GetAtoms():
        formalCharge += atom.GetFormalCharge()
    return formalCharge


def num_rings(mol):
    return rdMolDescriptors.CalcNumRings(mol)


def fraction_csp3(mol):
    return round(Descriptors.FractionCSP3(mol), 3)


def num_halogens(mol):
    return Chem.Fragments.fr_halogen(mol)


def get_mr(mol):
    '''Wildman-Crippen MR value
    Uses an atom-based scheme based on the values in the paper:
    Wildman and G. M. Crippen JCICS 39 868-873 (1999)
    '''
    return round(Crippen.MolMR(mol),3)



