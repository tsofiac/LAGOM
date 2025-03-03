from rdkit import Chem
from rdkit.Chem import Crippen
import pandas as pd

# --------------- LogP ------------------

def assign_logp_token(logp):
    """Assign a token based on predefined LogP intervals of size 1."""
    # Hard-coded intervals as tuples (min, max)
    intervals = [
        (-4.0, -3.0), (-3.0, -2.0), (-2.0, -1.0), (-1.0, 0.0), (0.0, 1.0),
        (1.0, 2.0), (2.0, 3.0), (3.0, 4.0), (4.0, 5.0), (5.0, 6.0),
        (6.0, 7.0), (7.0, 8.0), (8.0, 9.0), (9.0, 10.0), (10.0, 11.0),
        (11.0, 12.0), (12.0, 13.0)
    ]

    if logp < -4.0:
        return "[LogP_BELOW_MIN]"
    elif logp >= 13.0:
        return "[LogP_ABOVE_MAX]"
    else:
        for lower, upper in intervals:
            if lower <= logp < upper:
                return f"[LogP_[{lower}, {upper})]"

def annotate_smiles_with_logp(smiles):
    mol = Chem.MolFromSmiles(smiles)
    logp = Crippen.MolLogP(mol)
    token = assign_logp_token(logp)
    return f"{token}{smiles}"

def annotate_data_with_logp(datafile, annotated_datafile):
    dataset = pd.read_csv(datafile)

    dataset['parent_smiles'] = dataset['parent_smiles'].apply(annotate_smiles_with_logp)
    dataset.to_csv(annotated_datafile, index=False)

# --------------- CSP3 ------------------




# ----------------------------------------------------------------------------------------------

if __name__ == "__main__":

    logp_annotations = True
    csp3_annotations = True

    datafile = 'dataset/curated_data/combined_smiles_clean.csv'
    annotated_datafile = 'dataset/curated_data/logp_combined_smiles_clean.csv' 

    if logp_annotations:
        annotate_data_with_logp(datafile, annotated_datafile)

    if csp3_annotations:
        