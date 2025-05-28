from rdkit import Chem
from rdkit.Chem import Crippen, Descriptors
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
    return f"{token}"

# def annotate_data_with_logp(datafile, annotated_datafile):
#     dataset = pd.read_csv(datafile)

#     dataset['parent_smiles'] = dataset['parent_smiles'].apply(annotate_smiles_with_logp)
#     dataset.to_csv(annotated_datafile, index=False)

# --------------- CSP3 ------------------

def assign_csp3_token(csp3):
    """Assign a token based on predefined CSP3 fractions intervals of size 0.1."""
    # Hard-coded intervals as tuples (min, max)
    intervals = [
        (0.0, 0.1), (0.1, 0.2), (0.2, 0.3), (0.3, 0.4), (0.4, 0.5), (0.5, 0.6), (0.6, 0.7), (0.7, 0.8), (0.8, 0.9) #, (0.9, 1.0) is included below
    ]

    if 0.9 <= csp3 <= 1:
        return "[csp3_[0.9, 1.0]]"
    else:
        for lower, upper in intervals:
            if lower <= csp3 < upper:
                return f"[csp3_[{lower}, {upper})]"

def annotate_smiles_with_csp3(smiles):
    mol = Chem.MolFromSmiles(smiles)
    csp3 = Descriptors.FractionCSP3(mol)
    token = assign_csp3_token(csp3)
    return f"{token}"

# def annotate_data_with_csp3(datafile, annotated_datafile):
#     dataset = pd.read_csv(datafile)

#     dataset['parent_smiles'] = dataset['parent_smiles'].apply(annotate_smiles_with_csp3)
#     dataset.to_csv(annotated_datafile, index=False)




# --------------- LogP & CSP3 ------------------

def annotate_data(logp_annotations, csp3_annotations, datafile, annotated_datafile): # logp3_annotations and csp3_annotations are either set as True or False

    dataset = pd.read_csv(datafile) # has column parent_smiles

    for index, row in dataset.iterrows():
        smiles = row["parent_smiles"]

        if logp_annotations and csp3_annotations:
            token1 = annotate_smiles_with_logp(smiles)
            token2 = annotate_smiles_with_csp3(smiles)
            smiles_annotated = f"{token1}{token2}{smiles}"
        elif logp_annotations:
            token = annotate_smiles_with_logp(smiles)
            smiles_annotated = f"{token}{smiles}"
        elif csp3_annotations:
            token = annotate_smiles_with_csp3(smiles)
            smiles_annotated = f"{token}{smiles}"

        dataset.at[index, "parent_smiles"] = smiles_annotated

    dataset.to_csv(annotated_datafile, index=False)

def annotate_data_finetune(logp_annotations, csp3_annotations, datafile, annotated_datafile): # if file already ready for finetuning

    dataset = pd.read_csv(datafile, sep='\t') # has column parent_smiles

    for index, row in dataset.iterrows():
        smiles = row["reactants"]

        if logp_annotations and csp3_annotations:
            token1 = annotate_smiles_with_logp(smiles)
            token2 = annotate_smiles_with_csp3(smiles)
            smiles_annotated = f"{token1}{token2}{smiles}"
        elif logp_annotations:
            token = annotate_smiles_with_logp(smiles)
            smiles_annotated = f"{token}{smiles}"
        elif csp3_annotations:
            token = annotate_smiles_with_csp3(smiles)
            smiles_annotated = f"{token}{smiles}"

        dataset.at[index, "reactants"] = smiles_annotated

    dataset.to_csv(annotated_datafile, sep='\t', index=False)   



#------------annotate enzymes----------------

def assign_enzyme_family(enzyme):

    if enzyme is None or not isinstance(enzyme, str):
        return '[unknown]'

    enzyme = enzyme.lower()

    if 'cyp' in enzyme or 'cytochrome p450' in enzyme:
        return '[CYP]'
    
    elif 'ugt'in enzyme or 'udp-glucuronosyltransferase' in enzyme:
        return '[UGT]'
    
    elif 'sult' in enzyme or 'sulfotransferase' in enzyme:
        return '[SULT]'
    
    elif 'gst' in enzyme or 'glutathione s-transferase' in enzyme:
        return '[GST]'
    
    elif 'nat' in enzyme or 'n-acetyltransferase' in enzyme:
        return '[NAT]'
    
    elif 'methyltransferase' in enzyme:
        return '[MT]'
    
    elif 'carboxylesterase' in enzyme:
        return '[CES]'
    
    elif 'kinase' in enzyme:
        return '[kinase]'

    else:
        print(enzyme)
        return '[other]'

def annotate_data_enzymefamily(datafile, output_datafile):

    dataset = pd.read_csv(datafile) # has column parent_smiles

    for index, row in dataset.iterrows():
        smiles = row["parent_smiles"]
        enzyme = row["enzymes"]

        token = assign_enzyme_family(enzyme)
        smiles_annotated = f"{token}{smiles}"
     

        dataset.at[index, "parent_smiles"] = smiles_annotated

    dataset.to_csv(output_datafile, index=False)

def annotate_data_enzymefamily_finetune(datafile, output_datafile):

    dataset = pd.read_csv(datafile, sep='\t') # has column parent_smiles

    for index, row in dataset.iterrows():
        smiles = row["reactants"]
        enzyme = row["enzymes"]

        token = assign_enzyme_family(enzyme)
        smiles_annotated = f"{token}{smiles}"

        dataset.at[index, "reactants"] = smiles_annotated

    dataset.to_csv(output_datafile, sep='\t', index=False)  

# ----------------------------------------------------------------------------------------------

def reformat_for_chemformer(input_file, output_file):
    df = pd.read_csv(input_file)

    df = df.rename(columns={
        "child_smiles": "products",
        "parent_smiles": "reactants",
    })

    df.to_csv(output_file, sep='\t', index=False)


def add_possible_products(input_file):
    df = pd.read_csv(input_file, sep='\t') 

    def agg_order_preserved(series):
        return list(series)

    grouped = df.groupby(['reactants'], sort=False).agg({
        'parent_name': 'first',  # Get the first entry for 'parent_name' in each group
        'child_name': agg_order_preserved, 
        'products': agg_order_preserved 
    }).reset_index()

    grouped.rename(columns={'products': 'possible_products'}, inplace=True)

    df = df.merge(grouped[['reactants', 'possible_products']], on='reactants', how='left')

    df['possible_products'] = df['possible_products'].apply(lambda x: '.'.join(x))

    df.to_csv(input_file, sep='\t', index=False)

def combine_datasets(df1, df2, output_csv):
    selected_df1 = df1.copy()
    selected_df2 = df2.copy()

    combined_df = pd.concat([selected_df1, selected_df2], ignore_index=True)
    combined_df.to_csv(output_csv, index=False)



if __name__ == "__main__":

    enzyme_annotations = False # Can not be done together with the other annotations
    logp_annotations = False
    csp3_annotations = False
    using_finetune_ready_file = False # True for evaluation set
    
    augment_parent_grandchild = False
    augment_parent_parent = True
   
    #OBS: look over all file names

    name = 'PP' # 'annotations' 'logp' 'csp3' 'enzymes' 'PG' 'PP' 'PG_PP'

    annotated_datafile = f'dataset/curated_data/annotated_data/{name}_combined_smiles_clean.csv' 
    finetune_file = f'dataset/finetune/{name}_finetune.csv'

    if logp_annotations or csp3_annotations:

        if using_finetune_ready_file is False:
            datafile = 'dataset/curated_data/combined_smiles_clean.csv'

            if enzyme_annotations is True:
                annotate_data_enzymefamily(datafile, annotated_datafile)  

            else:
                annotate_data(logp_annotations, csp3_annotations, datafile, annotated_datafile)
                
            reformat_for_chemformer(annotated_datafile, finetune_file)
            add_possible_products(finetune_file)

        else:
            datafile = 'dataset/finetune/combined_evaluation_finetune.csv'
            finetune_file = f'dataset/finetune/{name}_evaluation_finetune.csv'

            if enzyme_annotations is True:
                annotate_data_enzymefamily_finetune(datafile, finetune_file)
                print('here')
            else:
                annotate_data_finetune(logp_annotations, csp3_annotations, datafile, finetune_file)


    if augment_parent_grandchild and not augment_parent_parent:
        combined_csv = 'dataset/curated_data/combined_smiles_clean.csv'
        parent_grandchild = 'dataset/curated_data/augmented_parent_grandchild.csv'


        combined_df = pd.read_csv(combined_csv)
        parent_grandchild_df = pd.read_csv(parent_grandchild)
        combine_datasets(combined_df, parent_grandchild_df, annotated_datafile)

        reformat_for_chemformer(annotated_datafile, finetune_file)
        add_possible_products(finetune_file)

    if augment_parent_parent and not augment_parent_grandchild:
        combined_csv = 'dataset/curated_data/combined_smiles_clean.csv'
        parent_parent = 'dataset/curated_data/augmented_parent_parent.csv'

        combined_df = pd.read_csv(combined_csv)
        parent_parent_df = pd.read_csv(parent_parent)
        combine_datasets(combined_df, parent_parent_df, annotated_datafile)

        reformat_for_chemformer(annotated_datafile, finetune_file)
        add_possible_products(finetune_file)

    if augment_parent_parent and augment_parent_grandchild:
        combined_csv = 'dataset/curated_data/combined_smiles_clean.csv'
        parent_grandchild = 'dataset/curated_data/augmented_parent_grandchild.csv'
        parent_parent = 'dataset/curated_data/augmented_parent_parent.csv'

        combined_df = pd.read_csv(combined_csv)
        parent_grandchild_df = pd.read_csv(parent_grandchild)
        combine_datasets(combined_df, parent_grandchild_df, annotated_datafile)

        annotated_df = pd.read_csv(annotated_datafile)
        parent_parent_df = pd.read_csv(parent_parent)
        combine_datasets(annotated_df, parent_parent_df, annotated_datafile)

        reformat_for_chemformer(annotated_datafile, finetune_file)
        add_possible_products(finetune_file)