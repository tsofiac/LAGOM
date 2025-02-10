import xml.etree.ElementTree as ET
import pandas as pd
from rdkit import Chem

def sdf_to_csv(sdf_file, csv_file, is_drug_information = False): # I buy this
    name_property = "GENERIC_NAME" if is_drug_information else "NAME"
    data = []
    suppl = Chem.SDMolSupplier(sdf_file)
    for mol in suppl:
        if mol is None:
            continue
        
        drugbank_id = mol.GetProp("DRUGBANK_ID") if mol.HasProp("DRUGBANK_ID") else "NaN"
        inchi = Chem.inchi.MolToInchi(mol) if mol else "NaN"
        smiles = Chem.MolToSmiles(mol) if mol else "NaN"
        name = mol.GetProp(name_property) if mol.HasProp(name_property) else "NaN"
        # products = mol.GetProp("PRODUCTS") if mol.HasProp("PRODUCTS") else "NaN"

        data.append({
            "dbid": drugbank_id,
            "name": name,
            "smiles": smiles,
            "inchi": inchi
        })

    df = pd.DataFrame(data)
    df.to_csv(csv_file, index=False)

def reformat_external_csv(input_file, output_file): # makes 100% sense
    db_external_df = pd.read_csv(input_file)
    db_external_df = db_external_df[["DrugBank ID","SMILES","InChI",'Name']]
    db_external_df = db_external_df.rename(columns={
        "DrugBank ID": "dbid",
        "Name": "name",
        "SMILES": "smiles",
        "InChI": "inchi"
    })
    db_external_df.to_csv(output_file,index=False)

def remove_bad_drug_metabolite_rows(input_file): # this makes 100% sense
    df = pd.read_csv(input_file)
    initial_row_count = len(df)
    # if the identifiers are empty we remove
    df.dropna(subset=['dbid', 'name'], how='all', inplace=True)
    # If the data we need is empty the 
    df.dropna(subset=['smiles', 'inchi'], how='all', inplace=True)
    final_row_count = len(df)
    dropped_count = initial_row_count - final_row_count
    print("Number of dropped rows: ", dropped_count)
    df.to_csv(input_file, index=False)


class Reaction:
    def __init__(self):
        self.left_side_name = ""
        self.left_side_id = ""
        self.right_side_name = ""
        self.right_side_id = ""
        self.enzymes = []

    def has_data(self):
        return any([getattr(self, attr) for attr in vars(self)])

def generate_reaction_pairs(input_file, output_file,num_iter = 20000): # possible to incorporate enzymes from prev master thesis code
    # Create an iterator for parsing the XML file incrementally
    context = ET.iterparse(input_file, events=('start', 'end'))
    parent_child_ids = []
    
    local_reaction = Reaction()
    is_in_local_reaction = False
    # is_in_enzymes = False
    right_or_left_reaction = "left"
    count = 0
    for i, (event, element) in enumerate(context):
        # Keep track on iteration and potenially break
        if (i % 1000000 == 0):
            print(f'iteration {i}')
        if(i == num_iter):
            break
        # Split the tag to remove garbage
        if(element.tag =="root"):
            tag = "root"
        else:
            tag = element.tag.split('}')[1]      
        if event == 'start':
            # We have entered a local reaction.
            if(tag == "reaction"):
                is_in_local_reaction = True
                # If local reaction has data then it is complete or it needs to be investigated
                if(local_reaction.has_data()):                   
                    # We need atleast an id or a name on both sides of the reaction
                    # To map data
                    if((local_reaction.left_side_id == "" and local_reaction.left_side_name == "") or (local_reaction.right_side_id == "" and local_reaction.right_side_name == "")):
                        print("Data is missing here. Remove and investigate this datapoint")
                        count += 1
                    else:
                        parent_child_ids.append([local_reaction.left_side_id, local_reaction.left_side_name, local_reaction.right_side_id, local_reaction.right_side_name])
                    # we reset here since we want to look for the new reaction
                    local_reaction = Reaction()

            # To keep track on what data we are looking at
            if(tag == "left-element" and is_in_local_reaction):
                right_or_left_reaction = "left"
            if(tag == "right-element" and is_in_local_reaction):
                right_or_left_reaction = "right"
            
            if(tag == "drugbank-id" and is_in_local_reaction):
                if(right_or_left_reaction == "left"):
                    local_reaction.left_side_id = element.text
                elif(right_or_left_reaction == "right"):
                    local_reaction.right_side_id = element.text
        
            if(tag == "name" and is_in_local_reaction):
                if(right_or_left_reaction == "left"):
                    local_reaction.left_side_name = element.text
                elif(right_or_left_reaction == "right"):
                    local_reaction.right_side_name = element.text  
            
            # We know its the end of the reaction if this tag shows
            if(tag == "snp-effects"):
                is_in_local_reaction = False

        elif event == 'end':
            # Clean up the element when the end tag is encountered
            element.clear()
    print('Number of rows with missing data: ', count)
    df = pd.DataFrame(parent_child_ids, columns=["parent_id",'parent_name',"child_id",'child_name'])
    df.to_csv(output_file, index=False)

def filter_endogenous_reaction(input_file):

    original_df = pd.read_csv(input_file)
    only_db_parents_df = original_df[~original_df['parent_id'].str.contains('DBMET').fillna(False)] # all parents that not contain DBMET
    original_df = original_df[original_df['parent_id'].str.contains('DBMET').fillna(False)] # all parents that contain DBMET

    num_new_rows = 1
    while( num_new_rows != 0 ):
        # Parent is DBMET. If parent is child in any reaction in only_db_parents_df. Then keep it
        
        child_ids = only_db_parents_df['child_id'].tolist()
        # store rows where parents are child in list
        metabolite_parent_with_drug_origin_df  = original_df[ original_df['parent_id'].isin(child_ids) ]
        # drop all childs that are not childs of drugs or drug metabolites
        original_df = original_df[~original_df['parent_id'].isin(child_ids) ]

        only_db_parents_df = pd.concat([only_db_parents_df,metabolite_parent_with_drug_origin_df])
        num_new_rows = len(metabolite_parent_with_drug_origin_df)
        print("Number of dropped rows: ", len(original_df))

    only_db_parents_df.to_csv(input_file, index=False)


def hmdb_metabolite_map(input_file, output_file, num_iter = 10000):

    # Create an iterator for parsing the XML file incrementally
    context = ET.iterparse(input_file, events=('start', 'end'))
    smiles_db_pair = []
    smiles_db_pairs = []
    prev_was_accession = False
    for i, (event, element) in enumerate(context):
        # Keep track on iteration and potenially break
        if (i % 1000000 == 0):
            print(f'iteration {i}')
        if(i == num_iter):
            break

        # Split the tag to remove garbage
        if(element.tag =="root"):
            tag = "root"
        else:
            tag = element.tag.split('}')[1]
                  

        if event == 'start':
            if(tag == "metabolite"):
                if(len(smiles_db_pair) == 3 ):
                    smiles_db_pairs.append(smiles_db_pair)
                    smiles_db_pair = []

                else:
                    smiles_db_pair = []
            if(tag == "smiles" or tag == "drugbank_id"):
                smiles_db_pair.append(element.text)
            if(tag == "accession"):
                prev_was_accession = True
                continue
            if(prev_was_accession and tag == "name"):
                smiles_db_pair.append(element.text)
                prev_was_accession = False

        elif event == 'end':
            # Clean up the element when the end tag is encountered
            element.clear()

    df = pd.DataFrame(smiles_db_pairs, columns=['name',"smiles","dbid"])
    df['inchi'] = "" # added by Sofia
    df.to_csv(output_file, index=False)


# combines two dataframes
def combine_datasets(first_file, second_file, output_file, remove_dupes=True):
    first_df = pd.read_csv(first_file, on_bad_lines='skip')
    print(first_df.columns)
    second_df = pd.read_csv(second_file, on_bad_lines='skip')
    print(second_df.columns)
    output_df = pd.concat([first_df, second_df], ignore_index=True)

    if remove_dupes:
        output_df = output_df.drop_duplicates()

    output_df.to_csv(output_file, index=False)


def enrich_reaction_pairs(reaction_pairs_file, structures_file, output_file):
    # Read the CSV files into DataFrames
    reaction_pairs_df = pd.read_csv(reaction_pairs_file)
    structures_df = pd.read_csv(structures_file, low_memory=False)

    # Merge to add parent structure data
    enriched_parents_df = reaction_pairs_df.merge(
        structures_df,
        left_on=['parent_id', 'parent_name'],
        right_on=['dbid', 'name'],
        how='left'
    ).rename(columns={
        'smiles': 'parent_smiles',
        'inchi': 'parent_inchi'
    }).drop(columns=['dbid', 'name'])

    # Merge to add child structure data
    enriched_children_df = enriched_parents_df.merge(
        structures_df,
        left_on=['child_id', 'child_name'],
        right_on=['dbid', 'name'],
        how='left'
    ).rename(columns={
        'smiles': 'child_smiles',
        'inchi': 'child_inchi'
    }).drop(columns=['dbid', 'name'])


    # Save the enriched DataFrame to the output file
    enriched_children_df.to_csv(output_file, index=False)


def clean_smiles(input_file): # this makes 100% sense
    df = pd.read_csv(input_file)
    initial_row_count = len(df)
    # Remove rows where either 'parent_smiles' or 'child_smiles' is missing
    df.dropna(subset=['parent_smiles', 'child_smiles'], how='any', inplace=True)
    final_row_count = len(df)
    dropped_count = initial_row_count - final_row_count
    print('Number of dropped rows: ', dropped_count)
    df.to_csv(input_file, index=False)


if __name__ == "__main__":

    get_drug_info = False
    get_metabolite_info = False
    get_external_drugs = False
    get_reaction_pairs = False
    get_hdmb_drugs = False
    
    combine_structures = False
    extend_dataset = True
   

    ## Drug structure - drugbank_drug_structures.sdf
    drugbank_drug_structures = "dataset/raw_data/drugbank_drug_structures.sdf"
    parsed_drug_structures = "dataset/processed_data/drugbank_drug_structures.csv"
    if get_drug_info: 
        sdf_to_csv(drugbank_drug_structures, parsed_drug_structures, is_drug_information=True)
        remove_bad_drug_metabolite_rows(parsed_drug_structures)
        # Number of rows dropped: 0

    ## Metabolite structure - drugbank_metabolite_structures.sdf
    drugbank_metabolite_structures = "dataset/raw_data/drugbank_metabolite_structures.sdf"
    parsed_metabolite_structures = "dataset/processed_data/drugbank_metabolite_structures.csv"
    if get_metabolite_info: 
        sdf_to_csv(drugbank_metabolite_structures, parsed_metabolite_structures, is_drug_information=False)
        remove_bad_drug_metabolite_rows(parsed_metabolite_structures)
        # Number of rows dropped: 0
            
    ## External drugs - drugbank_external_structures.csv
    drugbank_external_structures = "dataset/raw_data/drugbank_external_structures.csv"
    parsed_external_stuctures = "dataset/processed_data/drugbank_external_structures_cleaned.csv"
    if get_external_drugs:
        reformat_external_csv(drugbank_external_structures, parsed_external_stuctures)
        remove_bad_drug_metabolite_rows(parsed_external_stuctures)
        # Number of rows dropped: 850

    # Get the reaction pairs - drugbank_full_database.xml
    drugbank_full_database = "dataset/raw_data/drugbank_full_database.xml"
    drugbank_reaction_pairs = "dataset/processed_data/drugbank_reaction_pairs.csv"
    if get_reaction_pairs: 
        generate_reaction_pairs(drugbank_full_database, drugbank_reaction_pairs,-1)
        filter_endogenous_reaction(drugbank_reaction_pairs)
        # Missing data: 45
        # Number of reactions: from 4045 to 3194

    ## HMDB drug and metabolites - hmdb_metabolites.xml
    hmdb_metabolites = "dataset/raw_data/hmdb_metabolites.xml"
    hmdb_cleaned = "dataset/processed_data/hmdb_cleaned.csv"
    if get_hdmb_drugs:        
        hmdb_metabolite_map(hmdb_metabolites,hmdb_cleaned, -1)
        remove_bad_drug_metabolite_rows(hmdb_cleaned)
        # Number of rows dropped: 2119


    drugbank_full_structures = "dataset/processed_data/drugbank_full_structures.csv"
    if combine_structures:          
        combine_datasets(parsed_drug_structures, parsed_metabolite_structures, drugbank_full_structures)
        combine_datasets(drugbank_full_structures, parsed_external_stuctures, drugbank_full_structures)
        combine_datasets(drugbank_full_structures, hmdb_cleaned, drugbank_full_structures)


    # Extend dataset with SMILES or INCHI
    matched_smiles_to_drugs = "dataset/processed_data/extended_reaction_pairs.csv"
    if extend_dataset:
        print('Nice Try Mohahahha')
        # enrich_reaction_pairs(drugbank_reaction_pairs, drugbank_full_structures, matched_smiles_to_drugs) # alot of duplicates in id
        # clean_smiles(matched_smiles_to_drugs)