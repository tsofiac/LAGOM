import xml.etree.ElementTree as ET
import pandas as pd
from rdkit import Chem

def sdf_to_csv(sdf_file, csv_file, is_drug_information = False):
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
            "inchi": inchi,
            "smiles": smiles,
            "name": name
        })

    df = pd.DataFrame(data)
    df.to_csv(csv_file, index=False)

def reformat_external_csv(input_file, output_file):
    db_external_df = pd.read_csv(input_file)
    db_external_df = db_external_df[["DrugBank ID","SMILES","InChI",'Name']]
    db_external_df = db_external_df.rename(columns={
        "DrugBank ID": "dbid",
        "SMILES": "smiles",
        "InChI": "inchi",
        "Name": "name"
    })
    db_external_df.to_csv(output_file,index=False)

def remove_bad_drug_metabolite_rows(input_file): # removes row if both identifiers are missing or if both smiles and inchi_key are missing
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

# Method generate the source of the reaction pairs
def generate_reaction_pairs(input_file, output_file,num_iter = 20000):
    # Create an iterator for parsing the XML file incrementally
    context = ET.iterparse(input_file, events=('start', 'end'))
    parent_child_ids = []
    
    local_reaction = Reaction()
    is_in_local_reaction = False
    # is_in_enzymes = False
    right_or_left_reaction = "left"
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
                    if((local_reaction.left_side_id == "" and local_reaction.left_side_name == "") or (local_reaction.right_side_id == "" and local_reaction.right_side_name == "")): # right...
                    # if((local_reaction.left_side_id == "" or local_reaction.right_side_id == "") and (local_reaction.left_side_name == "" or local_reaction.right_side_name == "")): # wrong
                        print("Data is missing here. Remove and investigate this datapoint")
                    else:
                        # We join the list of enzymes to single string 
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
    df = pd.DataFrame(parent_child_ids, columns=["parent_id",'parent_name',"child_id",'child_name'])
    df.to_csv(output_file, index=False)

def filter_endogenous_reaction(input_file):

    original_df = pd.read_csv(input_file)
    only_db_parents_df = original_df[~original_df['parent_id'].str.contains('DBMET').fillna(False)]
    original_df = original_df[original_df['parent_id'].str.contains('DBMET').fillna(False)]

    num_new_rows = 1
    iter = 0
    while( num_new_rows != 0 ):
        # if reaction (parent,child). Parent is DBMET. If parent is child in any reaction in only_db_parents_df. Then keep it
        child_ids = only_db_parents_df['child_id'].tolist()

        # store rows where parents are child in list
        metabolite_parent_with_drug_origin_df  = original_df[ original_df['parent_id'].isin(child_ids) ]
        # Drop these rows to not get duplicates
        original_df = original_df[~original_df['parent_id'].isin(child_ids) ]

        only_db_parents_df = pd.concat([only_db_parents_df,metabolite_parent_with_drug_origin_df])
        num_new_rows = len(metabolite_parent_with_drug_origin_df)
        print(len(original_df))

        iter +=1
    only_db_parents_df.to_csv(input_file, index=False)

def extend_and_filter_reaction_pairs(reaction_pairs_file, structure_data_file, output_file):
    # Load data from CSV files into DataFrames
    reaction_pairs_df = pd.read_csv(reaction_pairs_file)
    structure_data_df = pd.read_csv(structure_data_file)

    # Create a mapping dictionary {key: (smiles, inchi)}
    structure_map = {}
    for _, row in structure_data_df.iterrows():
        dbid = row['dbid']
        name = row['name']
        smiles = row['smiles']
        inchi = row['inchi']

        # Map using both dbid and name if they exist
        if pd.notna(dbid):
            structure_map[dbid] = (smiles, inchi)
        if pd.notna(name):
            structure_map[name] = (smiles, inchi)

    # Function to find smiles or inchi
    def map_smiles_inchi(row):
        # Check with parent_id
        key = row['parent_id']
        if key in structure_map:
            smiles, inchi = structure_map[key]
            return smiles if pd.notna(smiles) else inchi

        # Check with parent_name if parent_id didn't match
        key = row['parent_name']
        if key in structure_map:
            smiles, inchi = structure_map[key]
            return smiles if pd.notna(smiles) else inchi

        return None

    # Apply the function to each row in the DataFrame
    reaction_pairs_df['mapped_info'] = reaction_pairs_df.apply(map_smiles_inchi, axis=1)

    # Count and remove rows with no matched data
    unmatched_count = reaction_pairs_df['mapped_info'].isna().sum()
    filtered_df = reaction_pairs_df.dropna(subset=['mapped_info'])

    # Inform the number of unmatched rows
    print(f"Number of unmatched entries removed: {unmatched_count}")

    # Save the filtered DataFrame to a new CSV file
    filtered_df.to_csv(output_file, index=False)


if __name__ == "__main__":

    get_drug_info = False
    get_metabolite_info = False
    get_external_drugs = False
    get_reaction_pairs = False
   

    ## Drug structure - drugbank_drug_structures.sdf
    drugbank_drug_structures = "dataset/raw_data/drugbank_drug_structures.sdf"
    parsed_drug_structures = "dataset/processed_data/drugbank_drug_structures.csv"
    if get_drug_info: 
        sdf_to_csv(drugbank_drug_structures, parsed_drug_structures, is_drug_information=True)
        remove_bad_drug_metabolite_rows(parsed_drug_structures)

    ## Metabolite structure - drugbank_metabolite_structures.sdf
    drugbank_metabolite_structures = "dataset/raw_data/drugbank_metabolite_structures.sdf"
    parsed_metabolite_structures = "dataset/processed_data/drugbank_metabolite_structures.csv"
    if get_metabolite_info: 
        sdf_to_csv(drugbank_metabolite_structures, parsed_metabolite_structures, is_drug_information=False)
        remove_bad_drug_metabolite_rows(parsed_metabolite_structures)
            
    ## External drugs - drugbank_external_structures.csv
    drugbank_external_structures = "dataset/raw_data/drugbank_external_structures.csv"
    drugbank_external_cleaned = "dataset/processed_data/drugbank_external_structures_cleaned.csv"
    if get_external_drugs:
        reformat_external_csv(drugbank_external_structures, drugbank_external_cleaned)
        remove_bad_drug_metabolite_rows(drugbank_external_cleaned)

    # Get the reaction pairs - drugbank_full_database.xml
    drugbank_full_database = "dataset/raw_data/drugbank_full_database.xml"
    drugbank_reaction_pairs = "dataset/processed_data/drugbank_reaction_pairs.csv"
    if get_reaction_pairs: 
        generate_reaction_pairs(drugbank_full_database, drugbank_reaction_pairs,-1)
        filter_endogenous_reaction(drugbank_reaction_pairs)

    # Example usage:
    matched_smiles_to_drugs = "dataset/processed_data/extended_reaction_pairs.csv"
    extend_and_filter_reaction_pairs(drugbank_reaction_pairs, parsed_drug_structures, matched_smiles_to_drugs)



 