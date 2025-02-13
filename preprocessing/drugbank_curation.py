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

        data.append({
            "dbid": drugbank_id,
            "name": name,
            "smiles": smiles,
            "inchi": inchi
        })

    df = pd.DataFrame(data, columns=["dbid", "name", "smiles", "inchi"])
    df.to_csv(csv_file, index=False)

def reformat_external_csv(input_file, output_file):
    db_external_df = pd.read_csv(input_file)
    db_external_df = db_external_df.rename(columns={
        "DrugBank ID": "dbid",
        "Name": "name",
        "SMILES": "smiles",
        "InChI": "inchi"
    })
    db_external_df = db_external_df[["dbid", "name", "smiles", "inchi"]]
    db_external_df.to_csv(output_file, index=False)

def remove_bad_drug_metabolite_rows(input_file):
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
    df['inchi'] = ""
    df.to_csv(output_file, index=False)

def count_nan_smiles(csv_file):
    df = pd.read_csv(csv_file)
    # Count the number of NaN values in the 'smiles' column
    nan_count = df['smiles'].isna().sum()
    print('Missing SMILES: ', nan_count)

def inchi_to_smiles(csv_file):

    df = pd.read_csv(csv_file)
    # Iterate over each row where `smiles` is NaN
    for index, row in df[df['smiles'].isna()].iterrows():
        try:
            # Convert InChI to molecule
            mol = Chem.inchi.MolFromInchi(row['inchi'])
            
            if mol is not None:
                # Convert molecule to SMILES
                smiles = Chem.MolToSmiles(mol)
                # Replace NaN value in `smiles` column
                df.at[index, 'smiles'] = smiles
            else:
                print(f"Conversion failed for entry {row['dbid']} with InChI: {row['inchi']}")
        
        except Exception as e:
            print(f"ERROR processing row {index}, DBID: {row['dbid']}, error: {e}")

    df.to_csv(csv_file, index=False)


# combines two dataframes
def combine_datasets(first_file,second_file,output_file):
    first_df = pd.read_csv(first_file,on_bad_lines='skip')
    print(first_df.columns)
    second_df = pd.read_csv(second_file,on_bad_lines='skip')
    print(second_df.columns)

    first_df.drop(columns=['inchi'], inplace=True, errors='ignore')
    second_df.drop(columns=['inchi'], inplace=True, errors='ignore')

    output_df = pd.concat([first_df,second_df],ignore_index=True)

    output_df.to_csv(output_file, index=False)


def remove_duplicates(csv_file):

    df = pd.read_csv(csv_file)

    print('Nr Rows 1: ', len(df))
    df=df[df['dbid'].isnull() | ~df[df['dbid'].notnull()].duplicated(subset='dbid',keep='first')]
    print('Nr Rows 2: ', len(df))
    df=df[df['name'].isnull() | ~df[df['name'].notnull()].duplicated(subset='name',keep='first')]
    print('Nr Rows 3: ', len(df))

    df.to_csv(csv_file, index=False)

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
    print('Number of rows: ', len(original_df))
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

    print('Number of rows: ', len(only_db_parents_df))
    unique_reactions_df = only_db_parents_df.drop_duplicates() # Keep unique rows only 
    print('Number of rows: ', len(unique_reactions_df))
    unique_reactions_df.to_csv(input_file, index=False)


def map_smiles_to_reaction_pairs(reaction_pairs_file, structures_file, output_file, is_drug=False):
    smiles_property = "parent_smiles" if is_drug else "child_smiles"
    id_property = "parent_id" if is_drug else "child_id"
    name_property = "parent_name" if is_drug else "child_name"

    reaction_pairs_df = pd.read_csv(reaction_pairs_file)
    structures_df = pd.read_csv(structures_file)

    reaction_pairs_df[smiles_property] = None

    # Fill the smiles property based on id_property and name_property
    for index, row in reaction_pairs_df.iterrows():
        # Check for id_property match
        if pd.notna(row[id_property]):
            # Lookup match in structures_df
            matched_row = structures_df[structures_df['dbid'] == row[id_property]]
            if not matched_row.empty:
                # Assign smiles from the matched structure
                reaction_pairs_df.at[index, smiles_property] = matched_row['smiles'].values[0]
        
        # Check for name_property match if smiles not found yet
        elif pd.notna(row[name_property]):
            # Lookup match in structures_df
            matched_row = structures_df[structures_df['name'] == row[name_property]]
            if not matched_row.empty:
                # Assign smiles from the matched structure
                reaction_pairs_df.at[index, smiles_property] = matched_row['smiles'].values[0]

    # Save the enriched DataFrame to CSV
    reaction_pairs_df.to_csv(output_file, index=False)


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

    get_drug_structures = False
    get_metabolite_structures = False
    get_external_structures = False
    get_hdmb_structures = False
    get_reaction_pairs = False
    
    combine_all_structures = True
    extend_dataset = True
   

    ## Drug structure - drugbank_drug_structures.sdf
    drugbank_drug_structures = "dataset/raw_data/drugbank_drug_structures.sdf"
    parsed_drug_structures = "dataset/processed_data/drugbank_drug_structures.csv"
    if get_drug_structures:
        sdf_to_csv(drugbank_drug_structures, parsed_drug_structures, is_drug_information=True)
        remove_bad_drug_metabolite_rows(parsed_drug_structures)
        count_nan_smiles(parsed_drug_structures)
        # Number of rows dropped: 0
        # Missing SMILES: 0

    ## Metabolite structure - drugbank_metabolite_structures.sdf
    drugbank_metabolite_structures = "dataset/raw_data/drugbank_metabolite_structures.sdf"
    parsed_metabolite_structures = "dataset/processed_data/drugbank_metabolite_structures.csv"
    if get_metabolite_structures: 
        sdf_to_csv(drugbank_metabolite_structures, parsed_metabolite_structures, is_drug_information=False)
        remove_bad_drug_metabolite_rows(parsed_metabolite_structures)
        count_nan_smiles(parsed_metabolite_structures)
        # Number of rows dropped: 0
        # Missing SMILES: 0
            
    ## External drugs - drugbank_external_structures.csv
    drugbank_external_structures = "dataset/raw_data/drugbank_external_structures.csv"
    parsed_external_stuctures = "dataset/processed_data/drugbank_external_structures_cleaned.csv"
    if get_external_structures:
        reformat_external_csv(drugbank_external_structures, parsed_external_stuctures)
        remove_bad_drug_metabolite_rows(parsed_external_stuctures)
        count_nan_smiles(parsed_external_stuctures)
        inchi_to_smiles(parsed_external_stuctures)
        count_nan_smiles(parsed_external_stuctures)
        # Number of rows dropped: 850
        # Missing SMILES: 3
        # Missing SMILES: 0

    ## HMDB drug and metabolites - hmdb_metabolites.xml
    hmdb_metabolites = "dataset/raw_data/hmdb_metabolites.xml"
    hmdb_cleaned = "dataset/processed_data/hmdb_cleaned.csv"
    if get_hdmb_structures:
        hmdb_metabolite_map(hmdb_metabolites,hmdb_cleaned, -1)
        remove_bad_drug_metabolite_rows(hmdb_cleaned)
        count_nan_smiles(hmdb_cleaned)
        # Number of rows dropped: 2119
        # Missing SMILES: 0

    # Get the reaction pairs - drugbank_full_database.xml
    drugbank_full_database = "dataset/raw_data/drugbank_full_database.xml"
    drugbank_reaction_pairs = "dataset/processed_data/drugbank_reaction_pairs.csv"
    if get_reaction_pairs: 
        generate_reaction_pairs(drugbank_full_database, drugbank_reaction_pairs,-1)
        filter_endogenous_reaction(drugbank_reaction_pairs)
        # Missing data: 45
        # Nr of rows: 4045 (original)
        # Nr of rows: 3194 (endogenous)
        # Nr of rows: 2610 (duplicates)

    # Combine all strucutres into one file
    drugbank_full_structures = "dataset/processed_data/drugbank_full_structures.csv"
    if combine_all_structures:          
        combine_datasets(parsed_drug_structures, parsed_metabolite_structures, drugbank_full_structures)
        combine_datasets(drugbank_full_structures, parsed_external_stuctures, drugbank_full_structures)
        # combine_datasets(drugbank_full_structures, hmdb_cleaned, drugbank_full_structures) 
        remove_duplicates(drugbank_full_structures)

    # Extend reaction pairs with SMILES
    mapped_smiles_to_reactions = "dataset/curated_data/drugbank_smiles.csv"
    if extend_dataset:
        map_smiles_to_reaction_pairs(drugbank_reaction_pairs, drugbank_full_structures, mapped_smiles_to_reactions, True)
        map_smiles_to_reaction_pairs(mapped_smiles_to_reactions, drugbank_full_structures, mapped_smiles_to_reactions, False)
        clean_smiles(mapped_smiles_to_reactions)
        # Number of dropped rows:  1312 (without HMDB)
        
        
    # Without hmdb
    # drugbank_full_structures: 225782 -> 15324
    # drugbank_smiles: 1300 -> 1299

    # Might be worth removing hmdb all in all... 