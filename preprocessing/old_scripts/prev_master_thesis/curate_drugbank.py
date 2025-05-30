import xml.etree.ElementTree as ET
import pandas as pd

class Reaction:
    def __init__(self):
        self.left_side_name = ""
        self.left_side_id = ""
        self.right_side_name = ""
        self.right_side_id = ""
        self.enzymes = []

    def has_data(self):
        return any([getattr(self, attr) for attr in vars(self)])


# method is used to parse a specific xml file to generate a df with id, name, smiles. 
# The source filename is drugbank_full_database.xml
# Method generate the source of the reaction pairs
def parse_and_save_metabolite(input_file, output_file,num_iter = 20000):
    # Create an iterator for parsing the XML file incrementally
    context = ET.iterparse(input_file, events=('start', 'end'))
    parent_child_ids = []
    
    local_reaction = Reaction()
    is_in_local_reaction = False
    is_in_enzymes = False
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
                    if((local_reaction.left_side_id == "" or local_reaction.right_side_id == "") and (local_reaction.left_side_name == "" or local_reaction.right_side_name == "")):
                        
                        print("Data is missing here. Remove and investigate this datapoint")
                       
                    # If all data exists then we add it to the list
                    else:
                      # We join the list of enzymes to single string 
                      enzyme_string = ", ".join(local_reaction.enzymes)
                    
                      parent_child_ids.append([local_reaction.left_side_id,local_reaction.left_side_name,local_reaction.right_side_id,local_reaction.right_side_name,enzyme_string])
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
            
            # get the names of the enzymes
            if(tag == "enzyme" and is_in_local_reaction):
                right_or_left_reaction = "none"
                is_in_enzymes = True
            if( is_in_enzymes):
                if(tag == "name"):
                    if(element.text != None):
                        local_reaction.enzymes.append(element.text)
                    is_in_enzymes = False
            
            # We know its the end of the reaction if this tag shows
            if(tag == "snp-effects"):
                is_in_enzymes = False
                is_in_local_reaction = False

        elif event == 'end':
            # Clean up the element when the end tag is encountered
            element.clear()
    df = pd.DataFrame(parent_child_ids, columns=["parent_id",'parent_name',"child_id",'child_name','enzymes'])
    df.to_csv(output_file, index=False)
    
# reads and generates a dataframe with relevant information such as name, id, smiles, inchi
def read_sdf_file(input_file, output_file,num_iter = 20000, skip_name=False, total_rows=514046, is_drug_information = False ):
  
  name_of_interest = "<GENERIC_NAME>" if is_drug_information else "<NAME>"

  id_smiles_pairs = []
  id_smiles_pair = []
  save_next_line = False
  f = open(input_file,'r',encoding="utf8")
  i = 0
  for line in f.readlines():
    i+=1
    # Keep track on iteration and potenially break
    if (i % 1000 == 0):
        print(f'iteration {i} out of {total_rows}, {(i+1)/total_rows:.2f} done')

    if(i == num_iter):
        break
    # Go to the next line. Not worth computing
    if( not('DATABASE_ID'  in line or
       'INCHI_KEY'  in line or
       'SMILES'   in line or
       'NAME'   in line or 
       'PRODUCT' in line or
       'M  END'  in line)):
        if(not save_next_line):
            continue
    # Means it is a new molecule
    if("M  END" in line):
        # If both contain data we append it
        if( id_smiles_pair != [] ):
            if(skip_name):
                id_smiles_pair.append('')
           
            id_smiles_pairs.append(id_smiles_pair)
        else:
            print("There was missing data here! Investigate")
        id_smiles_pair = []

    if('<DATABASE_ID>' in line or
       '<INCHI_KEY>'in line or
       '<SMILES>' in line or
       name_of_interest in line or 
       '<PRODUCTS>' in line):
        save_next_line = True

        if('<PRODUCTS>' in line):
            was_product = True
            continue
        was_product = False
        continue
    if(save_next_line):
        
        line = line.replace('\n','')
        if(was_product):
            line = "Had product(s)" 
        id_smiles_pair.append(line)
        save_next_line = False
    
    df = pd.DataFrame(id_smiles_pairs, columns=["dbid","smiles","inchi_key",'name','products'])

    df.to_csv(output_file, index=False)

# Method combines two files. 
# input_pairs_name is the file of reaction pairs
# input_id_smiles is file with the smiles, id, name information for any molecule
def combine_id_with_smiles_inchi(input_pairs_name, input_id_smiles, output_file,output_error_file ):
    pairs_name_df = pd.read_csv(input_pairs_name)
    id_smiles_df = pd.read_csv(input_id_smiles)

    pairs_with_no_match_df = pd.DataFrame([],columns=['who_failed',"parent_id",'child_id','parent_name','child_name','enzymes'])

    output_pairs_smiles = pd.DataFrame(columns=[
                    'parent_name',
                    'parent_smiles',
                    'parent_inchi',
                    'metabolite_name',
                    'metabolite_smiles',
                    'metabolite_inchi',
                    'enzymes'])
    for index, row in pairs_name_df.iterrows():
        if(index % 100 == 0):
            print(f"current iteration: {index}")

        # Check by ID
        parent_smiles = id_smiles_df.loc[id_smiles_df['dbid'] == row['parent_id']]['smiles']
        parent_inchi = id_smiles_df[id_smiles_df['dbid'] == row['parent_id']]['inchi_key']
        
        metabolite_smiles = id_smiles_df.loc[id_smiles_df['dbid'] == row['child_id']]['smiles']
        metabolite_inchi = id_smiles_df[id_smiles_df['dbid'] == row['child_id']]['inchi_key']
        
        # Clean the data
        parent_name = row['parent_name']
        parent_smiles = parent_smiles.iloc[0] if (len(parent_smiles) >= 1) else ""
        parent_inchi = parent_inchi.iloc[0] if(len(parent_inchi) >= 1) else ""
    
        metabolite_name = row['child_name']
        metabolite_smiles = metabolite_smiles.iloc[0] if(len(metabolite_smiles) >= 1 ) else ""
        metabolite_inchi = metabolite_inchi.iloc[0] if(len(metabolite_inchi) >= 1) else ""
        
 
        # get smiles by the name
        if(parent_smiles == "" ):
            parent_smiles = id_smiles_df.loc[id_smiles_df['name'] == row['parent_name']]['smiles']
            parent_smiles = parent_smiles.iloc[0] if (len(parent_smiles) >= 1) else ""

        if(metabolite_smiles == ""):
            metabolite_smiles = id_smiles_df.loc[id_smiles_df['name'] == row['child_name']]['smiles']
            metabolite_smiles = metabolite_smiles.iloc[0] if(len(metabolite_smiles) >= 1 ) else ""

        # get inchi by the name
        if(parent_inchi == "" ):
            parent_inchi = id_smiles_df.loc[id_smiles_df['name'] == row['parent_name']]['inchi_key']
            parent_inchi = parent_inchi.iloc[0] if(len(parent_inchi) >= 1) else ""

        if(metabolite_inchi == ""):
            metabolite_inchi = id_smiles_df.loc[id_smiles_df['name'] == row['child_name']]['inchi_key']
            metabolite_inchi = metabolite_inchi.iloc[0] if(len(metabolite_inchi) >= 1) else ""

        if(parent_smiles == "" and parent_inchi == "")  or (metabolite_smiles == "" and metabolite_inchi == ""):
            if(parent_smiles == "" and parent_inchi == ""):
                who_failed = "parent"
            elif(metabolite_smiles == "" and metabolite_inchi == ""):
                who_failed = "child"
            else:
                who_failed = "both"
            pair_failed_df = pd.DataFrame( [[who_failed,row['parent_id'], row['child_id'],parent_name ,metabolite_name, row['enzymes']]], columns=['who_failed',"parent_id",'child_id','parent_name','child_name','enzymes']) 

            pairs_with_no_match_df = pd.concat([pairs_with_no_match_df, pair_failed_df],ignore_index=True)
            continue
        else:
          enzymes = row['enzymes']
          pair_df = pd.DataFrame( 
              [[parent_name, parent_smiles, parent_inchi, metabolite_name, metabolite_smiles, metabolite_inchi, enzymes]],
                columns=[
                    'parent_name',
                    'parent_smiles',
                    'parent_inchi',
                    'metabolite_name',
                    'metabolite_smiles',
                    'metabolite_inchi',
                    'enzymes']) 
          
          output_pairs_smiles = pd.concat([output_pairs_smiles, pair_df],ignore_index=True)
    
    output_pairs_smiles.to_csv(output_file, index=False)
    pairs_with_no_match_df.to_csv(output_error_file, index=False)

#MDB data as a complementary source of molecule data such as name, inchi, id, smiles
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
    df.to_csv(output_file, index=False)

# combines two dataframes
def combine_datasets(first_file,second_file,output_file, remove_dupes=False):
    first_df = pd.read_csv(first_file,on_bad_lines='skip')
    print(first_df.columns)
    second_df = pd.read_csv(second_file,on_bad_lines='skip')
    print(second_df.columns)
    output_df = pd.concat([first_df,second_df],ignore_index=True)

    output_df.to_csv(output_file, index=False)


def remove_bad_drug_metabolite_rows(input_file):
    df = pd.read_csv(input_file)
    initial_row_count = len(df)
    # if the identifiers are empty we remove
    df.dropna(subset=['dbid', 'name'], how='all', inplace=True)
    # If the data we need is empty the 
    df.dropna(subset=['smiles', 'inchi_key'], how='all', inplace=True)
    final_row_count = len(df)
    dropped_count = initial_row_count - final_row_count
    print("Number of dropped rows: ", dropped_count)
    df.to_csv(input_file, index=False)


def filter_endogenous_reaction(input_file, output_file):

# parent_id, parent_name, child_id, child_name, enzymes
# DB00001, Lepirudin, DBMET03462, M1 (1-64),
# DBMET03462, M1 (1-64), DBMET03463, M2 (1-63),

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
    only_db_parents_df.to_csv(output_file, index=False)

def drop_products_from_drug_information(input_File, output_file):
    df = pd.read_csv(input_File, on_bad_lines='skip')
    count = len(df)
    print(df.columns)
    df = df.dropna(subset=['products'])
    print("dropped: ", count - len(df))
    df.to_csv(output_file, index=False)
        

if __name__ == "__main__":
    # Conditionals 
    get_reaction_pairs = False
    remove_endogenous_pairs = False
    get_drug_info = False
    get_metabolite_info = False
    get_external_drugs = False
    get_hdmb_drugs = False
    drop_drugs_with_no_products = True

    combine_the_data = True
    # Methods 
    # These can be run all together but easier if run one at a time.
    drugbank_filepath = "dataset/raw_data/drugbank_full_database.xml"
    parsed_pair_name_output = "dataset/raw_data/drugbank_pairs_names.csv"
    # Get the reaction pairs
    if get_reaction_pairs: 
        parse_and_save_metabolite(drugbank_filepath,parsed_pair_name_output,-1 )
    

    parsed_pair_name_output_exogenous = "dataset/raw_data/drugbank_pairs_names_exogenous.csv"
    if(remove_endogenous_pairs):
        filter_endogenous_reaction(parsed_pair_name_output, parsed_pair_name_output_exogenous)

    # get the drug & metabolite information
    ## Metabolites
    drugbank_metabolite_structures_filepath = "dataset/raw_data/drugbank_metabolite_structures.sdf"
    parsed_metabolite_structures_output = "dataset/raw_data/drugbank_metabolite_structures.csv"
    if get_metabolite_info: 
        read_sdf_file(drugbank_metabolite_structures_filepath, parsed_metabolite_structures_output, -1,is_drug_information=False)
        remove_bad_drug_metabolite_rows(parsed_metabolite_structures_output)


    ## Drugs
    drugbank_drug_structures = "dataset/raw_data/drugbank_drug_structures.sdf"
    parsed_drug_structures_output = "dataset/raw_data/drugbank_drug_structures.csv"
    if get_drug_info: 
        read_sdf_file(drugbank_drug_structures, parsed_drug_structures_output, -1, False, 2057177,is_drug_information=True)
        remove_bad_drug_metabolite_rows(parsed_drug_structures_output)
        
    ## external drugs
    drugbank_external_drugs = "dataset/raw_data/drugbank_external_structures.csv"
    drugbank_external_drugs_cleaned = "dataset/raw_data/drugbank_external_structures_cleaned.csv"
    if get_external_drugs:        
        db_external_df = pd.read_csv(drugbank_external_drugs)
        db_external_df = db_external_df[["dbid","smiles","inchi_key",'name']]
        db_external_df.to_csv(drugbank_external_drugs_cleaned,index=False)
        remove_bad_drug_metabolite_rows(drugbank_external_drugs_cleaned)
    
    ## HMDB drug and metabolites
    hmdb_metabolites = "dataset/raw_data/hmdb_metabolites.xml"
    hmdb_cleaned = "dataset/raw_data/hmdb_cleaned.csv"
    if get_hdmb_drugs:        
        hmdb_metabolite_map(hmdb_metabolites,hmdb_cleaned, -1 )
        remove_bad_drug_metabolite_rows(hmdb_cleaned)
    
    ## Combine the drug metabolite into one file
    ## The resulting information files we have are:
    ### Drug, metabolite, hmdb, external
    if combine_the_data:
        drugbank_drugs_only_products = 'dataset/raw_data/drugbank_drugs_only_products.csv'
        if(drop_drugs_with_no_products):
            drop_products_from_drug_information(parsed_drug_structures_output, drugbank_drugs_only_products)
        
        drugbank_full_structures = "dataset/raw_data/drugbank_full_structures.csv"
        combine_datasets(parsed_metabolite_structures_output, drugbank_drugs_only_products, drugbank_full_structures)
        combine_datasets(drugbank_full_structures, drugbank_external_drugs_cleaned, drugbank_full_structures)
        combine_datasets(drugbank_full_structures, hmdb_cleaned, drugbank_full_structures)

        # map reaction pairs to smiles and inchi
        parent_child_smiles_output = "dataset/raw_data/drugbank_pairs_smiles.csv"
        parent_child_error_map = "dataset/raw_data/drugbank_error_map.csv"
        combine_id_with_smiles_inchi(parsed_pair_name_output_exogenous, drugbank_full_structures, parent_child_smiles_output,parent_child_error_map)
    
        ## Check hmdb also
        parent_child_smiles_output_hmdb = "dataset/raw_data/drugbank_pairs_smiles_hmdb.csv"
        parent_child_error_map_hmdb = "dataset/raw_data/drugbank_error_map_hmdb.csv"
        parent_child_smiles_finalized = "dataset/raw_data/drugbank_pairs_smiles_finalized.csv"
        combine_id_with_smiles_inchi(parent_child_error_map, hmdb_cleaned, parent_child_smiles_output_hmdb,parent_child_error_map_hmdb)
        combine_datasets(parent_child_smiles_output, parent_child_smiles_output_hmdb, parent_child_smiles_finalized)
