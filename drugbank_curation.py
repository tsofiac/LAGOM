import xml.etree.ElementTree as ET
import pandas as pd



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


def remove_bad_drug_metabolite_rows(input_file): # removes row if both identifiers are missing or if both smiles and inchi_key are missing
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


if __name__ == "__main__":

    get_drug_info = False # funkar
    get_metabolite_info = False # funkar inte
    get_external_drugs = False # funkar


    ## Drug structure - drugbank_drug_structures.sdf
    drugbank_drug_structures = "dataset/raw_data/drugbank_drug_structures.sdf"
    parsed_drug_structures_output = "dataset/raw_data/drugbank_drug_structures.csv"
    if get_drug_info: 
        read_sdf_file(drugbank_drug_structures, parsed_drug_structures_output, -1, False, 2057177,is_drug_information=True)
        remove_bad_drug_metabolite_rows(parsed_drug_structures_output)


    ## Metabolite structure - drugbank_metabolite_structures.sdf
    drugbank_metabolite_structures_filepath = "dataset/raw_data/drugbank_metabolite_structures.sdf"
    parsed_metabolite_structures_output = "dataset/raw_data/drugbank_metabolite_structures.csv"
    if get_metabolite_info: 
        read_sdf_file(drugbank_metabolite_structures_filepath, parsed_metabolite_structures_output, -1, False, 514046, is_drug_information=False)
        remove_bad_drug_metabolite_rows(parsed_metabolite_structures_output)
            

    ## External drugs - drugbank_external_structures.csv
    drugbank_external_drugs = "dataset/raw_data/drugbank_external_structures.csv"
    drugbank_external_drugs_cleaned = "dataset/raw_data/drugbank_external_structures_cleaned.csv"
    if get_external_drugs:
        db_external_df = pd.read_csv(drugbank_external_drugs)
        db_external_df = db_external_df[["DrugBank ID","SMILES","InChIKey",'Name']]
        db_external_df = db_external_df.rename(columns={
            "DrugBank ID": "dbid",
            "SMILES": "smiles",
            "InChIKey": "inchi_key",
            "Name": "name"
        })
        db_external_df.to_csv(drugbank_external_drugs_cleaned,index=False)
        remove_bad_drug_metabolite_rows(drugbank_external_drugs_cleaned)
