from rdkit import Chem
import pandas as pd


def load_metxbiodb(file):
    try:
        metxbiodb_df = pd.read_csv(file)
        metxbiodb_df = metxbiodb_df[
            ["substrate_name", "substrate_inchi", "prod_name", "prod_inchi", "enzyme"]
        ]
        return metxbiodb_df

    except FileNotFoundError:
        print("Error: The file ", file, " was  not found. Please check the file path.")
    except Exception as e:
        print(f"An error occurred: {e}")


def metxbiodb_inchi_to_smiles(data, output_file):
    parent_child = []
    counter = 0
    for ind in data.index:
        try:
            parent_mol = Chem.inchi.MolFromInchi(data.loc[ind]["substrate_inchi"])
            child_mol = Chem.inchi.MolFromInchi(data.loc[ind]["prod_inchi"])
            parent_smiles = Chem.rdmolfiles.MolToSmiles(parent_mol)
            child_smiles = Chem.rdmolfiles.MolToSmiles(child_mol)

            parent_name = data.loc[ind]["substrate_name"]
            child_name = data.loc[ind]["prod_name"]
            enzyme = data.loc[ind]["enzyme"]

            parent_child.append(
                [parent_name, parent_smiles, child_name, child_smiles, enzyme]
            )
        except Exception as e:
            counter += 1
            print("ERROR MISSING DATA, number: ", counter)
            print("Exception:", e)
            print(data.loc[ind])

    smiles_metxbiodb_df = pd.DataFrame(parent_child)
    smiles_metxbiodb_df.columns = [
        "parent_name",
        "parent_smiles",
        "child_name",
        "child_smiles",
        "enzymes",
    ]

    smiles_metxbiodb_df.to_csv(output_file, index=False)


def find_drug_origin_metxbiodb(data_file):
    df = pd.read_csv(data_file)

    # Initialize origins as parent_name or 'metxbiodb unknown' if NaN
    df["origin"] = df["parent_name"].fillna("metxbiodb unknown")

    max_iterations = len(df)
    for _ in range(max_iterations):
        update_made = False
        for index, row in df.iterrows():
            current_origin = row["origin"]

            # Check if the current origin is a child elsewhere, update if necessary
            for _, pot_parent_row in df.iterrows():
                if current_origin == pot_parent_row["child_name"]:
                    new_origin = pot_parent_row["origin"]
                    if current_origin != new_origin:
                        df.at[index, "origin"] = new_origin
                        update_made = True

            if current_origin == "metxbiodb unknown":
                current_origin = row["parent_smiles"]
                for _, pot_parent_row in df.iterrows():
                    if current_origin == pot_parent_row["child_smiles"]:
                        new_origin = pot_parent_row["origin"]
                        if current_origin != new_origin:
                            df.at[index, "origin"] = new_origin
                            update_made = True
                        else:
                            df.at[index, "origin"] = current_origin

        if not update_made:
            break

    df["source"] = "MetXBioDB"
    df.to_csv(data_file, index=False)


if __name__ == "__main__":
    dataset = "dataset/raw_data/metxbiodb.csv"
    final = "dataset/extracted_data/metxbiodb_smiles.csv"

    metxbiodb = load_metxbiodb(dataset)
    metxbiodb_inchi_to_smiles(metxbiodb, final)
    find_drug_origin_metxbiodb(final)
