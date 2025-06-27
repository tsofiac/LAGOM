from rdkit import Chem, DataStructs
from rdkit.Chem import Descriptors, AllChem
import pandas as pd

from standardize_smiles import standardize_smiles_collection
from sklearn.model_selection import GroupShuffleSplit, train_test_split
import datetime


# --------------------------Standardising SMILES and removing duplicates in each data set -----------------------------------
def standardize_smiles(input_file):
    df = pd.read_csv(input_file)
    df["parent_smiles"] = standardize_smiles_collection(df["parent_smiles"], False)
    df["child_smiles"] = standardize_smiles_collection(df["child_smiles"], False)
    return df


def remove_duplicates(df, duplicates_data_file):
    len_before = len(df)
    duplicates_df = df[
        df.duplicated(subset=["parent_smiles", "child_smiles"], keep=False)
    ]
    df = df.drop_duplicates(subset=["parent_smiles", "child_smiles"], keep="first")
    len_after = len(df)
    print("Total data points removed with remove_duplicates:", len_before - len_after)

    if not duplicates_df.empty:
        duplicates_df.to_csv(duplicates_data_file, index=False)

    return df


def non_equal_smiles(parent_smiles, child_smiles):
    return parent_smiles != child_smiles


def remove_equal_parent_child(data, removed_data_file):
    total_removed = 0

    allowed_molecules = [
        non_equal_smiles(row["parent_smiles"], row["child_smiles"])
        for index, row in data.iterrows()
    ]
    filtered_data = data[allowed_molecules]
    removed_data = data[[not allowed for allowed in allowed_molecules]]
    total_removed += allowed_molecules.count(False)

    if not removed_data.empty:
        removed_data.to_csv(removed_data_file, index=False)

    print(f"Total data points removed with remove_equal_parent_child: {total_removed}")

    return filtered_data


# --------------------------Combining datasets and removing duplicates with each other and with test data set -----------------------------------
def add_source_column(df, source):
    if source == "drugbank":
        df["source"] = "DrugBank"
    elif source == "metxbiodb":
        df["source"] = "MetXBioDB"
    elif source == "gloryx":
        df["source"] = "GLORYx"
    else:
        print["Look over this."]

    return df


def combine_datasets(df1, df2, output_csv):
    selected_df1 = df1.copy()
    selected_df2 = df2.copy()

    combined_df = pd.concat([selected_df1, selected_df2], ignore_index=True)
    combined_df.to_csv(output_csv, index=False)


def remove_duplicates_combined(combined_csv, removed_duplicates_csv):
    combined_df = pd.read_csv(combined_csv)
    len_before = len(combined_df)

    drugbank_df = combined_df[
        combined_df["source"].str.contains("DrugBank", case=False, na=False)
    ]
    metxbiodb_df = combined_df[
        combined_df["source"].str.contains("MetXBioDB", case=False, na=False)
    ]

    dups_drugbank_df = drugbank_df[
        drugbank_df.duplicated(subset=["parent_smiles", "child_smiles"], keep=False)
    ]
    drugbank_df = drugbank_df.drop_duplicates(
        subset=["parent_smiles", "child_smiles"], keep="first"
    )

    dups_metxbiodb_df = metxbiodb_df[
        metxbiodb_df.duplicated(subset=["parent_smiles", "child_smiles"], keep=False)
    ]
    metxbiodb_df = metxbiodb_df.drop_duplicates(
        subset=["parent_smiles", "child_smiles"], keep="first"
    )

    combined_df = pd.concat([drugbank_df, metxbiodb_df])
    removed_df = pd.concat([dups_drugbank_df, dups_metxbiodb_df])

    combined_df["is_duplicate"] = combined_df.duplicated(
        subset=["parent_smiles", "child_smiles"], keep=False
    )
    duplicates_df = combined_df[combined_df["is_duplicate"]]
    combined_df.loc[combined_df["is_duplicate"], "source"] = "Both"
    combined_df = combined_df.drop(columns="is_duplicate")
    duplicates_df = duplicates_df.drop(columns="is_duplicate")
    combined_df = combined_df.drop_duplicates(
        subset=["parent_smiles", "child_smiles"], keep="first"
    )

    full_df = pd.concat([removed_df, duplicates_df])

    len_after = len(combined_df)
    print(
        "Total data points removed due to duplicates in both datasets:",
        len_before - len_after,
    )
    combined_df.to_csv(combined_csv, index=False)
    full_df.to_csv(removed_duplicates_csv, index=False)


def compare_datasets(combined_csv, testdata_csv, removed_file):
    df1 = pd.read_csv(combined_csv)
    df2 = pd.read_csv(testdata_csv)

    duplicate_parent_smiles = set(df1["parent_smiles"]).intersection(
        set(df2["parent_smiles"])
    )

    non_duplicates_in_df1 = df1[
        ~df1["parent_smiles"].isin(duplicate_parent_smiles)
    ].copy()
    duplicates_df = df1[df1["parent_smiles"].isin(duplicate_parent_smiles)].copy()

    print(
        "Number of overlapping reactions with test data removed: ",
        len(df1) - len(non_duplicates_in_df1),
    )
    non_duplicates_in_df1.to_csv(combined_csv, index=False)
    duplicates_df.to_csv(removed_file, index=False)


# -------------------------Splitting the data-------------------------------
def set_distribution(
    data_file, evaluation_csv, val_size, test_size=0, set_random_state=56
):
    df = pd.read_csv(data_file)
    df["set"] = None
    train_val_size = 1 - test_size

    if "origin" in df.columns:
        if test_size != 0:
            splitter = GroupShuffleSplit(
                test_size=test_size, n_splits=1, random_state=56
            )
            train_val_inds, test_inds = next(splitter.split(df, groups=df["origin"]))

            train_val_df = df.iloc[train_val_inds]
            test_df = df.iloc[test_inds]

            train_val_df.loc[:, "set"] = "train/val"
            test_df.loc[:, "set"] = "test"
        else:
            train_val_df = df
            test_df = pd.DataFrame()

        effective_val_size = val_size / train_val_size
        splitter = GroupShuffleSplit(
            test_size=effective_val_size, n_splits=1, random_state=set_random_state
        )
        train_inds, val_inds = next(
            splitter.split(train_val_df, groups=train_val_df["origin"])
        )

        train_df = train_val_df.iloc[train_inds]
        val_df = train_val_df.iloc[val_inds]

        train_df.loc[:, "set"] = "train"
        val_df.loc[:, "set"] = "val"

    else:  # if 'origin' does not exist
        if test_size != 0:
            train_val_df, test_df = train_test_split(
                df, test_size=test_size, random_state=set_random_state
            )
            train_val_df.loc[:, "set"] = "train/val"
            test_df.loc[:, "set"] = "test"
        else:
            train_val_df = df
            test_df = pd.DataFrame()

        effective_val_size = val_size / train_val_size
        train_df, val_df = train_test_split(
            train_val_df, test_size=effective_val_size, random_state=set_random_state
        )
        train_df.loc[:, "set"] = "train"
        val_df.loc[:, "set"] = "val"

    df = pd.concat([train_df, val_df, test_df]).reset_index(drop=True)

    df.to_csv(data_file, index=False)

    if not test_df.empty:
        test_df.to_csv(evaluation_csv, index=False)


# --------------------------Filtering data--------------------------------------
def filter_data_on_both_sides(
    data_file, filter_method, removed_data_file, save_removed=True
):
    data = pd.read_csv(data_file)

    data["child_smiles"] = data["child_smiles"].fillna("").astype(str)
    data["parent_smiles"] = data["parent_smiles"].fillna("").astype(str)

    total_removed = 0

    # Filtering based on child molecules
    allowed_children = [
        filter_method(metabolite) for metabolite in data["child_smiles"]
    ]
    filtered_data = data[allowed_children]
    removed_data_children = data[[not allowed for allowed in allowed_children]]
    total_removed += allowed_children.count(False)

    # Filtering based on parent molecules
    allowed_parents = [
        filter_method(molecule) for molecule in filtered_data["parent_smiles"]
    ]
    removed_data_parents = filtered_data[[not allowed for allowed in allowed_parents]]
    filtered_data = filtered_data[allowed_parents]
    total_removed += allowed_parents.count(False)

    removed_data = pd.concat([removed_data_children, removed_data_parents])

    filtered_data.to_csv(data_file, index=False)

    if save_removed:
        if not removed_data.empty:
            removed_data.to_csv(removed_data_file, index=False)

    print(f"Total data points removed with {filter_method.__name__}: {total_removed}")


def filter_data_on_one_side(
    data_file, filter_method, removed_data_file, if_parent=True, save_removed=True
):
    name_property = "parent_smiles" if if_parent else "child_smiles"
    data = pd.read_csv(data_file)

    data[name_property] = data[name_property].fillna("").astype(str)

    total_removed = 0

    # Filtering based on parent molecules
    allowed_molecules = [filter_method(molecule) for molecule in data[name_property]]
    filtered_data = data[allowed_molecules]
    removed_data = data[[not allowed for allowed in allowed_molecules]]
    total_removed += allowed_molecules.count(False)

    filtered_data.to_csv(data_file, index=False)

    if save_removed:
        if not removed_data.empty:
            removed_data.to_csv(removed_data_file, index=False)

    print(f"Total data points removed with {filter_method.__name__}: {total_removed}")


def valid_smiles(molecule):
    try:
        return Chem.MolFromSmiles(molecule) is not None
    except Exception:
        return False


def atoms_allowed_in_molecules(molecule):
    try:
        atoms_to_include = ["C", "N", "S", "O", "H", "F", "I", "P", "Cl", "Br"]
        mol = Chem.MolFromSmiles(molecule)
        atoms = [atom.GetSymbol() for atom in mol.GetAtoms()]
        return set(atoms).issubset(set(atoms_to_include))
    except Exception:
        return False


def molecule_allowed_based_on_weight(molecule, max_weight=750, min_weight=100):
    try:
        mol_weight = Descriptors.ExactMolWt(Chem.MolFromSmiles(molecule))
        return min_weight <= mol_weight <= max_weight
    except Exception:
        return False


def define_fingerprint_similarity(dataset):
    parent_smiles = dataset["parent_smiles"].tolist()
    child_smiles = dataset["child_smiles"].tolist()

    parent_mol = [Chem.MolFromSmiles(x) for x in parent_smiles]
    child_mol = [Chem.MolFromSmiles(x) for x in child_smiles]

    parent_fps = [
        AllChem.GetMorganFingerprintAsBitVect(x, radius=2, nBits=1024)
        for x in parent_mol
    ]
    child_fps = [
        AllChem.GetMorganFingerprintAsBitVect(x, radius=2, nBits=1024)
        for x in child_mol
    ]

    fingerprint_similarities = []
    for i in range(len(parent_smiles)):
        s = DataStructs.TanimotoSimilarity(parent_fps[i], child_fps[i])
        fingerprint_similarities.append(s)

    dataset["tanimoto"] = fingerprint_similarities

    return dataset


def fingerprints_allowed(similarity, min_similarity=0.20):
    if similarity < min_similarity:
        return False
    return True


def filter_fingerprint_similarity(data_file, removed_data_file, save_removed=True):
    data = pd.read_csv(data_file)
    data = define_fingerprint_similarity(data)

    total_removed = 0

    allowed_molecules = [fingerprints_allowed(value) for value in data["tanimoto"]]
    filtered_data = data[allowed_molecules]
    removed_data = data[[not allowed for allowed in allowed_molecules]]
    total_removed += allowed_molecules.count(False)

    filtered_data.to_csv(data_file, index=False)

    if save_removed:
        if not removed_data.empty:
            removed_data.to_csv(removed_data_file, index=False)

    print(
        f"Total data points removed with fingerprint_similarity_filter: {total_removed}"
    )


# ---------------------Reformat for Chemformer-------------------------------------


def get_unique_parents(input_file, output_file):
    df = pd.read_csv(input_file)
    if "parent_smiles" not in df.columns:
        raise ValueError("The DataFrame must contain a 'parent_smiles' column.")

    unique_parent_dataset = []

    for parent in df["parent_smiles"].unique():
        row = df[df["parent_smiles"] == parent].iloc[0]
        unique_parent_dataset.append(row)

    unique_parent_df = pd.DataFrame(unique_parent_dataset).reset_index(drop=True)
    unique_parent_df = unique_parent_df[df.columns]

    unique_parent_df.to_csv(output_file, index=False)


def reformat_for_chemformer(input_file, output_file):
    df = pd.read_csv(input_file)

    df = df.rename(
        columns={
            "child_smiles": "products",
            "parent_smiles": "reactants",
        }
    )

    df.to_csv(output_file, sep="\t", index=False)


# --------------------- Augmentation -------------------------------------
def augment_drugbank(input_file):
    df_drugbank = pd.read_csv(input_file)

    new_reactions = []
    for index, row in df_drugbank.iterrows():
        # Check if the reactant is different from the origin
        if row["parent_id"] != row["origin"]:
            # Create a new reaction with the origin as the reactant and keep the product same
            new_reaction = {
                "parent_id": row["origin"],
                "child_id": row["child_id"],
                "child_name": row["child_name"],
                "child_smiles": row["child_smiles"],
                "origin": row["origin"],
            }
            new_reactions.append(new_reaction)
    augmented_reactions = pd.DataFrame(new_reactions)

    unique_df_drugbank = df_drugbank[
        ["parent_id", "parent_name", "parent_smiles"]
    ].drop_duplicates(subset="parent_id")
    augmented_reactions = pd.merge(
        augmented_reactions,
        unique_df_drugbank[["parent_id", "parent_name", "parent_smiles"]],
        on="parent_id",
        how="left",
    )
    augmented_reactions["source"] = "DrugBank Augmented"
    desired_column_order = [
        "parent_id",
        "parent_name",
        "child_id",
        "child_name",
        "parent_smiles",
        "child_smiles",
        "origin",
        "source",
    ]
    augmented_reactions = augmented_reactions[desired_column_order]
    return augmented_reactions


def augment_metxbiodb(input_file):
    df_metxbiodb = pd.read_csv(input_file)

    new_reactions = []
    for index, row in df_metxbiodb.iterrows():
        # Check if the reactant is different from the origin
        if row["parent_name"] != row["origin"]:
            # Create a new reaction with the origin as the reactant and keep the product same
            new_reaction = {
                "parent_name": row["origin"],
                "child_name": row["child_name"],
                "child_smiles": row["child_smiles"],
                "origin": row["origin"],
            }
            new_reactions.append(new_reaction)
    augmented_reactions = pd.DataFrame(new_reactions)

    unique_df_metxbiodb = df_metxbiodb[
        ["parent_name", "parent_smiles"]
    ].drop_duplicates(subset="parent_name")
    augmented_reactions = pd.merge(
        augmented_reactions,
        unique_df_metxbiodb[["parent_name", "parent_smiles"]],
        on="parent_name",
        how="left",
    )
    augmented_reactions["source"] = "MetXBioDB Augmented"
    desired_column_order = [
        "parent_name",
        "child_name",
        "parent_smiles",
        "child_smiles",
        "origin",
        "source",
    ]
    augmented_reactions = augmented_reactions[desired_column_order]
    return augmented_reactions


def join(df1, df2, output_file):
    combined = pd.concat([df1, df2], ignore_index=True)
    combined.to_csv(output_file, index=False)


def select_reactions(input_file, comparison_file):
    df = pd.read_csv(input_file)
    comparison_df = pd.read_csv(comparison_file)

    train_df = comparison_df[comparison_df["set"] == "train"]
    filtered_df = df[df["origin"].isin(train_df["origin"])].copy()

    filtered_df.loc[:, "set"] = "train"

    filtered_df.to_csv(input_file, index=False)


def parent_to_parent(input_file, output_file):
    df = pd.read_csv(input_file)

    new_df = pd.DataFrame(
        {
            "parent_name": df["parent_name"],
            "parent_smiles": df["parent_smiles"],
            "child_name": df["parent_name"],
            "child_smiles": df["parent_smiles"],
            "origin": df["origin"],
            "set": df["set"],
        }
    )

    new_df["source"] = "Parent-Parent Augmented"

    new_df.to_csv(output_file, index=False)


def add_possible_products(input_file):
    df = pd.read_csv(input_file, sep="\t")

    def agg_order_preserved(series):
        return list(series)

    grouped = (
        df.groupby(["reactants"], sort=False)
        .agg({"products": agg_order_preserved})
        .reset_index()
    )

    grouped.rename(columns={"products": "possible_products"}, inplace=True)

    df = df.merge(
        grouped[["reactants", "possible_products"]], on="reactants", how="left"
    )

    df["possible_products"] = df["possible_products"].apply(lambda x: ".".join(x))

    df.to_csv(input_file, sep="\t", index=False)


# ------------Time log------------------
def log_time(message):
    # Get the current time
    current_time = datetime.datetime.now()
    # Format the time (optional, for easier reading)
    formatted_time = current_time.strftime("%Y-%m-%d %H:%M:%S")
    # Print or log the time with a message
    print(f"[{formatted_time}] {message}")


if __name__ == "__main__":
    """ rows for VA """
    start_row = 0
    end_row = 1101304

    # start_row = 1101304
    # end_row = 2202607

    # start_row = 2202607
    # end_row = 3303912

    # start_row = 3303912
    # end_row = 4405216

    # start_row = 4405216
    # end_row = 5506519

    # start_row = 5506519
    # end_row = 6607824

    # start_row = 6607824
    # end_row = 7709127

    # start_row = 7709127
    # end_row = 8810431

    # start_row = 8810431
    # end_row = 9911735

    # start_row = 9911735
    # end_row = None #11013037

    name = "LAGOM"  # [ 'LAGOM' 'VA_filter_part'(uncomment a start_row and corresponding end_row) 'VA_last_filtering' 'metatrans']

    clean_csv = f"dataset/curated_data/{name}_smiles_clean.csv"
    dataset_gloryx = "dataset/curated_data/gloryx_smiles_clean.csv"

    removed_duplicates = f"dataset/removed_data/{name}_removed_duplicates.csv"
    removed_equal = f"dataset/removed_data/{name}_removed_equal.csv"
    removed_valid_smiles = f"dataset/removed_data/{name}_removed_valid_smiles.csv"
    removed_atoms_allowed = f"dataset/removed_data/{name}_removed_atoms_allowed.csv"
    removed_weights_allowed = f"dataset/removed_data/{name}_removed_weights_allowed.csv"
    removed_fingerprints = f"dataset/removed_data/{name}_removed_fingerprints.csv"

    evaluation_csv = f"dataset/curated_data/{name}_evaluation.csv"
    evaluation_unique_csv = f"dataset/curated_data/{name}_evaluation_unique.csv"
    evaluation_finetune_csv = f"dataset/finetune/{name}_evaluation_finetune.csv"
    finetune_csv = f"dataset/finetune/{name}_finetune.csv"

    if name == "LAGOM":
        val_size = 0.1  # val
        eval_size = 0.05  # test

        print("Starting preprocessing of LAGOM dataset")
        metxbiodb_csv = "dataset/extracted_data/metxbiodb_smiles.csv"
        drugbank_csv = "dataset/extracted_data/drugbank_smiles.csv"
        LAGOM_csv = "dataset/curated_data/LAGOM_smiles_clean.csv"

        df_metx = pd.read_csv(metxbiodb_csv)
        df_drugbank = pd.read_csv(drugbank_csv)

        combine_datasets(df_drugbank, df_metx, LAGOM_csv)

        LAGOM_df = standardize_smiles(LAGOM_csv)
        LAGOM_df.to_csv(LAGOM_csv, index=False)
        remove_duplicates_combined(LAGOM_csv, removed_duplicates)
        LAGOM_df = pd.read_csv(LAGOM_csv)
        LAGOM_df = remove_equal_parent_child(LAGOM_df, removed_equal)
        LAGOM_df.to_csv(LAGOM_csv, index=False)

        compare_datasets(
            LAGOM_csv,
            dataset_gloryx,
            "dataset/removed_data/LAGOM_compare_gloryx_removed_duplicates.csv",
        )

        LAGOM_df = pd.read_csv(LAGOM_csv)
        LAGOM_df.to_csv("dataset/removed_data/LAGOM_smiles_before_filter.csv", index=False)

        filter_data_on_both_sides(LAGOM_csv, valid_smiles, removed_valid_smiles)
        filter_data_on_both_sides(
            LAGOM_csv, atoms_allowed_in_molecules, removed_atoms_allowed
        )
        filter_data_on_one_side(
            LAGOM_csv,
            molecule_allowed_based_on_weight,
            removed_weights_allowed,
            True,
        )
        filter_fingerprint_similarity(LAGOM_csv, removed_fingerprints, False)

        set_distribution(clean_csv, evaluation_csv, val_size, eval_size)
        get_unique_parents(evaluation_csv, evaluation_unique_csv)
        reformat_for_chemformer(evaluation_unique_csv, evaluation_finetune_csv)

        # ------- parent - grandchild --------------
        print("Starting preprocessing of parent-grandchild")
        parent_grandchild = "dataset/curated_data/augmented_parent_grandchild.csv"
        augmented_drugbank = augment_drugbank(drugbank_csv)
        augmented_metxbiodb = augment_metxbiodb(metxbiodb_csv)
        join(augmented_drugbank, augmented_metxbiodb, parent_grandchild)

        parent_grandchild_df = standardize_smiles(parent_grandchild)
        parent_grandchild_df = remove_duplicates(
            parent_grandchild_df,
            "dataset/removed_data/augmented_removed_duplicates.csv",
        )
        parent_grandchild_df = remove_equal_parent_child(
            parent_grandchild_df, "dataset/removed_data/augmented_removed_equal.csv"
        )

        parent_grandchild_df.to_csv(parent_grandchild, index=False)

        filter_data_on_both_sides(
            parent_grandchild,
            valid_smiles,
            "dataset/removed_data/augmented_removed_valid_smiles.csv",
        )
        filter_data_on_both_sides(
            parent_grandchild,
            atoms_allowed_in_molecules,
            "dataset/removed_data/augmented_removed_atoms_allowed.csv",
        )
        filter_data_on_one_side(
            parent_grandchild,
            molecule_allowed_based_on_weight,
            "dataset/removed_data/augmented_removed_weights_allowed.csv",
            True,
        )
        filter_fingerprint_similarity(
            parent_grandchild, "dataset/removed_data/augmented_removed_fingerprints.csv"
        )

        select_reactions(parent_grandchild, LAGOM_csv)
        print("Finishing preprocessing of parent-grandchild")
        # --------------------------------------

        # ------- parent - parent --------------
        print("Starting preprocessing of parent-parent")
        parent_parent = "dataset/curated_data/augmented_parent_parent.csv"
        LAGOM_df = pd.read_csv(LAGOM_csv)
        train_df = LAGOM_df[LAGOM_df["set"] == "train"]
        train_df.to_csv(parent_parent, index=False)
        get_unique_parents(parent_parent, parent_parent)
        parent_to_parent(parent_parent, parent_parent)
        print("Finishing preprocessing of parent-parent")
        # --------------------------------------

        reformat_for_chemformer(LAGOM_csv, finetune_csv)
        add_possible_products(finetune_csv)

    if "VA" in name:
        val_size = 0.005  # val
        eval_size = 0  # test

        print("VA")

        if name == "VA_last_filtering":
            clean_csv = "dataset/VA/VA_smiles_clean.csv"
            finetune_csv = "dataset/VA/VA_finetune.csv"

            print("Comparing with Gloryx")
            compare_datasets(
                clean_csv,
                dataset_gloryx,
                "dataset/VA/VA_removed/VA_compare_gloryx_removed_duplicates.csv",
            )
            log_time("Filtered on Gloryx")

            print("Comparing with Metabolic dataset")
            compare_datasets(
                clean_csv,
                "dataset/curated_data/LAGOM_smiles_clean.csv",
                "dataset/VA/VA_removed/VA_compare_LAGOM_removed_duplicates.csv",
            )
            log_time("Filtered on LAGOM dataset")

            print("Removing duplicates")
            df = pd.read_csv(clean_csv)

            log_time("Begin filtering on duplicates")
            df = remove_duplicates(df, removed_duplicates)
            log_time("Duplicates removed")
            df.to_csv(clean_csv, index=False)
            log_time("New csv created")
            log_time("Setting training and validation distributions")
            set_distribution(clean_csv, evaluation_csv, val_size, eval_size)
            reformat_for_chemformer(clean_csv, finetune_csv)
            log_time("Done.")

        elif name == "VA_filter_part":
            dataset = f"dataset/curated_data/VA_rows_{start_row}_to_{end_row}.csv"

            log_time("Begin filtering")
            df = standardize_smiles(dataset)
            log_time("Smiles are standardised")
            df = remove_equal_parent_child(df, removed_equal)
            log_time("Equal_parent_child removed")
            df.to_csv(clean_csv, index=False)

            filter_data_on_both_sides(
                clean_csv, valid_smiles, removed_valid_smiles, save_removed=False
            )
            log_time("Filtered valid smiles")
            filter_data_on_both_sides(
                clean_csv,
                atoms_allowed_in_molecules,
                removed_atoms_allowed,
                save_removed=False,
            )
            log_time("Filtered atoms allowed")
            filter_data_on_one_side(
                clean_csv,
                molecule_allowed_based_on_weight,
                removed_weights_allowed,
                True,
                save_removed=False,
            )
            log_time("Filtered on weight")
            filter_fingerprint_similarity(clean_csv, removed_fingerprints)
            log_time("Filtered on fingerprint similarity")

        else:
            print("Incorrect task name for VA")

    if name == "metatrans":
        val_size = 0.005  # val
        eval_size = 0  # test

        dataset = "dataset/extracted_data/metatrans_smiles.csv"

        log_time("Begin filtering metatrans")
        df = standardize_smiles(dataset)
        log_time("Smiles are standardised")
        df = remove_duplicates(df, removed_duplicates)
        log_time("Duplicates removed")
        df = remove_equal_parent_child(df, removed_equal)
        log_time("Equal_parent_child removed")
        df.to_csv(clean_csv, index=False)

        filter_data_on_both_sides(
            clean_csv, valid_smiles, removed_valid_smiles, save_removed=False
        )
        log_time("Filtered valid smiles")
        filter_data_on_both_sides(
            clean_csv,
            atoms_allowed_in_molecules,
            removed_atoms_allowed,
            save_removed=False,
        )
        log_time("Filtered atoms allowed")
        filter_data_on_one_side(
            clean_csv,
            molecule_allowed_based_on_weight,
            removed_weights_allowed,
            True,
            save_removed=False,
        )
        log_time("Filtered on weight")
        filter_fingerprint_similarity(clean_csv, removed_fingerprints)
        log_time("Filtered on fingerprint similarity")

        compare_datasets(
            clean_csv,
            dataset_gloryx,
            "dataset/removed_data/metatrans_compare_gloryx_removed_duplicates.csv",
        )
        log_time("Filtered on Gloryx")

        compare_datasets(
            clean_csv,
            "dataset/curated_data/LAGOM_smiles_clean.csv",
            "dataset/removed_data/metatrans_compare_LAGOM_removed_duplicates.csv",
        )
        log_time("Filtered on Metabolic dataset")

        set_distribution(clean_csv, evaluation_csv, val_size, eval_size)
        log_time("Distribution is set")

        reformat_for_chemformer(clean_csv, finetune_csv)
        log_time("Reformated for Chemformer")

        add_possible_products(finetune_csv)
        log_time("Added possible products")
