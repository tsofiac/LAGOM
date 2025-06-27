import pandas as pd
from rdkit import Chem
from rdkit import DataStructs
from rdkit.ML.Cluster import Butina
from rdkit.Chem import rdMolDescriptors as rdmd
from tqdm import tqdm


def butina_cluster(mol_list, cutoff=0.8):
    print("Calculating fingerprints...")
    fp_list = [
        rdmd.GetMorganFingerprintAsBitVect(Chem.MolFromSmiles(m), radius=2, nBits=1024)
        for m in tqdm(mol_list)
    ]
    dists = []
    nfps = len(fp_list)
    print("Calculating bulk similarity...")
    for i in tqdm(range(1, nfps)):
        sims = DataStructs.BulkTanimotoSimilarity(fp_list[i], fp_list[:i])
        dists.extend([1 - x for x in sims])
    mol_clusters = Butina.ClusterData(dists, nfps, cutoff, isDistData=True)
    cluster_id_list = [0] * nfps
    for idx, cluster in enumerate(mol_clusters, 1):
        for member in cluster:
            cluster_id_list[member] = idx

    # Calculates number of unique clusters
    num_clusters = max(cluster_id_list)
    print(f"Number of clusters formed: {num_clusters}")
    return cluster_id_list


def count_clusters(input_file):
    df = pd.read_csv(input_file)
    transformations = df["TB_Cluster"]

    transformation_counts = transformations.value_counts()
    num_unique_transformations = transformation_counts.size

    print(f"There are {num_unique_transformations} types of transformations.")

    for transformation, count in transformation_counts.items():
        print(f"Cluster: {transformation}, Count: {count}")


if __name__ == "__main__":
    # Choose if clustering based on children or parents:
    name = "child"  # 'child' 'parent'

    df = pd.read_csv("dataset/curated_data/LAGOM_smiles_clean.csv")
    output_csv = f"dataset/curated_data/tb_output_{name}.csv"

    cid_list = butina_cluster(list(df[f"{name}_smiles"]))
    df.loc[:, "TB_Cluster"] = cid_list
    df.to_csv(output_csv, index=False)

    count_clusters(output_csv)
