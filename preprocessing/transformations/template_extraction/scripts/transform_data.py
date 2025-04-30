import sys

import pandas as pd


data = pd.read_csv(sys.argv[1], sep="\t")
pd.DataFrame(
    {
        "id": [f"row{idx}" for idx in range(len(data))],
        "source": "chemformer",
        "date": "2023-05-29",
        "rsmi": data.mapped_rxn,
        "classification": data.NMC,
    }
).to_csv("imported_reactions.csv", sep="\t", index=False)