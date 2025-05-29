# Master's Thesis: Optimising the Chemformer Model for Metabolite Prediction in Drug Discovery

The aim of this Master’s thesis is to optimise the Chemformer model to reliably
predicting drug metabolites.

## Table of Contents

- [Installation](#installation)
- [Extracting Data](#extracting-data)
- [Preprocess Data](#preprocess-data)
- [Pre-training](#pre-training)
- [Augmentation and Annotation](#augmentation-and-annotation)
- [Splitting Data](#splitting-data)
- [Evaluation](#evaluation)
- [Contributing](#contributing)
- [License](#license)
- [Contact](#contact)

## Installation
Here are the instructions on how to install and set up this project.

First you need set up the Chemformer by following the instructions given the repository [Aizynthmodels](https://github.com/MolecularAI/aizynthmodels/tree/main).

To set up this project you need to clone this repository:

```bash
git clone https://github.com/yourusername/your-project.git
cd your-project
```

You also need to download an extra package:

```bash
conda activate aizynthmodels
conda install chembl_structure_pipeline
```

You also need the pre-trained Chemformer model

## Extracting Data

To extract the data needed for this project you need to download the raw files from the following webpages:

* [Drugbank](https://go.drugbank.com/releases/latest)
  * drugbank_drug_structures.sdf
  * drugbank_metabolite_structures.sdf
  * drugbank_external_structures.csv
  * drugbank_full_database.xml
* [MetXBioDB](https://zenodo.org/records/13235312)
  * metxbiodb.csv
* [GLORYx Test Dataset](https://github.com/christinadebruynkops/GLORYx/tree/master/datasets/test_dataset)
  * gloryx_test_dataset.json

To extract the data you need to run the following files:

```bash
load_drugbank.py
load_metxbiodb.py
load_gloryx.py
```

## Preprocess Data

To preprocess the data, you run the following file:

```bash
preprocess_data.py
```

## Running Chemformer

Here we need to explain how to fine-tune a model

```bash
run_fine_tuning.sh
run_inference_score.sh
```


## Pre-training

To extract the data needed for pre-training the model, you need to download the raw files from the following webpages:

* [Virtual Analogues Dataset](https://zenodo.org/records/45807)
  * 1297204_Virtual__Analogs.dat
* [ChEMBL](https://chembl.gitbook.io/chembl-interface-documentation/downloads)
  * chembl_35_chemreps.txt
* [MetaTrans Dataset](https://github.com/KavrakiLab/MetaTrans/tree/master/datasets)
  * source_train.txt
  * target_train.txt
  * source_train.txt
  * target_train.txt

To extract the data from the files, you run the following files:

```bash
load_mmp.py
load_metatrans.py
```

```bash
preprocess_data.py
```

## Augmentation and Annotation

To augment the training dataset with parent-parent or parent-granchild reactions, or to annotate the dataset with either csp3 fraction or logp the following file should be run:

```bash
augmentation_annotation.py
```

## Splitting Data

For splitting the data at random the following file should be run:

```bash
split_for_ensemble.py
```

To split the data based on the Burtina Clustering algorithm, an additional file need to be run beforehand:

```bash
tb_cluster.py
split_for_ensemble.py
```

## Evaluation

To evaluate the result files obtained from the inference score, the following file can be run:

```bash
prediction_reader.py
```

## Contributing

Miranda Carlsson

Sofia Larsson

Filip Miljkovic

Richard Beckmann

Rocío Mercado Oropeza

Previous Master's Thesis

Chemformer


