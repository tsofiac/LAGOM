# LAGOM: Language-model-Assisted Generation Of Metabolites 

This repository contains the code used to generate the results reported in the [LAGOM article](https://doi.org/10.26434/chemrxiv-2025-0mqf8).

LAGOM (Language-model Assisted Generation Of Metabolites) is a Transformer-based [Chemformer model](https://github.com/MolecularAI/chemformer) fine-tuned to predict likely metabolic transformations of drug candidates. As part of this work, we assembled a rigorously curated and standardised collection of publicly available datasets ([DrugBank](https://go.drugbank.com/releases/latest) and [MetXBioDB](https://zenodo.org/records/13235312)) specifically for metabolite prediction. To enhance model generalisation and predictive accuracy, we incorporated a range of data augmentation strategies. In addition, we explored the impact of different pre-training approaches, including general chemical pre-training using the [Virtual Analogs dataset](https://zenodo.org/records/45807) and metabolite-specific pre-training with [MetaTrans dataset](https://github.com/KavrakiLab/MetaTrans/tree/master/datasets).

## Table of Contents

- [Installation](#installation)
- [Data](#data)
- [Running Chemformer](#running-chemformer)
- [Augmentation and Annotation](#augmentation-and-annotation)
- [Splitting Data for an Ensemble Model](#splitting-data-for-an-ensemble-model)
- [Evaluation](#evaluation)
- [Visualisations](#visualisations)
- [Contributors](#contributors)
- [License](#license)

## Installation
Here are the instructions on how to set up this project.

First, set up Chemformer by following the instructions given in the repository [aizynthmodels](https://github.com/MolecularAI/aizynthmodels/). This includes cloning the repository using Git and creating the conda environment called ``aizynthmodels``. The release version used in this project was ``aizynthmodels 1.0.0`` and we used Python 3.10.16.

Then clone this repository using Git: 

```bash
git clone https://github.com/tsofiac/LAGOM.git 
cd LAGOM/
```

Download this additional package:

```bash
conda activate aizynthmodels
pip install chembl_structure_pipeline==1.2.2
```

The pre-trained Chemformer model used in this project is available via a link provided in the README file in the [chemformer](https://github.com/MolecularAI/chemformer) repository. We used the public model ``models/pre-trained/combined/step=1000000.ckpt`` with updated checkpoints.

## Data 
 Before proceeding with data extraction and curation, create a folder named `dataset` and, within it, the following subfolders:

- `raw_data`: Place the raw and unmodified data files here.
- `extracted_data`: This is where data files will be located after running the load scripts specified later.
- `curated_data`: Location for data files obtained after the curation process.
- `fine_tuning`: Location for data files that will be used for fine-tuning. These are also created after the curation process.
- `removed_data`: For the data that is removed during the curation process.
- `VA`: For processed files of the Virtual Analogs (VA) dataset. This dataset gets a separate folder, as these files are too large to be committed to GitHub.
  - `VA_removed`: Subfolder within the `VA` folder for removed data associated with the VA dataset.

### Extracting Data

To extract the data needed for this project, download the raw files from the webpages listed below. The versions stated in the list have been used for this project. Datasets from DrugBank and MetXBioDB are used for fine-tuning, the GLORYx test dataset is used for external evaluation, and the Virtual Analogs (VA) and MetaTrans datasets are used for additional pre-training. All datasets are publically available, but a license was required to obatin data from DrugBank. For this project, a free academic license was used. 

* [DrugBank](https://go.drugbank.com/releases/latest) (Version 5.1.13) 
  * drugbank_drug_structures.sdf
  * drugbank_metabolite_structures.sdf
  * drugbank_external_structures.csv
  * drugbank_full_database.xml
* [MetXBioDB](https://zenodo.org/records/13235312) (Version NORMAN-SLE-S73.0.1.7)
  * metxbiodb.csv
* [GLORYx Test Dataset](https://github.com/christinadebruynkops/GLORYx/tree/master/datasets/test_dataset)
  * gloryx_test_dataset.json
* [Virtual Analogs Dataset](https://zenodo.org/records/45807) (Version v1)
  * 1297204_Virtual_Analogs.dat
* [ChEMBL35](https://chembl.gitbook.io/chembl-interface-documentation/downloads)
  * chembl_35_chemreps.txt
* [MetaTrans Dataset](https://github.com/KavrakiLab/MetaTrans/tree/master/datasets)
  * source_train.txt
  * target_train.txt
  * source_valid.txt
  * target_valid.txt

To extract the data, run the following files in this repository for each corresponding data source. The extracted data files will be saved in `dataset/extracted_data`.

* `preprocessing/load_drugbank.py` for the DrugBank dataset
* `preprocessing/load_metxbiodb.py` for the MetXBioDB dataset
* `preprocessing/load_gloryx.py` for the GLORYx test dataset
* `optimisation/load_metatrans.py` for the MetaTrans dataset
* Due to the large size of the VA dataset, it is loaded in parts:
  * Open `preprocessing/load_VA_parts.py` and uncomment a `start_row` and corresponding `end_row`
  * Run the shell script `preprocessing/submit_load_mmp.sh`
  * Repeat this for all start- and end-rows.


### Preprocess Fine-Tuning Data

To preprocess the data, use the file `preprocessing/preprocess_data.py`.

To combine and preprocess the DrugBank and MetXBioDB datasets, set the `name` parameter to `'LAGOM'`. LAGOM is the name of the combined fine-tuning dataset. This will generate the file for finetuning (``LAGOM_finetune.csv``), along with the other files, including an evaluation file and augmentation files. 

(Observe that the GLORYx test dataset is not preprocessed.) 

### Preprocess the VA Dataset for Additional Pre-Training

For complete curation and filtering of the VA dataset, the following steps should be conducted.

#### Initial Filtering

The intital filtering is done separately for each split of the extracted dataset. 

* Go to file `preprocessing/preprocessing_data.py`
* Set `name = 'VA_filter_part'` and uncomment a `start_row` and corresponding `end_row` (an example is provided below):
```bash
if __name__ == "__main__":

    ''' rows for VA '''
    start_row = 0
    end_row = 1101304

    # start_row = 1101304
    # end_row = 2202607

    ...

    name = 'VA_filter_part'
```
* Run file `preprocessing/submit_preprocess.sh`
* Repeat the steps above for all splits
* Run the file `optimisation/combine_VA.py` to combine the different VA files into one

#### Final Filtering

* Go back to the file `preprocessing/preprocessing_data.py`
* Set `name = 'VA_last_filtering'`
* Run the shell script `preprocessing/submit_load_mmp.sh`

This outputs the file `VA/VA_finetune.csv` that can be used for additionally pre-training the base Chemformer model.

### Preprocess the MetaTrans Dataset for Additional Pre-Training

For filtering of the MetaTrans dataset, proceed with the following steps.

* Open file `preprocessing/preprocessing_data.py`
* Set `name = 'metatrans'`
* Run the script (no shell script needed)

This outputs the file `finetune/metatrans_finetune.csv` that can be used for additionally pre-training the ChemVA model.

## Running Chemformer

Chemformer model training involves three main steps: **pre-training**, **fine-tuning**, and **inference score**. Below, each step is described using an example. Note that the dataset used in either of these steps need to be tab-separated, with the columns ``reactants`` and ``products``, which corresponds to the datasets in ``dataset/finetune``. For more detailed information about Chemformer, see the [chemformer](https://github.com/MolecularAI/chemformer) repository.

### Pre-Training

Pre-training allows the model to learn chemical transformations from a large dataset. The following scripts are used:

* `pretrain_bart.yaml`: Hyperparameter settings for pre-training
* `pretrain.sh`: SBATCH script for submitting the job
* `bart_vocab_downstream.json`: Vocabulary file
* `run_pretrain.sh`: Script to execute pre-training


In order to additionally pre-train a model:

1. Edit `run_pretrain.sh` as needed, specifying:

    * Dataset file (e.g., `VA_finetune.csv` or `metatrans_finetune.csv`)
    * An already pre-trained model (ckpt-file), or `null` if starting from scratch
    * Batch size (e.g., ``128``)
    * Randomisation probability (e.g., ``0.5``; can be changed between 0 and 1)
    * Masking probability (e.g., ``0.1``; can be changed between 0 and 1)

2. Run the script:

    ```bash
    ./run_pretrain.sh
    ```

### Fine-Tuning

Fine-tuning adapts a pre-trained model to a specific task, such as metabolite prediction. The required scripts are analogous to pre-training:

* `fine_tune_bart.yaml`: Hyperparameters for fine-tuning
* `fine_tune.sh`: SBATCH script
* `bart_vocab_downstream.json`: Vocabulary file
* `run_fine_tuning.sh`: Script to run fine-tuning


In order to fine-tune the model for the task of metabolite prediction:

1. Edit `run_fine_tuning.sh`:

    * Fine-tunig dataset file (e.g., `LAGOM_finetune.csv`)
    * A pre-trained model (ckpt-file)
    * Batch size (e.g., ``64``)

2. Run the script:

    ```bash
    ./run_fine_tuning.sh
    ```

### Evaluation (Inference Score)

To evaluate a trained model, the following scripts are used:

* `inference_score.yaml`: Hyperparameters for inference scoring
* `inference_score.sh`: SBATCH script
* `bart_vocab_downstream.json`: Vocabulary file
* `run_inference_scoring.sh`: Script for evaluation


To evaluate the model, the `run_pretrain.sh` is adjusted accordingly:

1. Edit `run_pretrain.sh`:

    * Set the evaluation dataset (e.g., `LAGOM_evaluation_finetune.csv` or `gloryx_finetune.csv`)
    * A pre-trained model (ckpt-file)
    * Batch size (e.g., ``38``)
    * Beam search parameter for the number of predictions per input (e.g., ``20``)

2. Run the script:

    ```bash
    ./run_inference_scoring.sh
    ```

## Augmentation and Annotation  

To augment the training dataset with parent-parent or parent-granchild reactions, or to annotate the dataset with either Csp3 fraction or LogP, proceed with the steps below.

* Open file `optimisation/augmentation_annotation.py`
* Set the relevant annotation or augmentataion techniques to *True*. The example below shows how to augment the data with parent-parent reactions.
  * Both annotation techniques can be applied together, and both augmentation techniques (parent-parent and parent-grandchild reactions) can also be used together. However, note that augmentation and annotation are conducted separately.

```bash
if __name__ == "__main__":
    logp_annotations = False
    csp3_annotations = False

    augment_parent_grandchild = False
    augment_parent_parent = True
```
  
* Run the script

Note: Specific data files are required for the augmentation techniques (namely `curated_data/augmented_parent_grandchild.csv` and `curated_data/augmented_parent_parent.csv`). Please ensure these files are available before proceeding with the augmentation process. They are generated when running `preprocessing/preprocessing_data.py` when the `name` parameter is set to `'LAGOM'`.

For inference scoring a model that has been fine-tuned on annotated data, a corresponing annotated test set should to be used. By default, if you select an annotation technique in the script, an annotated version of the LAGOM test set will also be created, e.g. `dataset/finetune/logp_evaluation_finetune.csv`

## Splitting Data for an Ensemble Model

The data can be split at random (stratified split), or split based the Butina clustering algorithm (parent or child similarity). 

To split based on similarity (pre-requisite for child-/parent-based splitting) the following preparation steps should be conducted.

* Open the file `optimisation/tb_cluster.py`
* Set `name = 'child'` or `'parent'`, depending on if you want to cluster based on the metabolites or drugs, respectively
* Run the script

For all different splitting methods,

* Go to the file `optimisation/split_for_ensemble.py`
* Set `split_type` to the splitting method of choice, i.e., `'stratified'`, `'parents'` or `'children'`

The separate files with the different splits, e.g. `dataset/finetune/strat_split1_finetune.csv` can then be used for individual fine-tuning. 

## Evaluation

To read and evaluate the predictions obtained from the inference scoring, the file `evaluation/prediction_reader.py` can be used. This is where scores such as recall and precision are obtained.

All main choices are specified in these lines:

```bash
if __name__ == "__main__":
    benchmark = False  # True if GLORYx, False if Evaluation
    status = 'score'  # 'score' 'combine' 'new'
    name = 'chemVA-Met_base_05'

    specification = 0  # 0 (all) 1 (only_child) 2 (more than 1) 3 (more than 2)
    fingerprint = None  # 1 (similarity = 1) 0.8 (similarity >= 0.8) None (exact SMILES string)
```
Note: Per default, the file `results/evaluation/scores/predictions0.json` is saved when running inference scoring. Thus, if no adjustments are made, `evaluation/prediction_reader.py` should be used before running a new inference scoring.

If predictions have been made on the LAGOM test set (during infrence scoring), set `benchmark` to `False`. If they have been made on the GLORYx test set, set it to `True`.

For evaluating an inference scoring for the first time, set `status = 'new'` and choose a suitable name which will be used when generating the new files, i.e.`f"evaluation/scores/predictions_{name}.csv"` and `f"evaluation/scores/result_{name}.csv"`. Unless interested in specific predictions, these files will not have to be opened; all scores are printed. For printing the scores again at a later point, set `status = 'score'` and use the same name as previously. 

`specification` should be set to `0` to evaluate predictions for all drugs. However, if interested in how the model performed only on for example drugs with two or more children, `specification` can be set to `2`. 

Set `fingerprint` to `None` to calculate scores based on exact matches with the actual predictions. To instead count predictions as correct when their Morgan fingerprint Tanimoto similarity is 0.8 or higher, set `fingerprint` to `0.8`.

### Evaluation of an Ensemble Model

To evaluate an ensemble model, each part of th ensemble model has to be evaluated separately (`status = 'new'`). After running these evaluations, specify the result files in `ensemble_list` as shown below with an example.  Also note the setting `samples_per_model` which determines how many predictions per drug from each model are combined into the ensemble. However, keep in mind that if you use four models with five predictions per drug each, the total number of predictions per drug may be less than 20 (= 4 x 5), as duplicate predictions are removed. Before running the code to combine the models, set status to `combine`.

```bash
 ensemble_list = [
        "evaluation/scores/ensemble/result_ChemVA-Met_Ssplit1.csv",
        "evaluation/scores/ensemble/result_ChemVA-Met_Ssplit2.csv",
        "evaluation/scores/ensemble/result_ChemVA-Met_Ssplit3.csv",
        "evaluation/scores/ensemble/result_ChemVA-Met_Ssplit4.csv",
    ]
samples_per_model = 5
```


## Visualisations

The folder ``for_article`` contains scripts for generating figures and tables for the article.

* `for_article/data_analysis_plots.ipynb`: Data analysis of the datasets and training progress
* `for_article/boxplots_article.ipynb`: Results from evaluation



## Contributors

* Sofia Larsson: [@tsofiac](https://github.com/tsofiac)
* Miranda Carlsson: [@mirandacarlsson](https://github.com/mirandacarlsson)
* Richard Beckmann: [@corteswain](https://github.com/corteswain)
* Filip Miljković: [@filipm90](https://github.com/filipm90)
* Rocío Mercado Oropeza: [@rociomer](https://github.com/rociomer)

## License 

The software is licensed under the Apache 2.0 license (see LICENSE file), and is free and provided as-is.


