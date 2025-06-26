# LAGOM: Language-model-Assisted Generation Of Metabolites 

The aim of this work is to optimise the Chemformer model to reliably predict drug metabolites. <!-- A bit longer summary of this methodological framework would be appropriate. You can check other similar packages (e.g., MetaTrans) for inspiration. This will also be followed by the link to the pre-print article/manuscript (once out). -->

<!-- In addition, how to pre-train/fine-tune models and also run them is missing. This information may be contained in the original Chemformer repository, but it would be valuable to provide instructions there as well once a user extracts, preprocesses and splits the data for the model. -->

## Table of Contents

- [Installation](#installation)
- [Extracting Data](#extracting-data)
- [Preprocess Data](#preprocess-data)
- [Pre-training](#pre-training)
- [Augmentation and Annotation](#augmentation-and-annotation)
- [Splitting Data for Ensemble Model](#splitting-data-for-ensemble-model)
- [Evaluation](#evaluation)
- [Contributing](#contributing) <!-- Contributors. -->
- [License](#license)
- [Contact](#contact)

## Installation
Here are the instructions on how to install and set up this project.

First, set up the Chemformer by following the instructions given in the repository [Aizynthmodels](https://github.com/MolecularAI/aizynthmodels/tree/main).

To set up this project you need to clone this repository: <!-- For this, you can provide actual link to the repository to be cloned instead of your-project. -->

```bash
git clone https://github.com/yourusername/your-project.git 
cd your-project
```

Download this extra package: <!-- For reproducibility of your code, you should list exact versions of aizynth and chembl_structure_pipeline. Typically, all required Python dependencies can be provided in .yaml file but you can list packages as well (also to avoid any potential conflicts between versions of different packages). -->

```bash
conda activate aizynthmodels
conda install chembl_structure_pipeline
```

You also need the pre-trained Chemformer model. <!-- How to access them? -->

## Extracting Data

To extract the data needed for this project, download the raw files from the following webpages: 

* [Drugbank](https://go.drugbank.com/releases/latest) <!-- Please indicate that one would need a license to access these files. Also, include versions of the databases used where appropriate (DrugBank and MetXBioDb), and then the user will decide what would be best to do (you focus on reproducibility). -->
  * drugbank_drug_structures.sdf
  * drugbank_metabolite_structures.sdf
  * drugbank_external_structures.csv
  * drugbank_full_database.xml
* [MetXBioDB](https://zenodo.org/records/13235312)
  * metxbiodb.csv
* [GLORYx Test Dataset](https://github.com/christinadebruynkops/GLORYx/tree/master/datasets/test_dataset)
  * gloryx_test_dataset.json

To extract the data, run the following files in this repository for each corresponding data source:

* `load_drugbank.py`
* `load_metxbiodb.py`
* `load_gloryx.py`


## Preprocess Data

To preprocess the data, run the file `preprocess_data.py`.

To combine and preprocess the MetaTrans datasets (check section _Pre-training_ for more information how to download MetaTrans), set the name parameter under *if __name__ == "__main__"* to *'combined'. 
The GLORYx dataset does not have to be further preprocessed. 



## Running Chemformer

Here we need to explain how to fine-tune a model <!-- Some work in the section required. -->

```bash
run_fine_tuning.sh
run_inference_score.sh
```


## Pre-training

To extract the data needed for additional pre-training, download the raw files from the following webpages:

* [Virtual Analogues Dataset](https://zenodo.org/records/45807)
  * 1297204_Virtual__Analogs.dat
* [ChEMBL](https://chembl.gitbook.io/chembl-interface-documentation/downloads)
  * chembl_35_chemreps.txt
* [MetaTrans Dataset](https://github.com/KavrakiLab/MetaTrans/tree/master/datasets)
  * source_train.txt
  * target_train.txt
  * source_train.txt
  * target_train.txt

To extract the data from the files, run the following files:

* `preprocessing/load_mmp.py` for the Virtual Analogues Dataset
* `preprocessing/load_metatrans.py` for the MetaTrans Dataset

<!-- VA abbreviation should be introduced at the place of first mention. -->
### Virtual Analogues (VA) Dataset 

Due to the large size of the dataset, it was loaded in parts. For complete curation and filtering, the following steps should be conducted.

#### Loading
* Go to file `preprocessing/load_mmp_parts.py`
* Uncomment a *start_row* and corresponding *end_row*
* Run the sh script `preprocessing/submit_load_mmp.sh`

#### Initial Filtering

The intital filtering is done separately on each split of the extracted dataset. 

* Go to file `preprocessing/NEW_preprocessing_all.py`
* Under *if __name__ == "__main__"*, set name to *'mmp_filter_part'* and uncomment a *start_row* and corresponding *end_row* (an example is provided below):
```bash
if __name__ == "__main__":

    ''' rows for mmp '''
    start_row = 0
    end_row = 1101304

    # start_row = 1101304
    # end_row = 2202607

...

    name = 'mmmp_filter_part'
```
* Run file `preprocessing/submit_preprocess.sh`
* Repeat the steps for all splits
* Run the file `preprocessing/combine_mmp.py` to combine the different VA files into one

#### Final Filtering

* Go back to the file `preprocessing/NEW_preprocessing_all.py`
* Set name to *'mmp_last_filtering'*
* Run the sh script `preprocessing/submit_load_mmp.sh`

This outputs the file `mmp_finetune.csv` that can be used for additionally pre-training training the base Chemformer model.

### MetaTrans Dataset

For filtering of the MetaTrans dataset, proceed with the following steps.

* Go to file `preprocessing/NEW_preprocessing_all.py`
* Set name to *'metatrans'
* Run the script (no sh script needed)

## Augmentation and Annotation

To augment the training dataset with parent-parent or parent-granchild reactions, or to annotate the dataset with either Csp3 fraction or LogP, follow these steps:

* Go to the file `optimisation/augmentation_annotation.py`
* Set the relevant augmentation or augmentataion techniques to *True*
  * Note that both annotation techniques can be done together, and that both augmentation techniques (parent-parent and parent-granchild reactions) can be done together. However, augmentation and annotation are conducted separately
* Run the script

Note: Specific data files are required for the augmentation techniques. Please ensure these files are available before proceeding with the augmentation process. They are generated by running `preprocessing/NEW_preprocessing_all.py` with the name parameter set to *'combined'* , which should have been done in the preprocessing step already. <!-- The last part of the sentence seems redundant. Alsom NEW_preprocessing_all.py does not sound very intuitive. -->

## Splitting Data for Ensemble Model

The data can be split at random, or split based the Burtina Clustering algorithm (parent or child similarity). 

To split based on similarity (pre-requisite for child-/parent-based splitting) the following preparation steps should be conducted.

* Go to the file `optimisation/tb_cluster.py`
* Set name to *'child'* or *'parent'*, depending on if you want to cluster based on the metabolites or drugs, respectively
* Run the script

For all different splitting methods,

* Go to the file `optimisation/split_for_ensemble.py`
* Set *split_type* to the splitting method of choice, i.e., *'random'*, *'parents'* or *'children'*

## Evaluation

To evaluate the result files obtained from the inference score, the following file can be run:

```bash
prediction_reader.py
```

## Contributing

Miranda Carlsson

Sofia Larsson

Filip Miljković

Richard Beckmann

Rocío Mercado Oropeza

Previous Master's Thesis

Chemformer <!-- I think you can place Chemformer and previous Master's Thesis in Acknowledgements. Also, at the end, once we have a paper out we can provide a citation as a bibtex. -->


