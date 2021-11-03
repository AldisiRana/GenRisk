# GenRisk

GenRisk is a package that implements different gene-based scoring schemes to analyze and find significant genes 
within a phenotype in a population

## Requirements
* plink >= 1.9 https://www.cog-genomics.org/plink/
* R version >= 3.6.3

## Installation
Option 1: The latest release of ``GenRisk`` can be installed on python3+ with:

    $ pip install genrisk

Option2: you can also install the package with the latest updates directly from `GitHub <https://github.com/AldisiRana/GenRisk>`_ with:

    $ pip install git+https://github.com/AldisiRana/GenRisk.git

## Usage

### Score genes
This command calculate the gene-based scores for a given dataset.

It requires an annotated vcf (i.e: annotated with variant ID , ALT, Gene, and deleterious score, for more information check out the example in toy_example)

    $ genrisk score-genes --annotated-vcf annotated_vcf_toy.vcf --temp-dir test/ --output-file test.tsv --weight-func beta --maf-threshold 0.01 --alt-col ALT --variant-col ID --af-col AF --del-col CADD --gene-col Gene

* For further CLI options and parameters use --help

### Calculate p-values
This function calculates the p-values across the genes between two given groups
    
    $ genrisk find-association --scores-file toy_example/toy_dataset_scores --info-file toy_example/toy.pheno 
    --cases-column trait1 --samples-column IID --test betareg --output-file toy_dataset_betareg.tsv --covariates age,sex
    --adj-pval bonferroni
* For further CLI options and parameters use --help

### Visualize
Visualize manhatten plot and qqplot for the data.

    $ genrisk visualize --pvals-file toy_example/toy_dataset_scores --info-file annotated_toy_dataset.vcf
    --qq-output toy_example/toy_dataset_qqplot.jpg --manhattan-output toy_example/toy_dataset_manhattanplot.jpg 
* For further CLI options and parameters use --help

### Create model
Create a prediction model (classifier or regressor) with given dataset

    $ genrisk create-model --data-file toy_example_regressor_features.tsv --model-type regressor --output-folder toy_regressor 
    --test-size 0.25 --test --model-name toy_regressor --target-col trait1 --imbalanced --normalize
* For further CLI options and parameters use --help

### Test model
Evaluate a prediction model with a given dataset.

    $ genrisk test-model --model-path regressor_model.pkl --input-file testing_dataset.tsv --model-type regressor 
    --labels-col target --samples-col IID
* For further CLI options and parameters use --help

### Get PRS scores
This command aquires a PGS file (provided by the user or downloaded from pgscatalog) then calculates the PRS scores for dataset.
Note: This command is interactive.

    $ genrisk get-prs
* For further CLI options and parameters use --help

### Get GBRS
Calculate gene-based risk scores for individuals. 
If users do not have weights for calculation, they can provide a file with the phenotype and weights will be calculated.

    $genrisk get-gbrs --scores-file scores_file.tsv --weights-file weights_file.tsv --weights-col zscore --sum
* For further CLI options and parameters use --help

## Contact
If you have any questions or problems with the tool or its installation please feel free to create an issue in the repository or contact me via email:
aldisi.rana@gmail.com