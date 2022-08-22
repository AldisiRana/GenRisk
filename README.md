# GenRisk

GenRisk is a package that implements different gene-based scoring schemes to analyze and find significant genes 
within a phenotype in a population

## Citation
Rana Aldisi, Emadeldin Hassanin, Sugirthan Sivalingam, Andreas Buness, Hannah Klinkhammer, Andreas Mayr, Holger Fröhlich, Peter Krawitz, Carlo Maj, GenRisk: a tool for comprehensive genetic risk modeling, Bioinformatics, Volume 38, Issue 9, 1 May 2022, Pages 2651–2653, https://doi.org/10.1093/bioinformatics/btac152

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

    $ genrisk score-genes -a ../path/to/toy_vcf_data.vcf -o toy_genes_scores.tsv -t toy_vcf_scoring -v ID -f AF -g gene -l ALT -d RawScore

* For further CLI options and parameters use --help

### Calculate p-values
This function calculates the p-values across the genes between two given groups
    
    $ genrisk find-association -s toy_genes_scores.tsv -i info.pheno -o linear_assoc_quan.tsv -t linear -c quan -a fdr_bh -v sex,age,bmi 

* For further CLI options and parameters use --help

### Visualize
Visualize manhatten plot and qqplot for the data.

    $ genrisk visualize -p logit_assoc_binary.tsv -i genes_info_ref.txt -q logit_assoc_binary_qqplot.png -m logit_assoc_binary_manhattan.png --genescol-1 genes

* For further CLI options and parameters use --help

### Create model
Create a prediction model (classifier or regressor) with given dataset

    $ genrisk create-model -d toy_dataset_feats.tsv -o quan_regression_model -n quan_regression_model --model-type regressor -l quan --normalize

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