# GenRisk

GenRisk is a package that implements different gene-based scoring schemes to analyze and find significant genes 
within a phenotype in a population

## Requirements
* annotated vcf (should contain variant ID, ALT, Gene, and deleterious score, for more information check out the example in toy_example)
* plink >= 1.9 https://www.cog-genomics.org/plink/
* R version >= 3.6.3


## Installation
``GenRisk`` can be installed on python3+ from the latest code on [GitHub](https://github.com/AldisiRana/GenRisk) with:

    $ pip install git+https://github.com/AldisiRana/GenRisk.git

## Usage

### Score genes
Calculate the gene-based scores for a given dataset.

    $ genrisk score_genes --annotated-vcf toy_example/annotated_toy_dataset.vcf --temp-dir toy_example/toy_dataset/ 
    --output-file toy_example/toy_dataset_scores --weight-func beta --remove-dir
```
required arguments:
  -a, --annotated-vcf   a file containing variant IDs, genes, alt, deleteriousness scores, allele frequency and samples info.
  -t, --temp-dir        a temporary directory to save temporary files before merging.
  -o, --output-file     the final output scores matrix.

optional arguments:
  -h, --help         show this help message and exit
  --plink            the directory of plink, if not set in environment
  --weight-func      the weighting function used in score calculation.
                     Choices: [beta, log10]. Default=beta.
  --beta-param       the parameters from beta weight function. Default=(1,25)
  --variant-col      the column containing the variant IDs. Default=SNP
  --gene-col         the column containing gene names. Default=Gene.refGene
  --af-col           the column containing allele frequency. Default=MAF.
  --del-col          the column containing the deleteriousness score. Default=CADD_raw
  --alt-col          the column containing the alternate base. Default=Alt
  --maf-threshold    the threshold for minor allele frequency. Default=0.01
  --remove-temp      if flagged the temporary directory will be deleted after process completion.
```

### Calculate p-values
This function calculates the p-values across the genes between two given groups
    
    $ genrisk find-association --scores-file toy_example/toy_dataset_scores --info-file toy_example/toy.pheno 
    --cases-column trait1 --samples-column IID --test betareg --output-path toy_dataset_betareg.tsv
```
required arguments:
  -s, --scores-file         a file containing variant IDs, genes, alt, deleteriousness scores and allele frequency.
  -i, --info-file       text file for genotype.
  -o, --output-path         file with variant information
  -t, --test                statistical test for calculating P value. Choices: ttest_ind, mannwhitneyu, logit, glm, betareg
  -c, --cases-column        the name of the column that contains the case/control type.
  -m, --samples-column      the name of the column that contains the samples.
  
optional arguments:
  -h, --help                show this help message and exit
  -g, --genes               a list containing the genes to calculate. if not provided all genes will be used.
  -p, --pc-file             Principle components values for logistic regression. if not in genotype file
  --adj-pval                adjust pvalue using one of these methods: bonferroni, sidak, holm-sidak, holm, simes-hochberg, hommel, fdr_bh, fdr_by, fdr_tsbh, fdr_tsbky
  --covariates              the covariates used for calculation. Default=PC1,PC2
```

### Visualize
Visualize manhatten plot and qqplot for the data.

    $ genrisk visualize --pvals-file toy_example/toy_dataset_scores --info-file annotated_toy_dataset.vcf
    --qq-output toy_example/toy_dataset_qqplot.jpg --manhattan-output toy_example/toy_dataset_manhattanplot.jpg 

```
required arguments:
  --pvals-file              the file containing p-values.
  --info-file               file containing variant/gene info.
  --qq-output               the name of the qq plot file. should end with jpg.
  --manhattan-output        the name of the manhatten plot file. should end with jpg.

optional arguments:
  -h, --help                show this help message and exit
  --genescol-1              the name of the genes column in pvals file. Default=gene.
  --genescol-2              the name of the genes column in info file. default=Gene.refGene.
  --pvalcol                 the name of the pvalues column. Default=p_value.
```

### Create model
Create a machine learning model (classifier or regressor) with given dataset

    $ genrisk create-model --data-file toy_example_regressor_features.tsv --model-type regressor --output-folder toy_regressor 
    --test-size 0.25 --test --model-name toy_regressor --target-col trait1 --imbalanced --normalize

```
required arguments:
  --data-file               file containing features and target.
  --output-folder           a folder path to save all outputs.
  --model-name              the name of the model to be saved.
  --model-type              the type of model (regressor or classifier).
  --target-col              the name of the target column in data file.

optional arguments:
  -h, --help                show this help message and exit
  --test                    if flagged, a test set will be created for evaluating the final model.
  --test-size               test size for cross validation and evaluation. Default=0.25
  --imbalanced              if flagged methods will be used to account for the imbalance.
  --normalize               if flagged the data will be normalized before training.
  --folds                   number of cross-validation folds in training. Default=10
  --metric                  the metric used to choose best model after training. Regressor default=RMSE, classifier default=AUC
  
```

### Test model
Evaluate a machine learning model with a given dataset.

    $ genrisk test-model --model-path regressor_model.pkl --input-file testing_dataset.tsv --model-type regressor 
    --labels-col target --samples-col IID

```
required arguments:
  -t, --model-type           type of prediction model.
  -i, --input-file           the testing dataset.
  -l, --label-col            the target/phenotype/label column
  -m, --model-path           the path to the ML model.

optional arguments:
  -h, --help                 show this help message and exit
  -s, --samples-col          the samples column.
  
```


### Get PRS scores
Gets a PGS file (provided by the user or downloaded from pgscatalog) then calculates the PRS scores for dataset.
Note: This command is interactive.

    $ genrisk get-prs

```

optional arguments:
  -h, --help                 show this help message and exit
  --plink                    provide plink path if not default in environment.
```
