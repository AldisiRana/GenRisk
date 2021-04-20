# CoGenAssess

CoGenAssess is a package that implements different gene-based scoring schemes to analyze and find significant genes 
within a phenotype in a population

## Requirements
* annotated vcf file (should contain variant ID, ALT, Gene, and deleterious score, for more information check out the example in toy_example)
* plink 1.9 https://www.cog-genomics.org/plink/


## Installation
``CoGenAssess`` can be installed on python3+ from the latest code on [GitHub](https://github.com/AldisiRana/CoGenAssess) with:

    $ pip install git+https://github.com/AldisiRana/CoGenAssess.git

## Usage

### Score genes
Calculate the gene-based scores for a given dataset.

``` sh
$ cogenassess score_genes --annotated-vcf toy_example/annotated_toy_dataset.vcf --bed toy_example/toy_dataset.bed 
--bim toy_example/toy_dataset.bim --fam toy_example/toy_dataset.fam --temp-dir toy_example/toy_dataset/ 
--output-file toy_example/toy_dataset_scores --weight-func beta --remove-dir

required arguments:
  -a, --annotated-vcf   a file containing variant IDs, genes, alt, deleteriousness scores and allele frequency.
  --bed                 text file for genotype.
  --bim                 file with variant information
  --fam                 text file for pedigree information
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
``` sh
$ cogenassess calculate-pval --scores-file path/to/scores/file --gentype-file file/containing/the/cases 
    --cases-column the-name-of-cases-column --samples-column name-of-samples-column 
    --test [ttest-ind|mannwhitneyu|logit|glm|betareg] --output-path path/of/output/file --covariates PC1,PC2
```

### Visualize

### Prediction model

### Normalize scores
This function normalizes the scores by gene length

    $ cogenassess normalize --matrix-file path/to/score/file --output-path path/to/output/file 
    --samples-col name-of-samples-column


