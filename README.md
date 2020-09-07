# CoGenAssess

CoGenAssess is a package that implements different gene-based scoring schemes to analyze and find significant genes 
within a phenotype in a population

## Requirements
* annotated vcf file (check vcfanno)
* plink 1.9 https://www.cog-genomics.org/plink/


## Installation
``CoGenAssess`` can be installed on python3+ from the latest code on [GitHub](https://github.com/AldisiRana/CoGenAssess) with:

    $ pip install git+https://github.com/AldisiRana/CoGenAssess.git

## Usage

### Get scores
This function calculates the gene-based scores of a gene list for a group of samples.

    $ cogenassess get-scores --vcf vcf/file/path --bed bed/file/path --bim bim/file/path --fam fam/file/path --temp-dir /dir/for/temp/files --ourput-file /final/output/file 
    --weight-func [log10|beta] 


### Normalize scores
This function normalizes the scores by gene length

    $ cogenassess normalize --matrix-file path/to/score/file --output-path path/to/output/file --samples-col name-of-samples-column


### Calculate p-values
This function calculates the p-values across the genes between two given groups

    $ cogenassess calculate-pval --scores-file path/to/scores/file --gentype-file file/containing/the/cases --cases-column the-name-of-cases-column
    --samples-column name-of-samples-column --test [ttest-ind|mannwhitneyu|logit|glm] --output-path path/of/output/file

