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

    $ cogenassess score_genes --annotated-vcf toy_example/annotated_toy_dataset.vcf --bed toy_example/toy_dataset.bed 
    --bim toy_example/toy_dataset.bim --fam toy_example/toy_dataset.fam --temp-dir toy_example/toy_dataset/ 
    --output-file toy_example/toy_dataset_scores --weight-func beta --remove-dir
```
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
    
    $ cogenassess calculate-pval --scores-file toy_example/toy_dataset_scores --gentype-file toy_example/toy.pheno 
    --cases-column trait1 --samples-column IID --test betareg --output-path toy_dataset_betareg.tsv
```
required arguments:
  -s, --scores-file         a file containing variant IDs, genes, alt, deleteriousness scores and allele frequency.
  -i, --genotype-file       text file for genotype.
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

    $cogenassess visualize --pvals-file toy_example/toy_dataset_scores --info-file annotated_toy_dataset.vcf
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

### Prediction model
TODO
