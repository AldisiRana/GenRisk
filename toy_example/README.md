# Toy Example

This folder contains an example on how to run the GenRisk pipeline.
The toy dataset is not real data. It contains real SNPs and gene info, but no real individuals.

## Files:

### Annotations file
(filename: toy_annotations.tsv)

The annotated file contains information about the SNPS, gene, deleteriousness score and allele frequency.
It should also contain samples genotypes. This information is important for the calculation of the gene-based scores.

### Plink binary files
(filenames: toy_data.bed, toy_data.bim, toy_data.fam)

The binary files contains all genotype information for the cohort.

**The annotations and binary files [here](https://uni-bonn.sciebo.de/s/Ih1Hhrf5ek8u8sa).**

### Phenotype file
(filename: info.pheno)

This file contains covariates (columns: bmi, age and sex), and two traits, one quantitative (column: quan) and one binary (column: binary).
This file can be used for running the association analysis.

### Genes info file
(filename: genes_info_ref.txt)

This file contains the genes names and locations in the genome. It can be used for the visualization pipeline.

### Features file
(filename: toy_dataset_feats.tsv)

This file contains a subset of the genes and the columns from the phenotype file. This can be used as an input file to generate prediction models either regressor or classifier.

### results files
We have run the pipeline on the dataset mentioned above and generated multiple results:
- the scoring pipeline was performed as default on the annotated vcf.
- the association pipeline was performed on the quantitative trait (linear and betareg) and binary trait (logit).
- the qqplots and manhattan plots were generated for the association analysis results using the gene info file.
- two models were generated, a regression model using the quantitative trait as target, and a classification model using the binary trait.
