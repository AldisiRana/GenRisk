# Toy Example

This folder contains an example on how to run the CoGenAssess pipeline.
The toy dataset is not real data. It contains about 10,000 SNP and 20 genes, 
2 of which (gene5 and gene17) were synthetically made to produce signals for the included phenotypes.
The dataset the signals and phenotypes were simulated using [SIMU](https://github.com/precimed/simu).

## Files:
### Annotated vcf 
The annotated VCF can be downloaded from [here](https://uni-bonn.sciebo.de/s/XUgOnCKxDhV60YI).
It contains information about the SNPS, gene, deleteriousness score and allele frequency.
This information is important for the calculation of the gene-based scores.

### bim/bam/fam files
The bim, bam and fam files are important for the PLINK processing:
- bed: text file for genotype
- bim: file with variant information
- fam: text file for pedigree information

### snps_list
This file contains the toy significant variants that were used to create the toy phenotypes.

