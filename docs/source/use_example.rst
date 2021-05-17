Example use case
================
This toy dataset is not real data. It contains about 10,000 SNP and 20 genes,
2 of which (gene5 and gene17) were synthetically made to produce signals for the included phenotypes.
The signals and phenotypes were simulated using `SIMU <https://github.com/precimed/simu>`_.
It was created for testing and evaluating the pipeline only.
In this example use case, we use the toy dataset to generate gene-based scores, association analysis and machine learning models.

Annotated VCF
**************
(`Click here <https://uni-bonn.sciebo.de/s/WQroVFBQ8NXNnF1>`_ to download VCF)

The annotated vcf contains information about the SNPS, gene, deleteriousness score and allele frequency.
It should also contain samples genotypes. This information is important for the calculation of the gene-based scores.

Gene-based scores
*****************
The gene-based scores can be found `here <https://github.com/AldisiRana/CoGenAssess/blob/master/toy_example/toy_dataset_scores>`_.
These scores are used as input for the association analysis and as features for the machine learning models.

Association analysis
********************
The association analysis (`found here <https://github.com/AldisiRana/CoGenAssess/blob/master/toy_example/toy_dataset_betareg_pvals.tsv>`_)
done on the dataset showed significance in gene5 and geen17 as expected.

The QQ-plot was also produced using the pipeline:

.. image:: images/toy_dataset_qqplot.jpeg

Manhatten plot was not generated because the dataset is small and doesn't contain real data.

Machine learning models
***********************

TBA


