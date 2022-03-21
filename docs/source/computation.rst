.. _computation-info:

Computation information
########################

It should be noted that the aim of our work is to provide a novel framework more comprehensive in terms of genetic risk assessment and currently is not yet optimized for computational performance.
For all the computation below, we use a “standard” workstation (RAM=64GB with 6 CPU dual core).

Gene-based scoring and analysis
********************************
In the following table we show the computation time for the gene-core computation and association analysis with linear regression of the biggest chromosome (1,727,756 variants, MAF filtering = <1%) which includes also the higher number of genes (1,972 genes) by considering different numbers of individuals.

+----------------------------------------------+------------+-------------+-----------------+
|                                              | 1K samples | 10K samples | 100K samples    |
+==============================================+============+=============+=================+
| gene-scoring (in mins)                       |     22     |     25      |      48         |
+----------------------------------------------+------------+-------------+-----------------+
| Find-association, linear regression (in sec) |      8     |     28      |       134       |
+----------------------------------------------+------------+-------------+-----------------+

While for prediction models the complete input matrix (i.e., samples and genes plus covariates) should be loaded in RAM, for gene-scoring we use the efficient score function implemented in PLINK v2 ().
For gene-association GenRisk the memory usage depends on the size of the input matrix, the larger the matrix the more memory it uses.

+------------+------------+-------------+-----------------+
|            | 1K samples | 10K samples | 100K samples    |
+============+============+=============+=================+
| Mem (in Gb)|    3.1     |     3.4     |      9.4        |
+------------+------------+-------------+-----------------+


Prediction models generation
******************************
GenRisk has “per-se” no limit in the number of features that can be used. However, there could be computational issues
according to the dimensionality of the input, that is samples and features (genes, covariates, etc..).
The tables below present the total run time (in seconds) and maximum memory usage (in GB) given different sample sizes
with increasing number of features. Please note that
it might be wise to run big data size (e.g 100K x 1000feats) using an HPC infrastructure.
Another point to consider is that the time and memory usage also depends on the models included in the analysis and the
best model fine-tuning and finalization. Some models, such as gradient boosting, might take more time than simpler models,
like linear or lasso regression, to be finalized.

Total run time of prediction model generation in seconds
---------------------------------------------------------

+------------+------------+-------------+-----------------+
|            | 1K samples | 10K samples | 100K samples    |
+============+============+=============+=================+
| 10 feats   |     14     |     19      |    1690         |
+------------+------------+-------------+-----------------+
| 100 feats  |      24    |     678     |       41649     |
+------------+------------+-------------+-----------------+
| 1000 feats |     143    |    1034     | 432000(≈ 5days) |
+------------+------------+-------------+-----------------+

Maximum memory used in GB
--------------------------

+------------+------------+-------------+-----------------+
|            | 1K samples | 10K samples | 100K samples    |
+============+============+=============+=================+
| 10 feats   |    2.81    |     2.93    |      2.97       |
+------------+------------+-------------+-----------------+
| 100 feats  |     2.93   |     2.95    |       3.29      |
+------------+------------+-------------+-----------------+
| 1000 feats |     3.51   |    3.82     |       8.29      |
+------------+------------+-------------+-----------------+

Feature Selection
------------------
In general in the context of prediction models for big datasets we would suggest a feature selection using the “association” module and then generate prediction models.
This is also in line with the expected genetic architecture of the majority of the traits in which only a small proportion of the genes plays a pivotal role.
If instead we have a really highly polygenic phenotype the computation of genome-wide polygenic risk score (PRS) is probably the most appropriate approach, as PRS is a value per individual only a vector of scores would be generated and therefore the computational burden is limited.


