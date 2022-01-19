.. _computation_info:

Computation information
########################


Prediction models generation
******************************
GenRisk has “per-se” no limit in the number of features that can be used. However, there could be computational issues
according to the dimensionality of the input, that is samples and features (genes, covariates, etc..).
The tables below present the total run time (in seconds) and maximum memory usage (in GB) given different sample sizes
with increasing number of features obtained on a “standard” workstation (RAM=64GB with 6 CPU dual core). Please note that
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
.. list-table:: Maximum memory used of prediction model generation in GB
   :widths: 25 25 50
   :header-rows: 1

* -
     - 1K samples
     - 10K samples
     - 100K samples
   * - 10 feats
     - 2.81
     - 2.93
     - 2.97
   * - 100 feats
     - 2.93
     - 2.95
     - 3.29
   * - 1000 feats
     - 3.51
     - 3.82
     - 8.29

| | 1K samples | 10K samples | 100K samples |
| ------------- | ------------- | :-------------: | -----:|
| 10 feats | 2.81 | 2.93 | 2.97 |
| 100 feats | 2.93 | 2.95 | 3.29 |
| 1000 feats | 3.51 | 3.82 | 8.29 |


In general in the context of prediction models for big datasets we would suggest a feature selection using the “association” module and then generate prediction models.
This is also in line with the expected genetic architecture of the majority of the traits in which only a small proportion of the genes plays a pivotal role.
If instead we have a really highly polygenic phenotype the computation of genome-wide polygenic risk score (PRS) is probably the most appropriate approach, as PRS is a value per individual only a vector of scores would be generated and therefore the computational burden is limited.


