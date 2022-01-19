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

Total run time in seconds

| | 1K samples | 10K samples | 100K samples |
| ------------- | ------------- | :-------------: | -----:|
| 10 feats | 14 | 19 | 1690 |
| 100 feats | 24 | 678 | 41649 |
| 1000 feats | 143 | 1034 | 432000 (5 days) |

