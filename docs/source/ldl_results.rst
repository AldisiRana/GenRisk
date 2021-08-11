Real use case
###############
We ran the pipeline using about 200K samples from UKBiobank.
We filtered for British to remove outliers and calculated the scores using 1% MAF threshold, beta 1-25 weighting parameter and CADD raw scores.

Association analysis
*********************
We ran association analysis on the samples with LDL direct measurements as phenotype(quantitative).
We used linear regression for the analysis and age, sex, BMI and PC1-4 were used as covaraites.
The association analysis highlighted PCSK9 and LDLR as significant genes, both are known to be associated with LDL.
QQ-plot and Manhattan plot are presented below.

The QQ-plot:

.. image:: ../../ldl_example/linear_assoc_ldl_200k_qqplot.png
    :width: 400
    :align: center

The Manhattan plot:

.. image:: ../../ldl_example/linear_assoc_ldl_200k_manhattan.png
    :width: 400
    :align: center

Regression model
*****************
Feature importance plot:

.. image:: ../../ldl_example/Feature_Importance.png
    :width: 400
    :align: center

Actual vs Predicted:

.. image:: ../../ldl_example/ldl_testing_regressor_scatterplot.png
    :width: 400
    :align: center

Model residuals:

.. image:: ../../ldl_example/Residuals.png
    :width: 400
    :align: center

