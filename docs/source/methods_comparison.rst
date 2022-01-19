.. _methods_comparison:
Methods Comparisons
#####################
We have compared GenRisk scores association test results with different burden test methods.
In this example, we use LDL measurements (adjusted for statin) as phenotype and ≃160,000 samples from UKbiobank.

We compared the results with SKATO, using the same input as GenRisk for LDL phenotype,
and found that SKATO has detected many genes (58 genes after p-value adjustment) as significant,
including PCSK9 and LDLR. However, the lambda of the p-values is inflated (1.325),
as opposed to GenRisk which had a lambda of 1.056, which means that some of the detected significant genes could be false-positive.

We, also performed two burden tests from rvtest (https://github.com/zhanxw/rvtests), CMC and
Zeggini, and no significant genes were detected with CMC analysis, while Zeggini was able to detect PCSK9 but not LDLR.


