.. _methods-comparison:
Methods Comparisons
#####################
We have compared GenRisk scores association test results with different burden test methods.
In this example, we use LDL measurements (adjusted for statin) as phenotype and ≃160,000 samples from UKbiobank.

We compared the results with SKATO, using the same input as GenRisk for LDL phenotype,
and found that SKATO has detected many genes (58 genes after p-value adjustment) as significant,
including PCSK9 and LDLR. However, the lambda of the p-values is inflated (1.325),
as opposed to GenRisk which had a lambda of 1.056, which means there is a risk of having false-positives.

We, also performed two burden tests from rvtest (https://github.com/zhanxw/rvtests), CMC and
Zeggini, and no significant genes were detected with CMC analysis, while Zeggini was able to detect PCSK9 but not LDLR.

Burden tests usually use genotypes only to score the genes, and can sometimes use filters like functional annotations and allele frequency, but this is just filtering (no values) and it has to be done as a pre-step before the actual analysis.
An example of that is Genebass, where they applied SKATO on three different sets (Loss of function, missense and synonymous), and the results for LDL presented PCSK9 and LDLR on along with other genes, however the lambda for the p-values is inflated (e.g the lambda for LDL direct across all burden sets for the SKATO is 1.18 and for Burden test is 1.16), which means there is a risk of having false positives.
https://genebass.org/gene/undefined/phenotype/continuous-30780-both_sexes--irnt?burdenSet=pLoF&resultIndex=gene-manhattan&resultLayout=full
Another example is astrazeneca phewas portal (azphewas), here they have multiple models that filter for frequency, categorical annotations and/or deleteriousness score threshold.
https://azphewas.com/phenotypeView/f87604bb-7293-44e8-8e29-bf58d9872841/4b20a1ff-bded-4f1e-8301-f2922f0b8499/glr



