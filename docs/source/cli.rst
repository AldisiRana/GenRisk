.. _cli:
Command Line Interface
#######################
The genrisk command line interface includes multiple commands which can be used as follows:

.. click:: genrisk.cli:score_genes
    :prog: genrisk score-genes

.. rubric:: Gene-scoring equation

    The gene scores are derived by the weighted sum of the variants in a gene.

    .. math::
        G_{sg}= \sum_{\it i=1}^{\it k} (D_i \times A_i) C_i

    D\ :sub:`i` is the functional annotation (e.g CADD)

    A\ :sub:`i` is the weighted allele frequency

    C\ :sub:`i` is the allele count.

.. note:: Weighting functions

    :beta: this option uses two parameters α and β, to create beta distribution. Depending on the parameters chosen, the distribution can change its shape, giving more flexibilty for the user to chose how to weight the variables.
        The default for this function is [1,25] which are the same parameters used in SKAT-O.

    .. image::  https://upload.wikimedia.org/wikipedia/commons/thumb/f/f3/Beta_distribution_pdf.svg/1920px-Beta_distribution_pdf.svg.png
        :width: 300
        :alt: Beta distribution
    `image source here <https://en.wikipedia.org/wiki/Beta_distribution>`_

    :log10: this option uses -log distribution to upweight rare variants. This has been applied previously in another `gene-based score tool <https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-019-2877-3>`_

    .. image::  https://ljvmiranda921.github.io/assets/png/cs231n-ann/neg_log.png
        :width: 300
        :alt: -log distribution
    `image source here <https://ljvmiranda921.github.io/notebook/2017/08/13/softmax-and-the-negative-log-likelihood/>`_

|
|


.. click:: genrisk.cli:normalize
    :prog: genrisk normalize

.. rubric:: Normalization methods

    Multiple methods have been implemented to normalize a dataset. Below is a brief describtion of each function.

    :gene_length: This method divides each gene-based score by the length of the gene. The genes lengths can be provided by the user, or retrieved from ensembl database. The gene length from ensembl database is calculated as such: gene length = gene end (bp) - gene start (bp)

    :minmax: This method rescales the values of each column to [0,1] by using the following formula x`= x - min(x) / max(x) - min(x)

    :maxabs: In this method, the values are normalized by the maximum absolute to [-1,1] using the following formula x` = x / max(|x|)

    :zscore: This method uses the mean and standard deviation to normalize the values. Formula is x`= x - mean(x) / std

    :robust: Great choice for dataset with many outliers. In this method, the values are substracted by the median then divided by the interquantile range (difference between the third and the first quartile). Formula x`= x - median(x) / Q3(x) - Q1(x)

    Every normalization method has it's advantages and disadvantages, so choose the method that works best with your dataset. To learn more about the normalization methods, check out this helpful `article <https://towardsdatascience.com/data-normalization-with-pandas-and-scikit-learn-7c1cc6ed6475>`_

|
|


.. click:: genrisk.cli:find_association
    :prog: genrisk find-association

.. click:: genrisk.cli:visualize
    :prog: genrisk visualize


.. click:: genrisk.cli:create_model
    :prog: genrisk create-model


.. click:: genrisk.cli:test_model
    :prog: genrisk test-model


.. click:: genrisk.cli:get_prs
    :prog: genrisk get-prs

|
|


