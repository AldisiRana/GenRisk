.. _cli:
Command Line Interface
#######################
The genrisk command line interface includes multiple commands which can be used as follows::

    $ genrisk COMMAND [ARGS]

score-genes
***********
::

    $ genrisk score-genes [ARGS]

.. autofunction:: genrisk.cli.score_genes

Gene-scoring equation
======================
The gene scores are derived by the weighted sum of the variants in a gene.

.. math::
    G_{sg}= \sum_{\it i=1}^{\it k} (D_i \times A_i) C_i

D\ :sub:`i` is the functional annotation (e.g CADD)

A\ :sub:`i` is the weighted allele frequency

C\ :sub:`i` is the allele count.

Weight functions
-----------------
The weighting function is applied to the variant frequency. I can be used to up-weight the biological importance of rare variants.

**beta** : this option uses two parameters α and β, to create beta distribution.
Depending on the parameters chosen, the distribution can change its shape, giving more flexibilty for the user to chose how to weight the variables.
The default for this function is [1,25] which are the same parameters used in SKAT-O.

.. image::  https://upload.wikimedia.org/wikipedia/commons/thumb/f/f3/Beta_distribution_pdf.svg/1920px-Beta_distribution_pdf.svg.png
    :width: 300
    :alt: Beta distribution
`image source here <https://en.wikipedia.org/wiki/Beta_distribution>`_

**log10** : this option uses -log distribution to upweight rare variants. This has been applied previously in another
`gene-based score tool <https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-019-2877-3>`_

.. image::  https://ljvmiranda921.github.io/assets/png/cs231n-ann/neg_log.png
    :width: 300
    :alt: -log distribution
`image source here <https://ljvmiranda921.github.io/notebook/2017/08/13/softmax-and-the-negative-log-likelihood/>`_


normalize
**********
::

    $ genrisk normalize [ARGS]

.. autofunction:: genrisk.cli.normalize


find-association
*****************
::

    $ genrisk find-association [ARGS]

.. autofunction:: genrisk.cli.find_association

visualize
**********
::

    $ genrisk visualize [ARGS]

.. autofunction:: genrisk.cli.visualize


create-model
*************
::

    $ genrisk create-model [ARGS]

.. autofunction:: genrisk.cli.create_model

Available models
=================
The types of models available for training can be found :ref:`model_types`

test-model
***********
::

    $ genrisk test-model [ARGS]

.. autofunction:: genrisk.cli.test_model


get-prs
********
::

    $ genrisk get-prs

.. autofunction:: genrisk.cli.get_prs

