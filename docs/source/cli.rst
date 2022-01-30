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


Weight functions
-----------------
**beta** : this option uses two parameters α and β, to create beta distribution.
Depending on the parameters chosen, the distribution can change its shape, giving more flexibilty for the user to chose how to weight the variables.
The default for this function is [1,25] which are the same parameters used in SKAT-O.

.. image::  https://upload.wikimedia.org/wikipedia/commons/thumb/f/f3/Beta_distribution_pdf.svg/1920px-Beta_distribution_pdf.svg.png
    :alt: Beta distribution
image source: https://en.wikipedia.org/wiki/Beta_distribution

**log10** : this option uses -log distribution to upweight rare variants.

.. image::  https://ljvmiranda921.github.io/assets/png/cs231n-ann/neg_log.png
    :alt: -log distribution
image source: https://ljvmiranda921.github.io/notebook/2017/08/13/softmax-and-the-negative-log-likelihood/


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
The types of models available for training can be found :ref:`_model_types`

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

