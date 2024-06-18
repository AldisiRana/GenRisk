
.. genrisk documentation master file, created by
    sphinx-quickstart on Wed May 12 10:33:11 2021.
    You can adapt this file completely to your liking, but it should at least
    contain the root `toctree` directive.

Welcome to GenRisk's documentation!
=======================================

GenRisk is a package that implements different gene-based scoring schemes to analyze and find significant genes
within a phenotype in a population

.. toctree::
    :maxdepth: 1
    :hidden:
    :caption: Content

   cli
   pipeline
   utils
   model_types
   use_example
   real_cases
   methods_comparison
   computation


Installation
----------------

Requirements
~~~~~~~~~~~~~~~
* `plink <https://www.cog-genomics.org/plink>`_ >= 1.9
* R version >= 3.6.3
* python >= 3.

Package installation
~~~~~~~~~~~~~~~~~~~~~~~~
Option 1: The latest release of ``GenRisk`` can be installed on python3+ with:

.. code-block:: bash

    pip install genrisk

Option2: you can also install the package with the latest updates directly from `GitHub <https://github.com/AldisiRana/GenRisk>`_ with:

.. code-block:: bash

    pip install git+https://github.com/AldisiRana/GenRisk.git


Citation
---------
Rana Aldisi, Emadeldin Hassanin, Sugirthan Sivalingam, Andreas Buness, Hannah Klinkhammer, Andreas Mayr, Holger Fröhlich, Peter Krawitz, Carlo Maj, GenRisk: a tool for comprehensive genetic risk modeling, Bioinformatics, Volume 38, Issue 9, 1 May 2022, Pages 2651–2653, https://doi.org/10.1093/bioinformatics/btac152


###############################


Indices and tables
==================
* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
