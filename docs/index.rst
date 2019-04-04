.. baghera_tool documentation master file, created by Viola Fanfani
   sphinx-quickstart
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.


Welcome to Bayesian Gene Heritability Analysis documentation!
=====================================================================
By Viola Fanfani (v.fanfani@sms.ed.ac.uk)

StracquadanioLab, University of Edinburgh

BAGHERA is a software able to estimate the single-gene contribution to the genome-wide heritability. It only requires
the summary statistics of a GWAS, ld-scores and a set of annotated variants.

To help the user BAGHERA comes with the ld-scores for the European population of 1000 Genomes
and a GENCODE file already available inside the `data/` folder.

You  can check here to see how the package is evolving

.. toctree::
   :maxdepth: 2

   changelog.rst


File parsing
------------

Before the real analysis the user needs to create an input dataset suitable for the analysis tools.

.. toctree::
   :maxdepth: 2
   :caption: File Parsing:

   createdata

Principal Features
------------------

The software allows two different analyses: the Bayesian estimation of genome-wide heritability and the gene-level heritability analysis.

.. toctree::
   :maxdepth: 2
   :caption: Principal Features:

   analysis

Examples
++++++++

Example pipeline::

   $ baghera-tool create-files <output-folder>
   $ baghera-tool generate-SNPs-file <UKBB_summary_filename> --input_type ukbb --output-file <output-filename>.csv
   $ baghera-tool regression <output-filename>.csv <SUFFIX> <output-folder>

This pipeline uses the default options and files of baghera-tool. Look at the documentation for each function
to change the parameters and files that are used.

Help and logs
-------------
.. toctree::
   :maxdepth: 1
   :caption: Contents:

   logging


Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
