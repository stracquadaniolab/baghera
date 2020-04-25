.. baghera_tool documentation master file, created by Viola Fanfani
   sphinx-quickstart
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Bayesian Gene Heritability Analysis
===================================
By Viola Fanfani (v.fanfani@sms.ed.ac.uk)

The Bayesian Gene Heritability Analysis software (BAGHERA) estimates the contribution
to the heritability of a trait/disease of all the SNPs in the genome (genome-wide heritability)
and those nearby protein-coding genes (gene-level heritability).

BAGHERA requires only summary statistics from a Genome-wide Association Study (GWAS),
LD scores calculated from a population matching the ethnicity of the GWAS study and
a gene annotation file in GTF format.

Workflow
---------

Alongside BAGHERA, we are providing a snakemake workflow repository with sample data.
`workflow-baghera <https://github.com/stracquadaniolab/workflow-baghera>`_

Getting Started
------------------

.. toctree::
   :maxdepth: 2

   gettingstarted.rst

Usage
---------
.. toctree::
   :maxdepth: 2

   createdata.rst
   analysis.rst

API
-----

.. toctree::
   :maxdepth: 2

   api.rst

Issues
------
BAGHERA is still under development, and a major refactoring will happen at some point.
If you find a bug, please report it on GitHub.

Changelog
---------------


.. toctree::
   :maxdepth: 2

   changelog.rst



Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
