.. baghera_tool documentation master file, created by Viola Fanfani
   sphinx-quickstart
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Bayesian Gene Heritability Analysis
===================================
By Viola Fanfani (v.fanfani@sms.ed.ac.uk)

StracquadanioLab, School of Biological Sciences, University of Edinburgh.

BAGHERA, the Bayesian Gene Heritability Analysis, is a software to estimate the contribution
to the heritability of a trait/disease of all the SNPs in the genome (genome-wide heritability)
and those nearby protein-coding genes (gene-level heritability).

BAGHERA requires only summary statistics from a Genome-wide Association Study (GWAS),
LD scores calculated from a population matching the ethnicity of the GWAS study and
a gene annotation file in GTF format.

Installation
------------

The easiest and fastest way to install BAGHERA using conda::

$ conda install -c stracquadaniolab -c bioconda -c conda-forge


Getting started
---------------

A typical BAGHERA analysis consists of 3 steps:

1. Build a SNP annotation file, where SNPs are annotated to genes and LD scores are assigned.
2. Annotate summary statistics with the SNP annotation built in step 2.
3. Run the regression.

Example
+++++++
Running BAGHERA on the UK Biobank summary statistics for breast cancer, using EUR LD scores
and the Gencode annotation. ::

  $ baghera-tool create-files baghera/primary_data/ -l baghera/eur_w_ld_chr/ -a baghera/gencode.v27lift37.basic.annotation.gtf
  $ baghera-tool generate-snp-file baghera/primary_data/100032.assoc.tsv --input_type ukbb --output-file baghera/primary_data/SNPs_breast.csv
  $ baghera-tool regression baghera/primary_data/SNPs_breast.csv breast baghera/results/  baghera/primary_data/genesTable.csv

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
