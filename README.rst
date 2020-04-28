Bayesian Gene Heritability Analysis
===================================
.. image:: https://img.shields.io/github/v/tag/stracquadaniolab/baghera
.. image:: https://github.com/stracquadaniolab/baghera/workflows/Release%20package/badge.svg
.. image:: https://anaconda.org/stracquadaniolab/baghera/badges/platforms.svg
.. image:: https://anaconda.org/stracquadaniolab/baghera/badges/version.svg

The Bayesian Gene Heritability Analysis software (BAGHERA) estimates the contribution
to the heritability of a trait/disease of all the SNPs in the genome (genome-wide heritability)
and those nearby protein-coding genes (gene-level heritability).to the heritability of
a trait/disease of all the SNPs in the genome (genome-wide heritability)
and those nearby protein-coding genes (gene-level heritability).

BAGHERA requires only summary statistics from a Genome-wide Association Study (GWAS),
LD scores calculated from a population matching the ethnicity of the GWAS study and
a gene annotation file in GTF format.

Installation
------------

The easiest and fastest way to install BAGHERA using conda::

$ conda install -c stracquadaniolab -c bioconda -c conda-forge


Tutorial
---------------

A typical BAGHERA analysis consists of 3 steps:

1. Build a SNP annotation file, where SNPs are annotated to genes and LD scores
are assigned. We use `precomputed ld-score <https://github.com/bulik/ldsc>`_ ,
from the set of variants for the European population of 1000 Genomes, and  the
genes in the `Gencode v31 annotations
<https://www.gencodegenes.org/releases/current.html>`_ , using only the protein coding terms.
To cope with overlapping genes, we clustered them, obtaining a dataset of
15000 non-overlapping genes. For the annotation, we use a 50 kb window. ::

    $ baghera-tool create-files -l <ldscore_folder> -a <annotation.gtf> -s <ld_annotated_snps> -g <genes_table>


2. Annotate summary statistics with the SNP annotation built in step 2. We used the summary statistics `here <http://www.nealelab.is/uk-biobank>`_::

    $ baghera-tool generate-snp-file -s <stats file> -i <input_type> -o <snps_file> -a <ld_annotated_snps>

3. Run the regression::

    $ baghera-tool gene-heritability <snps_file> <results_table> <summary_table> <log_file> --sweeps <samples> --burnin <tuning> --n-chains <chains> --n-cores <cores> -m <models>



Example
+++++++

Running BAGHERA on the UK Biobank summary statistics for breast cancer, using EUR LD scores
and the Gencode annotation. ::

  $ baghera-tool create-files -l data/eur_w_ld_chr/ -a data/gencode.v31lift37.basic.annotation.gtf -s data/ld_annotated_gencode_v31.csv -g data/genes_gencode_v31.csv
  $ baghera-tool generate-snp-file -s data/C50.gwas.imputed_v3.both_sexes.tsv -i position_ukbb -o data/c50.snps.csv -a data/ld_annotated_gencode_v31.csv
  $ baghera-tool gene-heritability data/c50.snps.csv data/results_normal_c50.csv data/summary_normal_c50.csv data/log_normal_c50.txt --sweeps 10000 --burnin 2500 --n-chains 4 --n-cores 4 -m normal


Workflow
++++++++

Alongside BAGHERA, we are providing a snakemake workflow `repository <https://github.com/stracquadaniolab/workflow-baghera>`_ with sample data.



Authors
-------
- Viola Fanfani (v.fanfani@sms.ed.ac.uk): mantainer.
- Giovanni Stracquadanio (giovanni.stracquadanio@ed.ac.uk)

Citation
--------
Gene-level heritability analysis explains the polygenic architecture of cancer.
Viola Fanfani, Luca Citi, Adrian L. Harris, Francesco Pezzella, Giovanni Stracquadanio
bioRxiv 599753; doi: https://doi.org/10.1101/599753

Issues
------

We just released a major upgrade of the code, please report any issue.
