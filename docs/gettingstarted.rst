
Installation
------------

The easiest and fastest way to install BAGHERA using conda::

$ conda install -c stracquadaniolab -c bioconda -c conda-forge



Tutorial
---------------


A typical BAGHERA analysis consists of 3 steps, we briefly explain them here,
more details can be found in the documentation and a practical example is
in the snakemake workflow.

1- Build a SNP annotation file
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Build a SNP annotation file, where SNPs are annotated to genes and LD scores
are assigned. We use `precomputed ld-score <https://github.com/bulik/ldsc>`_ ,
from the set of variants for the European population of 1000 Genomes, and  the
genes in the `Gencode v31 annotations
<https://www.gencodegenes.org/releases/current.html>`_ , using only the protein coding terms.
To cope with overlapping genes, we clustered them, obtaining a dataset of
15000 non-overlapping genes. For the annotation, we use a 50 kb window. ::

    $ baghera-tool create-files -l <ldscore_folder> -a <annotation.gtf> -s <ld_annotated_snps> -g <genes_table>

2- Annotate summary statistics
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Annotate summary statistics with the SNP annotation built in step 1::

    $ baghera-tool generate-snp-file -s <stats file> -i <input_type> -o <snps_file> -a <ld_annotated_snps>

3- Run the regression
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Run the regression::

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


