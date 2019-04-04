Creating the annotated files
============================

Annotated variants
------------------

To run the GWAS analysis the variants in the study need to be annotated with the ld-score and to a gene.

We use `precomputed ld-score <https://github.com/bulik/ldsc>`_ , from the set of variants for the European population of 1000 Genomes, and  the genes in the `Gencode v27 annotations <https://www.gencodegenes.org/releases/current.html>`_ , only the protein coding ones. To cope with overlapping genes, we clustered them, obtaining a dataset of 15000 non-overlapping genes. For the annotation, we use a 50 kb window.
The resulting dataset of annotated variants has around 1.3 millions SNPs, 55% of which are annotated with a gene.

*Please note that this file has already been created, to process the data skip to the next section*

It is possible to annotate a different set of variants, for example another reference panel, using the `create-files` function.
For the moment it only supports .gtf files for the genes annotation and the LD-score folder with the structure in
`<https://github.com/bulik/ldsc>`_

To create the SNPs dataset use the `baghera-tool create-files` command

.. autofunction:: baghera_tool.preprocess.create_files



Create the dataset
------------------

There are only two accepted inputs for the variants statistics:

- a **sumstats** file, .sumstats.txt, a tab-delimited text file which, for our code needs to have at least these three columns: SNP, Z, N, representing the SNP id, its effect size and the number of patiens in the study.

- an **UKBB** file, .assoc.tsv, a tab-delimited text file which, for our code needs to have at least these three columns: rsid, nCompleteSamples,tstat, representing the SNP id, the number of patiens in the study and the variant's effect size.

The actual dataset used for the analysis is created merging the SNPs in the study and those that are the Outputof the annotated LD file.

To create the SNPs dataset use the `baghera-tool generate-SNPs-file` command

.. autofunction:: baghera_tool.preprocess.generate_snp_file
