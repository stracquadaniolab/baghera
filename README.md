# Bayesian Gene Heritability Analysis

![alt-ci](https://github.com/stracquadaniolab/baghera/workflows/Release%20package/badge.svg)
![alt-version](https://img.shields.io/github/v/tag/stracquadaniolab/baghera)
![alt-version](https://anaconda.org/stracquadaniolab/baghera/badges/version.svg)
![alt-platforms](https://anaconda.org/stracquadaniolab/baghera/badges/platforms.svg)


The Bayesian Gene Heritability Analysis software (BAGHERA) estimates the contribution
to the heritability of a trait/disease of all the SNPs in the genome (genome-wide heritability)
and those nearby protein-coding genes (gene-level heritability).to the heritability of
a trait/disease of all the SNPs in the genome (genome-wide heritability)
and those nearby protein-coding genes (gene-level heritability).

BAGHERA requires only summary statistics from a Genome-wide Association Study (GWAS),
LD scores calculated from a population matching the ethnicity of the GWAS study and
a gene annotation file in GTF format.

## Installation


The easiest and fastest way to install BAGHERA using conda

```
$ conda install -c stracquadaniolab -c bioconda -c conda-forge baghera
```

Tutorial
---------------

A typical BAGHERA analysis consists of 3 steps:

1. Build a SNP annotation file, where SNPs are annotated to genes and are
    assigned an LD score. We used precomputed LD scores
    (https://github.com/bulik/ldsc), from the set of variants for the European
    population in 1000 Genomes, and protein coding genes as annotated in Gencode
    v31 (https://www.gencodegenes.org/releases/current.html). Overlapping genes
    within 50Kb were considered together, obtaining a dataset of 15,000
    non-overlapping genes. To build your own annotation files, you should run
    the following command:

```
    $ baghera-tool create-files -l <ldscore_folder> -a <annotation.gtf> -s <ld_annotated_snps> -g <genes_table>
```

2. Annotate summary statistics with the SNP annotation built in step 2. We used summary statistics available at http://www.nealelab.is/uk-biobank, followed by the command below:

```
    $ baghera-tool generate-snp-file -s <stats file> -i <input_type> -o <snps_file> -a <ld_annotated_snps>
```

3. Run the regression.

```
    $ baghera-tool gene-heritability <snps_file> <results_table> <summary_table> <log_file> --sweeps <samples> --burnin <tuning> --n-chains <chains> --n-cores <cores> -m <models>
```


## Example


Running BAGHERA on the UK Biobank summary statistics for breast cancer, using
EUR LD scores and the Gencode annotation.
```
  $ baghera-tool create-files -l data/eur_w_ld_chr/ -a data/gencode.v31lift37.basic.annotation.gtf -s data/ld_annotated_gencode_v31.csv -g data/genes_gencode_v31.csv
  $ baghera-tool generate-snp-file -s data/C50.gwas.imputed_v3.both_sexes.tsv -i position_ukbb -o data/c50.snps.csv -a data/ld_annotated_gencode_v31.csv
  $ baghera-tool gene-heritability data/c50.snps.csv data/results_normal_c50.csv data/summary_normal_c50.csv data/log_normal_c50.txt --sweeps 10000 --burnin 2500 --n-chains 4 --n-cores 4 -m normal
```

## Workflow

Alongside BAGHERA, we are providing a Snakemake workflow https://github.com/stracquadaniolab/workflow-baghera, including sample data to test our method.

## Authors

- Viola Fanfani (v.fanfani@sms.ed.ac.uk): mantainer.
- Giovanni Stracquadanio (giovanni.stracquadanio@ed.ac.uk)

## Citation

_Gene-level heritability analysis explains the polygenic architecture of cancer_.
Viola Fanfani, Luca Citi, Adrian L. Harris, Francesco Pezzella, Giovanni Stracquadanio
bioRxiv 599753; doi: https://doi.org/10.1101/599753

## Issues
We just released a major upgrade of the code, please report any issue.
