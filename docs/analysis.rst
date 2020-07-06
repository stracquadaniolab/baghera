Running the analysis
=============================

To date, the available code for the analysis allows to perform:

- a **gene-level Bayesian analysis**, that returns the probability of the single gene to affect the trait heritability with an higher probability than by chance

- a **genome-wide Bayesian regression**, able to estimate the heritability of the trait.


Gene-level analysis
+++++++++++++++++++++

The algorithm to perform a gene-level analysis is runnable as `baghera-tool gene-regression -h`
here is the complete list of options of the function:

.. code-block:: bash

    usage: baghera-tool gene-heritability [-h] [--sweeps SWEEPS] [-b BURNIN]
                                        [--n-chains N_CHAINS]
                                        [--n-cores N_CORES] [-N N_1KG]
                                        [-c CHROMOSOME] [--snp-thr SNP_THR]
                                        [--sep SEP] [-m MODEL] [-f]
                                        input-snp-filename output-genes-filename
                                        output-summary-filename logger-filename


        Performs bayesian gene-level heritability analysis.
        As input it needs the annotated snps file created with generate-snps-file.

        From command line one can specify all the parameters for the sampler
        (sweeps, burnin, chains and cores) and the parameters for the SNPs
        and genes filtering.

        Specify the gamma model by passing --model gamma


    positional arguments:
    input-snp-filename    Data Input, use the SNPs file from dataParse
    output-genes-filename
                            output file for gene-level results, use .csv
    output-summary-filename
                            output file for the genomewide results summary, use
                            .csv
    logger-filename       file for the logger, use a txt

    optional arguments:
    -h, --help            show this help message and exit
    --sweeps SWEEPS       number of samples for each chain (default: 1000)
    -b BURNIN, --burnin BURNIN
                            number of burnin samples (default: 1000)
    --n-chains N_CHAINS   number of chains of the sampler (default: 4)
    --n-cores N_CORES     number of parallel cores to use (default: 4)
    -N N_1KG, --N-1kG N_1KG
                            number of SNPs onwhich the LD-score is calculated
                            (default: 1290028)
    -c CHROMOSOME, --chromosome CHROMOSOME
                            chromosome on which the analysis is run (default:
                            'all')
    --snp-thr SNP_THR     threshold for the minimum number of SNPs in a gene
                            (default: 10)
    --sep SEP             separator for the input files, use t for tab separated
                            (not ) (default: ',')
    -m MODEL, --model MODEL
                            specify the model for the regression, one betwenn
                            normal/gamma (default: 'normal')
    -f, --fix-intercept   False





From the input annotated file

1. **gene_level_analysis.log** : is a log file that reports the parameters of
the analysis and the output results. *This is just a synthesis of the analysis,
not the run-time log file which can be found in the logs folder*

2. **summary.csv** : output summary statistics of the traces. *e* is the
intercept and *mi* is the regression slope (heritability), *W* is the variance
term and *herTOT* is the weighted summation of the signle gene term. For each of
them the principal stats are reported.

3. **results.csv**: is the output file, where the stats for the single gene are
reported.

Results Files format
---------------------

The output file ``results.csv`` is a csv file, which contains all the results for each gene. Here is the a description of the format:

- *chrom* , *start* , *stop* , *name* : we are considering all the genes after the clustering, some of them have a composite name, and that are included in the analysis.
- *LDvariance*, variance of the ld-score within the gene
- *SNPs*, number of SNPs within the gene
- *StatsMax,StatsMean,StatsMin,StatsVariance*, for the chi-squared statistics within the gene max, min, average and variance values
- *P*, output probability
- *bg_mean, bg_var, bg_median, bg_5perc, bg_95perc*, are the stats for the single gene heritability weights.
- *mi_mean, mi_median*, genome-wide heritability posterior
- *h2g*, weighted sum of the single gene heritability


Genome-wide Heritability
+++++++++++++++++++++++++

The algorithm to get the heritability estimate, performing a genomewide regression on the SNPs

.. code-block:: bash

    usage: baghera-tool gw-heritability [-h] [--sweeps SWEEPS] [-b BURNIN]
                                        [--n-chains N_CHAINS] [--n-cores N_CORES]
                                        [-N N_1KG] [-c CHROMOSOME] [--sep SEP]
                                        [-m MODEL] [-f]
                                        input-snp-filename output-summary-filename
                                        logger-filename

        Computes the genome-wide estimate heritability using Bayesian regression.


    positional arguments:
    input-snp-filename    Data Input, use the SNPs file from dataParse
    output-summary-filename
                            output file for the genomewide results summary, use
                            .csv
    logger-filename       file for the logger, use a txt

    optional arguments:
    -h, --help            show this help message and exit
    --sweeps SWEEPS       number of samples for each chain (default: 1000)
    -b BURNIN, --burnin BURNIN
                            number of burnin samples (default: 1000)
    --n-chains N_CHAINS   number of chains of the sampler (default: 4)
    --n-cores N_CORES     number of parallel cores to use (default: 4)
    -N N_1KG, --N-1kG N_1KG
                            number of SNPs onwhich the LD-score is calculates
                            (default: 1290028)
    -c CHROMOSOME, --chromosome CHROMOSOME
                            chromosome on which the analysis is run (default:
                            'all')
    --sep SEP             separator for the input files, use t for tab separated
                            (not ) (default: ',')
    -m MODEL, --model MODEL
                            regression model (default: 'normal')
    -f, --fix-intercept   False


As output, two files are created.

1. **Log file.txt** : is a log file that reports the parameters of the analysis
and the output results. *This is just a synthesis of the analysis, not the
run-time log file which can be found in the logs folder*

2. **Results_file.csv** : output summary statistics of the traces. *e* is the
intercept and *mi* is the regression slope (heritability). For each of them the
principal stats are reported.

This analysis is the bayesian version of the LDSC hertiability model.
