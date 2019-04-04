Running the analysis
=============================

To date, the available code for the analysis allows to perform:

- a genome-wide Bayesian regression, able to estimate the heritability of the trait
- a gene-level Bayesian analysis, that returns the probability of the single gene to affect the trait heritability with an higher probability than by chance.


Heritability
++++++++++++

The algorithm to get the heritability estimate, performing a genomewide regression on the SNPs, is heritability.py and is runnable as:

.. autofunction:: baghera_tool.heritability.heritability

The output files are written in the output folder specified by the user and
have the suffix specified in the name.

As output, two files are created.

1. **heritability_suffix_.log** : is a log file that reports the parameters of the analysis and the output results. *This is just a synthesis of the analysis, not the run-time log file which can be found in the logs folder*

2. **heritability_$suffix.csv** : output summary statistics of the traces. *e* is the intercept and *mi* is the regression slope (heritability). For each of them the principal stats are reported.

Gene-based analysis
+++++++++++++++++++

The algorithm to perform a gene-level analysis is runnable as `baghera-tool regression`
here is the complete list of options of the function.

.. autofunction:: baghera_tool.regression.regression

*Please notice* that the regression with the BAGHERA-gamma model can be
run using the --gamma flag.

*Please notice* to make sure the mapping to the genes is correct, use the flag -m
which checks that there is a matching in the name, regardless of the order of the
clustering names.

The output files are written in the output folder specified by the user and
have the suffix specified in the name.

As output, two files are created.

1. **regression_$suffix.log** : is a log file that reports the parameters of the analysis and the output results. *This is just a synthesis of the analysis, not the run-time log file which can be found in the logs folder*

2. **regression_$suffix.csv** : output summary statistics of the traces. *e* is the intercept and *mi* is the regression slope (heritability), *W* is the variance term and *herTOT* is the weighted summation of the signle gene term. For each of them the principal stats are reported.

3. **genesResults_$suffix.csv**: is the output file, where the stats for the single gene are reported.

Output Files format
-------------------

The output file ``genesResults_$SUFFIX.csv`` is a csv file, which contains all the results for each gene. Here is the a description of the format:

- *chrom* , *start* , *stop* , *name* : we are considering all the genes after the clustering, some of them have a composite name, and that are included in the analysis.
- *LDvariance*, variance of the ld-score within the gene
- *SNPs*, number of SNPs within the gene
- *StatsMax,StatsMean,StatsMin,StatsVariance*, for the chi-squared statistics within the gene max, min, average and variance values
- *P*, output probability
- *bg_mean, bg_var, bg_median, bg_5perc, bg_95perc*, are the stats for the single gene heritability weights.
- *mi_mean, mi_median*, genome-wide heritability posterior
- *h2g*, weighted sum of the single gene heritability
