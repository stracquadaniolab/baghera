import sys
import logging
import numpy as np
import theano.tensor as tt
import pymc3 as pm

import pandas as pd
import HTSeq
import os

from baghera_tool import gene_regression as gr
from baghera_tool import heritability
import baghera_tool.logging as log
from baghera_tool import snps as s
from baghera_tool import genes as g


#####################################################################
############### REGRESSION ANALYSIS #################################
#####################################################################


def gene_heritability(
    input_snp_filename: "Data Input, use the SNPs file from dataParse",
    output_genes_filename: 'output file for gene-level results, use .csv',
    output_summary_filename: 'output file for the genomewide results summary, use .csv',
    logger_filename: 'file for the logger, use a txt',
    sweeps: "number of samples for each chain" = 1000,
    burnin: "number of burnin samples" = 1000,
    n_chains: "number of chains of the sampler" = 4,
    n_cores: "number of parallel cores to use" = 4,
    N_1kG: "number of SNPs onwhich the LD-score is calculated" = 1290028,
    chromosome: "chromosome on which the analysis is run" = "all",
    snp_thr: "threshold for the minimum number of SNPs in a gene" = 10,
    sep: "separator for the input files, use t for tab separated (not \t)" = ",",
    model: 'specify the model for the regression, one betwenn normal/gamma' = 'normal',
    fix_intercept = False,
    ):

    """
    Performs bayesian gene-level heritability analysis.
    As input it needs the annotated snps file created with generate-snps-file.

    From command line one can specify all the parameters for the sampler
    (sweeps, burnin, chains and cores) and the parameters for the SNPs
    and genes filtering.

    Specify the gamma model by passing --model gamma
    """

    # Initialisation of the logger
    output_logger = log.setup_logger("output_logger", logger_filename)
    log.initialise_log(output_logger,
                    'gene level regression, model: %s' %model,
                    [input_snp_filename],
                    [output_genes_filename,output_summary_filename],
                    sweeps,
                    burnin,
                    chromosome = str(chromosome),
                    other_params_diz = {'chains': n_chains, 'cores': n_cores, 'SNP threshold': snp_thr})

    # Initialisation function, it reads the summary stats file, filters the SNPs,
    # creates the output files

    logging.info("Start Analysis")

    snps = s.Snps()
    # read table
    logging.info("Reading SNP file: %s,\n\t with %s delimiter"%(input_snp_filename, sep))
    snps.read_table(input_snp_filename, separator=sep)
    # generate chi squared stats
    snps.generate_stats()
    # update the summary stats
    snps.update_summary()
    output_logger.info(" Sample size " + str(snps.n_patients) + "\n")



    snps.apply_filter_table(s.baghera_filter)
    snps.update_summary()
    output_logger.info("After baghera init filter.\n\t Number of SNPs: %s\n\t Number of genes: %s\n" \
        %(str(snps.n_snps), str(snps.n_genes)) )

    # Non coding SNPs are assigned to a dummy gene, such that the regression is done on the entire SNPs' set
    snps.rename_non_annotated(name='NonCoding')

    if chromosome != "all":
        snps.apply_filter_table(snps.cut_single_chrom, **{'chromosome': chromosome})
        output_logger.info(
            "Analysis restricted to chr %s" %str(chromosome) )

        snps.update_summary()
        output_logger.info("Analysis. Number of SNPs: %s\n,  Number of genes: %s\n" \
            %(str(snps.n_snps), str(snps.n_genes)) )

    # Creates the genes table with the number of SNPs for each gene and the basic stats values
    genes=g.Genes()
    genes.initialise_genes(snps.table.copy(), snps_thr=snp_thr)

    output_logger.info("Output gene table initialised:\nNumber of genes: %s\n" \
        %(str(genes.n_genes)) )

    snps.set_non_annotated(genes.cut_genes, 'NonCoding')

    if model == 'gamma':
        result = gr.analyse_gamma(snps, output_summary_filename, output_logger,
                                 sweeps, burnin, n_chains, n_cores, N_1kG, fix_intercept,
                                 )
    else:
        result = gr.analyse_normal(snps, output_summary_filename, output_logger,
                sweeps, burnin, n_chains, n_cores, N_1kG, fix_intercept,
        )

    logging.info("Saving genes table")
    genes.table = genes.table.merge(
        result, left_index=False, left_on="name", right_on="name")

    k = genes.table.n_snps / float(N_1kG)
    genes.table["h2g"] = genes.table.bg_mean.astype("float") * k

    genes.table = genes.table.sort_values(by=["P", "bg_median"])

    genes.save_table(output_genes_filename)

    non_coding = genes.table[genes.table.name == "NonCoding"]
    h2g_tot = np.sum(genes.table["h2g"].values) - non_coding["h2g"].values

    output_logger.info(" Non coding heritability : " +
                       str(non_coding["h2g"].values) + "\n")
    output_logger.info(" Coding heritability : " + str(h2g_tot) + "\n")


##################################################################
############# GENOME-WIDE HERITABILITY ###########################
##################################################################



def gw_heritability(
    input_snp_filename: "Data Input, use the SNPs file from dataParse",
    output_summary_filename: 'output file for the genomewide results summary, use .csv',
    logger_filename: 'file for the logger, use a txt',
    sweeps: "number of samples for each chain" = 1000,
    burnin: "number of burnin samples" = 1000,
    n_chains: "number of chains of the sampler" = 4,
    n_cores: "number of parallel cores to use" = 4,
    N_1kG: "number of SNPs onwhich the LD-score is calculates" = 1290028,
    chromosome: "chromosome on which the analysis is run" = "all",
    sep: "separator for the input files, use t for tab separated (not \t)" = ",",
    model: 'regression model'='normal',
    fix_intercept = False,
    ):

    """
    Computes the genome-wide estimate heritability using Bayesian regression.
    """

    # Initialisation of the logger
    output_logger = log.setup_logger("output_logger", logger_filename)
    log.initialise_log(output_logger,
                    'genome-wide regression, model: %s' %model,
                    [input_snp_filename],
                    [output_summary_filename],
                    sweeps,
                    burnin,
                    chromosome = str(chromosome),
                    other_params_diz = {'chains': n_chains, 'cores': n_cores})

    # Initialisation function, it reads the summary stats file, filters the SNPs,
    # creates the output files

    logging.info("Start Analysis")

    snps = s.Snps()
    # read table
    snps.read_table(input_snp_filename, separator=sep)
    # generate chi squared stats
    snps.generate_stats()
    # update the summary stats
    snps.update_summary()
    output_logger.info(" Sample size " + str(snps.n_patients) + "\n")


    snps.apply_filter_table(s.baghera_filter)
    snps.update_summary()
    output_logger.info("After baghera init filter.\nNumber of SNPs: %s\nNumber of genes: %s\n" \
        %(str(snps.n_snps), str(snps.n_genes)) )

    # Non coding SNPs are assigned to a dummy gene, such that the regression is done on the entire SNPs' set
    snps.rename_non_annotated(name='NonCoding')

    if chromosome != "all":
        snps.apply_filter_table(snps.cut_single_chrom, **{'chromosome': chromosome})
        output_logger.info(
            "Analysis restricted to chr %s" %str(chromosome) )

        snps.update_summary()
        output_logger.info("Analysis. Number of SNPs: %s\n,  Number of genes: %s\n" \
            %(str(snps.n_snps), str(snps.n_genes)) )


    if model =='normal':
        [intercept, slope] = heritability.gw_normal(snps, output_summary_filename, output_logger,
                                            sweeps, burnin, n_chains, n_cores, N_1kG, fix_intercept)
    elif model=='gamma':
        [intercept, slope] = heritability.gw_normal(snps, output_summary_filename, output_logger,
                                            sweeps, burnin, n_chains, n_cores, N_1kG, fix_intercept)
    else:
        logging.info('Normal model by default')
        [intercept, slope] = heritability.gw_normal(snps, output_summary_filename, output_logger,
                                            sweeps, burnin, n_chains, n_cores, N_1kG, fix_intercept)
    logging.info("Analysis complete")
