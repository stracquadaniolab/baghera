#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jan 19 15:54:38 2018

@author: viola
"""
import sys
import logging
import numpy as np
import theano.tensor as tt
import pymc3 as pm
import datetime
import pandas as pd
import HTSeq
import os


def setup_logger(name, log_file, level=logging.INFO):
    """Function setup as many loggers as you want"""

    handler = logging.FileHandler(log_file)
    formatter = logging.Formatter("[%(levelname)s]: %(message)s")
    handler.setFormatter(formatter)

    logger = logging.getLogger(name)
    logger.setLevel(level)
    logger.addHandler(handler)

    return logger


def trace_sd(x):
    return pd.Series(np.std(x, 0), name="sd")


def trace_quantiles(x):
    return pd.DataFrame(pm.quantiles(x, [5, 50, 95]))


def trace_median(x):
    return pd.Series(np.median(x, 0), name="median")


def gw_heritability(snp_dataset, folder, output_logger,
                    SWEEPS, TUNE, CHAINS, CORES, N_1kG, SUFFIX,):
    """ Bayesian heritability analysis: requires a dataFrame with the SNPs as input. The file can be created with the dataParse function.
        """

    logging.info("Bayesian analysis started")
    snp_dataset = snp_dataset.reset_index(drop=True)
    nPATIENTS = snp_dataset["sample_size"][0]

    with pm.Model() as model:

        e = pm.Normal("e", 1, 1)
        mi = pm.Beta("mi", 1, 1)
        fixed_variable = pm.Normal(
            "fxd",
            mu=(nPATIENTS / N_1kG) * mi * snp_dataset["l"] + e,
            sd=np.sqrt(np.asarray(snp_dataset["l"])),
            observed=snp_dataset["z"],
        )
        trace = pm.sample(
            SWEEPS,
            tune=TUNE,
            chains=CHAINS,
            cores=CORES,
            nuts_kwargs=dict(target_accept=0.90),
        )  # ,chains=2, ,trace=db

        if CHAINS > 1:
            output_logger.info("evaluating Gelman-Rubin")
            GR = pm.diagnostics.gelman_rubin(trace, varnames=["mi"])
            output_logger.info(
                "DIAGNOSTIC (gelman-rubin) "
                + str(GR)
                + "\n"
                + "(If this number is >> 1 the method has some convergence problem, \n try increasing the number of s and b)"
            )

    su = pm.summary(
        trace,
        extend=True,
        stat_funcs=[trace_median, trace_quantiles],
    )

    logging.info("Writing output")
    su.to_csv(
        folder + "heritability_" + SUFFIX + ".csv", sep=",", mode="w"
    )

    mi_mean = np.mean(trace["mi"])
    mi_median = np.mean(trace["mi"])
    mi_std = np.std(trace["mi"])
    mi_95perc = np.percentile(trace["mi"], 95)
    mi_5perc = np.percentile(trace["mi"], 5)

    intercept = np.mean(trace["e"])

    output_logger.info(" heritability mean: " + str(mi_mean) + "\n")
    output_logger.info(" heritability median: " + str(mi_median) + "\n")
    output_logger.info(" heritability std: " + str(mi_std) + "\n")
    output_logger.info(" heritability 95perc: " + str(mi_95perc) + "\n")
    output_logger.info(" heritability 5perc: " + str(mi_5perc) + "\n")

    return [intercept, mi_mean]


def heritability(
    snp_file: "Data Input, use the SNPs file from dataParse",
    SUFFIX: "suffix for the output file",
    output_folder: "folder where to put the results",
    SWEEPS: "number of samples for each chain" = 1000,
    TUNE: "number of burnin samples" = 1000,
    CHAINS: "number of chains of the sampler" = 4,
    CORES: "number of parallel cores to use" = 4,
    N_1kG: "number of SNPs onwhich the LD-score is calculates" = 1290028,
    CHR: "chromosome on which the analysis is run" = "all",
    sep: "separator for the input files, use t for tab separated (not \t)" = ",",
):
    """
    This function computes the genome wide estimate for observed heritability
    in a Bayesian fashion. The output files are going to be saved in the specified output folder
    with the given suffix. A step by step output logger is saved as well.
    """

    folder = output_folder
    logging.info('Output files are in %s' % folder)
    # create output folder
    now = datetime.datetime.now()
    logging.info("Folder with the results: %s" % folder)

    # outputText = folder + "regression_" + SUFFIX + ".log"  # input file with SNPs and LD-scores

    output_logger = setup_logger(
        "output_logger", folder + "heritability_" + SUFFIX + ".log")

    output_logger.info(
        "Regression, bayesian gene-level regression analysis results\n "
        + "Current date & time "
        + now.strftime("%Y-%m-%d %H:%M")
    )
    output_logger.info("File: " + snp_file)
    output_logger.info(" Analysis on chr: " + CHR + "\n")
    output_logger.info(" Sweeps: " + str(SWEEPS) +
                       " , Burn: " + str(TUNE) + "\n")

    # Initialisation function, it reads the summary stats file, filters the SNPs,
    # creates the output files

    logging.info("Start Analysis")

    if sep == ',':
        with open(snp_file) as f:
            try:
                snp_dataset = pd.read_csv(f, sep=',')
            except ValueError:
                logging.exception("Wrong format of the input file")
    elif 't' == sep:
        with open(snp_file) as f:
            try:
                snp_dataset = pd.read_csv(f, sep='\t')
            except ValueError:
                logging.exception("Wrong format of the input file")

    else:
        with open(snp_file) as f:
            try:
                snp_dataset = pd.read_csv(f, sep=sep)
            except ValueError:
                logging.exception("Wrong format of the input file")

    print(snp_dataset.head())
    n_patients = snp_dataset["sample_size"][0]

    output_logger.info(" Sample size " + str(n_patients) + "\n")
    output_logger.info(" Initial Number of SNPs: " +
                       str(len(snp_dataset)) + "\n")

    snp_dataset["l"] = 1 + snp_dataset["l"] * \
        (snp_dataset["l"] > 0)  # ld-score [1,+inf)
    snp_dataset["z"] = snp_dataset["z"] ** 2  # chi-square

    # MAF filtering
    snp_dataset[snp_dataset["maf"] > 0.01]
    output_logger.info(
        "Number of SNPs after MAF>0.01 filter:" + str(len(snp_dataset)) + "\n")

    # Filter chromosome 6
    snp_dataset = snp_dataset[
        (snp_dataset.chr != 6) | ((snp_dataset.position >=
                                   34000000) | (snp_dataset.position <= 26000000))
    ]

    output_logger.info(" Number of SNP safter chr6 filter: " +
                       str(len(snp_dataset)) + "\n")

    # Non coding SNPs are assigned to a dummy gene, such that the regression is done on the entire SNPs' set
    snp_dataset = snp_dataset.replace(np.nan, "NonCoding", regex=True)

    if CHR != "all":
        snp_dataset = snp_dataset[
            snp_dataset.chr == int(CHR)
        ]  # comment to run the analysis on the whole genome
        output_logger.info(
            "Analysis restricted to a single chromosome: chr" + CHR)

    output_logger.info(
        " Number of SNPs which the analysis is conducted on : " +
        str(len(snp_dataset)) + "\n"
    )

    [intercept, slope] = gw_heritability(snp_dataset, folder, output_logger,
                                         SWEEPS, TUNE, CHAINS, CORES, N_1kG, SUFFIX)

    logging.info("Analysis complete")

# def runAnalysis():

#     # open file with the datasets
#     with open(fileSNP) as f:
#         SNP = pd.read_table(f)

#     logger.info("Filtering Data")

#     nPATIENTS = SNP["sample_size"][0]

#     output_logger.info(" Sample size " + str(nPATIENTS) + "\n")
#     output_logger.info(
#         "Initial Number of SNPs: " + str(len(SNP)) + "\n"
#     )

#     SNP[SNP["maf"] > 0.01]
#     output_logger.info(
#         " Number of SNPs after MAF>0.01 filter: " + str(len(SNP)) + "\n"
#     )

#     SNP = SNP[
#         (SNP.chr != 6)
#         | ((SNP.position >= 34000000) | (SNP.position <= 26000000))
#     ]
#     output_logger.info(
#         " Number of SNPsafter chr6 filter: " + str(len(SNP)) + "\n"
#     )

#     if CHR != "all":
#         SNP = SNP[SNP.chr == int(CHR)]

#     SNP = SNP[SNP["l"] > 0]

#     SNP["l"] = 1 + SNP["l"] * (SNP["l"] > 0)  # cleans LD-scores
#     SNP["z"] = SNP["z"] ** 2  # chi-squared

#     output_logger.info(
#         " Number of SNPs which the analysis is conducted on : "
#         + str(len(SNP))
#         + "\n"
#     )

#     logger.info("Bayesian analysis started")
#     [intercept, slope] = bayesianAnalysis(SNP)

#     logger.info("Analysis complete")


# def parserFunc():
#     parser = argparse.ArgumentParser(
#         description="Parse command line options."
#     )
#     parser.add_argument(
#         "--fileInput",
#         "-f",
#         type=str,
#         required=True,
#         help="Data Input, use the SNP_suffix file from dataParse",
#     )
#     parser.add_argument(
#         "--sweeps",
#         "-s",
#         type=int,
#         required=True,
#         help="Number of sweeps for the sampler",
#     )
#     parser.add_argument(
#         "--burn",
#         "-b",
#         type=int,
#         required=True,
#         help="Number of burnt sweeps during sampling",
#     )
#     parser.add_argument(
#         "--chromosome",
#         "-c",
#         type=str,
#         required=False,
#         default="all",
#         help="Chromosome to run the analysis on",
#     )
#     parser.add_argument(
#         "--NSNP",
#         "-N",
#         type=int,
#         default=1290028,
#         help="Number of SNPs on which the LD-score is calculated",
#     )
#     parser.add_argument(
#         "--chains",
#         "-ch",
#         type=int,
#         default=4,
#         help="Number of chains, default=4",
#     )
#     parser.add_argument(
#         "--cores",
#         "-co",
#         type=int,
#         default=1,
#         help="Number of parallel chains, by default they are sequential",
#     )
#     parser.add_argument(
#         "--outputSuffix",
#         "-o",
#         type=str,
#         required=True,
#         help="Output suffix",
#     )

#     parser.add_argument(
#         "-v",
#         "--verbose",
#         dest="verbose_count",
#         action="count",
#         default=0,
#         help="increases log verbosity for each occurence.",
#     )
#     return parser
