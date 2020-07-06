import sys
import logging
import numpy as np
import theano.tensor as tt
import pymc3 as pm
import datetime
import pandas as pd
import HTSeq
import os

from baghera_tool.logging import setup_logger


def trace_sd(x):
    return pd.Series(np.std(x, 0), name="sd")


def trace_quantiles(x):
    return pd.DataFrame(pm.quantiles(x, [5, 50, 95]))


def trace_median(x):
    return pd.Series(np.median(x, 0), name="median")


def gw_normal(
    snps_object,
    output_summary_filename,
    output_logger,
    SWEEPS,
    TUNE,
    CHAINS,
    CORES,
    N_1kG,
    fix_intercept=False,
):
    """ Bayesian heritability analysis: requires a dataFrame with the SNPs as input. The file can be created with the dataParse function.
        """

    logging.info("Bayesian analysis started")
    snps_object.table = snps_object.table.reset_index(drop=True)
    nPATIENTS = snps_object.table["sample_size"][0]

    with pm.Model() as model:
        if fix_intercept:
            e = pm.Uniform("e", 0.9999, 1.0000001)
        else:
            e = pm.Normal("e", 1, 1)
        mi = pm.Beta("mi", 1, 2)
        fixed_variable = pm.Normal(
            "fxd",
            mu=(nPATIENTS / N_1kG) * mi * snps_object.table["l"] + e,
            sd=np.sqrt(np.asarray(snps_object.table["l"])),
            observed=snps_object.table["stats"],
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

    su = pm.summary(trace, extend=True, stat_funcs=[trace_median, trace_quantiles],)

    logging.info("Writing output")
    su.to_csv(output_summary_filename, sep=",", mode="w")

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

def gw_gamma(
    snps_object,
    output_summary_filename,
    output_logger,
    SWEEPS,
    TUNE,
    CHAINS,
    CORES,
    N_1kG,
    fix_intercept=False,
):
    """ Bayesian heritability analysis: requires a dataFrame with the SNPs as input. The file can be created with the dataParse function.
        """

    logging.info("Bayesian analysis started")
    snps_object.table = snps_object.table.reset_index(drop=True)
    n_patients = snps_object.table["sample_size"][0]

    with pm.Model() as model:

        e = pm.Normal("e", 1, 1)
        mi = pm.Beta("mi", 1, 1)

        k = pm.Gamma("k", alpha=1, beta=1)
        fixed_variable = pm.Gamma(
            "fxd",
            alpha=k * ((n_patients) * mi * (snps_object.table["l"]) + 1) ** 2,
            beta=k * ((n_patients) * mi * (snps_object.table["l"]) + 1),
            observed=snps_object.table["stats"],
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

    su = pm.summary(trace, extend=True, stat_funcs=[trace_median, trace_quantiles],)

    logging.info("Writing output")
    su.to_csv(output_summary_filename, sep=",", mode="w")

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
