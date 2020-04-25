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

def subtract(x, y):
    return x - y


def analyse_normal(snps_object, output_summary_filename, output_logger,
            SWEEPS, TUNE, CHAINS, CORES, N_1kG,fix_intercept = False,
            ):

    """
    Bayesian hierarchical regression on the dataset with the normal model.

    :params snps_object: snps instance
    :params output_summary_filename: output summary table
    :params output_logger: logger instance
    :params SWEEPS: samples for each chain
    :params TUNE: burn-in samples
    :params CHAINS: number of chains
    :params CORES: number of cores
    :params N1kG: number of SNPs
    :params fix_intercept: if True the model fixes the intercept.
    """


    snp_dataset = snps_object.table.copy().reset_index(drop=True)
    n_patients = snps_object.n_patients
    nSNP = snps_object.n_snps

    # to run the regression as a mixed effect model, I need a vector (cat) to assign each SNP to its gene
    idx = 0
    cat = np.zeros(nSNP)
    Mg = []
    genes = []
    # g2 are all the SNPs inside the gene
    for k2, g2 in snp_dataset.groupby("gene"):
        cat[g2.index] = int(idx)
        idx += 1
        genes.append(k2)
        Mg.append(len(g2))
    cat = cat.astype(int)

    logging.info("Model Evaluation Started")
    logging.info("Model Evaluation Started")
    print('Average stats: %f' %np.mean(snps_object.table["stats"].values))

    with pm.Model() as model:
        W = pm.InverseGamma("W", alpha=1.0, beta=1.0)
        e = pm.Normal("e", mu=1, sd=1)
        mi = pm.Beta("mi", 1, 1)
        beta = pm.Normal("beta", mu=mi, sd=W, shape=idx)
        diff = pm.Deterministic("diff", subtract(beta, mi))
        herTOT = pm.Deterministic("herTOT", tt.sum(beta * Mg / N_1kG))
        if fix_intercept:
            fixed_variable = pm.Normal(
                "fxd",
                mu=(n_patients / N_1kG) * beta[cat] * (snp_dataset["l"]) + 1,
                sd=np.sqrt(np.asarray(snp_dataset["l"])),
                observed=snp_dataset["stats"],
            )  #

        else:
            fixed_variable = pm.Normal(
                "fxd",
                mu=(n_patients / N_1kG) * beta[cat] * (snp_dataset["l"]) + e,
                sd=np.sqrt(np.asarray(snp_dataset["l"])),
                observed=snp_dataset["stats"],
            )  #

        trace = pm.sample(
            SWEEPS,
            tune=TUNE,
            chains=CHAINS,
            cores=CORES,
            nuts_kwargs=dict(target_accept=0.90),
        )

        if CHAINS > 1:
            logging.info("evaluating Gelman-Rubin")
            GR = pm.diagnostics.gelman_rubin(trace, varnames=["mi"])
            output_logger.info(
                "DIAGNOSTIC (gelman-rubin) "
                + str(GR)
                + "\n"
                + "(If this number is >> 1 the method has some convergence problem, \n try increasing the number of s and b)"
            )

    logging.info("Writing output")
    # save general stats to summary file
    su = pm.summary(
        trace,
        varnames=["mi", "herTOT", "e", "W"],
        extend=True,
        stat_funcs=[trace_median, trace_quantiles],
    )

    su.to_csv(output_summary_filename, sep=",", mode="w")

    d = dict()
    d["beta"] = trace["beta"]

    e_GW = np.mean(trace["e"])
    e_GW_sd = np.std(trace["e"])
    output_logger.info(" Intercept: " + str(e_GW) +
                       " (sd= " + str(e_GW_sd) + ")\n")

    W = np.mean(trace["W"])
    W_sd = np.std(trace["W"])
    output_logger.info(" W: " + str(W) + " (sd= " + str(W_sd) + ")\n")

    herTOT = np.median(trace["herTOT"])
    herTOT_sd = np.std(trace["herTOT"])
    output_logger.info(
        " heritability from genes: " +
        str(herTOT) + " (sd= " + str(herTOT_sd) + ")\n"
    )

    mi_mean = np.mean(trace["mi"], axis=0)
    mi_median = np.median(trace["mi"], axis=0)
    mi_std = np.std(trace["mi"], axis=0)
    mi_5perc = np.percentile(trace["mi"], 5, axis=0)
    mi_95perc = np.percentile(trace["mi"], 95, axis=0)
    output_logger.info(
        " Heritability: "
        + str(mi_mean)
        + " (std= "
        + str(mi_std)
        + ")\n"
        + "[ 5th perc= "
        + str(mi_5perc)
        + ","
        + " 95 perc= "
        + str(mi_95perc)
        + "]\n"
    )

    Prob = (np.sum(trace["diff"] > 0, axis=0) /
            len(trace["diff"]))[:, np.newaxis]

    data = np.hstack((np.asarray(genes)[:, np.newaxis], Prob))
    df = pd.DataFrame(data, columns=("name", "P"))

    df["bg_mean"] = np.mean(d["beta"], axis=0)[:, np.newaxis]
    df["bg_median"] = np.median(d["beta"], axis=0)[:, np.newaxis]
    df["bg_var"] = np.var(d["beta"], axis=0)[:, np.newaxis]
    df["bg_5perc"] = np.percentile(d["beta"], 5, axis=0)[:, np.newaxis]
    df["bg_95perc"] = np.percentile(d["beta"], 95, axis=0)[:, np.newaxis]

    df["mi_mean"] = mi_mean
    df["mi_median"] = mi_median

    logging.info("analysis done")
    return df

def analyse_gamma(snps_object, output_summary_filename, output_logger,
            SWEEPS, TUNE, CHAINS, CORES, N_1kG,fix_intercept = False,):

    """
    Bayesian hierarchical regression on the dataset with the gamma model.

    :params snps_object: snps instance
    :params output_summary_filename: output summary table
    :params output_logger: logger instance
    :params SWEEPS: samples for each chain
    :params TUNE: burn-in samples
    :params CHAINS: number of chains
    :params CORES: number of cores
    :params N1kG: number of SNPs
    :params fix_intercept: if True the model fixes the intercept.
    """

    snp_dataset = snps_object.table.copy().reset_index(drop=True)
    n_patients = snps_object.n_patients
    nSNP = snps_object.n_snps

    # to run the regression as a mixed effect model, I need a vector (cat) to assign each SNP to its gene
    idx = 0
    cat = np.zeros(nSNP)
    Mg = []
    genes = []
    # g2 are all the SNPs inside the gene
    for k2, g2 in snp_dataset.groupby("gene"):
        cat[g2.index] = int(idx)
        idx += 1
        genes.append(k2)
        Mg.append(len(g2))
    cat = cat.astype(int)

    logging.info("Model Evaluation Started")
    print('Average stats: %f' %np.mean(snps_object.table["stats"].values))

    with pm.Model() as model:
        e = pm.Normal("e", mu=1, sd=0.001)
        mi = pm.Beta("mi", 1, 1)
        beta = pm.Gamma("beta", alpha=mi, beta=N_1kG, shape=idx)
        diff = pm.Deterministic("diff", subtract(beta, mi / N_1kG))
        herTOT = pm.Deterministic("herTOT", tt.sum(beta * Mg))
        if fix_intercept:
            fixed_variable = pm.Normal(
                "fxd",
                mu=(n_patients ) * beta[cat] * (snps_object.table["l"]) + 1,
                sd=np.sqrt(np.asarray(snps_object.table["l"])),
                observed=snps_object.table["stats"],
            )
        else:
            fixed_variable = pm.Normal(
                "fxd",
                mu=(n_patients ) * beta[cat] * (snps_object.table["l"]) + e,
                sd=np.sqrt(np.asarray(snps_object.table["l"])),
                observed=snps_object.table["stats"],
            )

        # step = pm.Metropolis()
        trace = pm.sample(
            SWEEPS,
            tune=TUNE,
            chains=CHAINS,
            cores=CORES,
            nuts_kwargs=dict(target_accept=0.90),
        )

        if CHAINS > 1:
            logging.info("evaluating Gelman-Rubin")
            GR = pm.diagnostics.gelman_rubin(trace, varnames=["mi"])
            output_logger.info(
                "DIAGNOSTIC (gelman-rubin) "
                + str(GR)
                + "\n"
                + "(If this number is >> 1 the method has some convergence problem, \n try increasing the number of s and b)"
            )

    logging.info("Writing output")
    # save general stats to summary file
    su = pm.summary(
        trace,
        varnames=["mi", "herTOT", "e"],
        extend=True,
        stat_funcs=[trace_median, trace_quantiles],
    )
    su.to_csv(
        output_summary_filename, sep=",", mode="w"
    )

    d={}
    d["beta"] = N_1kG * trace["beta"]

    e_GW = np.mean(trace["e"])
    e_GW_sd = np.std(trace["e"])
    output_logger.info(
        " Intercept: " + str(e_GW) + " (sd= " + str(e_GW_sd) + ")\n"
    )

    herTOT = np.median(trace["herTOT"])
    herTOT_sd = np.std(trace["herTOT"])
    output_logger.info(
        " heritability from genes: "
        + str(herTOT)
        + " (sd= "
        + str(herTOT_sd)
        + ")\n"
    )

    mi_mean = np.mean(trace["mi"], axis=0)
    mi_median = np.median(trace["mi"], axis=0)
    mi_std = np.std(trace["mi"], axis=0)
    mi_5perc = np.percentile(trace["mi"], 5, axis=0)
    mi_95perc = np.percentile(trace["mi"], 95, axis=0)
    output_logger.info(
        " Heritability: "
        + str(mi_mean)
        + " (std= "
        + str(mi_std)
        + ")\n"
        + "[ 5th perc= "
        + str(mi_5perc)
        + ","
        + " 95 perc= "
        + str(mi_95perc)
        + "]\n"
    )

    Prob = (np.sum(trace["diff"] > 0, axis=0) / len(trace["diff"]))[
        :, np.newaxis
    ]

    data = np.hstack((np.asarray(genes)[:, np.newaxis], Prob))
    df = pd.DataFrame(data, columns=("name", "P"))


    df["bg_mean"] = np.mean(d["beta"], axis=0)[:, np.newaxis]
    df["bg_median"] = np.median(d["beta"], axis=0)[:, np.newaxis]
    df["bg_var"] = np.var(d["beta"], axis=0)[:, np.newaxis]
    df["bg_5perc"] = np.percentile(d["beta"], 5, axis=0)[:, np.newaxis]
    df["bg_95perc"] = np.percentile(d["beta"], 95, axis=0)[
        :, np.newaxis
    ]

    df["mi_mean"] = mi_mean
    df["mi_median"] = mi_median

    return df
