import inspect
import sys
import logging
import numpy as np
import theano.tensor as tt
import pymc3 as pm
import datetime
import glob as g
import pandas as pd
import HTSeq
import csv
import os
import argparse


def trace_sd(x):
    return pd.Series(np.std(x, 0), name="sd")


def trace_quantiles(x):
    return pd.DataFrame(pm.quantiles(x, [5, 50, 95]))


def trace_median(x):
    return pd.Series(np.median(x, 0), name="median")

def subtract(x, y):
    return x - y

def analyse(snp_dataset,genes_final_file,folder, output_logger, SWEEPS, TUNE, CHAINS, CORES, N_1kG, SUFFIX):

    """ Gamma Bayesian hierarchical regression on the dataset: it takes two dataframes as input,
    one with the SNPs summary stats and one with the genes ( already initialised with "genesInitialise" ) """


    snp_dataset = snp_dataset.reset_index(drop=True)
    nPATIENTS = snp_dataset["sample_size"][0]

    d = dict()

    idx = 0
    nSNP = len(snp_dataset)
    cat = np.zeros(nSNP)
    Mg = []
    genes = []
    for k2, g2 in snp_dataset.groupby(
        "gene"
    ):  # g2 are all the SNPs inside the gene
        cat[g2.index] = int(idx)
        idx += 1
        genes.append(k2)
        Mg.append(len(g2))
    cat = cat.astype(int)

    logging.info("Model Evaluation Started")

    with pm.Model() as model:
        e = pm.Beta("e", alpha=1, beta=1000)
        b = pm.Gamma('b',1,1)
        mi = pm.Beta("mi", 1, 1)
        beta = pm.Gamma("beta", alpha=mi, beta=b, shape=idx)
        diff = pm.Deterministic("diff", subtract(beta, mi))
        herTOT = pm.Deterministic("herTOT", tt.sum(beta/N_1kG * Mg))
        fixed_variable = pm.ChiSquared(
            "fxd",
            nu=(nPATIENTS)/N_1kG * beta[cat] * (snp_dataset["l"]) + 1 +e,
            observed=snp_dataset["z"]+0.0000001,
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
        folder + "regression_" + SUFFIX + ".csv", sep=",", mode="w"
    )

    d["beta"] = trace["beta"]

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

    # import genes csv file, where the results are going to be written
    with open(genes_final_file, "r+") as f:
        genes_table = pd.read_csv(f, sep=',')

    genes_table = genes_table.merge(
        df, left_index=False, left_on="name", right_on="name"
    )

    k = genes_table.SNPs  # /float(N_1kG)
    genes_table["h2g"] = genes_table.bg_mean.astype("float") * k

    genes_table = genes_table.sort_values(by=["P", "bg_median"])

    genes_table.to_csv(genes_final_file, index=False, sep=",", mode="w")

    non_coding = genes_table[genes_table.name == "NonCoding"]
    h2g_tot = np.sum(genes_table["h2g"]) - non_coding["h2g"]

    output_logger.info(
        " Non coding heritability : " + str(non_coding["h2g"]) + "\n"
    )
    output_logger.info(" Coding heritability : " + str(h2g_tot) + "\n")


