import sys
import logging
import numpy as np
import theano.tensor as tt
import pymc3 as pm
import datetime
import pandas as pd
import HTSeq
import os

from baghera_tool import regression_gamma
from baghera_tool.logging import setup_logger

def trace_sd(x):
    return pd.Series(np.std(x, 0), name="sd")


def trace_quantiles(x):
    return pd.DataFrame(pm.quantiles(x, [5, 50, 95]))


def trace_median(x):
    return pd.Series(np.median(x, 0), name="median")

def subtract(x, y):
    return x - y


def initialise_genes(genes_file, genes_final_file, snp_dataset, SNPthr, output_logger, merge_flag):
    """ Function to create the results csv file, """
    logging.info("Initialisation of genes table")

    with open(genes_file) as f:
        try:
            genes_table = pd.read_csv(
                f, usecols=["chrom", "start", "stop", "name"], sep=",")
        except ValueError:
            logging.error(
                "Wrong format of the genes file, check the file has chrom, start, stop, name columns"
            )

    genes_table.loc[genes_table.index.max() + 1] = [None, None,
                                                    None, "NonCoding"]

    if merge_flag:
        # new
        names = genes_table.name.str.replace(" ", "").str.split(";")
        genes_table["name_set"] = [set(i) for i in names]

        output_logger.info(" Number of genes (overall): " +
                           str(len(genes_table)) + "\n")

        vTot = np.var(snp_dataset["l"])
        variances = []
        g = []
        d = {}
        n_snp = []
        chiVar = []
        chiMean = []
        chiMin = []
        chiMax = []

        chrom = []
        start = []
        stop = []

        snp_dataset.head()

        counter = 0
        # drop all the genes with less than 10 SNPs
        for k1, g1 in snp_dataset.groupby("gene"):
            if k1 != "NonCoding":
                if len(g1) < SNPthr:
                    # if a gene contains less than SNPthr SNPs, those are assigned to the NonCoding gene
                    snp_dataset["gene"] = snp_dataset["gene"].replace(
                        k1, "NonCoding")
                    counter += 1
                else:
                    g.append(str(k1))
                    n_snp.append(float(len(g1)))
                    variances.append(np.var(g1["l"]) / vTot)
                    chiVar.append(np.var(g1["z"]))
                    chiMean.append(np.mean(g1["z"]))
                    chiMin.append(np.min(g1["z"]))
                    chiMax.append(np.max(g1["z"]))
                    gtemp = genes_table[
                        genes_table["name_set"] == set(
                            k1.replace(" ", "").split(";"))
                    ]
                    chrom.append(gtemp.chrom.values[0])
                    start.append(gtemp.start.values[0])
                    stop.append(gtemp.stop.values[0])
                    counter += 1

        # Add row for non coding with respective statistics
        non_coding = snp_dataset[snp_dataset.gene == "NonCoding"]
        g.append("NonCoding")
        n_snp.append(float(len(non_coding)))
        variances.append(np.var(non_coding["l"]) / vTot)
        chiVar.append(np.var(non_coding["z"]))
        chiMean.append(np.mean(non_coding["z"]))
        chiMin.append(np.min(non_coding["z"]))
        chiMax.append(np.max(non_coding["z"]))
        chrom.append(None)
        start.append(None)
        stop.append(None)

        d["name"] = g
        d["chrom"] = chrom
        d["start"] = start
        d["stop"] = stop
        d["LDvariance"] = variances
        d["StatsVariance"] = chiVar
        d["StatsMean"] = chiMean
        d["StatsMax"] = chiMax
        d["StatsMin"] = chiMin
        d["SNPs"] = n_snp

        df = pd.DataFrame(d)

        output_logger.info(
            " Number of genes (with more than 10 SNPs): " + str(len(df)) + "\n"
        )

        df.to_csv(genes_final_file, index=False, sep=",", mode="w")
        logging.info("output gene file created")
        return g

    else:

        output_logger.info(" Number of genes (overall): " +
                           str(len(genes_table)) + "\n")

        vTot = np.var(snp_dataset["l"])
        variances = []
        g = []
        d = {}
        n_snp = []
        chiVar = []
        chiMean = []
        chiMin = []
        chiMax = []

        snp_dataset.head()

        counter = 0
        # drop all the genes with less than 10 SNPs
        for k1, g1 in snp_dataset.groupby("gene"):
            if k1 != "NonCoding":
                if len(g1) < SNPthr:
                    # if a gene contains less than SNPthr SNPs, those are assigned to the NonCoding gene
                    snp_dataset["gene"] = snp_dataset["gene"].replace(
                        k1, "NonCoding")
                    counter += 1
                else:
                    g.append(str(k1))
                    n_snp.append(float(len(g1)))
                    variances.append(np.var(g1["l"]) / vTot)
                    chiVar.append(np.var(g1["z"]))
                    chiMean.append(np.mean(g1["z"]))
                    chiMin.append(np.min(g1["z"]))
                    chiMax.append(np.max(g1["z"]))
                    counter += 1

        # Add row for non coding with respective statistics
        non_coding = snp_dataset[snp_dataset.gene == "NonCoding"]
        g.append("NonCoding")
        n_snp.append(float(len(non_coding)))
        variances.append(np.var(non_coding["l"]) / vTot)
        chiVar.append(np.var(non_coding["z"]))
        chiMean.append(np.mean(non_coding["z"]))
        chiMin.append(np.min(non_coding["z"]))
        chiMax.append(np.max(non_coding["z"]))

        d["name"] = g
        d["LDvariance"] = variances
        d["StatsVariance"] = chiVar
        d["StatsMean"] = chiMean
        d["StatsMax"] = chiMax
        d["StatsMin"] = chiMin
        d["SNPs"] = n_snp

        df = pd.DataFrame(d)

        output_logger.info(
            " Number of genes (with more than 10 SNPs): " + str(len(df)) + "\n"
        )

        check_g = len(genes_table)
        genes_table = genes_table.merge(
            df, left_index=False, left_on="name", right_on="name")
        if len(genes_table) != check_g:
            logging.warning(
                "something might have gone wrong with the merging process, check gene names compatibility"
            )
        logging.info("generating output genes file")
        genes_table.to_csv(genes_final_file, index=False, sep=",", mode="w")
        logging.info("output genes file generated")
        return g


def analyse(snp_dataset, genes_final_file, folder, output_logger,
            SWEEPS, TUNE, CHAINS, CORES, N_1kG, SUFFIX,
            ):
    """ Bayesian hierarchical regression on the dataset: it takes two dataframes as input,
    one with the SNPs summary stats and one with the genes ( already initialised with "initialise_genes" ) """

    snp_dataset = snp_dataset.reset_index(drop=True)
    n_patients = snp_dataset["sample_size"][0]

    d = dict()

    # to run the regression as a mixed effect model, I need a vector (cat) to assign each SNP to its gene
    idx = 0
    nSNP = len(snp_dataset)
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

    with pm.Model() as model:
        W = pm.InverseGamma("W", alpha=1.0, beta=1.0)
        e = pm.Normal("e", mu=1, sd=1)
        mi = pm.Beta("mi", 1, 1)
        beta = pm.Normal("beta", mu=mi, sd=W, shape=idx)
        diff = pm.Deterministic("diff", subtract(beta, mi))
        herTOT = pm.Deterministic("herTOT", tt.sum(beta * Mg / N_1kG))
        fixed_variable = pm.Normal(
            "fxd",
            mu=(n_patients / N_1kG) * beta[cat] * (snp_dataset["l"]) + e,
            sd=np.sqrt(np.asarray(snp_dataset["l"])),
            observed=snp_dataset["z"],
        )  #

        trace = pm.sample(
            SWEEPS,
            tune=TUNE,
            chains=CHAINS,
            cores=CORES,
            nuts_kwargs=dict(target_accept=0.90),
        )  # njobs=2,,init="advi+adapt_diag"
        # trace = pm.sample(SWEEPS,step=step, tune=TUNE,chains=CHAINS,cores=CORES) #njobs=2,,init="advi+adapt_diag"

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

    su.to_csv(folder + "regression_" + SUFFIX + ".csv", sep=",", mode="w")

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

    # import genes csv file, where the results are going to be written
    with open(genes_final_file, "r+") as f:
        genes_table = pd.read_csv(f, sep=",")

    genes_table = genes_table.merge(
        df, left_index=False, left_on="name", right_on="name")

    k = genes_table.SNPs / float(N_1kG)
    genes_table["h2g"] = genes_table.bg_mean.astype("float") * k

    genes_table = genes_table.sort_values(by=["P", "bg_median"])

    genes_table.to_csv(genes_final_file, index=False, sep=",", mode="w")

    non_coding = genes_table[genes_table.name == "NonCoding"]
    h2g_tot = np.sum(genes_table["h2g"].values) - non_coding["h2g"].values

    output_logger.info(" Non coding heritability : " +
                       str(non_coding["h2g"].values) + "\n")
    output_logger.info(" Coding heritability : " + str(h2g_tot) + "\n")

    logging.info("analysis done")

#####################################################################
############### REGRESSION ANALYSIS #################################
#####################################################################


def regression(
    snp_file: "Data Input, use the SNPs file from dataParse",
    SUFFIX: "suffix for the output files",
    output_folder: "folder where to put the results",
    genes_file: "file with the genes table" ,
    SWEEPS: "number of samples for each chain" = 1000,
    TUNE: "number of burnin samples" = 1000,
    CHAINS: "number of chains of the sampler" = 4,
    CORES: "number of parallel cores to use" = 4,
    N_1kG: "number of SNPs onwhich the LD-score is calculates" = 1290028,
    CHR: "chromosome on which the analysis is run" = "all",
    SNPthr: "threshold for the minimum number of SNPs in a gene" = 10,
    merge_flag: "flag for compatibility of previous version" = False,
    sep: "separator for the input files, use t for tab separated (not \t)" = ",",
    gamma: "use BAGHERA-gamma regression" = False,
):

    """
    Performs bayesian gene-level heritability analysis.
    As input it needs the annotated snps file created with generate-SNPs-file
    and the gene table created by create-files.
    The output folder and the suffix for the output names are used to save the
    output files as follows: output_folder/<filetype>_<suffix>.<fmt>

    From command line one can specify all the parameters for the sampler
    (sweeps, burnin, chains and cores) and the parameters for the SNPs
    and genes filtering.

    We recommend the use of the -m flag for the gene parsing,
    however if the whole pipeline is run, there is no need
    to check for annotations compatibility.

    With the -gamma flag, the BAGHERA-gamma regression is run.
    """

    folder = output_folder
    logging.info('Output files are in %s' % folder)
    # create output folder
    now = datetime.datetime.now()
    logging.info("Folder with the results: %s" % folder)

    # output files
    genes_final_file = folder + "genesResults_" + SUFFIX + ".csv"

    # outputText = folder + "regression_" + SUFFIX + ".log"  # input file with SNPs and LD-scores

    output_logger = setup_logger(
        "output_logger", folder + "regression_" + SUFFIX + ".log")

    output_logger.info(
        "Regression, bayesian gene-level regression analysis results\n "
        + "Current date & time "
        + now.strftime("%Y-%m-%d %H:%M")
    )
    output_logger.info("File: " + snp_file)
    output_logger.info(" Analysis on chr: " + CHR + "\n")
    output_logger.info(" Sweeps: " + str(SWEEPS) +
                       " , Burn: " + str(TUNE) + "\n")
    output_logger.info(
        " SNPs threshold: " + str(SNPthr) + " , Output: " + str(SUFFIX) + "\n"
    )

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
                logging.exception(
                    "Wrong format of the input file, can't recognise the \t separator")

    else:
        with open(snp_file) as f:
            try:
                snp_dataset = pd.read_csv(f, sep=sep)
            except ValueError:
                logging.exception("Wrong format of the input file")

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

    # Creates the genes table with the number of SNPs for each gene and the basic stats values
    gene_list = initialise_genes(
        genes_file, genes_final_file, snp_dataset, SNPthr, output_logger, merge_flag
    )

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

    if gamma:
        regression_gamma.analyse(snp_dataset, genes_final_file, folder, output_logger,
                                 SWEEPS, TUNE, CHAINS, CORES, N_1kG, SUFFIX,
                                 )
    else:
        analyse(snp_dataset, genes_final_file, folder, output_logger,
                SWEEPS, TUNE, CHAINS, CORES, N_1kG, SUFFIX,
                )
