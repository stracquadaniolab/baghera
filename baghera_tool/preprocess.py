import pandas as pd
import HTSeq
import os
import glob as g
import inspect
import sys
import logging
import datetime


def import_position_ukbb(fileInput):
    """
    Imports UK Biobank files.
    Three fields are mandatory in the assoc.tsv file:
    "variant", "nCompleteSamples", "tstat"
    and the merging is done on the position.

    This function directly splits the variant column
    if those are not found an exception is raised
    """

    with open(fileInput) as f:
        try:
            SNP = pd.read_csv(
                f, usecols=["variant", "n_complete_samples", "tstat"], sep="\t"
            )
        except ValueError:
            logging.exception(
                "Wrong format of the input file, check the file has ['rsid','n_complete_samples','tstat'] columns"
            )
    SNP = SNP.dropna(subset=["variant"])

    SNP["chr"], SNP["position"], _, _ = SNP["variant"].str.split(":").str
    # SNP['chr'],SNP['position'],_
    del SNP["variant"]

    print(SNP.head())

    SNP = SNP.rename(columns={"tstat": "z", "n_complete_samples": "sample_size",})
    return SNP


def import_position(fileInput):
    """
    Imports stats files.
    Mandatory fields are "chrom", "pos", "nCompleteSamples", "tstat"
    and the merging is done on the position

    if those are not found an exception is raised
    """

    with open(fileInput) as f:
        try:
            SNP = pd.read_csv(
                f, usecols=["chrom", "pos", "nCompleteSamples", "tstat"], sep="\t"
            )
        except ValueError:
            logging.exception(
                "Wrong format of the input file, check the file has ['rsid','nCompleteSamples','tstat'] columns"
            )

    SNP = SNP.rename(
        columns={
            "chrom": "chr",
            "pos": "position",
            "tstat": "z",
            "nCompleteSamples": "sample_size",
        }
    )
    return SNP


def import_ukbb(fileInput):
    """
    Imports UK Biobank files.
    Three fields are mandatory in the assoc.tsv file:
    "rsid", "nCompleteSamples", "tstat"
    and the merging is done on the rs_id

    if those are not found an exception is raised
    """

    with open(fileInput) as f:
        try:
            SNP = pd.read_csv(
                f, usecols=["rsid", "nCompleteSamples", "tstat"], sep="\t"
            )
        except ValueError:
            logging.exception(
                "Wrong format of the input file, check the file has ['rsid','nCompleteSamples','tstat'] columns"
            )

    SNP = SNP.rename(
        columns={"rsid": "rs_id", "tstat": "z", "nCompleteSamples": "sample_size",}
    )
    return SNP


def import_ldsc(fileInput):
    """
    Imports sumstats files.
    The file must have the fields "SNP", "Z", "N".
    and the merging is done on the rs_id

    If those aren't in the file an exceprion is raised.
    """

    with open(fileInput) as f:
        try:
            SNP = pd.read_csv(f, usecols=["SNP", "Z", "N"], sep="\t")
        except ValueError:
            logging.exception(
                "Wrong format of the input file, check the file has SNP, Z, N columns"
            )

    SNP = SNP.rename(columns={"SNP": "rs_id", "Z": "z", "N": "sample_size"})

    return SNP


# create a GenomicArray object with all the features in a GTF file
def load_gtf_annotation(
    filename, chrom_list, feature="gene", transcript="protein_coding"
):
    # GenomicArray will store feature objects because it will be needed
    # when computing distance from boundaries
    gene_annotation = HTSeq.GenomicArrayOfSets(chrom_list, stranded=False)

    # loop through all the features and retain just the one of type "feature"
    # by default it retains features of type "gene"

    counter = 0
    for feat in HTSeq.GFF_Reader(filename):
        if feat.iv.chrom in chrom_list:
            if (feat.type == feature) & ("gene_type" in feat.attr.keys()):
                if feat.attr["gene_type"] == transcript:
                    gene_annotation[feat.iv] += feat.attr["gene_name"]

                    counter += 1
    logging.info("Number of genes in annotation file: " + str(counter))

    # return the annotation
    return gene_annotation


# assign SNPs to genes, based on nearest distance criteria
# the overlap can be extended by an x-bp using the window parameter
# the function takes in input a GenomicArray, a GenomicLocation and an optional
# window size (10kb by default)


def annotate_snp(annotation, snp, window):

    # create a query SNP to lookup features +/- window from it
    query_snp = HTSeq.GenomicInterval(
        snp.chrom, snp.pos - window * (snp.pos - window > 0), snp.pos + window,
    )

    # keep track of the closest gene and distance
    closest_gene = None
    closest_distance = float("inf")

    # query the genomic array for features overlapping the SNP
    region = annotation[query_snp].steps()

    # for each gene overlapping the SNP, calculate the distance
    for iv, gene in region:
        curr_dist = float("inf")

        if gene != set():
            # the SNP is within the genes

            if snp.pos >= iv.start and snp.pos <= iv.end:
                curr_dist = 0
            # otherwise compute the minimum distance from the gene boundaries
            else:
                curr_dist = min(abs(snp.pos - iv.start), abs(snp.pos - iv.end))
            # if current_dist is closer than the best one, update
            if curr_dist < closest_distance:
                closest_distance = curr_dist
                closest_gene = gene

    # return the gene_id
    if closest_gene:
        return closest_gene
    else:
        return None


def cluster_genes(genes, chrom_list):
    """ cleans overlapping regions, all partially or completely
        overlapping genes are clustered into a single gene
    """

    genes2 = HTSeq.GenomicArrayOfSets(chrom_list, stranded=False)
    region = genes.steps()
    last = set()
    num = 0
    FLAG = False
    iv0 = HTSeq.GenomicInterval("chr1", 0, 1)

    for iv, gene in region:

        if len(gene) == 0:
            if FLAG == False:
                last = set([])
                num = 0
            else:
                genes2[iv0] = last

                last = set()
                num = 0
        else:

            FLAG = True
            last = set.union(last, gene)
            num += 1
            if num > 1:
                iv0.extend_to_include(iv)
            else:
                iv0 = iv

    return genes2


####################################################################
############# CREATE ANNOTATED LD ##################################
####################################################################


def create_files(
    ldscore_folder: "folder with LD score as in 1kG" = "data/eur_w_ld_chr/",
    annotation_file: "gtf file for the annotation" = "data/gencode.v31lift37.basic.annotation.gtf",
    snps_output: "annotated ld-score snps file" = "data/ld_annotated_gencode_v31.csv",
    genes_output: "genes table with clustered genes" = "data/genes_gencode_v31.csv",
    chrom_list: "list of chromosomes used by HTSeq, if None chr<no> is used " = None,
):
    """
    Builds the annotated set of SNPs required for downstream analysis.
    It requires a genome annotation in GTF format (preferably Gencode)
    and the LD score folders with $CHR.l2.ldscore files.

    :param ldscore_folder: snps table filename
    :param annotation_file: separator of the table
    :param snps_output: ldscore snps annotated to genes
    :param genes_output: genes table, the same used to annotate the snps
    :param chrom_list: list of chromosomes used by HTSeq, if None chr<no> is used
    """

    if chrom_list == None:
        chrom_list = [
            "chr1",
            "chr2",
            "chr3",
            "chr4",
            "chr5",
            "chr6",
            "chr7",
            "chr8",
            "chr9",
            "chr10",
            "chr11",
            "chr12",
            "chr13",
            "chr14",
            "chr15",
            "chr16",
            "chr17",
            "chr18",
            "chr19",
            "chr20",
            "chr21",
            "chr22",
        ]

    folder = g.glob(ldscore_folder + "*.l2.ldscore")
    if len(folder) > 0:
        logging.info("LD score folders found")
    else:
        logging.error("LD score folder not found")
        sys.exit()

    file_genes = annotation_file

    LD = []
    for f in folder:
        with open(f) as file:
            LD.append(
                pd.read_csv(
                    file, usecols=["CHR", "SNP", "BP", "CM", "MAF", "L2"], sep="\t",
                )
            )

    ld = pd.concat(LD)
    ld = ld.rename(columns={"CHR": "chr", "SNP": "rs_id", "BP": "position", "L2": "l",})

    genes = load_gtf_annotation(file_genes, chrom_list)  # loads annotation
    logging.info("Annotation file loaded")

    genes = cluster_genes(genes, chrom_list)  # non overlapping set of genes
    logging.info("Genes clustering successfull")

    # creating table of genes
    genes_table = pd.DataFrame(columns=["chrom", "start", "stop", "name"])
    for iv, gene in genes.steps():
        if len(gene) > 0:
            df2 = pd.DataFrame(
                [
                    [
                        str(iv.chrom)[3:],
                        int(iv.start),
                        int(iv.end),
                        "; ".join(str(s) for s in gene),
                    ]
                ],
                columns=["chrom", "start", "stop", "name"],
            )
            genes_table = genes_table.append(df2, ignore_index=True)

    # Saving genes table
    genes_table.to_csv(
        genes_output,
        # directory + "genesTable.csv",
        sep=",",
        index=False,
        float_format="%.3f",
    )

    logging.info("genes table created")

    # merging SNPs and annotations
    last_col = []
    for index, v in ld.iterrows():

        s = HTSeq.GenomicPosition("chr" + str(v.chr), int(v.position))
        count = 0
        gList = genes[s]

        if len(gList) > 0:
            last_col.append("; ".join(str(s) for s in gList))
        else:
            gAnn = annotate_snp(genes, s, window=50000)
            if gAnn != None:
                last_col.append("; ".join(str(s) for s in gAnn))
            else:
                last_col.append("")
        count += 1

    ld["gene"] = last_col

    header_list = ["chr", "rs_id", "position", "cm", "maf", "l", "gene"]

    # write data to the output file
    ld.to_csv(
        snps_output,
        # directory + "annotatedLD.csv",
        sep=",",
        header=header_list,
        index=False,
        float_format="%.3f",
    )

    logging.info("SNPs annotated")


####################################################################
############# CREATE FILE FOR REGRESSION ANALYSIS ##################
####################################################################


def generate_snp_file(
    stats_input_file: "SNPs file" = "data/c50_breast_snps.csv",
    input_type: "ldsc or ukbb, position or position_ukbb" = "position_ukbb",
    annotated_ld_file: "previously generated annotated LD file" = "data/ld_annotated_gencode_v31.csv",
    output_file: "output filename (.csv)" = "data/c50_snps.csv",
):
    """
    Annotate input summary statistics file using the the annotation
        generated by BAGHERA create-files.



    :param snps_input_file: snps table with summary stats
    :param input_type: type of the input, use one between ldsc, ukbb, position, position_ukbb
    :param annotated_ld_file: ldscore snps annotated to genes can be generated with create-files
    :param output_file: output snps file with stats ldscore and gene annotation
    """

    if input_type == "ldsc":
        SNP = import_ldsc(stats_input_file)
    elif input_type == "ukbb":
        SNP = import_ukbb(stats_input_file)
    elif input_type == "position":
        SNP = import_position(stats_input_file)
    elif input_type == "position_ukbb":
        SNP = import_position_ukbb(stats_input_file)
    else:
        logging.error("Unknown input-type parameter")

    logging.info("The input file has %d SNPs" % len(SNP))

    try:
        with open(annotated_ld_file) as f:
            annotatedLD = pd.read_csv(f, sep=",")
    except ValueError:
        logging.exception("Missing or wrong file with annotated variants")

    logging.info("There are %d annotate variants" % len(annotatedLD))

    annotatedLD["chr"] = annotatedLD["chr"].astype(str)
    annotatedLD["position"] = annotatedLD["position"].astype(str)
    SNP["chr"] = SNP["chr"].astype(str)
    SNP["position"] = SNP["position"].astype(str)

    if (input_type == "position") | (input_type == "position_ukbb"):
        SNP = SNP.merge(
            annotatedLD,
            left_index=True,
            left_on=["chr", "position"],
            right_on=["chr", "position"],
        )
        logging.info(
            "We were able to merge %d SNPs between the input file and the annotations"
            % len(SNP)
        )
    else:
        SNP = SNP.merge(annotatedLD, left_index=True, left_on="rs_id", right_on="rs_id")
        logging.info(
            "We were able to merge %d SNPs between the input file and the annotations"
            % len(SNP)
        )

    header_list = list(SNP)

    # write data to the output file
    SNP.to_csv(output_file, index=False, header=header_list, sep=",", mode="w")
    logging.info("File created")
