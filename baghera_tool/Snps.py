import sys
import logging
import numpy as np
import theano.tensor as tt
import pymc3 as pm

import pandas as pd
import HTSeq
import os
import baghera_tool.logging as log

def square(a):
    return a**2

class Snps(object):

    def __init__(self):
        self.table = pd.DataFrame()
        self.n_patients = None
        self.n_snps = None    
        self.n_genes = None    
        self.min_maf = None
        self.max_maf = None
        self.min_ld = None
        self.max_ld = None
        self.avg_stats = None

    def read_table(self, input_snp_filename, separator =','):

        if separator == ',':
            with open(input_snp_filename) as f:
                try:
                    self.table = pd.read_csv(f, sep=',')
                    self.table['chr'] =self.table['chr'].astype(str)                  
                except ValueError:
                    logging.exception("Wrong format of the input file")

        elif 't' == separator:
            with open(input_snp_filename) as f:
                try:
                    self.table = pd.read_csv(f, sep='\t')
                    self.table['chr'] =self.table['chr'].astype(str)
                except ValueError:
                    logging.exception(
                        "Wrong format of the input file, can't recognise the \t separator")

        else:
            with open(input_snp_filename) as f:
                try:
                    self.table = pd.read_csv(f, sep=separator)
                    self.table['chr'] =self.table['chr'].astype(str)
                except ValueError:
                    logging.exception("Wrong format of the input file")

    def generate_stats(self, from_column = 'z', apply = square):
        ''' 
        generates the stats, specify the origin column
        and the function applied to the column values
        '''
        if 'stats' in self.table.columns:
            logging.info('stats column already exists and is going to be replaced')

        self.table['stats']=apply(self.table[from_column])


    def update_summary(self):

        self.n_patients = self.table["sample_size"][0]
        self.n_snps = len(self.table)
        self.min_maf = np.min(self.table["maf"].values)
        self.max_maf = np.max(self.table["maf"].values)
        self.min_ld = np.min(self.table["l"].values)
        self.max_ld = np.max(self.table["l"].values)
        self.avg_stats = np.average(self.table["stats"].values)
        self.n_genes = len(set(self.table['gene'].values.tolist()))

    def apply_filter_table(self, fltr, **args):
        self.table = fltr(self.table, **args)

    def rename_non_annotated(self, name='NonCoding'):
        # Non coding SNPs are assigned to a dummy gene, such that the regression is done on the entire SNPs' set
        self.table['gene'] = self.table['gene'].replace(np.nan, name, regex=True)
    
    def set_non_annotated(self, names, non_annotated_name):
        self.table['gene']=self.table['gene'].replace(to_replace=names, value=non_annotated_name)

########################################################
############# FILTERS ##################################
########################################################

def baghera_filter(table):

    table["l"] = 1 + table["l"] * \
        (table["l"] > 0)  # ld-score [1,+inf)

    logging.info('table length: %d' %len(table))
    # MAF filtering
    table = table[table["maf"] > 0.01]

    logging.info('table length maf > 0.01: %d' %len(table))
    # Filter chromosome 6
    table = table[
        (table.chr != 6) | ((table.position >=
                                   34000000) | (table.position <= 26000000))
    ]

    logging.info('table no HPC chrom 6: %d' %len(table))

    return table

def cut_single_chrom(table, chromosome = 1):
    print(table.head())
    print(chromosome)
    print(str(chromosome))
    table= table[table['chr'] == str(chromosome)]
    print(table.head())
    return table

