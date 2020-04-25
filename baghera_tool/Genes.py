import sys
import logging
import numpy as np
import theano.tensor as tt
import pymc3 as pm

import pandas as pd
import HTSeq
import os
import baghera_tool.logging as log

def get_stats(group):
    
    return pd.Series({'chrom': int(group['chr'].values[0]),
    'n_snps': len(group),
    'ld_min': group['l'].min(),
    'ld_max': group['l'].max(),
    'ld_var': np.var(group['l'].values),
    'stats_min': group['stats'].min(),
    'stats_max': group['stats'].max(),
    'stats_avg': group['stats'].mean(),
    })

class Genes(object):

    def __init__(self):
        self.filename = None
        self.n_genes = None    
        self.table =pd.DataFrame()
        self.cut_genes = []

    def update_summary(self):
        pass

    def initialise_genes(self, snps_table, snps_thr, non_annotated_name='NonCoding'):

        """ Function to create the results csv file, """

        logging.info("Initialisation of genes table")

        vTot = np.var(snps_table["l"])

        self.table = snps_table[snps_table['gene']!=non_annotated_name].groupby(["gene"]).apply(get_stats)
        self.table.reset_index(inplace=True)
        self.table = self.table[self.table['n_snps']>=snps_thr]
        self.cut_genes = self.table[self.table['n_snps']<snps_thr].gene.values.tolist()
        self.gene_names = self.table['gene'].values.tolist()

        
        non_coding = snps_table[~(snps_table['gene'].isin(self.gene_names))]
        self.table=self.table.append({"gene":non_annotated_name, 
                    "chrom": None, 
                    'n_snps': float(len(non_coding)),
                    "ld_var": np.var(non_coding["l"].values), 
                    'ld_min':np.min(non_coding["l"].values),
                    'ld_max':np.max(non_coding["l"].values),
                    'stats_avg':np.average(non_coding["stats"]),
                    'stats_min':np.min(non_coding["stats"]),
                    'stats_max':np.max(non_coding["stats"])}, ignore_index=True)
        self.table['ld_var']=self.table['ld_var']/vTot
        self.table.rename(columns={'gene':'name'}, inplace=True)
        self.n_genes = len(self.table)
        
    def save_table(self, output_filename, separator=','):
        self.table.to_csv(output_filename, sep=',', index=False)