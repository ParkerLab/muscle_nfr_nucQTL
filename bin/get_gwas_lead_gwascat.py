#!/usr/bin/env python
# coding: utf-8
# adpated from Arushi's code



import numpy
import sys
import os
import glob
import time
import pybedtools
import argparse
from string import Template
import yaml
import argparse
import pandas as pd

#def getOpts():
parser =  argparse.ArgumentParser(description='Start with GWAS summary stats, sort GW significant variants by p value and make xbp regions centered on the lead variant until all significant variants fall within a region.')
parser.add_argument('--summary', required=True, help="""directory with trait1 summary stats and harmonized .h.tsv.gz""")
parser.add_argument('--window', type=int, default=250000, required=True, help="""Flank window on trait1 lead SNP""")
parser.add_argument('--traitname', required=True, help="""trait name""")    
parser.add_argument('--output', required=True, help="""Output file name""")
args = parser.parse_args()
    #return args

def get_signal_regions(gfile, flank, traitname):
    
    # Load the header row to check for the presence of column names
    with open(gfile, 'r') as f:
        first_line = f.readline().strip()

    # Split the header into column names
    header_columns = first_line.split('\t')

    # Determine the column names to use
    column_mapping = {
        'hm_pos': 'hm_pos' if 'hm_pos' in header_columns else 'hg38_pos'
    }

    required_columns = ['hm_chrom', column_mapping['hm_pos'], 'hm_rsid', 'p_value']

    # Read the file with the required columns
    tf = pd.read_csv(
        gfile, 
        sep='\t', 
        usecols=required_columns, 
        dtype={
            'hm_chrom': 'object', 
            column_mapping['hm_pos']: 'Int64'
        }
    )

    # To handle the dynamic column name more uniformly, rename column if needed
    if 'hm_pos' not in tf.columns and 'hg38_pos' in tf.columns:
        tf.rename(columns={'hg38_pos': 'hm_pos'}, inplace=True)

    usecols = ['#snp_chrom', 'snp_start', 'snp_end', 'SNP', 'pval']
    #trait = os.path.basename(gfile).replace(".bed.gz", "")
    #tf = pd.read_csv(gfile, sep='\t', usecols=['hm_chrom', 'hm_pos', 'hm_rsid', 'p_value'], dtype={'hm_chrom': 'object','hm_pos':'Int64'})
    tf['#snp_chrom'] = tf['hm_chrom']
    tf['snp_start'] = tf['hm_pos']
    tf['snp_end'] = tf['hm_pos']+1
    tf['SNP'] = tf['hm_rsid']
    tf['pval'] = tf['p_value']
    tf = tf[tf['pval'] < 5e-8 ]
    tf = tf[['#snp_chrom','snp_start', 'snp_end', 'SNP', 'pval']]
    if tf.empty:
        print("No significant variants")
        return pandas.DataFrame()
    
    else:
        rlist = []
        while not tf.empty:
            tf.sort_values('pval', inplace=True)
            chrom = tf.iloc[0]['#snp_chrom']
            pos = tf.iloc[0]['snp_end']
            snp = tf.iloc[0]['SNP']
            trait1_locus =  f"{traitname}__{snp}__P"
            start = tf.iloc[0]['snp_start'] - flank
            start = 0 if start < 0 else start
            end = tf.iloc[0]['snp_end'] + flank
            nstart = len(tf.index)
            tbed = pybedtools.BedTool.from_dataframe(tf)
            rbed = pybedtools.BedTool(f"{chrom} {start} {end}", from_string=True)
            tf = tbed.intersect(rbed, wa=True, v=True).to_dataframe(names = usecols)
            nend = len(tf.index)
            
            if nend < nstart:
                rlist.append([chrom, pos-1, pos, trait1_locus])
            
        rdf = pd.DataFrame(rlist, columns = ['chrom', 'start', 'end', 'trait1_locus'])
        pybedtools.helpers.cleanup()
        
        return rdf

#if __name__ == '__main__':
    
#args = getOpts()

tdf = get_signal_regions(args.summary, args.window, args.traitname)
tdf.to_csv(f"{args.output}", sep='\t', index=False, header=False)

pybedtools.helpers.cleanup()