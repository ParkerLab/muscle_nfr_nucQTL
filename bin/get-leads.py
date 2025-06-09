#!/usr/bin/env python
# coding: utf-8

# In[28]:


import pandas
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

def getOpts():
    parser =  argparse.ArgumentParser(description='Start with GWAS summary stats, sort GW significant variants by p value and make xbp regions centered on the lead variant until all significant variants fall within a region.')
    parser.add_argument('--summary', required=True, help="""directory with trait1 summary stats indexed .bed.gz files""")
    parser.add_argument('--window', type=int, default=250000, required=True, help="""Flank window on trait1 lead SNP""")
    parser.add_argument('--traitname', required=True, help="""trait name""")    
    parser.add_argument('--output', required=True, help="""Output file name""")
    args = parser.parse_args()
    return args

    
def get_signal_regions(gfile, flank, traitname):
    usecols = ['#snp_chrom', 'snp_start', 'snp_end', 'SNP', 'pval']
    trait = os.path.basename(gfile).replace(".bed.gz", "")
    tf = pandas.read_csv(gfile, sep='\t', usecols=usecols, dtype={'snp_start':int, 'snp_end':int})
    tf = tf[tf['pval'] < 5e-8 ]
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
            
        rdf = pandas.DataFrame(rlist, columns = ['chrom', 'start', 'end', 'trait1_locus'])
        pybedtools.helpers.cleanup()
        
        return rdf
    
if __name__ == '__main__':
    
    args = getOpts()
    
    tdf = get_signal_regions(args.summary, args.window, args.traitname)
    tdf.to_csv(f"{args.output}", sep='\t', index=False, header=False)
    
    pybedtools.helpers.cleanup()
