#!/usr/bin/env python
# coding: utf-8

# In[6]:


import pandas
import os
import argparse

def getOpts():
    parser = argparse.ArgumentParser(description='Compile LDSC outputs')
    parser.add_argument('--results', required=True, nargs='+', help ="""ldsc outputs """)
    parser.add_argument('--output', required=True, help ="""outputfile """)

    args = parser.parse_args()
    return args


def fixdf(filename):
    d = pandas.read_csv(filename, sep='\t')
    d['trait'] = os.path.basename(filename).replace(".baseline.results", "")
    return d
    
if __name__ == '__main__':
    args = getOpts()

    df = pandas.concat([fixdf(filename) for filename in args.results])
    df.to_csv(args.output, sep='\t', index=False, na_rep="NA")
