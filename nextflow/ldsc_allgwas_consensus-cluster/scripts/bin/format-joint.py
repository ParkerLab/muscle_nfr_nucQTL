#!/usr/bin/env python
# coding: utf-8

# In[6]:


import pandas
import os
import argparse

def getOpts():
    parser = argparse.ArgumentParser(description='Add annot names to the joint model')
    parser.add_argument('--input', required=True, help ="""ldsc joint output """)
    parser.add_argument('--annotfile', required=True, help ="""annotations file """)
    parser.add_argument('--baseline', action="store_true", help ="""annotations file """)
    parser.add_argument('--trait', required=True, help ="""trait name """)
    parser.add_argument('--model-type', required=True, help ="""modeltype """)
    parser.add_argument('--output', required=True, help ="""outputfile """)

    args = parser.parse_args()
    return args

if __name__ == '__main__':
    args = getOpts()
    # all annots and baseline
    f = pandas.read_csv(args.input, sep='\t')
    f['trait'] = args.trait
    f['modeltype'] = args.model_type
    f = f[f['Category'].str.startswith("L2_")]

    d = {}
    al = open(args.annotfile, 'r').readline().split(',')
    a = [os.path.basename(x).replace('.annot.@.ld', '') for x in al] 

    
    for i, x in enumerate(a):
        if "baseline" in args.model_type or args.baseline:
            d[f"L2_{i+1}"] = x
        else:
            d[f"L2_{i}"] = x


    f['cluster'] = f['Category'].map(d)

    f.to_csv(args.output, sep='\t', index=False)
