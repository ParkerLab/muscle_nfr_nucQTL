#!/usr/bin/env python
# coding: utf-8

# In[6]:


import pandas
import os
import argparse

def getOpts():
    parser = argparse.ArgumentParser(description='Compile LDSC outputs')
    parser.add_argument('--annotlists', required=True, nargs='+', help ="""file with one column of annotations to include in each joint model set""")
    parser.add_argument('--prefix', required=True, help ="""outputfile prefix """)
    parser.add_argument('--dir', required=True, help ="""directory with ld files """)
    parser.add_argument('--baseline', help ="""baseline path """)

    args = parser.parse_args()
    return args

    
if __name__ == '__main__':
    args = getOpts()

    for annotfile in args.annotlists:
        annotlist = pandas.read_csv(annotfile, sep='\t', header=None, squeeze=True).tolist()
        annotset = os.path.basename(annotfile).replace(".txt", "")

        filelist = []
        with open(f"{args.prefix}.{annotset}.txt", 'w+') as out:
            if "baseline" in annotlist:
                filelist.append(args.baseline)
                annotlist.remove("baseline")
            filelist += [f"{args.dir}/{annot}.annot.@.ld" for annot in annotlist]
            
            out.write(','.join(filelist))
