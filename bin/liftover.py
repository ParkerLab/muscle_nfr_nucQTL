#!/usr/bin/env python
# coding: utf-8
# adpated from Arushi's code

import pandas as pd
import numpy
import sys
import os
import glob
import time
import argparse
from string import Template
import yaml
import argparse
import gzip
from pyliftover import LiftOver
import pysam

#def getOpts():
parser =  argparse.ArgumentParser(description='Start with GWAS summary stats, sort GW significant variants by p value and make xbp regions centered on the lead variant until all significant variants fall within a region.')
parser.add_argument('--summary', required=True, help="""directory with trait1 summary stats and harmonized .tsv.bgz""")
parser.add_argument('--vcf', required=True, help="""Which vcf file to rely on for rsid""")
parser.add_argument('--output', required=True, help="""Output file name""")
args = parser.parse_args()

# read in raw sumstats
file_path = args.summary 
with gzip.open(file_path, 'rt') as f:
    hg19_df = pd.read_csv(f, sep='\t')

print("Done reading the raw sumstats!")
# Initialize the LiftOver object
lo = LiftOver('hg19', 'hg38')  # Converts from hg19 to hg38

# Example DataFrame with 'chrom' and 'pos' columns
chr_num_l = hg19_df['chr'].tolist()
common_string = "chr"
chr_l = [f"{common_string}{elem}" for elem in chr_num_l]
data = {
    'chrom': chr_l,
    'pos': hg19_df['pos'].tolist(),
    'chrom_num': chr_num_l
}
df = pd.DataFrame(data)

# Function to apply liftover to each row
def liftover_row(row):
    chrom, pos = row['chrom'], row['pos']
    result = lo.convert_coordinate(chrom, pos)
    if result:
        # Extract the first conversion result (there could be more than one)
        new_chrom, new_pos, _, _ = result[0]
        return pd.Series([new_chrom, new_pos])
    else:
        # Return None or NaN if conversion fails
        return pd.Series([None, None])

# Apply liftover to each row and store the results in new columns
df[['hg38_chrom', 'hg38_pos']] = df.apply(liftover_row, axis=1)
print("Done lifting over!")
# merge with original
merged_df = hg19_df.merge(df, how='left', left_on=['chr', 'pos'], right_on = ['chrom_num', 'pos'])

# Path to the indexed VCF file with rsIDs
vcf_file = args.vcf

# Open the VCF file using pysam
vcf = pysam.VariantFile(vcf_file)

def get_rsid_from_vcf(chrom, pos):
    # Fetch records for the given chromosome and position
    for record in vcf.fetch(chrom, pos-1, pos):
        return record.id  # Assuming first match is desired
    return None

# Example usage
#rsid = get_rsid_from_vcf('chr1', 800909)
#print(rsid)
# tidy hg38_chrom col
# drop na
print(len(merged_df))
# Drop rows with NaN in the 'chromosome' column
merged_df = merged_df.dropna(subset=['hg38_chrom'])
print(len(merged_df))
# filter for clean chr pattern
pattern = r'^chr([1-9]|1[0-9]|2[0-2])$'
merged_df = merged_df[merged_df['hg38_chrom'].str.match(pattern)]
# get rsid
merged_df['hm_rsid'] = merged_df.apply(lambda row: get_rsid_from_vcf(row['hg38_chrom'], row['hg38_pos']), axis=1)
print("Done adding rsid from vcf!")
# tidy
merged_df['hm_chrom'] = merged_df['hg38_chrom'].str.replace('chr', '')
merged_df = merged_df.rename(columns={"hm_pos": "hg38_pos"})
merged_df['p_value'] = 10**((-1)*merged_df['neglog10_pval_EUR'])

merged_df.to_csv(f"{args.output}", sep='\t', index=False, header=True)
print("Saved!")