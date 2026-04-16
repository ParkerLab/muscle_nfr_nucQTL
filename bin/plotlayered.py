import pandas
import numpy
import math
import matplotlib
import matplotlib.pyplot as plt
import seaborn as sns
import sys
import os
import glob
import time
import pybedtools
import subprocess as sp
#import font_manager
#%load_ext autoreload
#%autoreload 2
#%matplotlib inline
dpi = 150
import matplotlib.font_manager as font_manager
font_path = sys.argv[8] if len(sys.argv) > 8 else ""
if font_path and os.path.exists(font_path):
    prop = font_manager.FontProperties(fname=font_path)
    plt.rcParams['font.family'] = prop.get_name()
else:
    plt.rcParams['font.family'] = 'DejaVu Sans'
matplotlib.rcParams['figure.dpi']= dpi

matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams.update({'font.size': 8})
#import pandas_extra
import plot_utils
#from scipy import stats
#from statsmodels.stats.multitest import multipletests
#from scipy.stats import mannwhitneyu
#import random
import pyBigWig
import pygenometracks.tracks as pygtk
import matplotlib.ticker as ticker
import humanize
#import genome_tracks
from genome_tracks import plot_bigwigs, plot_beds, layer_by_gt
import vcfpy
"""Fix the directory where all figures will be saved"""
rename = {"Mesenchymal_Stem_Cell": "Fibro-adipogenic\nProgenitors",
         "_": " "}
def renamec(x):
    for key in rename.keys():
        x = x.replace(key, rename[key])
    return x
#vcf="../../data/sample_info/example.filtered-vcf.vcf.gz"
#gene_bed= "/lab/work/arushiv/muscle-sn/analyses_hg38/coloc/loci/data/in.sorted.gene_name.bed12.gz"
#cluster_bigwig = "/lab/work/arushiv/muscle-sn/analyses_hg38/coloc/loci/data/bigwigs/{cluster}.bw"
colors = {
    "Type_1" : "#a6cee3",
    "Type_2a": "#1f78b4",
    "Type_2x": "#8dd3c7"
    #"Muscle_Fiber_Mixed": "#33a02c",
    #"Neuromuscular_junction": "#b2df8a",
    #"Mesenchymal_Stem_Cell": "#6a3d9a",
    #"Endothelial":"#fb9a99",
    #"Smooth_Muscle":"#fdbf6f",
    #"Satellite_Cell": "#cab2d6",
    #"T_cell": "#e31a1c",
    #"Macrophage": "#b15928",
    #"Adipocyte": "#999999",
    #"Neuronal": "#ffff33"
}
#bw_files = [cluster_bigwig.format(cluster=key) for key in colors.keys()]
bw_colors = [colors[n] for n in colors.keys()]
bw_names = [renamec(k) for k in colors.keys()]
plot_df = pandas.read_csv(sys.argv[1], sep='\t')
outdir = sys.argv[2]
nominal_path = sys.argv[3] if len(sys.argv) > 3 else '../../results/results_05102024/QTL_opt_results/'
vcf = sys.argv[4] if len(sys.argv) > 4 else "../../data/sample_info/fusion.filtered-vcf.maf0.05-hwe1e6.vcf.gz"
bigwig_sample_map_nfr = sys.argv[5] if len(sys.argv) > 5 else "../../data/bw_maps/{cluster}_nfr.map"
bigwig_sample_map_nuc = sys.argv[6] if len(sys.argv) > 6 else "../../data/bw_maps/{cluster}_nuc.map"
bigwig_sample_map_nuc_smooth = sys.argv[7] if len(sys.argv) > 7 else "../../data/bw_maps/{cluster}_nuc_smooth.map"
os.makedirs(outdir, exist_ok=True)
for i in range(len(plot_df)):
    print(i)
    print(plot_df.iloc[i])

    celltype = plot_df.iloc[i]['Celltype']
    # chr
    chr = plot_df.iloc[i]['QTL_peak'].split('.')[0]
    # nominal_pass result
    path_nfr = nominal_path + celltype +'_nfr/nominal_' + chr + '.txt'
    path_nuc = nominal_path + celltype +'_nuc/nominal_' + chr + '.txt'
    # nfrPeak
    nfrPeak = plot_df.iloc[i]['QTL_peak']
    hit_nfr = plot_df.iloc[i]['snp']
    # read in nominal result
    df_nfr = pandas.read_csv(path_nfr, sep=' ', names = ['phe_id', 'phe_chr', 'phe_from', 'phe_to', 'phe_strd', 'n_var_in_cis', 'dist_phe_var', 'var_id', 
                                      'var_chr', 'var_from', 'var_to', 'nom_pval', 'r_squared', 'slope', 'slope_se', 'best_hit'])
    # select var - nfr
    var_df = df_nfr[(df_nfr['phe_id'] == nfrPeak) & (df_nfr['var_id'] == hit_nfr)]
    #var_df = df_nfr[(df_nfr['phe_id'] == nfrPeak)]
    #var_df=var_df[(var_df['dist_phe_var'] < 100) & (var_df['dist_phe_var'] > -100)]
    # process info
    chr = chr
    var_start = str(var_df.iloc[0]['var_from'])
    var_end = str(var_df.iloc[0]['var_to'])
    hit_nfr = hit_nfr
    phe_name = nfrPeak.replace(".", ":")
    pval = str(var_df.iloc[0]['nom_pval'] )
    slope = str(var_df.iloc[0]['slope'])
    start = str(var_df.iloc[0]['phe_from'])
    end = str(var_df.iloc[0]['phe_to'])
    pos_list = [chr, int(var_start), int(var_end), str(var_df.iloc[0]['var_id']), phe_name, float(pval), float(slope), chr, int(var_start), int(var_end), nfrPeak, celltype]
    # nucpeak info
    nucPeak1 = plot_df.iloc[i]["Nucleosomal_peak"].replace(".", ":")
    # plot
    pos = {'snp_chrom': pos_list[0], 'snp_end': pos_list[1], 'snp': pos_list[3],
        'feature': pos_list[4], 'y_max_signal': 5, 'palette': "volcano_earth", 'slope': pos_list[6]}
    pos_nuc = {'snp_chrom': pos_list[0], 'snp_end': pos_list[1], 'snp': pos_list[3],
        'feature': nucPeak1, 'y_max_signal': 2, 'palette': "volcano_earth", 'slope': pos_list[6]}
    #pos_nuc1 = {'snp_chrom': pos_list[0], 'snp_end': pos_list[1], 'snp': pos_list[3],
    #    'feature': nucPeak2, 'y_max_signal': 3, 'slope': pos_list[6]}
    #pos_nuc2 = {'snp_chrom': pos_list[0], 'snp_end': pos_list[1], 'snp': pos_list[3],
    #    'feature': nucPeak3, 'y_max_signal': 3, 'palette': "volcano"}
    #pos_nuc3 = {'snp_chrom': pos_list[0], 'snp_end': pos_list[1], 'snp': pos_list[3],
    #       'feature': pos_list[9], 'y_max_signal': 3, 'palette': "volcano"}

    pos['x_range'] = (int(pos['snp_end']-5e2), int(pos['snp_end']+5e2))

    nrows = 4
    #ratios = [4,4] + [0.5] 
    ratios=[2, 2, 0.2, 0.2]
    fig, axs = plt.subplots(nrows, 1, figsize=(4, 5), squeeze=False, gridspec_kw={'height_ratios': ratios})

    font = {'family' : plt.rcParams['font.family'],
            'weight' : 'normal',
            'size'   : 12}
    matplotlib.rc('font', **font)

    layer_by_gt(vcf, bigwig_sample_map_nfr.format(cluster = celltype), pos, axs[0,0], showx=False, palette="colorbrewer2", ylabel=celltype+"_nfr")
    #layer_by_gt(vcf, bigwig_sample_map_nuc_old.format(cluster = "Type_2a"), pos, axs[1,0], showx=False, palette="orangesr", ylabel="Type_2a_nuc_nosmoothing")
    layer_by_gt(vcf, bigwig_sample_map_nuc_smooth.format(cluster = celltype), pos, axs[1,0], showx=False, palette="colorbrewer2", ylabel=celltype+"_nuc")
    plot_beds(pos['snp_chrom'], pos['x_range'], features = [f"{pos['feature']}:nfrPeak"], axslist=[axs[2,0]], labelx=False)
    plot_beds(pos['snp_chrom'], pos['x_range'], features = [f"{pos_nuc['feature']}:nucPeak"], axslist=[axs[3,0]], labelx=False)
    #plot_beds(pos['snp_chrom'], pos['x_range'], features = [f"{pos_nuc1['feature']}:nucPeak"], axslist=[axs[4,0]], labelx=False)
    #plot_beds(pos['snp_chrom'], pos['x_range'], features = [f"{pos_nuc2['feature']}:nucPeak"], axslist=[axs[5,0]], labelx=True)
    plot_stub = os.path.join(outdir, f"{celltype}+{hit_nfr}_{i}.png")
    if i==len(plot_df)-1:
        fig.savefig(plot_stub, dpi=300, bbox_inches='tight')
        fig.savefig(os.path.join(outdir, "last.png"), dpi=100, bbox_inches='tight')
    else:
        fig.savefig(plot_stub, dpi=300, bbox_inches='tight')
