#! /usr/bin/env python3

import sys
import re
from pathlib import Path
from collections import OrderedDict
import warnings
import os
import glob
import humanize
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import numpy as np
import pyBigWig
import vcfpy
from tqdm import tqdm
from cycler import cycler
import pygenometracks.tracks as pygtk
import pandas
import argparse
import tabix
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import pybedtools
from adjustText import adjust_text
import tempfile
mpl.rcParams['pdf.fonttype'] = 42
mpl.rcParams.update({'font.size': 12})

palettes = {
    "vcustom": ["#440154", "#AED4E5", "#FDE725"],  # viridis
    "viridis": ["#440154", "#21908C", "#FDE725"],  # viridis
    "custom": ["#4F443F", "#AED4E5", "#27357E"],  # custom
    "inferno": ["#991B1E", "#FCA50A", "#27357E"],  # inferno (viridis)
    "magma": ["#000004", "#b73779", "#fcfdbf"], # magma (viridis)
    "oranges": ["#F7941D", "#8C6B46", "#0C0C0C"],
    "orangesr": ["#0C0C0C", "#8C6B46", "#F7941D"],
    "volcano": ["#111111","#FCB930","#D82909"],
    "volcano_earth": ["#D82909","#FCB930","#111111"],
    "colorbrewer1": ["#a6cee3", "#1f78b4", "#b2df8a"],
    "colorbrewer2": ["#66c2a5", "#fc8d62", "#8da0cb"]
    
}
mpl.rcParams["savefig.pad_inches"] = 0
mpl.rcParams["legend.framealpha"] = 0
mpl.rcParams['pdf.fonttype'] = 42
mpl.rc("font", size=8)
mpl.rc("axes", titlesize=8)
mpl.rc("legend", fontsize=7)
mpl.rc('ytick',labelsize=8)


@ticker.FuncFormatter
def x_tick_formatter(x, pos):
    # return humanize.filesize.naturalsize(x, format="%.3f").replace("B", "b") # 
    return humanize.filesize.naturalsize(x, format="%.3f").split(" ")[0]

class BigwigObj:
    def __init__(self, url):
        myurl = url
        self.bw = pyBigWig.open(myurl)

    def get_scores(self, pos):
        return self.bw.values(*pos)


def get_value_from_pos(bw_object, pos, xval, yval):
    """Fetch values from BigWig object. See class defintion above."""
    scores = None
    try:
        scores = bw_object.get_scores(
            [pos['snp_chrom'], xval, yval]
        )
    except Exception as e:
        print("error: was trying to get scores: {0}".format(e))

    if scores is not None:
        if len(scores) != 0:
            return [pos, np.mean(scores), scores]

    return None


def remove_spines(axs, poslist=['top','left','bottom','right']):
    for spinepos in poslist:
        axs.spines[spinepos].set_visible(False)

def get_mid(x):
    s = x.split(":")
    st = int(s[1])
    en = int(s[2])
    mid = st + int((en - st)/2)
    return mid

def checkfile(path):
    if not Path(path).exists:
        print(f"{path} not found")
        sys.exit(1)

def plot_bigwigs(chrom, x_range, bwfiles, names, axslist=None, colors=None, labely=True, labelx=True, show_data_range=False, max_val=None, optimize_max_y=False):
    """plot a set of bigwig files in a region"""
    if colors is not None:
        assert len(colors) == len(bwfiles)
    else:
        colors = ["#F7941D" for f in bwfiles]
        
    assert len(bwfiles) == len(names)
    if axslist is None:
        fig, axslist = plt.subplots(len(bwfiles), 1)
        if len(bwfiles) == 1:
            #axslist is just 1 ax. make it a list to iterate 
            axslist = [axslist]            
    maxy = 0
    for i, (axs, bw, name, color) in enumerate(zip(axslist, bwfiles, names, colors)):
        track_config = dict(
            file = bw,
            height = 0.1,
            title = name,
            color = color,
            nans_to_zeros = True,
            max_value = max_val,
            show_data_range = show_data_range
        )
        tk = pygtk.BigWigTrack(track_config)

        tk.plot(axs, chrom, *x_range)
        axs.set_xlim(*x_range)
        if i==len(bwfiles)-1 and labelx:
            axs.tick_params(labelsize=7)
            axs.xaxis.set_major_locator(ticker.LinearLocator(4))
            axs.xaxis.set_major_formatter(x_tick_formatter)
            axs.set_xlabel(f"{chrom} (Mb)")
        else:
            axs.get_xaxis().set_visible(False)

        remove_spines(axs)
        axs.get_yaxis().set_ticks([])#set_visible(False)
        if labely:
            axs.set_ylabel(name, rotation=0, size=8, horizontalalignment='right', verticalalignment="center")
        else:
            axs.set_ylabel("")
        if optimize_max_y:
            maxy = max(maxy, axs.get_ylim()[1])
    if optimize_max_y:
        for axs in axslist:
            axs.set_aspect('auto')
            axs.set_ylim(0, round(maxy,0))

    return axslist

def plot_gene_track(chrom, x_range, gene_bed_file, ax=None, labels=True, labels_in_margin=False, all_labels_inside=True, arrowhead_included=False,
                    max_labels=0, gene_rows=1, style="flybase", merge_transcripts=True, fontsize=8, **kwargs):
    """plot gene track using a gene bed file"""

    track_config = dict(
        file=gene_bed_file,
        file_type="bed",
        labels=labels,
        labels_in_margin=labels_in_margin,
        arrowhead_included=arrowhead_included,
        max_labels=max_labels,
        gene_rows=gene_rows,
        style=style,
        merge_transcripts=merge_transcripts,
        fontsize=fontsize,
        all_labels_inside=all_labels_inside,
        section_name='test'
    )

    track_config.update(kwargs)
        

    if ax is None:
        fig, ax = plt.subplots(1, 1)
        
    tk = pygtk.BedTrack(track_config)
    try:
        tk.plot(ax,  chrom, *x_range)
    except:
        ax.axis('off')
    ax.set_xlim(*x_range)
    remove_spines(ax)
    ax.get_yaxis().set_visible(False)
    ax.get_xaxis().set_visible(False)

    

def plot_beds(chrom, x_range, features=[], beds=[], axslist=None, color="black", labelx=True):
    """plot a set of regions on a chromosome"""
    if features:
        nrows = 1
        df = pandas.DataFrame({'n': features})
        tmp = tempfile.NamedTemporaryFile(delete=False)
        with open(tmp.name, 'w+') as f:
            df['n'].str.split(":", expand=True).to_csv(f, sep='\t', index=False, header=False)
        beds = [tmp.name]
    elif beds:
        nrows = len(beds)
    else:
        print("Provide list of features or beds to plot")
        sys.exit(1)
        
    if axslist is None:
        fig, axslist = plt.subplots(nrows, 1)      
        if nrows == 1:
            #axslist is just 1 ax. make it a list to iterate 
            axslist = [axslist]            

    for i, (axs, bed) in enumerate(zip(axslist, beds)):
        track_config = dict(
            file=bed,
            file_type="bed",
            labels_in_margin=True,
            merge_transcripts=True,
            fontsize=6,
            arrowhead_included=True,
            color=color,
            # arrow_interval=5,
            # max_labels=5,
            # gene_rows=1,
            style="UCSC",
            section_name='test'
        )
        tk = pygtk.BedTrack(track_config)

        tk.plot(axs, chrom, *x_range)
        axs.set_xlim(*x_range)
        if i==len(beds)-1 and labelx:
            axs.tick_params(labelsize=7)
            axs.xaxis.set_major_locator(ticker.LinearLocator(4))
            axs.xaxis.set_major_formatter(x_tick_formatter)
            axs.set_xlabel(f"{chrom} (Mb)")
        else:
            axs.get_xaxis().set_visible(False)

        remove_spines(axs)
        axs.get_yaxis().set_visible(False)

    if features:
        for f in beds:
            os.remove(f)
    return axslist


def layer_by_gt(vcfpath, bwmap, pos, axs, palette="custom", ylabel=None,
               debug=False, showx=False, progress=False):
    """Given a set for bigwig files, aggregate into layers by genotype.
    Plotted on the supplied axs.
    pos is a dictionary containing plotting info. snp_chrom, snp_end, snp, slope. optionally y_max_signal for the y axis max value"""
    
    vcf = Path(vcfpath).absolute()
    bigwig_sample_map = Path(bwmap).absolute()
    checkfile(vcf)
    checkfile(bwmap)

    required_keys = ["snp_chrom","snp_end","snp","x_range"]
    for k in required_keys:
        if k not in pos.keys():
            print(f"key {k} not found")
            sys.exit(1)
            
    if palette not in palettes.keys():
        print("info: palette not found. using custom.")
        palette = 'custom'
    # mpl.rcParams["axes.prop_cycle"] = cycler(color=palettes[palette])
            
    x_range = pos['x_range']
    xmin = x_range[0]
    xmax = x_range[1]
    print(f"plotting {(xmax-xmin)/1000} kb region:  {x_range}")

    b = pandas.read_csv(bigwig_sample_map, sep='\t', header=None, names=['path', 'sample'], dtype={'sample':str})
    bigwigfiles = b['path'].tolist()
    sample_ids = b['sample'].tolist()

    for k in bigwigfiles:
        if not Path(k).exists():
            print(f"error: {k} not found.", file=sys.stderr)
            status = 1
    
    if not debug:
        warnings.filterwarnings("ignore")

    try:
        vcf_file = vcfpy.Reader.from_path(vcf, parsed_samples=sample_ids)
    except AssertionError:
        vcf_file = vcfpy.Reader.from_path(vcf)
        samples_vcf = vcf_file.header.samples.names
        not_found = [x for x in sample_ids if x not in samples_vcf]
        print(f"error: not all samples present in the vcf file.")
        print(f"error: following not found - {not_found}")
        sys.exit(1)

    # NOTE: BED is 0-based
    data = {}  # "sample": { "pos1": {...}, "pos2": {...} }

    if progress:
        bar = tqdm(total=len(sample_ids) + 1)
    for _, sample in enumerate(bigwigfiles):
        sample_id = b[b['path']==sample].iloc[0]['sample']
        if debug:
            print(f"info: processing sample {sample_id}")
        bw = BigwigObj(str(sample))

        slot = {}  # "Position": { "gt": .., }
        # for i, pos in positions.iterrows():
        key = f"{pos['snp_chrom']}:{pos['snp_end']}"

        # get bigwig values
        bw_values = get_value_from_pos(
            bw, pos, *x_range
        )

        # get genotype status
        if debug:
            print(f"info: querying {pos}..")
            # NOTE: 0-based query (but VCF is 1-based)
            # TODO: Check consistency of chr/no-chr for the chromosome field
        query = vcf_file.fetch(pos['snp_chrom'], pos['snp_end']-1, pos['snp_end'])
        for record in query:
            # if not record.is_snv():
            #     print("info: not a SNV!")
            #     continue
            info = record.INFO
            sample_call = record.call_for_sample.get(sample_id)
            alt = record.ALT
            if len(alt) > 1:
                print("warning: skipping multi-allelic variant!")
                continue

            if debug:
                print(sample_call)
                print(
                    f"info: found: {record.CHROM}:{record.POS},"
                    f" ref/alt: {record.REF}/{record.ALT[0].value},"
                    f" maf: {info.get('MAF')},"
                    f" af: {info.get('AF')}"
                )
            slot[key] = {
                "gt": "/".join(sample_call.gt_bases),
                "ref": record.REF,
                "alt": alt[0].value,
                "gt_type": sample_call.gt_type,
                "het": sample_call.is_het,
                "pos": bw_values[0],
                "signal_mean": bw_values[1],
                "signal_values": bw_values[2],
            }

        data[sample_id] = slot
        if progress:
            bar.update(1)

    #########
    #### Start plot axs
    #########
            
    key = f"{pos['snp_chrom']}:{pos['snp_end']}"
    genotype_signal = {}

    # Collect genotype signals from all samples
    for _, v in data.items():
        if key in v:
            vals = v[key]
            if vals["gt_type"] == 0:
                gt_key = vals["ref"] + vals["ref"]
            elif vals["gt_type"] == 1:
                gt_key = vals["ref"] + vals["alt"]
            elif vals["gt_type"] == 2:
                gt_key = vals["alt"] + vals["alt"]
            else:
                print("error: weird error -- check script plot logic!")
                sys.exit(1)

            if gt_key not in genotype_signal:
                genotype_signal[gt_key] = []
            genotype_signal[gt_key].append(vals["signal_values"])
            
    # Sort so that lowest signal track is in the front of the plot
    # actually, just use the caQTL slope if available because averaging across entire window sometimes isn't correct
    
    if 'slope' in pos:
        order = [f"{vals['alt']}{vals['alt']}", f"{vals['ref']}{vals['alt']}", f"{vals['ref']}{vals['ref']}"]
        o = [gtype for gtype in order if gtype in genotype_signal.keys()]
        
        if pos['slope'] < 0:
            # alt-alt is lower signal and is at the front
            o.reverse()
        genotype_signal = {k: genotype_signal[k] for k in o} 

    else:
        # Sort so that lowest signal track is in the front of the plot
        # For this, we will use the sum of average signal across the region.
        print(len(genotype_signal))
        genotype_signal = OrderedDict(
            sorted(
                genotype_signal.items(), key=lambda item: np.sum(np.sum(item[1], axis=0)/len(item[1])),
                reverse=True
            )
        )
        
        
    for fillcolor, (k, v) in zip(palettes[palette], genotype_signal.items()):
        # Add signal across samples for each genotype and average
        genotype_signal[k] = np.sum(v, axis=0) / len(v)
        
        # Make an area plot with genotype as legend
        axs.fill_between(
            range(*x_range),
            genotype_signal[k],
            label=f"{k} ({len(v)})",
            **{"edgecolor": "black", "linewidth": 0.1, "color": fillcolor}
        )
        
    # Prettify 
    axs.legend(loc=0, frameon=False, fontsize=7)
    remove_spines(axs, ['top','right','bottom'])
    axs.set_xlim(*x_range)
    axs.set_ylim(ymin=0)
    if 'y_max_signal' in pos.keys():
        axs.set_ylim(ymax=pos['y_max_signal'])

    axs.yaxis.set_major_locator(ticker.LinearLocator(2))
    yticks = axs.get_yticks()
    axs.set_yticks([yticks[-1]])

    if not ylabel:
        axs.set_ylabel("")
    else:
        axs.set_ylabel(ylabel, rotation=0, size=8, horizontalalignment='right', verticalalignment="center")

    # axs.set_ylabel(ylabel, rotation=0, size=8, horizontalalignment='right', verticalalignment="center", x=1.0)
    # axs.xaxis.set_major_locator(ticker.LinearLocator(4))
    # axs.xaxis.set_major_formatter(x_tick_formatter)
    if pos['snp_end'] >= xmin and pos['snp_end'] <= xmax: 
        axs.axvline(x=pos['snp_end'], color='k', linestyle='--')
        # axs.text(pos['snp_end'] - 100, 0, pos['snp'], rotation=90, **{'fontsize': 7})
    
    leg = axs.legend(handlelength=0, handleheight=0, title=pos['snp'], title_fontsize=6)
    for handle, text in zip(leg.legendHandles, leg.get_texts()):
        text.set_color(handle.get_facecolor())
        handle.set_visible(False)    
    
    if not showx:
        axs.get_xaxis().set_visible(False)
    else:
        axs.tick_params(labelsize=7)
        axs.xaxis.set_major_locator(ticker.LinearLocator(2))
        axs.xaxis.set_major_formatter(x_tick_formatter)
        axs.set_xlabel(f"{pos['snp_chrom']} (Mb)")

    if progress:
        bar.update(1)
        bar.close()
    return axs
