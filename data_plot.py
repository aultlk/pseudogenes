"""
   This script is for generating test plots on all merged data from data_process.py output.
   Finalized plots will be incorporated into main.py program.  

"""

from pathlib import Path
from scipy import stats
from itertools import product, combinations
from matplotlib.offsetbox import AnchoredText

import argparse
import pandas as pd

import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib.ticker as tkr
import statannot

ROOT_DIR = Path(__file__).resolve().parent.parent
APPENDED_META = f'{ROOT_DIR}/appended_meta.csv'
LOCI = ('IGH', 'IGK', 'IGL')
GENE_REARRS = ('productive', 'passenger', 'pseudogene')
LOCUS_BINS = ('productive', 'passenger', 'pseudogene')

def parse_args():

    parser = argparse.ArgumentParser()
    parser.add_argument('-n', '--graph-name', required=True, type=str, help=
                        'Name formmatted for saved graph pdf')
    parser.add_argument('-t', '--graph-type', required=True, type=str,  help=
                        'Indicate graph type, data set for type is already set')
    parser.add_argument('-x', default='locus', type=str, help=
                        'Defaults to locus if boxplot selected, other input required for regplot')
    parser.add_argument('-y', required=True, type=str, help=
                        'No default for y variable, must be selected')
    parser.add_argument('--row', type=str, help=
                        'For boxplots comparing distributions subsetted by cell features')
    parser.add_argument('--stat-annot', action='store_true', help=
                        'If provided, stat annotations will be provided for boxplot')
    parser.add_argument('--sharey', action='store_true', help=
                        'If provided, sharey will be true for boxplot facetgrid')
    

    args = parser.parse_args()

    return args

def main(args):

    sns.set(style='ticks')        
    (meta, meta_mixed) = process_data(APPENDED_META)

    if args.graph_type == 'boxplot':
        g = boxplot(meta, y=args.y)
    if args.graph_type == 'regplot':
        g = regplot(meta_mixed, y=args.y, x=args.x)

    plt.savefig(f'{ROOT_DIR}/figures/{args.graph_name}.pdf')
    plt.show()

def process_data(path):

    meta = pd.read_csv(path)
    meta.drop(meta[meta.locus_bin == 'undetected'].index, inplace=True)

    def is_mixed(x):
        has_prod = (x.rearrangement == 'productive').any()
        has_nonprod = (x.rearrangement != 'productive').any()
        return has_prod and has_nonprod

    # Filter for cells containing both productive and nonproductive rearrangements
    meta_mixed = meta.groupby('cell_id').filter(is_mixed)

    # Aggregate nonproductive and productive values into meta mixed data set
    meta_mixed.rearrangement.replace(['passenger', 'pseudogene'], 'nonproductive', inplace=True)
    meta_mixed = (meta_mixed
                  .groupby(['cell_id', 'rearrangement'])
                  .agg({'duplicate_count': 'sum',
                        'read_count': 'sum',
                        'gene_count': 'sum',
                        'RPM': 'sum', 
                        'SHM': 'mean'})
                  .unstack('rearrangement')
                  .swaplevel(0, 1, axis=1))

    # Format meta mixed for modular plotting e.g. nonproductive_RPM
    meta_mixed.columns = meta_mixed.columns.map('_'.join)

    # Set up locus bin and gene rearrangement categorical orders for order consistency
    meta.rearrangement = pd.Categorical(meta.rearrangement, GENE_REARRS)
    meta.locus = pd.Categorical(meta.locus, LOCI)
    meta = meta.sort_values(by=['rearrangement', 'locus'])

    return (meta, meta_mixed)


def boxplot(data, y, x='locus', row=None, stat_annot=False, sharey=True):

    boxplot_kw = dict(x=x, y=y, hue='locus_bin', palette='husl', fliersize=0.5)
    facetgrid_kw = dict(col='rearrangement', row=args.row, height=6, aspect=1.9, 
                        sharey=args.sharey, dropna=False, margin_titles=True) 

    g = sns.FacetGrid(data=data, **facetgrid_kw)

    # Set up y-axis based on y parameter
    if y == 'RPM':
        g.set(yscale='log', ylim=(10, 1e8))
        ylabel = 'RPM (log10)'
    if y == 'gene_count':
        ylabel = 'Number of Genes Detected'

    g.map_dataframe(annotate_boxplot, stat_annot=args.stat_annot, **boxplot_kw)
    g.set_xlabels('Locus')
    g.set_ylabels(ylabel)

    # Add figure legend and adjust layout
    plt.legend(title='Cell Bins', loc='center left', bbox_to_anchor=(1.1, 0.5))
    plt.subplots_adjust(left=0.05, right=0.85, bottom=0.1, top=0.9, wspace=0.2, hspace=0.2)

    return g


def annotate_boxplot(*args, **kwargs):

    stat_annot = kwargs.pop('stat_annot')
    ax = sns.boxplot(*args, **kwargs)

    # Generate box pairs for stats annotation 
    if stat_annot==True:
        box_pairs = []
        for locus in LOCI:
            pairs = list(product([locus], GENE_REARRS))
            paired_pairs = list(combinations(pairs, 2))
            box_pairs.extend(paired_pairs)
            
        statannot.add_stat_annotation(
            ax, plot='boxplot',
            data=kwargs['data'], x=kwargs['x'], y=kwargs['y'], 
            hue=kwargs['hue'], box_pairs=box_pairs,
            test='Mann-Whitney'
            )

def regplot(data, y, x):

    regplot_kw = dict(y=y, x=x, color='k', ci=None, 
                      scatter_kws={'s': 10, 'alpha': 0.5, 'color': 'b', 'edgecolor': 'k'}, 
                      line_kws={'linewidth': 2, 'color': 'k'}) 

    ax = sns.regplot(data=data, **regplot_kw)

    # Generate Pearson correlation coefficient and p value
    (r2, pval) = stats.pearsonr(data[x], data[y])
    stats_vals = r'$R^2$:'+ f'{r2:.2f} p:{pval:.2f}'

    # Annotate stats values
    stats_annot = AnchoredText(stats_vals, loc='upper right', frameon=False)
    ax.add_artist(stats_annot)

    def get_axis_label(var):
        label1 = 'Productive' if var.startswith('productive') else 'Nonproductive'
        label2 = (
            'SONAR Duplicate Count' if var.endswith('duplicate_count')
            else 'Number of Genes Detected' if var.endswith('gene_count')
            else 'Raw Read Count' if var.endswith('read_count')
            else var.rpartition('_')[2]
            )
        label = ' '.join([label1, label2])

        return label

    # Label axes
    ax.set(xlabel=get_axis_label(x), ylabel=get_axis_label(y))

    # Format x- and y-axis tick labels
    tick_formatter = lambda x, pos: '{:,.0f}'.format(x/1000) + 'K'
    if not (x.endswith('gene_count') or x.endswith('SHM')):
        ax.xaxis.set_major_formatter(tkr.FuncFormatter(tick_formatter))
    if not (y.endswith('gene_count') or y.endswith('SHM')):
        ax.yaxis.set_major_formatter(tkr.FuncFormatter(tick_formatter))

    plt.tight_layout()

    return ax


if __name__ == '__main__':
    args = parse_args()
    main(args)








    


    






