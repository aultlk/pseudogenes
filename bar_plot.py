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
LOCUS_BINS = ('productive', 'passenger', 'pseudogene', 'undetected')
FLATUI = ['g', 'darkorange', 'deepskyblue', 'deeppink']  # productive, passenger, pseudogene, undetected



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
    parser.add_argument('--hue', type=str, help=
                        'For barplots comparing distributions subsetted by cell features')
    parser.add_argument('--stat-annot', action='store_true', help=
                        'If provided, stat annotations will be provided for boxplot')
    parser.add_argument('--sharey', action='store_true', help=
                        'If provided, sharey will be true for boxplot facetgrid')
    parser.add_argument('--log-yscale', action='store_true', help=
                        'If provided, y axis will be logscale')
    parser.add_argument('--log-xscale', action='store_true', help=
                        'If provided, x axis will be logscale')
    

    args = parser.parse_args()

    return args

def main(args):

    sns.set(style='ticks')        
    (meta, meta_mixed) = process_data(APPENDED_META)

    if args.graph_type == 'boxplot':
        g = boxplot(meta, y=args.y)
    if args.graph_type == 'regplot':
        g = regplot(meta_mixed, y=args.y, x=args.x)
    if args.graph_type == 'barplot':
        g = barplot(meta, hue=args.hue)

    plt.savefig(f'{ROOT_DIR}/figures/{args.graph_name}.pdf')
    plt.show()

def process_data(path):

    meta = pd.read_csv(path)
    # meta.drop(meta[meta.locus_bin == 'undetected'].index, inplace=True)

    # Aggregate nonproductive and productive values into meta mixed data set
    meta_mixed = meta.copy()


    meta_mixed.rearrangement_type.replace(['passenger', 'pseudogene'], 'nonproductive', inplace=True)

    # Filter out cells with no baldr bin assigned ***** COME BACK TO LATER *******
    meta_mixed = meta_mixed[~meta_mixed.baldr_bin.isnull()]

    import ipdb; ipdb.set_trace()

    meta_mixed = (meta_mixed
                  .groupby(['cell_id', 'rearrangement_type', 'locus', 'baldr_bin'])
                  .agg({'read_count': 'sum',
                        'gene_count': 'sum',
                        'RPM': 'sum'})) 
    import ipdb; ipdb.set_trace()


    meta_mixed = meta_mixed.reset_index().set_index('cell_id')
    

    # # Format meta mixed for modular plotting e.g. nonproductive_RPM
    # meta_mixed.columns = meta_mixed.columns.map('_'.join)
    # meta_mixed.reset_index(inplace=True)

    # Set up locus bin and gene rearrangement categorical orders for order consistency
    meta.rearrangement_type = pd.Categorical(meta.rearrangement_type, GENE_REARRS)
    meta.locus = pd.Categorical(meta.locus, LOCI)
    meta = meta.sort_values(by=['rearrangement_type', 'locus'])


    return (meta, meta_mixed)


def boxplot(data, y, x='locus', row=None, stat_annot=False, sharey=False, log_yscale=False, log_xscale=False):

    boxplot_kw = dict(x=x, y=y, hue='baldr_bin', palette=FLATUI, fliersize=0.5)
    facetgrid_kw = dict(col='rearrangement_type', row=args.row, height=6, aspect=1.9, 
                        sharey=args.sharey, dropna=False, margin_titles=True) 

    g = sns.FacetGrid(data=data, **facetgrid_kw)

    # Set up y-axis based on y parameter
    if y == 'RPM':
        ylabel = 'RPM'
    if y == 'gene_count':
        ylabel = 'Number of Genes Detected'
    if args.log_yscale == True:
        g.set(yscale='log', ylim=(1e0, 1e6))
        ylabel += ' (log10)'

    g.map_dataframe(annotate_boxplot, stat_annot=args.stat_annot, **boxplot_kw)
    g.set_xlabels('Locus')
    g.set_ylabels(ylabel)

    # Add figure legend and adjust layout
    plt.legend(title='Cell Bins', loc='center left', bbox_to_anchor=(1.1, 0.5))
    plt.subplots_adjust(left=0.05, right=0.85, bottom=0.1, top=0.9, wspace=0.1, hspace=0.2)
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

def regplot(data, y, x, log_yscale=False, log_xscale=False):

    regplot_kw = dict(edgecolors='black', linewidths=0.1, s=20)
    facetgrid_kw = dict(col='locus', dropna=False, hue='locus_bin', palette=FLATUI, 
                        margin_titles=True, height=6, aspect=1.9) 


    g = sns.FacetGrid(data=data, **facetgrid_kw)
    g.map(plt.scatter, data[data.rearrangement == 'nonproductive'].RPM, data[data.rearrangement == 'nonproductive'].SHM, **regplot_kw)



    # Set axis based on y parameter
    if args.log_yscale == True:
        g.set(yscale='log')
        g.set(ylabel=x + ' (log10)', ylim=(1e0, 1e6))       
    if args.log_xscale == True:
        g.set(xscale='log')
        g.set(xlabel=x + ' (log10)', xlim=(1e0, 1e6))

    plt.legend(title='Cell Bins', loc='center left', bbox_to_anchor=(1.1, 0.5))
    plt.subplots_adjust(left=0.05, right=0.85, bottom=0.1, top=0.8, wspace=0.1, hspace=0.2)

    # Annotate number of data points per rearrangement and plot
    for coords,data in g.facet_data():
        ax = g.facet_axis(*coords[:2])
        label = data.locus_bin.unique()[0]
        title = f'{ax.get_title()}\n<{label}={data.shape[0]}>'
        ax.set_title(title)

    # # Getting error which might be due to version of numpy
    # data = data.loc[:, [x, y]]
    # data = data.replace([np.inf, -np.inf], np.nan).dropna(subset=[x, y], how="all")

    # (r2, pval) = stats.pearsonr(data[x], data[y])
    # stats_vals = r'$R^2$:'+ f'{r2:.2f} p:{pval:.2f}'
    # stats_annot = AnchoredText(stats_vals, loc='upper right', frameon=False)    
    # g.map_dataframe(stats_annot)

    # # Generate Pearson correlation coefficient and p value
    # data = data.dropna()
    # (r2, pval) = stats.pearsonr(data[x], data[y])
    # stats_vals = r'$R^2$:'+ f'{r2:.2f} p:{pval:.2f}'

    # g.add_artist(stats_annot)

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
    g.set(xlabel=get_axis_label(x), ylabel=get_axis_label(y))

    return g

# def bar_plot(meta, x='locus_bin', hue=args.hue, stat_annot=False, sharey=False, log_yscale=False, log_xscale=False):

#     barplot_kw = dict(x=x, y=len(meta), hue=hue)
  
#     ax = sns.barplot(data=meta, **barplot_kw)

#     ax.set_xlabels('Locus')
#     ax.set_ylabels('Number of Cells')

#     # Add figure legend and adjust layout
#     plt.legend(title=hue, loc='center left', bbox_to_anchor=(1.1, 0.5))
#     plt.subplots_adjust(left=0.05, right=0.85, bottom=0.1, top=0.9, wspace=0.1, hspace=0.2)
#     return g




if __name__ == '__main__':
    args = parse_args()
    main(args)








    


    






