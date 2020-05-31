"""
   This script is for generating test plots on all merged data from data_process.py output.
   Finalized plots will be incorporated into main.py program.  

"""

from pathlib import Path
from scipy import stats
from itertools import product, combinations

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

FIG1 = f'{ROOT_DIR}/figures/RPM_rearrangements.pdf'
FIG2 = f'{ROOT_DIR}/figures/gene_counts_rearrangements.pdf'

sns.set(style='ticks')    


def load_data(path):

    data = pd.read_csv(path)

    #Set up locus bin and gene rearrangement categorical orders
    data.rearrangement = pd.Categorical(data.rearrangement, GENE_REARRS)
    data.locus = pd.Categorical(data.locus, LOCI)
    data = data.sort_values(by=['rearrangement', 'locus'])

    return data


def boxplot(data, y, row=None, stat_annote=True):
    
    titles = ['V Productive Reads\n', 'V Passenger Reads\n', 'V Pseudogene Reads\n']
    boxplot_kw = dict(x='locus', y=y, hue='locus_bin', palette='husl', fliersize=0.5)
    facetgrid_kw = dict(col='rearrangement', row=row, height=6, aspect=1.9, dropna=False) 

    if row is not None:
        facetgrid_kw.update(margin_titles=True)

    g = sns.FacetGrid(data, **facetgrid_kw)

    # Set up y-axis based on y parameter
    if y == 'RPM':
        set_kw = dict(yscale='log', ylim=(10, 1e8)) 
        ylabel='RPM (log10)'

    if y == 'gene_count':
        set_kw = dict(ylim=(0, 50)) 
        ylabel='Number of Genes'
    
    g.set(**set_kw)

    # Add boxplot and statistical significance annotations (optional)
    g.map_dataframe(annotate_boxplot, stat_annote=stat_annote, **boxplot_kw)

    # Figure formatting
    plt.legend(title='Cell Bins', loc='center left', bbox_to_anchor=(1.1, 0.5))
    plt.subplots_adjust(left=0.05, right=0.85, bottom=0.1, top=0.9, wspace=0.2, hspace=0.2)
    for (ax, title) in zip(g.axes.flat, titles):
        ax.set_title(title)
        ax.set_ylabel(ylabel)
        ax.set_xlabel('Locus')

    return g


def annotate_boxplot(*args, **kwargs):

    stat_annote = kwargs.pop('stat_annote')
    ax = sns.boxplot(*args, **kwargs)

    # Generate box pairs for stats annotation 
    if stat_annote==True:
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

#### edit for modularity for other correlation checks

def regplot(data):

    regplot_kw = dict(x='duplicate_count', y='read_count', color='k', ci=None, 
                      scatter_kws={'s': 5, 'color': 'k'}) 

    (r2, pval) = stats.pearsonr(data.duplicate_count, data.read_count)
    text = r'$R^2$:'+ f'{r2:.2f} p:{pval:.2f}'

    g = sns.regplot(data=data, **regplot_kw)
    g.text(2.5e7, 10e4, text)
    g.set_ylabel('Raw Read Count')
    g.set_xlabel('SONAR Duplicate Count')

    return g




if __name__ == "__main__":


    meta = load_data(APPENDED_META)

    RPM_boxplot = boxplot(meta, y='RPM')
    gene_count_boxplot = boxplot(meta, y='gene_count')
    specificity_boxplot = boxplot(meta, y='RPM', row='phenotype', stat_annote=False)
    phenotype_boxplot = boxplot(meta, y='RPM', hue='phenotype')
    dup_count_corrs_plot = regplot(meta[meta.rearrangement == 'productive'])

    plt.show()
#    plt.savefig()
