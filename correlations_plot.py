"""
   New script to plot correlations without modularity
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



def process_data(path):
    meta = pd.read_csv(path)

    # Label passenger and pseudogene rearrangement types as nonproductive
    meta_mixed = meta.copy()
    meta_mixed.rearrangement_type.replace(['passenger', 'pseudogene'], 'nonproductive', inplace=True)
    
    meta_mixed = (meta_mixed
                  .groupby(['cell_id', 'rearrangement_type', 'locus', 'baldr_bin'])
                  .agg({'read_count': 'sum',
                        'n_genes': 'sum',
                        'RPM': 'sum', 
                        'SHM': 'mean'}))
    meta_mixed = meta_mixed.reset_index().set_index('cell_id')

    return (meta, meta_mixed)

def regplot(data):

    regplot_kw = dict(edgecolors='black', linewidths=0.1, s=20)
    facetgrid_kw = dict(col='locus', dropna=False, hue='baldr_bin', palette=FLATUI, 
                        margin_titles=True, height=6, aspect=1.9)
    
    # Uncomment for nonproductive RPM vs productive RPM
    # data = data.set_index(['rearrangement_type', 'baldr_bin', 'locus'], append=True).RPM.unstack('rearrangement_type').reset_index()
    # g = sns.FacetGrid(data, **facetgrid_kw)
    # g.map(plt.scatter, 'nonproductive', 'productive')
    # g.set(yscale='log')
    # g.set(ylabel='productive RPM (log10)', ylim=(1e0, 1e6)) 
    
    # Uncomment for nonproductive RPM versus SHM
    data = data.drop(columns=['read_count', 'n_genes']).set_index(['locus', 'baldr_bin', 'rearrangement_type'], append=True).unstack()
    data.columns = ['_'.join(x[::-1]) for x in data.columns]
    data.reset_index(inplace=True)
    g = sns.FacetGrid(data, **facetgrid_kw)

    g.map(plt.scatter, 'nonproductive_RPM', 'productive_SHM', **regplot_kw)

    g.set(xscale='log')
    g.set(xlabel='nonproductive RPM (log10)', xlim=(1e0, 1e6))

    plt.legend(title='BALDR Locus Bins', loc='center left', bbox_to_anchor=(1.1, 0.5))
    plt.subplots_adjust(left=0.05, right=0.85, bottom=0.1, top=0.8, wspace=0.1, hspace=0.2)

    # Annotate number of data points per rearrangement and plot
    # for coords,data in g.facet_data():
    #     ax = g.facet_axis(*coords[:2])
    #     label = data.locus_bin.unique()[0]
    #     title = f'{ax.get_title()}\n<{label}={data.shape[0]}>'
    #     ax.set_title(title)


    plt.show()
     
    return g

(meta, meta_mixed) = process_data(APPENDED_META)
regplot(meta_mixed)

