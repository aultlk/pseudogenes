"""
This script calculates chromosomal IG gene distances between functionally rearranged genes 
and nonproductively rearranged genes (passenger/pseudogene). Test correlations and a plot
is returned on processed data.  

"""

from pathlib import Path
import pandas as pd
from sklearn.metrics import pairwise_distances
from itertools import product
import numpy as np

import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib.ticker as tkr

idx = pd.IndexSlice
ROOT_DIR = Path(__file__).resolve().parent.parent
BED = f'{ROOT_DIR}/genes.sorted.bed'
SC_DATA = f'{ROOT_DIR}/final_merged_raw_reads.csv' 




def shaper(data):

    """
    data: single cell data with V gene assignments and characteristics

    Returns: single cell DataFrame reshaped for analyses  
    """

    data = pd.read_csv(data, header=[0, 1, 2], index_col=[0])
    
    shaped_df = data.loc[:, idx[:, :, ['gene_name', 'SHM']]].stack(['locus', 'rearrangement'])
    shaped_df.loc[:, 'gene_name'] = shaped_df.loc[:, 'gene_name'].map(eval, na_action='ignore').values
    shaped_df = shaped_df.explode('gene_name')                   

    return shaped_df

sc_df = shaper(SC_DATA)


def paired_dis(bed):

    """
    bed: Immunoglobulin mapped V genes bed file  

    Returns: Paired distances matrix where distance is bp length from median 
             position of gene 1 to gene 2 for a given gene pair
             
    """
    
    IG_distances = pd.read_csv(BED, sep='\t', header=None, usecols=[0, 1, 2, 3],
                           names=['chromosome', 'position_1', 'position_2', 'gene'])

    IG_distances['center'] = IG_distances.loc[:, ['position_1', 'position_2']].mean(axis=1)


    # Calculate pairwise distances between genes
    paired_dis = pairwise_distances(IG_distances.center.values.reshape(-1, 1))
    paired_dis = pd.DataFrame(paired_dis, index=IG_distances.gene, columns=IG_distances.gene)

    # Remove pairwise distances between genes on separate chromosomes
    for chromosome in ['chr2', 'chr14', 'chr22']:
        gene_names = IG_distances.loc[IG_distances.chromosome == chromosome, 'gene']
        other_genes = IG_distances.loc[IG_distances.chromosome != chromosome, 'gene']
        paired_dis.loc[gene_names, other_genes] = float('nan')

    return paired_dis
    
paired_dis = paired_dis(BED)


def compute_dis(data):

    """
    data = Subsetted single cell data for each cell and locus
     
    Returns = DataFrame with average pseudogene/passenger distance to functional gene
               for a given locus and single cell
    """

    data = data.groupby('rearrangement').agg(list)
    
    result = {}

    if any(data.index == 'functional') & any(data.index == 'passenger'):
        gene1 = data.functional    
        gene2 = data.passenger
        
        pass_dis = [paired_dis.loc[x, y] for (x, y) in product(gene1, gene2)]
        pass_dis = np.mean(pass_dis)

        result['passenger'] = pass_dis

    if any(data.index == 'functional') & any(data.index == 'pseudogene'):
        gene1 = data.functional
        gene3 = data.pseudogene

        pseud_dis = [paired_dis.loc[x, y] for (x, y) in product(gene1, gene3)]
        pseud_dis = np.mean(pseud_dis)

        result['pseudogene'] = pseud_dis

    return result


# Subset locus genes for single cells and retrieve average distance
result = {}

cells = sc_df.index.get_level_values('cell_id')
loci = ['VH', 'VL']

for i, (cell, locus) in enumerate(product(cells, loci)): 
    print("{:,}/{:,}".format(i, len(sc_df)*2), end='\r') # loader bar

    data = sc_df.loc[cell, 'gene_name']

    if locus in pd.MultiIndex.get_level_values(data.index, 'locus'):
        result[(cell, locus)] = compute_dis(data.loc[locus])

distance_df = pd.DataFrame(result)


def merge(df1, df2):
    """ 
    df1: single cell DataFrame 
    df2: single cell calculated V gene distances
    
    Returns: Merged DataFrame
    """
    
    df2 = df2.T
    df2.index.names = ['cell_id', 'locus']

    df1 = df1.xs('functional', level='rearrangement').SHM
    
    merged_df = (df2.merge(df1, left_index=True, right_index=True)
                 .reset_index()
                 .melt(id_vars=['cell_id', 'locus', 'SHM'], 
                       var_name='rearrangement', value_name='distance')
    )

    return merged_df

merged_df = merge(sc_df, distance_df)



def plot(data):
    """
     data: DataFrame with merged SHM and locus distances
       
     Returns: Plots (1) avg locus distance versus SHM %
    """
    
    # Set up panels for each locus and rearrangement
    sns.set(style="whitegrid")
    plot = sns.FacetGrid(data, col='rearrangement', row='locus', hue='locus')
    
    # Map kdeplot and scatterplot on panels
    plot.map(sns.kdeplot, 'distance', 'SHM', shade_lowest=False, shade=True)
    plot.map(sns.scatterplot, 'distance', 'SHM', 
             alpha=0.4, edgecolor=None, s=1, color='black'
    )
    
    # Format numbering scale on x-axis (bp length)
    for ax in plot.axes[0]:
        ax.xaxis.set_major_formatter(
            tkr.FuncFormatter(lambda x, p: "{:,}k".format(int(x/1000)))
        )
        plt.subplots_adjust(hspace = 0.4, wspace = 0.3)
    

    # Format graph 
    plot.set(
        xlim=(0, 8e5), ylim=(0, 20), 
        xlabel='Average Interlocus Distance (bp)', 
        ylabel='Average SHM (%)'
    )
    plot.set_titles(template='{col_name}: {row_name}', 
                    fontweight='bold', fontsize=18
    )
               
    return plot


plot(merged_df)

plt.show()














        



