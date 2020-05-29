"""
This script calculates chromosomal IG gene distances between functionally rearranged genes 
and nonproductively rearranged genes (passenger/pseudogene) using data generated from 
bedtools mapped raw reads (main/loader.py pipeline). Output is saved as a csv for generation
of plots in data_plot.py. 

"""

from pathlib import Path
import pandas as pd
from sklearn.metrics import pairwise_distances
from itertools import product
import numpy as np

import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib.ticker as tkr


ROOT_DIR = Path(__file__).resolve().parent.parent
BED = f'{ROOT_DIR}/genes.sorted.bed'
SC_DATA = f'{ROOT_DIR}/final_merged_raw_reads.csv' 
FIG1 = f'{ROOT_DIR}/figures/locus_distances.pdf' 

def shaper(data):

    """
    data: single cell data with V gene assignments and characteristics

    Returns: single cell DataFrame reshaped for analyses  
    """

    data = pd.read_csv(data, index_col=[0])
    
    # Explode gene name lists
    data.loc[:, 'gene_name'] = data.loc[:, 'gene_name'].map(eval, na_action='ignore').values
    shaped_df = data.explode('gene_name')                   

    return shaped_df


sc_df = shaper(SC_DATA)


def paired_dis(bed):

    """
    bed: Immunoglobulin mapped V genes bed file  

    Returns: Paired distances matrix where distance is bp length from median 
             position of gene 1 to gene 2 for a given gene pair
             
    """
    
    V_distances = pd.read_csv(bed, sep='\t', header=None, usecols=[0, 1, 2, 3],
                           names=['chromosome', 'position_1', 'position_2', 'gene'])

    V_distances['center'] = V_distances.loc[:, ['position_1', 'position_2']].mean(axis=1)


    # Calculate pairwise distances between genes
    paired_dis = pairwise_distances(V_distances.center.values.reshape(-1, 1))
    paired_dis = pd.DataFrame(paired_dis, index=V_distances.gene, columns=V_distances.gene)

    # Remove pairwise distances between genes on separate chromosomes
    for chromosome in ['chr2', 'chr14', 'chr22']:
        gene_names = V_distances.loc[V_distances.chromosome == chromosome, 'gene']
        other_genes = V_distances.loc[V_distances.chromosome != chromosome, 'gene']
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
        gene1 = data.loc['functional', 'gene_name']    
        gene2 = data.loc['passenger', 'gene_name']

        pass_dis = [paired_dis.loc[x, y] for (x, y) in product(gene1, gene2)]
        pass_dis = np.nanmean(pass_dis)

        result['passenger'] = pass_dis

    if any(data.index == 'functional') & any(data.index == 'pseudogene'):
        gene1 = data.loc['functional', 'gene_name']
        gene3 = data.loc['pseudogene', 'gene_name']

        pseud_dis = [paired_dis.loc[x, y] for (x, y) in product(gene1, gene3)]
        pseud_dis = np.nanmean(pseud_dis)

        result['pseudogene'] = pseud_dis


    return result


# Subset locus genes for single cells and retrieve average distance
result = {}
loci = ['VH', 'VK', 'VL']

for i, (cell, locus) in enumerate(product(sc_df.index, loci)):
    print("{:,}/{:,}".format(i, len(sc_df.locus)), end='\r') # loader bar                       

    if locus in set(sc_df.loc[cell, 'locus']):
        data = sc_df[(sc_df.index==cell) & (sc_df.locus==locus)]
        result[(cell, locus)] = compute_dis(data)
        
distance_df = pd.DataFrame(result)



def merge(df1, df2):
    """ 
    df1: single cell DataFrame 
    df2: single cell calculated V gene distances
    
    Returns: Merged DataFrame
    """
    # Format for merge
    df1 = df1[df1.rearrangement=='functional']
    df1 = (df1.reset_index().set_index(['cell_id', 'locus'])
           .drop(columns=['rearrangement', 'bedtools_RPM', 'duplicate_count'])
           )

    df2 = df2.T
    df2.index.names = ['cell_id', 'locus']
    df2 = (df2.reset_index()
           .melt(id_vars=['cell_id', 'locus'], var_name='rearrangement', 
                 value_name='distance')
           .set_index(['cell_id', 'locus'])
           )
           
    # Format output
    merged_df = df2.merge(df1, left_index=True, right_index=True)    
    merged_df = merged_df.rename(columns={'gene_name': 'functional_gene'})
    merged_df = merged_df.reset_index().set_index('cell_id')

    return merged_df


merged_df = merge(sc_df, distance_df)



def plot1(merged_df):
    """
    Returns: Plot 1 = Average Interlocus Distance vs. Average SHM %
    """
      
    sns.set_style("whitegrid") 

    # Plot data
    plot1 = sns.FacetGrid(data=locus_df, col='rearrangement', row='locus', hue='rearrangement')
    plot1.map(sns.kdeplot, 'distance', 'SHM', shade_lowest=False, shade=True)
    plot1.map(sns.scatterplot, 'distance', 'SHM', alpha=0.4, edgecolor=None, s=1, color='black')

    # Format graph
    plot1.set(xlim=(0, 1e6), ylim=(0, 20), 
              xlabel='Average Interlocus Distance (bp)', ylabel='Average SHM (%)')     
    plot1.set_titles(template='{col_name}: {row_name}', fontweight='bold', fontsize=18)
    plt.subplots_adjust(hspace = 0.4, wspace = 0.3)

    for ax in plot1.axes.flat:
        ax.xaxis.set_major_formatter(
            tkr.FuncFormatter(lambda x, p: "{:,}k".format(int(x/1000))))

    
    plt.savefig(FIG1)
    
    plt.show()














        



