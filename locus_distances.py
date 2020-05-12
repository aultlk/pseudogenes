"""
This script calculates chromosomal IG gene distances between functionally rearranged genes and nonproductively 
rearranged genes (passenger/pseudogene) for a given single cell. Correlations are then tested against the average
chromosomal distance between productive/nonproductively rearranged genes and average SHM on functional rearrangement.
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



# Load IG mapped genes bed file and calculate 'center' of position 1 and 2 for each gene
IG_distances = pd.read_csv(f'{ROOT_DIR}/genes.sorted.bed', sep='\t', header=None, usecols=[0, 1, 2, 3],
                           names=['chromosome', 'position_1', 'position_2', 'gene'])
IG_distances['center'] = IG_distances.loc[:, ['position_1', 'position_2']].mean(axis=1)

# Calculate pairwise distances between genes
paired_distances = pairwise_distances(IG_distances.center.values.reshape(-1, 1))
paired_distances_df = pd.DataFrame(paired_distances, index=IG_distances.gene, columns=IG_distances.gene)

# Remove pairwise distances between genes on separate chromosomes
for chromosome in ['chr2', 'chr14', 'chr22']:
    gene_names = IG_distances.loc[IG_distances.chromosome == chromosome, 'gene']
    other_genes = IG_distances.loc[IG_distances.chromosome != chromosome, 'gene']
    paired_distances_df.loc[gene_names, other_genes] = float('nan')

    

# Retrieve single cell data and reformat DataFrame to distinguish IGH, IGL, IGK loci (heavy, lambda, kappa)
single_cell_data = (pd.read_csv(f'{ROOT_DIR}/final_merged_raw_reads.csv', header=[0,1], index_col=[0]).swaplevel(0, 1, axis=1)
                    .rename_axis(index='cell_id'))

single_cell_data.loc[:,'gene_name'] = single_cell_data.loc[:,'gene_name'].fillna('float("nan")').applymap(eval).values
single_cell_data = (single_cell_data.drop(columns=['bedtools_read_count', 'duplicate_count', 'gene_count'])
                    .swaplevel(0, 1, axis=1).drop(columns=['VHL_functional', 'VHL_passenger', 'VHL_pseudogene'])
                    .drop(columns=('VL_pseudogene', 'SHM')))

single_cell_data.columns = pd.MultiIndex.from_tuples([x.split('_')+[y] for (x,y) in single_cell_data.columns])
single_cell_data = single_cell_data.stack(level=[0, 1]).explode('gene_name')
single_cell_data.index.names = ['cell_id', 'chain', 'rearrangement']
single_cell_data['locus'] = single_cell_data.gene_name.map(lambda x: x[:3])
single_cell_data = single_cell_data.reset_index().set_index(['cell_id', 'locus', 'rearrangement'])

single_cell_data = single_cell_data.dropna(subset=['gene_name'], how='any')


# Retrieve gene distances and calculate average
def compute_distance(data):
    """
    data = Subsetted locus genes for each individual cell
       
    Returns = DataFrame with average pseudogene/passenger distance to functional gene
              for a given locus and single cell
    """

    data = data.groupby(level=0).gene_name.agg(list)

    result = {}
    
    if any(data.index == 'functional') & any(data.index == 'passenger'):
        gene1 = data.functional    
        gene2 = data.passenger
        
        passenger_distance = [paired_distances_df.loc[x, y] for (x, y) in product(gene1, gene2)]
        passenger_distance = np.mean(passenger_distance)

        result['passenger'] = passenger_distance

    if any(data.index == 'functional') & any(data.index == 'pseudogene'):
        gene1 = data.functional
        gene3 = data.pseudogene

        pseudogene_distance = [paired_distances_df.loc[x, y] for (x, y) in product(gene1, gene3)]
        pseudogene_distance = np.mean(pseudogene_distance)

        result['pseudogene'] = pseudogene_distance

    return result


result = {}
for i, (cell_id, locus) in enumerate(single_cell_data.reset_index()[['cell_id', 'locus']].values):
    print("{:,}/{:,}".format(i, len(single_cell_data)), end='\r') # loader bar
    if single_cell_data.loc[(cell_id, locus)].size > 0:
        result[(cell_id, locus)] = compute_distance(single_cell_data.loc[(cell_id, locus)])

result_df = pd.DataFrame(result)


# Map to distances to functional SHM for plotting
result_df = result_df.T
result_df.index.names = ['cell_id', 'locus']

SHM_data = single_cell_data.xs('functional', level='rearrangement').SHM
result_df = result_df.merge(SHM_data, left_index=True, right_index=True)
result_df = result_df.reset_index().melt(id_vars=['cell_id', 'locus', 'SHM'], var_name='rearrangement', value_name='distance')

sns.set(style="whitegrid")
g = sns.FacetGrid(data=result_df, col='rearrangement', row='locus', hue='locus')
g.map(sns.kdeplot, 'distance', 'SHM', shade_lowest=False, shade=True)
g.map(sns.scatterplot, 'distance', 'SHM', alpha=0.4, edgecolor=None, s=1, color='black')
g.set(xlim=(0, 8e5), ylim=(0, 20), xlabel='Average Interlocus Distance (bp)', ylabel='Average SHM (%)')
g.set_titles(template='{col_name}: {row_name}', fontweight='bold')

for ax in g.axes[0]:
    ax.xaxis.set_major_formatter(
        tkr.FuncFormatter(lambda x, p: "{:,}k".format(int(x/1000))))
plt.subplots_adjust(hspace = 0.4, wspace = 0.1)

plt.show()











        



