"""
This script calculates chromosomal IG gene distances between functionally rearranged genes and nonproductively 
rearranged genes (passenger/pseudogene) for a given single cell. Correlations are then tested against the average
chromosomal distance between productive/nonproductively rearranged genes and average SHM on functional rearrangement.
"""

from pathlib import Path
import pandas as pd
from sklearn.metrics import pairwise_distances
import seaborn as sns
import matplotlib.pyplot as plt

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
    
# Retrieve single cell rearrangements (functional/pseudogene/passenger) and functional gene SHM percentage
single_cell_data = (pd.read_csv(f'{ROOT_DIR}/final_merged_raw_reads.csv', header=[0,1], index_col=[0]).swaplevel(0, 1, axis=1)
                    .rename_axis(index='cell_id'))

# Reformat string of gene names to list of gene names, drop columns not needed
single_cell_data.loc[:,'gene_name'] = single_cell_data.loc[:,'gene_name'].fillna('float("nan")').applymap(eval).values
single_cell_data = (single_cell_data.drop(columns=['bedtools_read_count', 'duplicate_count', 'gene_count'])
                    .swaplevel(0, 1, axis=1)
                    .drop(columns=['VHL_functional', 'VHL_passenger', 'VHL_pseudogene'])
                    .drop(columns=('VL_pseudogene', 'SHM'))
                    .swaplevel(0, 1, axis=1))

# Filter cells with any NaN values
cells_w_nan = single_cell_data[single_cell_data.isnull().any(axis=1)].index
single_cell_data = single_cell_data.drop(index=cells_w_nan)

# Create columns for extracted distances 
single_cell_data.reindex(list(single_cell_data)+
                        [('avg_distance', 'VH_passenger'), ('avg_distance', 'VH_pseudogene'),
                         ('avg_distance', 'VL_passenger'), ('avg_distance', 'VL_pseudogene')], axis=1)

for cell in single_cell_data.index:
    VH_functional_genes = single_cell_data.loc[cell, ('gene_name', 'VH_functional')]
    VH_passenger_genes = single_cell_data.loc[cell, ('gene_name', 'VH_passenger')]
    VH_pseudogenes = single_cell_data.loc[cell, ('gene_name', 'VH_pseudogene')]
    VL_functional_genes = single_cell_data.loc[cell, ('gene_name', 'VL_functional')]
    VL_passenger_genes = single_cell_data.loc[cell, ('gene_name', 'VL_passenger')]
    VL_pseudogenes = single_cell_data.loc[cell, ('gene_name', 'VL_pseudogene')]
    
    avg_distance = paired_distances_df.loc[VH_functional_genes, VH_passenger_genes].mean(axis=1)
    single_cell_data.loc[cell, ('avg_distance', 'VH_passenger')] = avg_distance.iloc[0].round(1)

    avg_distance = paired_distances_df.loc[VH_functional_genes, VH_pseudogenes].mean(axis=1)
    single_cell_data.loc[cell, ('avg_distance', 'VH_pseudogene')] = avg_distance.iloc[0].round(1)

    avg_distance = paired_distances_df.loc[VL_functional_genes, VL_passenger_genes].mean(axis=1)
    single_cell_data.loc[cell, ('avg_distance', 'VL_passenger')] = avg_distance.iloc[0].round(1)
   
    avg_distance = paired_distances_df.loc[VL_functional_genes, VL_pseudogenes].mean(axis=1)
    single_cell_data.loc[cell, ('avg_distance', 'VL_pseudogene')] = avg_distance.iloc[0].round(1)
   
# Reshape single cell data for plotting
single_cell_data.columns = pd.MultiIndex.from_tuples([[lvl1] + lvl2.split('_') for (lvl1, lvl2) in single_cell_data.columns])
single_cell_data.columns.names=['variable', 'locus', 'rearrangement']
df1 = single_cell_data['SHM'].reset_index().melt(id_vars='cell_id', value_name='SHM')
df2 = single_cell_data['avg_distance'].reset_index().melt(id_vars='cell_id', value_name='avg_distance')
single_cell_data = df1.drop(columns='rearrangement').merge(df2, on=['cell_id', 'locus'])

# Plot average locus distance to functional (VH/VL passenger/pseudogene) versus functional SHM % 
sns.set()
# g = sns.FacetGrid(col='locus', hue='rearrangement', data=single_cell_data, sharey=False)
# g.map(sns.jointplot, 'avg_distance', 'SHM')#, scatter_kws={'s': 3, 'alpha': 0.5})
# g.add_legend()
# g = sns.lmplot(x='avg_distance', y='SHM', hue='rearrangement', data=single_cell_data)
# g.set_axis_labels("average locus distance", "average functional SHM%")
# sns.jointplot('avg_distance', 'SHM', data=single_cell_data[(single_cell_data.rearrangement == 'passenger') & (single_cell_data.locus == 'VH')], kind='kde')

plt.show()










        



