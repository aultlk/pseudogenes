""" 
This script plots preliminary data on collected bedtools mapped read counts and average SHM for rearrangements 
(productive/nonproductive) at VH/VL locus or all loci: VHL.
"""

from pathlib import Path
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

ROOT_DIR = Path(__file__).resolve().parent.parent

# Load all single cell data for analysis control checks
mapped_rearrangements = (pd.read_csv(f'{ROOT_DIR}/final_merged_raw_reads.csv', 
                                     header=[0,1], 
                                     index_col=0)
                         .swaplevel(0,1, axis=1))

# Subset productive rearrangements for correlation check
productive_counts = (mapped_rearrangements
                     .loc[:, ['bedtools_read_count', 'duplicate_count']]
                     .swaplevel(0,1, axis=1)
                     .drop(columns=['VH_pseudogene', 'VH_passenger', 'VL_pseudogene', 'VL_passenger', 'VHL_pseudogene', 'VHL_passenger']))

# Test correlation between bedtools read counts and sonar duplicate counts for productive rearrangements 
productive_counts_corr = (productive_counts
                          .xs('bedtools_read_count', level=1, axis=1)
                          .corrwith(productive_counts.xs('duplicate_count', level=1, axis=1)))

print(f'Correlations: \n{productive_counts_corr}')

# Subset counts for all rearrangement types and reformat for plotting
rearrangement_counts = mapped_rearrangements.loc[:, 'bedtools_read_count'].rename_axis(index='cell_id')
rearrangement_counts.columns = pd.MultiIndex.from_tuples(
    [col.split('_') for col in rearrangement_counts.columns])

rearrangement_counts = rearrangement_counts.reset_index().melt(id_vars='cell_id', var_name=['locus', 'rearrangement'], value_name='bedtools_read_count')

collective_sums = (mapped_rearrangements.loc[:, 'bedtools_read_count'].sum().map(lambda x: '{:,}'.format(x)))
print(f'Collective sums: \n{collective_sums}')

# Plot number of reads per rearrangement type at loci
ax = sns.boxplot(x="locus", y="bedtools_read_count", hue='rearrangement', data=rearrangement_counts)
ax = sns.stripplot(x="locus", y="bedtools_read_count", hue='rearrangement', dodge=True, data=rearrangement_counts)

plt.show()



