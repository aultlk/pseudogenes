""" 
This script plots preliminary data on collected bedtools mapped read counts and average SHM for rearrangements 
(productive/nonproductive) at VH/VL locus or all loci: VHL.
"""

from pathlib import Path
import pandas as pd
from scipy.stats import mannwhitneyu
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
rearrangement_counts = rearrangement_counts.set_index('cell_id')

# Bin each rearrangement as productive or nonproductive
rearrangement_counts['bin'] = (rearrangement_counts.rearrangement
                               .replace({'functional': 'productive', 'passenger': 'nonproductive', 'pseudogene': 'nonproductive'}))

# Plot number of reads per rearrangement type at loci
g1 = sns.boxplot(x="locus", y="bedtools_read_count", hue='rearrangement', data=rearrangement_counts, fliersize=1)
g1 = sns.stripplot(x="locus", y="bedtools_read_count", hue='rearrangement', s=2, dodge=True, data=rearrangement_counts)

plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
g1.set_yscale("log")
g1.set(ylabel="Mapped Raw Read Counts (log10)")
g1.set(xlabel="Locus")

plt.show()

g2 = sns.violinplot(x="locus", y="bedtools_read_count", hue="bin", data=rearrangement_counts, palette="Set2", split=True, scale="count")
g2.set(ylabel="Mapped Raw Read Counts")
g2.set(xlabel="Locus")

plt.show()

#collective_sums = (mapped_rearrangements.loc[:, 'bedtools_read_count'].sum().map(lambda x: '{:,}'.format(x)))
#print(f'Collective sums: \n{collective_sums}')


# Calculate p values to compare nonproductive versus productive groups at VHL loci, and individual VH and VL locus
# def prod_nonprod_mannwhitneyu(loci="all_case"):
#    """
#    loci: By default, compares nonproductive and productive reads across cells for all three cases (all_case) - VH only, VL only, and VHL
#    Returns: Mann Whitney U p-value and corresponding case (loci) 
#    """

#    functional_reads = (rearrangement_counts['rearrangement'] == 'functional')
#    passenger_reads = (rearrangement_counts['rearrangement'] == 'passenger')
#    pseudogene_reads = (rearrangement_counts['rearrangement'] == 'pseudogene')
#    VH_loci = (rearrangement_counts['locus'] == 'VH')
#    VL_loci = (rearrangement_counts['locus'] == 'VL')
#    VHL_loci = (rearrangement_counts['locus'] == 'VHL')
   
#    # VH loci: passenger reads
#    VH_p_value1 = mannwhitneyu(
#       x=rearrangement_counts.loc
#       [functional_reads & VH_loci]['bedtools_read_count'],
#       y=rearrangement_counts.loc[passenger_reads & VH_loci]['bedtools_read_count'])

#    # VH loci: pseudogene reads
#    VH_p_value2 = mannwhitneyu(
#       x=rearrangement_counts.loc
#       [functional_reads & VH_loci]['bedtools_read_count'],
#       y=rearrangement_counts.loc[pseudogene_reads & VH_loci]['bedtools_read_count'])

#    # VL loci: passenger reads
#    VL_p_value1 = mannwhitneyu(x=rearrangement_counts.loc
#                 [functional_reads & VL_loci]['bedtools_read_count'],
#                 y=rearrangement_counts.loc
#                 [passenger_reads & VL_loci]['bedtools_read_count'])

#    # VL loci: pseudogene reads
#    VL_p_value2 = mannwhitneyu(x=rearrangement_counts.loc
#                 [functional_reads & VL_loci]['bedtools_read_count'],
#                 y=rearrangement_counts.loc
#                 [pseudogene_reads & VL_loci]['bedtools_read_count'])

   
#    # VHL loci: passenger reads
#    VHL_p_value1 = mannwhitneyu(x=rearrangement_counts.loc
#                 [functional_reads & VHL_loci]['bedtools_read_count'],
#                 y=rearrangement_counts.loc
#                 [passenger_reads & VHL_loci]['bedtools_read_count'])

#    # VHL loci: pseudogene reads
#    VHL_p_value2 = mannwhitneyu(x=rearrangement_counts.loc
#                 [functional_reads & VHL_loci]['bedtools_read_count'],
#                 y=rearrangement_counts.loc
#                 [pseudogene_reads & VHL_loci]['bedtools_read_count'])

#    stats_display = "Mann Whitney U Test Results:"
#    VH_stats = f"\nVH locus: \n\t Passenger: {VH_p_value1} \n\t Pseudogene: {VH_p_value2}"
#    VL_stats = f"\nVL locus: \n\t Passenger: {VL_p_value1} \n\t Pseudogene: {VL_p_value2}"
#    VHL_stats = f"\nVHL locus: \n\t Passenger: {VHL_p_value1} \n\t Pseudogene: {VHL_p_value2}"
   
#    return [print(f"{stats_display} \n\t{VH_stats} \n\t{VL_stats} \n\t{VHL_stats}")]
   

#prod_nonprod_mannwhitneyu()



