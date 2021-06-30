"""
   This script loads the merged data set with read counts from Judah's pipeline and sonar-baldr
   gene assignments and rearrangement types. Then it plots the top 2 highest read counts as well
   as ratios. Purpose is to use these plots for developing a criteria for manually assigning 
   cell bins for cells that clearly have a nonproductive rearrangement in addition to productive.
"""

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

# # Choose mode: conditions to subset cells on
# # NoProd: No productive rearrangement assigned by baldrsonar
# # HasProd: Productive rearrangement assigned by baldrsonar 
# has_prod = dict(heavy=True, light=True)

def second_largest(l):
    """
    Computes 2nd largest read count.
    If len(read_counts) == 1, returns 0
    """
    if len(l) == 1:
        return 0
    return sorted(l)[-2]

df = pd.read_csv('troubleshooting/merged_data.csv')
cols = ['cell_name', 'mapped_reads']
meta = pd.read_csv('../cellInfo_allMaster.csv', sep='\t', index_col='cell_name', usecols=cols)
meta.index.name = 'cell_id'

# Filter out cells with single gene read counts
df.locus = pd.Categorical(df.locus, ['IGH', 'IGK', 'IGL'])
groups = df.groupby(['cell_id', 'locus'])

# Create new dataframe with  NP/P indication
rearr_info = groups.rearrangement_type.agg(
    has_P=lambda s: any(s=='productive'), 
    has_NP=lambda s: any((s=='pseudogene') | (s=='passenger'))
)

# Add mean SHM for productive rearrangements per cell, locus
rearr_info['mean_prod_SHM'] = groups.apply(lambda df: df[df.rearrangement_type == 'productive'].SHM.mean())

# Add column for light chain productive (kappa or lambda)
cells = rearr_info.index.get_level_values('cell_id').unique()
kappa_prod = set(cells[rearr_info.xs('IGK', level='locus').has_P.fillna(False)])
lambda_prod = set(cells[rearr_info.xs('IGL', level='locus').has_P.fillna(False)])

rearr_info['light_chain_prod'] = ['both' if x in kappa_prod and x in lambda_prod 
                                  else 'kappa' if x in kappa_prod 
                                  else 'lambda' if x in lambda_prod
                                  else 'neither' for x in rearr_info.index.get_level_values('cell_id')]

# Add column with highest read count, 2nd highest read count, and ratio between the two
rearr_info['largest_read_count'] = groups.read_count.max()
rearr_info['2nd_largest_read_count'] = df.groupby(['cell_id', 'locus']).read_count.agg(second_largest)
#rearr_info['ratio'] = rearr_info['largest_read_count'] / rearr_info['2nd_largest_read_count']

import ipdb; ipdb.set_trace()

# Add read counts and mapped reads from meta file
# rearr_info = (rearr_info
#               .merge(meta, on='cell_id', how='left')
#               .merge(df.set_index('cell_id').read_count, how='outer', right_index=True, left_index=True))



# Calculate RPM
# rearr_info['RPM'] = rearr_info['read_count']/(rearr_info['mapped_reads'] / 1e6)




# Add column with highest RPM, 2nd highest RPM, and ratio between the two
#rearr_info['largest_RPM'] = groups.RPM.max()
#rearr_info['2nd_largest_RPM'] = df.groupby(['cell_id', 'locus']).RPM.agg(second_largest)



# For each cell, tell us if there is a productive assignment for the heavy and the light chain
prod_heavy = rearr_info.xs('IGH', level='locus').has_P.fillna(False)
prod_light = ((rearr_info.xs('IGK', level='locus').has_P).fillna(False) 
              | (rearr_info.xs('IGL', level='locus').has_P).fillna(False))

# # Subset the data based on the `mode` values (L15: condition for subsetting cells)
# if not has_prod['heavy'] and has_prod['light']: 
#     cells = cells[~prod_heavy & prod_light]
# elif has_prod['heavy'] and not has_prod['light']:
#     cells = cells[prod_heavy & ~prod_light]
# elif has_prod['heavy'] and has_prod['light']:
#     cells = cells[prod_heavy & prod_light]
# elif not has_prod['heavy'] and not has_prod['light']:
#     cells = cells[~prod_heavy & ~prod_light]
# else:
#     raise ValueError('unknown parameter for has_prod')

rearr_info = rearr_info.loc[cells].reset_index()

#temp
# bad_cells = {c.strip() for f in ['../cellstoremove.txt', '../IGHignore.txt', '../IGKLignore.txt'] for c in open(f)}
# print(len(rearr_info.cell_id.unique()))
# rearr_info = rearr_info[~rearr_info.cell_id.isin(bad_cells)]
# print(len(rearr_info.cell_id.unique()))


# Plot the largest 2 read counts
#g = sns.relplot(data=rearr_info, x='largest_read_count', y='2nd_largest_read_count', hue='has_NP', col='locus', facet_kws=dict(sharey=False, sharex=False), palette= 'Set2').set(xscale="log").set(yscale="log")

tmp=rearr_info[rearr_info.has_NP==True]

g = sns.relplot(data=rearr_info, x='largest_read_count', y='mean_prod_SHM', col='locus', facet_kws=dict(sharey=False, sharex=False), palette='Set2', edgecolor=None, s=20).set(xscale="log")
#g.set(xlim=(1, 1e7))
#g.map(plt.plot, 'largest_read_count', 'largest_read_count', color='k', linewidth=1)
#g.set(xlabel='largest read count', ylabel='2nd largest read count')
#plt.subplots_adjust(right=0.9)


# # Save largest read counts plots as figures
# plt.savefig('troubleshooting/figures/{}_top2_read_counts.pdf'.format(fig_name))
# plt.savefig('troubleshooting/figures/{}_top2_read_counts_log.pdf'.format(fig_name))

# # Generate ratios for ratio plot
# tmp = tmp.sort_values(by=['locus', 'ratio'])
# counts = tmp.locus.value_counts()
# tmp['index'] = [i for l in tmp.locus.unique() for i in range(counts[l])]

# # Plot the ratios between largest 2 read counts & save plots as figures 
# g = sns.FacetGrid(data=tmp, col='locus', sharey=False, sharex=False).map(plt.scatter, 'index', 'ratio')
# g.set(yscale='log', ylabel='log ratio of largest 2 read counts', xlabel='sorted cell index')
# plt.savefig('troubleshooting/figures/{}_top2_read_counts_ratio_log.pdf'.format(fig_name))

plt.show()
import ipdb; ipdb.set_trace()


