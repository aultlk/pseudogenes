"""
   This script loads the merged data set with read counts from Judah's pipeline and sonar-baldr
   gene assignments and rearrangement types. Then it plots the top 2 highest read counts as well
   as ratios. Purpose is to use these plots for developing a criteria for manually assigning 
   cell bins for cells that clearly have a nonproductive rearrangement in addition to productive.
"""

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
from scipy import stats

# Choose mode: conditions to subset cells on
# NoProd: No productive rearrangement assigned by baldrsonar
# HasProd: Productive rearrangement assigned by baldrsonar 
#has_prod = dict(heavy=True, light=True)

def second_largest(l):
    """
    Computes 2nd largest read count.
    If len(read_counts) == 1, returns 0
    """
    if len(l) == 1:
        return 0
    return sorted(l)[-2]

def third_largest(l):
    """
    Computes 3rd largest read count.
    If len(read_counts) == 1 or 2, returns 0
    """
    if len(l) == 1 or len(l) == 2:
        return 0
    return sorted(l)[-3]

def remain_sum(l):
    """
    Computes sum of RPM/read counts
    after substracting the largest read count/RPM
    """
    if len(l) == 1:
        return 0
    return l.dropna().sort_values()[:-1].sum()

def annotate(data, x=None, y=None, logx=False, **kws):
    """
    Annotates graph with statistical values for 
    correlation and significance. Also, plots
    a trendline.
    """
    df = data.dropna(subset=[x, y], how='any').sort_values(by=x)

    if logx:
        df = df[df[x] > 0]
        x_ = np.log(df[x])

    else:
        x_ = df[x]
    y_ = df[y]

    slope, intercept, r_value, p_value, std_err = stats.linregress(x_, y_)

    # Display correlation coefficient
    ax = plt.gca()
    ax.text(.05, .8, f'y={slope:.2g}x + {intercept:.2g} r={r_value:.2f}, p={p_value:.2g}',
            transform=ax.transAxes)

    # plot trendline
    if logx:
        ax.plot(np.exp(x_), intercept + slope*x_, 'r--')
    else:
        ax.plot(x_, intercept + slope*x_, 'r--')

def filter_RPM(df, col, n):
    """
    df = df to filter
    col = column to subset in df
    n = minimum value to subset on
    
    returns filtered dataframe
    """
    filtered_df = df[df[col] > n]

    return filtered_df[col]


df = pd.read_csv('troubleshooting/merged_data.csv')
df.locus = pd.Categorical(df.locus, ['IGH', 'IGK', 'IGL'])

# Filter out RPM < 100 unless it's productive
df = df[(df.RPM > 100) | (df.rearrangement_type == 'productive')]

# Generate df with info about cells to be filtered with > 3 gene cluster bins
many_bins = df.groupby(['cell_id', 'locus']).filter(lambda x: len(x) > 3)
bin_counts = many_bins.groupby(['cell_id', 'locus']).clustered_vcall.agg(lambda x: len(x)).fillna(0).astype(int)
many_bins['bins'] = bin_counts.loc[[(cell_id, locus) for (cell_id, locus) in zip(many_bins.cell_id, many_bins.locus)]].values

# Filter out cells with > 3 gene cluster bin
df = df.groupby(['cell_id', 'locus']).filter(lambda x: len(x) <= 3)

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

# Add column with largest read count, 2nd largest read count, and ratio between the two
rearr_info['largest_read_count'] = groups.read_count.max()
rearr_info['2nd_largest_read_count'] = df.groupby(['cell_id', 'locus']).read_count.agg(second_largest)
# rearr_info['ratio'] = rearr_info['largest_read_count'] / rearr_info['2nd_largest_read_count']

# Add column with largest RPM and 2nd largest RPM
rearr_info['largest_RPM'] = groups.RPM.max()
rearr_info['2nd_largest_RPM'] = groups.RPM.agg(second_largest)

# Add column for total np-RPM, which sums the remaining read counts/ RPM (top read count = productive)
rearr_info['remaining_RPM'] = groups.RPM.agg(remain_sum)

rearr_info['3rd_largest_RPM'] = df.groupby(['cell_id', 'locus']).RPM.agg(third_largest)

# For each cell, tell us if there is a productive assignment for the heavy and the light chain
prod_heavy = rearr_info.xs('IGH', level='locus').has_P.fillna(False)
prod_light = ((rearr_info.xs('IGK', level='locus').has_P).fillna(False) 
              | (rearr_info.xs('IGL', level='locus').has_P).fillna(False))

# Subset the data based on the `mode` values (L15: condition for subsetting cells)
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
# rearr_info = rearr_info.loc[cells].reset_index()

# Filter out dropped cells
# bad_cells = {c.strip() for f in ['../cellstoremove.txt', '../IGHignore.txt', '../IGKLignore.txt'] for c in open(f)}
# rearr_info = rearr_info[~rearr_info.cell_id.isin(bad_cells)]

# Plot scatter plot correlations 
tmp = rearr_info.reset_index()

g = sns.relplot(
    data=tmp, x='remaining_RPM', y='mean_prod_SHM', col='locus', palette='Set2', kind='scatter', 
    edgecolor=None, color='k', s=15, facet_kws=dict(sharey=False, sharex=False)
).set(xscale="log")
#g.set(xlim=(1,1e5))
g.map_dataframe(annotate, data=tmp, x='remaining_RPM', y='mean_prod_SHM', logx=True)
g.set(xlabel='remaining RPM', ylabel='mean productive SHM')


# # Plot 2 largest read counts/RPM; hue based on NP rearr detected or productive light chain assignment 
# tmp=rearr_info[rearr_info.light_chain_prod != 'neither'].reset_index()
# g = sns.relplot(
#     data=tmp, x='largest_RPM', y='remaining_RPM', hue='has_NP', col='locus', palette='Set2', edgecolor=None, s=10, 
#     facet_kws=dict(sharey=False, sharex=False)
# ).set(xscale="log", yscale="log")
# g.map(plt.plot, 'largest_RPM', 'largest_RPM', color='k', linewidth=1)
# g.set(xlabel='largest RPM', ylabel='remaining RPM')
# plt.subplots_adjust(right=0.9)




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


