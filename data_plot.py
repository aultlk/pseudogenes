"""
   This script is for generating test plots on all merged data from data_process.py output.
   Finalized plots will be incorporated into main.py program.  

"""

from pathlib import Path
from statannot import add_stat_annotation
from scipy.stats import mannwhitneyu
from itertools import product, combinations

import pandas as pd

import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib.ticker as tkr


ROOT_DIR = Path(__file__).resolve().parent.parent
APPENDED_META = f'{ROOT_DIR}/appended_meta.csv'
#BINNED_DATA = f'{ROOT_DIR}/final_mapped_read_gene_rearrangements.csv'

FIG1 = f'{ROOT_DIR}/figures/cell_bin_counts.pdf' 
FIG2 = f'{ROOT_DIR}/figures/reads_all.pdf' 
FIG3 = f'{ROOT_DIR}/figures/reads_cell_bins.pdf' 
FIG3B = f'{ROOT_DIR}/figures/gene_counts_cell_bins1.pdf' 
FIG4 = f'{ROOT_DIR}/figures/gene_counts_cell_bins2.pdf' 
FIG5 = f'{ROOT_DIR}/figures/reads_subsetted.pdf' 
FIG7 = f'{ROOT_DIR}/figures/distance_correlations.pdf'

meta_df = pd.read_csv(APPENDED_META)


def count_cell_bins(df, *groups):
        
    counts = df.groupby(list(groups)).cell_id.agg(lambda x: len(set(x)))
    data = pd.DataFrame(counts).rename(columns={'cell_id': 'count'})

    if 'productive' not in data.index:
        data.loc['productive'] = 0
        
    return data
 
   

def plot1(df1, df2, var1):
    """
        df1 = locus bins
        df2 = cell bins
        var1 = sample size

        Returns: Plot 1 = Nested pie chart with locus bins and bar plot with cell bins
                          Bin types assignments based on bedtools assignments
    """
    
    # Figure set up
    sns.set(style='ticks', font_scale=1.5) 
    fig, axes = plt.subplots(1, 2, figsize=(15, 8))

    # Set up nested pie chart and outer/inner wedge groups    
    cmap = plt.get_cmap("tab20c")
    outer_colors = cmap(np.arange(3)*4)
    inner_colors = cmap(np.array([1, 2, 5, 6, 9, 10]))
    size=0.3
    outer_groups = df1.groupby('locus').sum()
    
    ############################################################################################
    # Drop small passenger group for now because percentage negligible for pie chart - fix later
    inner_groups = df1.loc[(slice(None), ['passenger', 'pseudogene']), :].droplevel('locus')

    axes[0].pie(x=outer_groups['count'], labels=outer_groups.index, radius=1, colors=outer_colors)
    axes[0].pie(x=inner_groups['count'], labels=inner_groups.index, autopct='%.0f%%',
                colors=inner_colors, radius=1-size, pctdistance=0.8, labeldistance=None)

    axes[0].set(aspect="equal", title='Locus Bins')


    # Plot bar plot
    colors = sns.color_palette(['salmon', 'limegreen', 'dodgerblue'])

    df2 = df2.reset_index()

    sns.barplot(x='bedtools_cell_bin', y='count', ax=axes[1], palette=colors, data=df2)

    # Format figure
    axes[1].set(ylabel='Number of Cells', xlabel='Category')
    plt.text(2.75, 500, f'n = {sample_size}', fontsize=14)
    plt.subplots_adjust(right=0.8, wspace = 0.4, left = 0.025)

    plt.savefig(FIG1)

    plt.show()    


def plot2(meta_df):
    """
    Returns: Plot 2 = Rearrangement vs. RPM for all cells
    """
       

    sns.set(style="ticks", font_scale=2)    

 
    # Plot Figure 2 
    fig, ax = plt.subplots(figsize=(15, 8))
    g = sns.boxplot(x='locus', y='RPM', hue='rearrangement', ax=ax, 
                        palette="husl", fliersize=1, data=meta_df)

    # Format
    g.set(yscale='log', ylim=(0, 10**7), ylabel='RPM (log10)', xlabel='Locus')
    g.set_title(label='Rearrangement vs. Total Mapped Raw Reads\n', fontweight='bold')
    plt.subplots_adjust(right=0.8)
    plt.legend(loc='center left', bbox_to_anchor=(1.01, 0.5), fontsize=15)


    plt.savefig(FIG2)

    plt.show()


def get_stats(meta_df):

    # Group RPM values and drop unassigned baldr bins labeled 'none'
    grouped = (meta_df
               .groupby(['cellBin', 'locus', 'rearrangement', 'cell_id'])
               .agg({'RPM': 'sum'})
               .drop('none', level=0)
               .reset_index()
               .drop(columns='cell_id'))

    # Filter out Nan RPM values and get log10 values
    grouped.drop(grouped[grouped.RPM.isnull()].index, inplace=True)
    grouped['log_RPM'] = np.log10(grouped['RPM'])

    # Collect stats for each gene rearr type per locus
    stats = []
    loci = grouped.locus.unique()
    gene_rearrs = grouped.rearrangement.unique()
    cellbins = grouped.cellBin.unique()
    
    for (locus, gene_rearr) in product(loci, gene_rearrs):
        df = grouped[(grouped.locus == locus) & (grouped.rearrangement == gene_rearr)]

        for (cellbin1, cellbin2) in combinations(cellbins, 2):
            stat = mannwhitneyu(df[df.cellBin == cellbin1].log_RPM, 
                                df[df.cellBin == cellbin2].log_RPM)
            stats.append([gene_rearr, locus, cellbin1, cellbin2, stat.pvalue])
            
    stats = pd.DataFrame(stats, columns=['rearrangement', 'locus', 'cellbin1', 'cellbin2', 'pvalue'])
    stats.set_index(stats.columns[: -1].tolist(), inplace=True)

    return stats




def plot3(meta_df, stats):
    """
    Returns: Plot 3 = RPM vs. Cell Bin Rearrangement Types (based on baldr assignments)
    """

    meta_df.drop(index=meta_df[meta_df.cellBin == 'none'].index, inplace=True)

    sns.set(style='ticks')    


    # Plot figure
    g = sns.catplot(kind='box', x='locus', y='RPM', col='rearrangement', hue='cellBin', 
                    palette="husl", fliersize=0.5, legend=False, data=meta_df)

    # Format
    g.set(yscale='log', ylim=(0, 10**7)).set_axis_labels('Locus', 'RPM (log10)')

    titles = ['Productive V Reads', 'Passenger V Reads', 'Pseudogene V Reads']
    axes = g.axes.flatten()
    for (i, title) in enumerate(titles):
        axes[i].set_title(title)

    plt.subplots_adjust(wspace=0.2, top=0.9, hspace=0.3, right=0.85)
    plt.legend(title='BALDR Cell Bins', loc='center left', bbox_to_anchor=(1.2, 0.5))
    

    # Annotate statistical significance (Mann Whitney U test)
    box_pairs = [
        (('IGH', 'passenger'), ('IGH', 'productive')),
        (('IGH', 'passenger'), ('IGH', 'pseudogene')),
        (('IGH', 'productive'), ('IGH', 'pseudogene')),
        (('IGK', 'passenger'), ('IGK', 'productive')),
        (('IGK', 'passenger'), ('IGK', 'pseudogene')),
        (('IGK', 'productive'), ('IGK', 'pseudogene')),
        (('IGL', 'passenger'), ('IGL', 'productive')),
        (('IGL', 'passenger'), ('IGL', 'pseudogene')),
        (('IGL', 'productive'), ('IGL', 'pseudogene'))
        ]


    gene_rearrs = stats.index.get_level_values(0).unique()
    for (ax, gene_rearr) in zip(axes, gene_rearrs):
        add_stat_annotation(ax=ax, data=meta_df, x='locus', y='RPM', hue='cellBin', 
                            perform_stat_test=False, box_pairs=box_pairs, 
                            pvalues=stats.loc[gene_rearr].pvalue.tolist(),
                            )
                   

    plt.savefig(FIG3)

    
    plt.show()



def plot3b(meta_df):
    """
         Returns: Plot 3 = Number of detected genes vs. Cell Bin Rearrangement Types 
                           (based on baldr assignments)
    """

    sns.set(style="ticks")    
    
    meta_df.drop(index=meta_df[meta_df.cellBin == 'none'].index, inplace=True)

    # Plot figure
    g = sns.catplot(kind='box', x='locus', y='gene_count', col='rearrangement', hue='cellBin',
                    palette="husl", fliersize=0.5, legend=False, data=meta_df)
    
    # Format
    g.set_axis_labels('Locus', 'Number of Genes')
    g.fig.suptitle('Number of Genes Detected from Mapped Raw Reads in Baldr Cell Bins\n', 
                   fontweight='bold')
    plt.subplots_adjust(top=0.8, right=0.85)
    plt.legend(title='BALDR Cell Bins', loc='center left', bbox_to_anchor=(1.1, 0.5))
    
    titles = ['Productive V Genes', 'Passenger V Genes', 'Pseudogene V Genes']
    for (i, title) in enumerate(titles):
        g.axes[i].set_title(title)


    plt.savefig(FIG3B)
    
    plt.show()




def plot4(meta_df):
    """
    Returns: Plot 4 = Gene Count vs. Cell Bin Rearrangement Types (baldr assignment)
    """
       
    sns.set(style="ticks")    

    # Drop undetected genes and cell bins with none
    drop_rows = meta_df[(meta_df.rearrangement == 'not detected') | (meta_df.cellBin == 'none')]
    data = meta_df.drop(index=drop_rows.index)
    
    g = sns.catplot(kind='violin', x='locus', y='gene_count', hue='cellBin',
                        palette="husl", fliersize=1, legend=False, height=6, aspect=3, data=data)
    
    # Format
    g.set_axis_labels('Locus', 'Number of Mapped Genes')
    g.fig.suptitle('Number of Mapped Genes vs. Locus Bin\n', fontweight='bold')
    plt.subplots_adjust(top=0.8, right=0.70)
    plt.legend(title='Locus Bin', loc='center left', bbox_to_anchor=(1.1, 0.5))
    
    
    plt.savefig(FIG4)
    
    plt.show()


def plot5(meta_df):
    """
    Returns: Plot 5 = Rearrangement vs. RPM for subsetted cells (epitope-specificity and
                      phenotype)
    """

    sns.set(style='ticks', font_scale=1.08)    

    # Plot Figure 3 
    g = sns.catplot(kind='box', x='locus', y='RPM', col='phenotype', 
                         row='specificity', margin_titles=True, hue='rearrangement',
                         palette="husl", fliersize=1, height=4, legend=False, data=meta_df)

    # Format
    g.set(yscale='log', ylim=(0, 10**7)).set_axis_labels('Locus', 'RPM (log10)')
    g.fig.suptitle('Rearrangement vs. Reads per Million (RPM)', fontweight='bold')

    plt.subplots_adjust(wspace=0.2, top=0.9, hspace=0.3, right=0.85)
    plt.legend(loc='center left', bbox_to_anchor=(1.2, 0.5))


    plt.savefig(FIG5)

    plt.show()


def plot6(meta_df):
    """
    Returns: Plot 6 = Plotted correlation between productive reads and duplicate counts
    """

    # Subset productive rearrangements
    productive = meta_df[meta_df.rearrangement != 'pseudogene']

    ax = sns.scatterplot(x='RPM', y='duplicate_count', data=productive)       
    
#    plt.savefig(FIG6)
  
    plt.show()




# Load meta data and sort rearrangement types for plots


rearr_order = ['productive', 'passenger', 'pseudogene']
meta_df.rearrangement = pd.Categorical(meta_df.rearrangement, categories=rearr_order)
meta_df = meta_df.sort_values(by='rearrangement')


# Load binned data and generate counts for plot 1, sort rearrangement types
# binned_df = pd.read_csv(BINNED_DATA)

# locus_bins = count_cell_bins(binned_df, 'bedtools_locus_bin', 'locus').swaplevel(0, 1).sort_index()
# cell_bins = count_cell_bins(binned_df, 'bedtools_cell_bin').reset_index()
# cell_bins.bedtools_cell_bin = pd.Categorical(cell_bins.bedtools_cell_bin, categories=rearr_order)
# cell_bins = cell_bins.sort_values(by='bedtools_cell_bin')

# sample_size = '{:,}'.format(len(binned_df.cell_id.unique()))

stats = get_stats(meta_df)
# # Generate plots
# plot1 = plot1(locus_bins, cell_bins, sample_size)
# plot2 = plot2(meta_df)
plot3 = plot3(meta_df, stats)
plot3b = plot3b(meta_df)
plot4 = plot4(meta_df)
# plot5 = plot5(meta_df)        
# plot6 = plot6(meta_df)








