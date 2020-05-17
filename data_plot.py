"""
This script plots all merged data from data_process.py output.  

"""

from pathlib import Path
import pandas as pd

import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib.ticker as tkr

ROOT_DIR = Path(__file__).resolve().parent.parent
LOCUS_DISTANCES = f'{ROOT_DIR}/locus_distances.csv'
APPENDED_META = f'{ROOT_DIR}/appended_meta.csv'

FIG1 = f'{ROOT_DIR}/figures/locus_distances.pdf' 
FIG2 = f'{ROOT_DIR}/figures/reads_all.pdf' 
FIG3 = f'{ROOT_DIR}/figures/reads_subsetted.pdf' 
FIG4 = f'{ROOT_DIR}/figures/distance_correlations.pdf'


def plot1(locus_df):
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
    


def plot2(meta_df):
    """
    Returns: Plot 2 = Rearrangement vs. RPM for all cells
    """
       
    sns.set(style="ticks", font_scale=2)    

    # Plot Figure 2 
    fig, ax = plt.subplots(figsize=(15, 8))
    plot2 = sns.boxplot(x='locus', y='bedtools_RPM', hue='rearrangement', ax=ax, 
                        palette="husl", fliersize=1, data=meta_df)

    # Format
    plot2.set(yscale='log', ylim=(0, 10**7), ylabel='RPM (log10)', xlabel='Locus')
    plot2.set_title(label='Rearrangement vs. Total Mapped Raw Reads\n', fontweight='bold')
    plt.subplots_adjust(right=0.8)
    plt.legend(loc='center left', bbox_to_anchor=(1.01, 0.5), fontsize=15)


    plt.savefig(FIG2)

    plt.show()



def plot3(meta_df):
    """
    Returns: Plot 3 = Rearrangement vs. RPM for subsetted cells (epitope-specificity and
                      phenotype)
    """

    sns.set(style='ticks')    

    # Plot Figure 3 
    plot3 = sns.catplot(kind='box', x='locus', y='bedtools_RPM', col='phenotype', 
                         row='specificity', margin_titles=True, hue='rearrangement',
                         palette="husl", fliersize=1, height=4, legend=False, data=meta_df)

    # Format
    plot3.set(yscale='log', ylim=(0, 10**7)).set_axis_labels('Locus', 'RPM (log10)')
    plot3.fig.suptitle('Rearrangement vs. Reads per Million (RPM)', fontweight='bold')

    plt.subplots_adjust(wspace=0.2, top=0.9, hspace=0.3, right=0.85)
    plt.legend(loc='center left', bbox_to_anchor=(1.2, 0.5))


    plt.savefig(FIG3)

    plt.show()


def plot4(locus_df):
    """
    Returns: Plot 4 = Interlocus Distance & SHM Correlations of Rearrangements
    """

    # Create correlation matrices
    corrs = locus_df.loc[:, :].groupby('rearrangement').apply(lambda x: x.corr())
    
    sns.set(style='ticks')
    
    fig, axs = plt.subplots(2)

    ax1 = sns.heatmap(cmap="YlGnBu", vmin=-1, annot=True, data=corrs.loc['passenger'])
    ax2 = sns.heatmap(cmap="YlGnBu", vmin=-1, annot=True, data=corrs.loc['pseudogene'])

    axs[0].plot4a
    axs[1].plot4b

#    plt.savefig(FIG4)
#    fig.suptitle('Vertically stacked subplots')
  
    plt.show()

    
def plot5(meta_df):


    # Plot data
    plot4 = sns.heatmap(cmap="YlGnBu", vmin=-1, annot=True, data=corrs.loc['functional'])


#    plt.savefig(f'{ROOT_DIR}/figures/diag_correlations.pdf')

    plt.show()

def plot6(meta_df):
    return




cols = ['bedtools_RPM', 'gene_count', 'duplicate_count', 'SHM']
# Load data for plotting
meta_df = pd.read_csv(APPENDED_META, index_col='cell_id')
locus_df = pd.read_csv(LOCUS_DISTANCES, index_col='cell_id')


#plot1 = plot1(locus_df)
#plot2 = plot2(meta_df)
#plot3 = plot3(meta_df)
plot4 = plot4(locus_df)









