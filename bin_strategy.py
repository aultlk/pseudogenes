"""
   This script contains code that can be used to bin loci and cells based on detected
   gene rearrangements from mapped raw reads. Can be integrated into loader.py if 
   we decide to compare binning strategies using BALDR and mapped raw reads. Additional
   code provided for plotting to compare these cell bins to BALDR-derived cell bins.

"""


FIG1 = f'{ROOT_DIR}/figures/cell_bin_counts.pdf'



def bin_cells(df):    


    def get_bin(rearr_serie):
        bin_assignment = 'not_detected'
        
        for bin_label in ['pseudogene', 'passenger', 'productive']:
            if (rearr_serie == bin_label).any():
                bin_assignment = bin_label
                break # if condition met, exit the loop
        return pd.Series(bin_assignment, index=rearr_serie.index)

    df_rpm_pos = df[df.RPM > 0]

    # Bin per locus 
    df_rpm_pos['bedtools_locus_bin'] = (
        df_rpm_pos.reset_index()
        .groupby(['cell_id', 'locus'])
        .rearrangement
        .transform(get_bin).tolist())

    # Bin per cell
    df_rpm_pos['bedtools_cell_bin'] = (df_rpm_pos.reset_index()
                              .groupby('cell_id')
                              .rearrangement
                              .transform(get_bin)
                              .tolist())

    return df_rpm_pos


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



# Load binned data and generate counts for plot 1, sort rearrangement types
binned_df = pd.read_csv(BINNED_DATA)

locus_bins = count_cell_bins(binned_df, 'bedtools_locus_bin', 'locus').swaplevel(0, 1).sort_index()
cell_bins = count_cell_bins(binned_df, 'bedtools_cell_bin').reset_index()
cell_bins.bedtools_cell_bin = pd.Categorical(cell_bins.bedtools_cell_bin, categories=rearr_order)
cell_bins = cell_bins.sort_values(by='bedtools_cell_bin')

sample_size = '{:,}'.format(len(binned_df.cell_id.unique()))
