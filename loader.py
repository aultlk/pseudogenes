"""
Script with functions called from Main script (main.py)
"""

from pathlib import Path
import sys
from itertools import product

import pandas as pd



def load_cluster_map(*paths):
    """
    Loads gene --> gene cluster map (Judah)
    Args:
    - paths: paths to IGHV, IGKV, IGLV gene --> gene cluster maps
    Returns:
    - concatenated data frame for all of these files 
    """

    cluster_map = []
    for p in paths:
        df = pd.read_csv(p, sep='\t', names=['allele', 'gene_cluster'], 
                         index_col='allele', header=None)
        cluster_map.append(df)
    cluster_map = pd.concat(cluster_map, axis=0)
    cluster_map.index = cluster_map.index.str.strip()
    cluster_map.gene_cluster = cluster_map.gene_cluster.str.strip()

    return cluster_map


def process_baldr_sonar(path, cluster_map):
    """
    Args:
    - path: path to SONAR+BALDR output table
    - cluster_map: gene cluster assignment (Judah)
    Returns:
    (sonar, baldr) separate dataframes. SONAR is indexed by (cell, gene) and BALDR by cell
    """

    baldr_sonar_data = pd.read_csv(path, sep='\t', index_col='cell_id') 

    # Add duplicate gene identifier to distinguish from ambiguous calls
    baldr_sonar_data.insert(0, 'duplicate_id', baldr_sonar_data.groupby(level=0).cumcount())
    # Process baldr-sonar data to extrat SHM for all clustered genes
    baldr_sonar_data['SHM'] = 100*(1-baldr_sonar_data.v_identity)
    # Explode v_call gene list
    baldr_sonar_data.v_call = baldr_sonar_data.v_call.str.split(',')
    baldr_sonar_data = baldr_sonar_data.explode('v_call')

    # Map v calls to cluster gene assignments
    baldr_sonar_data['clustered_vcall'] = cluster_map.loc[baldr_sonar_data.v_call, 'gene_cluster'].values

    # Drop cells with multiple rearrangement types for a single gene cluster
    ## Drop the cells with contradicting 'productive' and 'rearrangement_type' columns
    ## Trust the 'rearrangement_type' or the 'productive' column?? For now, assuming 'rearr_type' column
    bad_cells = (baldr_sonar_data.reset_index()
                 .groupby(['cell_id', 'clustered_vcall']).rearrangement_type
                 .filter(lambda x: len(set(x)) > 1))
    # Drop the 16 cells with duplicated clustered_vcalls and differing rearrangement type (at any locus)
    baldr_sonar_data.drop(baldr_sonar_data.index[bad_cells.index], inplace=True)




    # # Extract BALDR bins
    # baldr_bins = (baldr_sonar_data[['cellBin.IGH', 'cellBin.IGK', 'cellBin.IGL']]
    #               .reset_index()
    #               .drop_duplicates()
    #               .set_index('cell_id'))

    # baldr_bins.columns = baldr_bins.columns.str.split('.').str[1]
    # baldr_bins = baldr_bins.stack().reset_index()
    # baldr_bins.columns = ['cell_id', 'locus', 'baldr_bin']

    # Extract BALDR bins
    baldr_bins = (baldr_sonar_data[['cellBin.final']]
                  .reset_index()
                  .drop_duplicates()
                  .set_index('cell_id'))

    baldr_bins.columns = baldr_bins.columns.str.split('.').str[1]
    baldr_bins = baldr_bins.stack().reset_index()
    baldr_bins.columns = ['cell_id', 'bin_type', 'baldr_bin']



    
    ## Extract SONAR data
    sonar_data = baldr_sonar_data[['v_call', 'clustered_vcall', 'duplicate_id', 'locus', 'rearrangement_type', 'SHM']].reset_index()
    # Compute aggregated SHM for each cell, locus and rearrangement_type
    shm_values = sonar_data.groupby(['cell_id', 'locus', 'rearrangement_type']).apply(lambda x: x.SHM.drop_duplicates().mean())
    # Map SHM values to original dataframe
    sonar_data.SHM = shm_values[[tuple(info) for info in 
                                 zip(sonar_data.cell_id, sonar_data.locus, sonar_data.rearrangement_type)]].values
    # Drop duplicated entries for cell and clustered_vcall
    sonar_data = sonar_data[~sonar_data[['cell_id', 'clustered_vcall']].duplicated()]


    return (sonar_data, baldr_bins)



def combine_data(path, sonar_judah_data, baldr_bins):
    """
    Combines all data sets
    Args:
       path: path to single cell meta file
       sonar_judah_data: processed sonar data with judah gene clusters and read counts
       baldr_bins: baldr assigned locus bins per cell
    Returns:
       Combined pandas dataframe
    """

    cols = ['cell_name', 'phenotype', 'epitope', 'specificity', 'mapped_reads']
    meta = pd.read_csv(path, sep='\t', index_col='cell_name', usecols=cols)
    meta.index.name = 'cell_id'

    appended_meta = (
        sonar_judah_data.merge(baldr_bins, how='outer', on=['cell_id']) #'locus'
        .merge(meta, on='cell_id', how='left')
    )

    # Drop cells with no rearrangement type (undetected in baldr bins)
    appended_meta = appended_meta[appended_meta.rearrangement_type.notnull()]

    # Compute RPM
    appended_meta['RPM'] = appended_meta['read_count']/(appended_meta['mapped_reads'] / 1e6)

    # # Generate Per Million (PM) scaling factor for single cell total read count
    # pm_factor = appended_meta.groupby('cell_id').mapped_reads.first()/1e6
    # appended_meta['pm_factor'] = pm_factor.loc[appended_meta.index]

    # # Compute RPM
    # appended_meta['RPM'] = appended_meta['read_count']/appended_meta['pm_factor']
                  
    return appended_meta 







