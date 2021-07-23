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
        df = pd.read_csv(p, sep='\t', names=['v_call', 'clustered_vcall'], 
                         index_col='v_call', header=None)
        cluster_map.append(df)
    cluster_map = pd.concat(cluster_map, axis=0)
    cluster_map.index = cluster_map.index.str.strip()
    cluster_map.clustered_vcall = cluster_map.clustered_vcall.str.strip()

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
    baldr_sonar_data['clustered_vcall'] = cluster_map.loc[baldr_sonar_data.v_call, 'clustered_vcall'].values

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
    # Compute aggregated SHM per cell, locus and rearrangement type
    shm_values = sonar_data.groupby(['cell_id', 'locus', 'rearrangement_type']).apply(lambda x: x.SHM.drop_duplicates().mean())

    # Map SHM values to original dataframe
    sonar_data.SHM = shm_values[[tuple(info) for info in 
                                 zip(sonar_data.cell_id, sonar_data.locus, sonar_data.rearrangement_type)]].values
    # Drop duplicated entries for cell and clustered_vcall
    sonar_data = sonar_data[~sonar_data[['cell_id', 'clustered_vcall']].duplicated()]


    return (sonar_data, baldr_bins)





