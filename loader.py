"""
Loads single cell bedtools mapped IG bed files; filters dropped cells and pulls in 
data from SONAR output file for single cells 
"""

from pathlib import Path
import sys
from collections import Counter
from itertools import product

import pandas as pd


idx = pd.IndexSlice


def parse_bedtools_output(path):
    """
    path: Path to bedtools output

    Returns: Collects Immunoglobulin gene names and read counts > 0 for a 
             given single cell. Normalizes reads for sequencing depth by dividing
             by per million scaling factor to return reads per million (RPM). 
    """
    mapped_IGs = pd.read_csv(path, sep='\t', header=None, usecols=[3, 5])
    
    # Extract gene name and read counts
    mapped_IGs.columns = ['gene_name', 'read_count']

    # Filter out read_counts < 0
    mapped_IGs = mapped_IGs.loc[mapped_IGs.read_count > 0]
    
    # Normalize reads
    PM_factor = mapped_IGs.read_count.sum()/1e6
    mapped_IGs.read_count/=PM_factor
    mapped_IGs.rename({'read_count': 'RPM'}, axis=1, inplace=True)

    return mapped_IGs 


def parse_sonar_output(sonar_output_file):
    """
    Loads data set from SONAR BCR seq output
    Returns gene name, productive/nonproductive, duplicate count, and SHM %
    """

    db = pd.read_csv(sonar_output_file, sep='\t', 
                     usecols=['cell_id', 'source_id', 'v_call', 'duplicate_count', 
                              'v_identity', 'productive', 'status'])

    # remove the byte order mark in cell_id column ??
    db.cell_id = db.cell_id.str.replace('\ufeff', '')

    # only use rows that have a "good" status
    sonar_good = db.loc[:, 'status'] == 'good'
    db = db.loc[sonar_good]

    db.v_call = db.v_call.str.split(',')
    db = db.explode('v_call')
    db = db.set_index('v_call')
    db['SHM'] = 100*(1-db.v_identity)
    
    return db.drop('v_identity', axis=1)


def parse_pseudogene_list(filename):
    pseudogene_list = open(filename).read().split('\n')
    
    return pseudogene_list


def parse_dropped_cells(*files):
    all_dropped_cells = set()
    for f in files:
        dropped_cells = {line.strip() for line in open(f)}
        all_dropped_cells |= dropped_cells

    return all_dropped_cells


def transform_all_data(cell_id, mapped_IGs, sonar_info, pseudogene_info):
    """    
    Returns: Data table with Productive vs. others read count, total read count
    """

    # Remove allele info
    sonar_info['gene'] = sonar_info.index.str.replace('\*\d+', '')

    # Extract cell_id only fron sonar info
    sonar_info = sonar_info.loc[sonar_info.cell_id == cell_id]

    # Subset the genes that are detected on the bed output
    ok_genes = set(mapped_IGs['gene_name']).intersection(sonar_info['gene'])

    if sonar_info.size == 0:
        with open('../cells_not_in_sonar.csv', 'a') as handle:
            handle.write('{}\n'.format(cell_id))
        return 

    if not ok_genes:
        with open('../cells_in_sonar_wo_gene-hits.csv', 'a') as handle:
            gene_names_bedtools = ','.join(mapped_IGs['gene_name'].tolist())
            gene_names_sonar = ','.join(sonar_info['gene'].tolist())
            handle.write('{}\t{}\t{}\n'.format(cell_id, gene_names_sonar, gene_names_bedtools))
        return
    
    # Map gene data from sonar with single cell ids 
    sonar_info = (sonar_info
                  .set_index('gene')
                  .reindex(index=mapped_IGs.gene_name)
                  .fillna(value={'cell_id': cell_id, 'source_id': sonar_info.source_id.iloc[0]}))

    # Create 'rearrangement' 
    sonar_info['rearrangement'] = ['pseudogene' if gene in pseudogene_info
                                   else 'functional' if sonar_info.loc[gene, 'productive']=='T'
                                   else 'passenger' for gene in sonar_info.index]

    # Get bedtools read counts for genes
    sonar_info['bedtools_RPM'] = mapped_IGs.set_index('gene_name').loc[:, 'RPM']

    # Group by rearrangement, report # of assigned genes and list of gene names
    grouped_VH = (sonar_info
                  .loc[sonar_info.index.str.contains('IGH'), :] 
                  .reset_index()
                  .groupby('rearrangement')
                  .agg({'bedtools_RPM': 'sum', 
                        'duplicate_count': 'sum', 
                        'SHM': 'mean', 
                        'rearrangement': len, 
                        'gene_name': list})
                  .rename(columns={'rearrangement': 'gene_count'}))

    grouped_VK = (sonar_info
                  .loc[sonar_info.index.str.match('^IGKV'), :]
                  .reset_index()
                  .groupby('rearrangement')
                  .agg({'bedtools_RPM': 'sum', 
                        'duplicate_count': 'sum', 
                        'SHM': 'mean', 
                        'rearrangement': len, 
                        'gene_name': list})
                  .rename(columns={'rearrangement': 'gene_count'}))

    grouped_VL = (sonar_info
                  .loc[sonar_info.index.str.match('^IGLV'), :]
                  .reset_index()
                  .groupby('rearrangement')
                  .agg({'bedtools_RPM': 'sum', 
                        'duplicate_count': 'sum', 
                        'SHM': 'mean', 
                        'rearrangement': len, 
                        'gene_name': list})
                  .rename(columns={'rearrangement': 'gene_count'}))

    merged_data = pd.concat({'VH': grouped_VH.stack(), 
                             'VK': grouped_VK.stack(),
                             'VL': grouped_VL.stack()})

    return merged_data



def clean_data(all_sc_info):    

    # Create DataFrame with all merged data
    all_sc_info = pd.DataFrame(all_sc_info).T
    all_sc_info.index.name = 'cell_id'
    all_sc_info.columns.names = ['locus', 'rearrangement', 'variable']
    all_sc_info = (all_sc_info.stack(level=['locus', 'rearrangement'])
                   .reset_index()
                   .set_index('cell_id')
                   )

    # Set data types and round values
    all_sc_info.SHM = all_sc_info.SHM.astype(float)
    
    int_cols = ['bedtools_RPM', 'duplicate_count', 'gene_count']
    all_sc_info.loc[:, int_cols] = all_sc_info.loc[:, int_cols].round(0).astype(int)

    return all_sc_info



def append_meta(meta_file, cleaned_sc_info):
    """
       Subsets data from mapped raw reads dataframe and appends to meta table
    """

    meta_df = pd.read_csv(meta_file, sep='\t', index_col='cell_id')
    
    # Adjust column order
    col_order = ['rearrangement', 'locus', 'bedtools_RPM', 
                 'gene_name', 'gene_count', 'duplicate_count', 'SHM']
    cleaned_sc_info = cleaned_sc_info.loc[:, col_order]

    appended_meta = meta_df.merge(cleaned_sc_info, left_index=True, right_index=True)

    return appended_meta 

