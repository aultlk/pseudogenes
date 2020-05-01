"""
Loads single cell bedtools mapped IG bed files; filters dropped cells and pulls in 
data from SONAR output file for single cells 
"""

from pathlib import Path
import sys

import pandas as pd

def parse_bedtools_output(path, gene_loci="all_IGV"):
    """
    path: Path to bedtools output
    gene_loci: Subsets bedtools output with only heavy IGV ("IGHV_only"), light 
               IGV ("IGLV_only"), or all loci ("all_IGV"); all loci is default
    Returns: IG gene names and read counts > 0 for all single cells; additionally,
             bins genes as productive/nonproductive and filters out single cells 
             that were removed for various reasons
    """
    mapped_IGs = pd.read_csv(path, sep='\t', header=None, usecols=[3, 5])
    
    # extract gene name and read counts
    mapped_IGs.columns = ['gene_name', 'read_count']
    # filter out read_counts < 0
    mapped_IGs = mapped_IGs.loc[mapped_IGs.read_count > 0]

    if gene_loci == "all_IGV":
        # take V genes from all loci
        condition = mapped_IGs.gene_name.str.match('^IG[HKL]V')
    
    elif gene_loci == "IGHV_only":
        # only take heavy V genes 
        condition = mapped_IGs.gene_name.str.startswith('IGHV')

    elif gene_loci == "IGLV_only":
        # only take light chain V genes
        condition = mapped_IGs.gene_name.str.match('^IG[KL]V')

    else:
        sys.exit(f'Unknown gene loci: {gene_loci}')
    # take gene loci where condition == True 
    mapped_IGs = mapped_IGs.loc[condition]

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

def transform_all_data(sc_id, sc_info, sonar_info, pseudogene_info):
    """    
    Returns: Data table with Productive vs. others read count, total read count
    """

    # Remove allele info
    sonar_info['gene'] = sonar_info.index.str.replace('\*\d+', '')
    # Extract sc_id only fron sonar info
    sonar_info = sonar_info.loc[sonar_info.cell_id == sc_id]

    # Subset the genes that are detected on the bed output
    ok_genes = set(sc_info['gene_name']).intersection(sonar_info['gene'])

    if sonar_info.size == 0:
        with open('../cells_not_in_sonar.csv', 'a') as handle:
            handle.write('{}\n'.format(sc_id))
        return 

    if not ok_genes:
        with open('../cells_in_sonar_wo_gene-hits.csv', 'a') as handle:
            gene_names_bedtools = ','.join(sc_info['gene_name'].tolist())
            gene_names_sonar = ','.join(sonar_info['gene'].tolist())
            handle.write('{}\t{}\t{}\n'.format(sc_id, gene_names_sonar, gene_names_bedtools))
        return
    
    # Set gene name as index, intersect gene data from sonar, and fill in NaN values for
    # cell_id and source_id for the given single cell
    sonar_info = (sonar_info
                  .set_index('gene')
                  .reindex(index=sc_info.gene_name)
                  .fillna(value={'cell_id': sc_id, 'source_id': sonar_info.source_id.iloc[0]}))

    # Create 'rearrangement' column and assign each gene (x) as productive, passenger, or pseudogene
    sonar_info['rearrangement'] = ['pseudogene' if gene in pseudogene_info
                                   else 'productive' if sonar_info.loc[gene, 'productive']=='T'
                                   else 'total_nonproductive' for gene in sonar_info.index]

    # Get bedtools read counts for genes
    sonar_info['bedtools_read_count'] = sc_info.set_index('gene_name').loc[:, 'read_count']

    # Group all loci (V), only heavy (VH), and only light (VL) by rearrangement  
##### Ask Chaim about 2 VH genes, should I still sum and take average??#####
    grouped_VH = (sonar_info
                  .loc[sonar_info.index.str.contains('IGH'), :] 
                  .groupby('rearrangement')
                  .agg({'bedtools_read_count': 'sum', 'duplicate_count': 'sum', 'SHM': 'mean'}))

    grouped_VL = (sonar_info
                  .loc[sonar_info.index.str.match('^IG[KL]V'), :]
                  .groupby('rearrangement')
                  .agg({'bedtools_read_count': 'sum', 'duplicate_count': 'sum', 'SHM': 'mean'}))
    
    grouped_V = (sonar_info
                 .groupby('rearrangement')
                 .agg({'bedtools_read_count': 'sum', 'duplicate_count': 'sum', 'SHM': 'mean'}))
    
    # Merged all grouped V data sets
    merged_data = pd.concat({'V': grouped_V.stack(),
                             'VL': grouped_VL.stack(),
                             'VH': grouped_VH.stack()})

    return merged_data

    

