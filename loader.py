"""
   Script with functions called from Main script (main.py)

"""

from pathlib import Path
import sys
from collections import Counter
from itertools import product

import pandas as pd



def load_dropped_cells(*paths):

    all_dropped_cells = set()
    for p in paths:
        dropped_cells = {line.strip() for line in open(p)}
        all_dropped_cells |= dropped_cells

    return all_dropped_cells


def process_raw_reads(path):
    
    mapped_IGs = pd.read_csv(path, sep='\t', header=None, usecols=[3, 5], 
                             names=['gene', 'read_count'], index_col='gene')

    # Extract V genes with read counts > 0
    mapped_Vs = mapped_IGs[(mapped_IGs.index.str.match('^IG[HLK][V]')) & 
                           (mapped_IGs.read_count > 0)]

    # Normalize reads to RPM
    PM_factor = mapped_Vs.read_count.sum()/1e6
    mapped_Vs.read_count /= PM_factor
    mapped_Vs.rename({'read_count': 'RPM'}, axis=1, inplace=True)

    # Round RPMs to integers
    mapped_Vs['RPM'] = mapped_Vs['RPM'].astype('int')

    return mapped_Vs 


def process_baldr_sonar(path):

    usecols=['cell_id', 'v_call', 'v_identity', 'duplicate_count', 'productive', 'locus',
             'cellBin.IGH', 'cellBin.IGK', 'cellBin.IGL']
    sonar_db = pd.read_csv(path, sep='\t', usecols=usecols)
    
    # Explode v_call gene list
    sonar_db.v_call = sonar_db.v_call.str.split(',')
    sonar_db = sonar_db.explode('v_call')

    # Reformat cell bin columns 
    sonar_db.set_index(sonar_db.columns[:-3].tolist(), inplace=True)    
    sonar_db.columns = sonar_db.columns.str.split('.', expand=True)
    sonar_db = sonar_db.stack(level=1).reset_index()
    sonar_db = sonar_db[sonar_db.locus == sonar_db.level_6].drop(columns='level_6')

    # Convert v_identity to SHM
    sonar_db['SHM'] = 100*(1-sonar_db.v_identity)
    
    # Convert v_call to gene by removing allele info
    sonar_db['gene'] = sonar_db.v_call.str.replace('\*\d+', '')

    return sonar_db.set_index('gene')


def agg_gene_duplicates(cell_id, sc_sonar): 

    sc_sonar = (sc_sonar
                .sort_values('productive')
                .groupby(['gene', 'cellBin', 'locus'])
                .agg({'productive': 'last', 
                      'duplicate_count': 'sum', 
                      'cell_id': 'first', 
                      'SHM': 'mean'})
                .reset_index(level=[1, 2]))

    return sc_sonar


def load_pseudogenes(path):

    pseudogenes = open(path).read().split('\n')
    
    return pseudogenes


def get_gene_rearrs(cell_id, mapped_Vs, sc_sonar, pseudogenes):

    # Merge mapped Vs for cell with corresponding cell sonar output
    merged_db = mapped_Vs.join(sc_sonar).fillna(value={'cell_id': cell_id})

    # Assign each gene a rearrangement type
    merged_db['rearrangement'] = [
        'pseudogene' if gene in pseudogenes
        else 'productive' if merged_db.loc[gene, 'productive']=='T'
        else 'passenger' for gene in merged_db.index
        ]
  
    # Group by rearrangement, report # of genes, total RPM and gene names
    grouped_IGH = (merged_db
                  .loc[merged_db.index.str.contains('IGH'), :] 
                  .reset_index()
                  .groupby('rearrangement')
                  .agg({'RPM': 'sum', 
                        'duplicate_count': 'sum', 
                        'SHM': 'mean', 
                        'rearrangement': len}) 
                  .rename(columns={'rearrangement': 'gene_count'}))

    grouped_IGK = (merged_db
                  .loc[merged_db.index.str.match('^IGKV'), :]
                  .reset_index()
                  .groupby('rearrangement')
                  .agg({'RPM': 'sum', 
                        'duplicate_count': 'sum', 
                        'SHM': 'mean', 
                        'rearrangement': len}) 
                  .rename(columns={'rearrangement': 'gene_count'}))

    grouped_IGL = (merged_db
                  .loc[merged_db.index.str.match('^IGLV'), :]
                  .reset_index()
                  .groupby('rearrangement')
                  .agg({'RPM': 'sum', 
                        'duplicate_count': 'sum', 
                        'SHM': 'mean', 
                        'rearrangement': len})
                  .rename(columns={'rearrangement': 'gene_count'}))

    gene_rearrs = pd.concat({'IGH': grouped_IGH.stack(), 
                             'IGK': grouped_IGK.stack(),
                             'IGL': grouped_IGL.stack()})

    return gene_rearrs


def clean_data(all_gene_rearrs):    

    # Asthetics
    cleaned_gene_rearrs = pd.DataFrame(all_gene_rearrs).T
    cleaned_gene_rearrs.index.name = 'cell_id'
    cleaned_gene_rearrs.columns.names = ['locus', 'rearrangement', 'variable']
    cleaned_gene_rearrs.stack('locus')
    
    return cleaned_gene_rearrs 


def append_meta(path, cleaned_gene_rearrs, sonar_db):

    meta = pd.read_csv(path, sep='\t', index_col='Cell')
    meta.index.name = 'cell_id'

    # Incorporate locus cell bins from SONAR-BALDR output
    sonar_db = sonar_db.reset_index().set_index(['cell_id', 'locus'])
    cleaned_gene_rearrs = (cleaned_gene_rearrs
                           .merge(sonar_db.cellBin, left_index=True, right_index=True))
    
    # Append merged gene rearrangements data to original meta
    appended_meta = meta.merge(cleaned_gene_rearrs, right_index=True, left_index=True)    

    return appended_meta 







