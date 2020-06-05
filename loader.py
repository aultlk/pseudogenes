"""
   Script with functions called from Main script (main.py)

"""

from pathlib import Path
import sys
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
                             names=['gene', 'read_count'])

    # Extract V genes with read counts > 0
    mapped_Vs = mapped_IGs[(mapped_IGs.gene.str.match('^IG[HLK][V]')) & (mapped_IGs.read_count > 0)]

    # Designate gene locus
    mapped_Vs['locus'] = ['IGH' if gene.startswith('IGH')
                          else 'IGK' if gene.startswith('IGK')
                          else 'IGL' for gene in mapped_Vs.gene]

    return mapped_Vs


def process_baldr_sonar(path):

    sonar_db = pd.read_csv(path, sep='\t')
    
    # Explode v_call gene list
    sonar_db.v_call = sonar_db.v_call.str.split(',')
    sonar_db = sonar_db.explode('v_call')

    # Reformat cell bin columns 
    sonar_db.set_index(sonar_db.columns[:-4].tolist(), inplace=True)
    sonar_db.columns = sonar_db.columns.str.split('.', expand=True).rename('bin_locus', level=1)
    sonar_db = sonar_db.stack(level='bin_locus').reset_index()
    sonar_db = (sonar_db[sonar_db.locus == sonar_db.bin_locus]
                .drop(columns='bin_locus')
                .rename(columns={'cellBin': 'locus_bin'}))

    # Convert v_identity to SHM
    sonar_db['SHM'] = 100*(1-sonar_db.v_identity)
    
    # Convert v_call to gene by removing allele info
    sonar_db['gene'] = sonar_db.v_call.str.replace('\*\d+', '')

    return sonar_db


def agg_gene_duplicates(cell_id, sc_sonar): 

    sc_sonar = (sc_sonar
                .sort_values('productive')
                .groupby(['gene', 'locus', 'locus_bin'])
                .agg({'productive': 'last',
                      'duplicate_count': 'sum', 
                      'cell_id': 'first', 
                      'SHM': 'mean'})
                .reset_index())

    return sc_sonar


def load_pseudogenes(path):

    pseudogenes = open(path).read().split('\n')
    
    return pseudogenes


def get_gene_rearrs(cell_id, mapped_Vs, sc_sonar, pseudogenes):

    unmapped_sonar_genes = []

    # Fill in BALDR-SONAR locus bins for mapped Vs
    locus_bins = sc_sonar.loc[:, ['locus', 'locus_bin']].drop_duplicates().set_index('locus')
    mapped_Vs = mapped_Vs.join(locus_bins, on='locus').fillna('undetected').set_index('gene')

    # Merge mapped Vs with intersecting BALDR-SONAR assigned gene and traits
    sonar_genes =  (sc_sonar.loc[:, ['gene', 'productive', 'duplicate_count', 'SHM']]
                    .set_index('gene'))
    merged_db = mapped_Vs.join(sonar_genes)

    # Troubleshooting cells with BALDR-SONAR genes not detected by bedtools
    genes = []
    for gene in sonar_genes.index:
        if gene not in mapped_Vs.index:
            genes.append(gene)
    if genes:
        unmapped_sonar_genes.extend(genes)
        unmapped_sonar_genes.append(len(genes))

    # Assign each gene a rearrangement type
    merged_db['rearrangement'] = [
        'pseudogene' if gene in pseudogenes
        else 'productive' if merged_db.loc[gene, 'productive'] == 'T'
        else 'passenger' for gene in merged_db.index]

    # Group by rearrangement and collect counts
    gene_rearrs = (merged_db
                   .groupby(['rearrangement', 'locus', 'locus_bin'])
                   .agg({'rearrangement': len,
                         'duplicate_count': 'sum',
                         'read_count': 'sum',
                         'SHM': 'mean'})
                   .rename(columns={'rearrangement': 'gene_count'})
                   .stack())

    return (gene_rearrs, unmapped_sonar_genes)


def clean_data(all_gene_rearrs):    

    # Asthetics
    cleaned_gene_rearrs = pd.DataFrame(all_gene_rearrs).T
    cleaned_gene_rearrs.index.name = 'cell_id'
    cleaned_gene_rearrs = (cleaned_gene_rearrs
                           .stack(['locus', 'rearrangement', 'locus_bin'])
                           .reset_index()
                           .set_index('cell_id'))

    return cleaned_gene_rearrs 


def append_meta(path, cleaned_gene_rearrs):

    meta = pd.read_csv(path, sep='\t', index_col='Cell')
    meta.index.name = 'cell_id'

    # Append merged gene rearrangements data to original meta
    appended_meta = meta.merge(cleaned_gene_rearrs, right_index=True, left_index=True)    

    # Generate Per Million (PM) scaling factor from single cell mapped reads
    pm_factors = (appended_meta
                  .groupby('cell_id')
                  .agg({'mapped_reads': 'sum'})
                  .applymap(lambda x: x/1e6)
                  .rename(columns={'mapped_reads': 'PM_factor'}))
                  
    # Convert read counts to RPM by dividing by single cell PM factor
    for cell in appended_meta.index.unique():
        pm_factor = pm_factors.loc[cell, 'PM_factor']
        appended_meta.loc[cell, 'RPM'] = appended_meta.loc[cell, 'read_count']/pm_factor.astype(int)
         
    return appended_meta 







