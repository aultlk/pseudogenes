# Loads data from 3 sources: BCR seq output from bedtools and SONAR, and pseudogenes list

from pathlib import Path
import sys

import pandas as pd

# reference to the path for this file, resolve it by converting to absolute path, and take the parent

def parse_bedtools_output(path, gene_loci="IGHV_only"):
    """
    path: Path to bedtools output
    gene_loci: Subsets bedtools output with either heavy IGV ("IGHV_only") or light IGV ("IGLV_only") 
    
    Returns: Productive vs. others read count, total read count
    """
    ig_gene_counts = pd.read_csv(path, sep='\t', header=None, usecols=[3, 5])
    # Should we remove duplicate genes (same chr, same position, same strand, same read count)?
    ig_gene_counts.drop_duplicates(inplace=True)
    ig_gene_counts.columns = ['gene_name', 'read_count']
    ig_gene_counts = ig_gene_counts.loc[ig_gene_counts.read_count > 0]

    if gene_loci == "IGHV_only":
        # only take heavy V genes 
        condition = ig_gene_counts.gene_name.str.startswith('IGHV')

    elif gene_loci == "IGLV_only":
        # only take light chain V genes
        condition = ig_gene_counts.gene_name.str.match('^IG[KL]V', regex=True)

    else:
        sys.exit(f'Unknown gene loci: {gene_loci}')

    # take gene_loci where condition == True 
    ig_gene_counts = ig_gene_counts.loc[condition]

    return ig_gene_counts 

def parse_sonar_output(sonar_output_file):
    """
    Loads data set from SONAR BCR seq output
    Returns gene name, productive/nonproductive, duplicate count, and SHM %
    """

    db = pd.read_csv(sonar_output_file, 
                     sep='\t', 
                     usecols=['productive', 'duplicate_count', 'v_identity', 'v_call'])
    db.v_call = db.v_call.str.split(",")
    db = db.explode('v_call')
    db = db.set_index('v_call')
    db['SHM'] = 100*(1-db.v_identity)
    
    return db.drop('v_identity', axis=1)

def parse_pseudogene_list():
    return

def transform_table_into_cool_info(my_necessary_arguments):
    cool_info = 0
    # do stuff
    return cool_info

