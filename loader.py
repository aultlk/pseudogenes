# Loads data from 3 sources: BCR seq output from bedtools and SONAR, and pseudogenes list

from pathlib import Path
import sys

import pandas as pd

# reference to the path for this file, resolve it by converting to absolute path, and take the parent

def parse_bedtools_output(path, gene_loci="all_IGV"):
    """
    path: Path to bedtools output
    gene_loci: Subsets bedtools output with only heavy IGV ("IGHV_only"), light IGV ("IGLV_only"), or all loci ("all_IGV"); 
               all loci is default
    Returns: Productive vs. others read count, total read count
    """
    ig_gene_counts = pd.read_csv(path, sep='\t', header=None, usecols=[3, 5])

    # ig_gene_counts.drop_duplicates(inplace=True)
    ig_gene_counts.columns = ['gene_name', 'read_count']
    ig_gene_counts = ig_gene_counts.loc[ig_gene_counts.read_count > 0]

    if gene_loci == "all_IGV":
        # take V genes from all loci
        condition = ig_gene_counts.gene_name.str.match('^IG[HKL]V')
    
    elif gene_loci == "IGHV_only":
        # only take heavy V genes 
        condition = ig_gene_counts.gene_name.str.startswith('IGHV')

    elif gene_loci == "IGLV_only":
        # only take light chain V genes
        condition = ig_gene_counts.gene_name.str.match('^IG[KL]V')

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

    db = pd.read_csv(sonar_output_file, sep='\t')

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

def transform_all_data(sc_id, sc_info, sonar_info, pseudogene_info):
    """    
    Returns: Data table with Productive vs. others read count, total read count
    """

    # to remove
    sonar_info['gene'] = sonar_info.index.str.replace('\*\d+', '')

    (index, lane, flow_id, _) = sc_id

    sonar_glob_ids = sonar_info.source_id.str.split('-', expand=True)
    sonar_index = sonar_glob_ids[3].str.replace('\.IG.*', '').str.replace('.*index_', '')
    sonar_lane = sonar_glob_ids[2]
    sonar_flowid = sonar_glob_ids[1]

    sonar_info = sonar_info[
        (sonar_flowid == flow_id) 
        & (sonar_lane == lane)
        & (sonar_index == index)
    ]

    ok_genes = set(sc_info['gene_name']).intersection(sonar_info['gene'])
    print(sonar_info.shape)

    if ok_genes:
        for gene in ok_genes:
            if sum(sonar_info['gene'] == gene) > 5:
                sonar_scid = sonar_info.source_id
                print(sonar_scid, sum(sonar_info['gene'] == gene))
                import ipdb; ipdb.set_trace()
    # import ipdb; ipdb.set_trace()

    # do stuff
    return []

