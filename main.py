"""
   Main script that calls functions from loader (loader.py) for getting productive, passenger,
   and pseudogene gene rearrangements using single cell mapped raw reads. It also merges 
   BALDR-SONAR derived single cell data: cell bins (productive, passenger, pseudogene), 
   immunoglobulin, V gene assignments as well as merges final outputs with meta single
   cell data. 

"""

from pathlib import Path
from tqdm import tqdm

import sys 
import argparse
import pandas as pd

from loader import load_cluster_map
from loader import process_baldr_sonar
from loader import combine_data

ROOT_DIR = Path(__file__).resolve().parent.parent
SONAR_TABLE = Path(f"{ROOT_DIR}/sonarOutput_withCellBins_20210712.tsv")
PSEUDOGENE_FILE = Path(f"{ROOT_DIR}/pseud.HKL")
META_TABLE = Path(f"{ROOT_DIR}/cellInfo_allMaster.csv")
CLUSTERED_READS = Path(f"{ROOT_DIR}/clustered_gene_counts_202107.tsv") 
IGHV_ALLELES = Path(f"{ROOT_DIR}/IGHV_clusterLookup.txt") 
IGKV_ALLELES = Path(f"{ROOT_DIR}/IGKV_clusterLookup.txt") 
IGLV_ALLELES = Path(f"{ROOT_DIR}/IGLV_clusterLookup.txt")  


def main():

    # Load clustered reads from Judah pipeline
    clustered_reads = pd.read_csv(CLUSTERED_READS, sep='\t', header=None, index_col=False,
                                  names=['cell_id', 'locus', 'clustered_vcall', 'read_count'])

    # Load map of v genes to clustered v genes
    cluster_map = load_cluster_map(IGHV_ALLELES, IGKV_ALLELES, IGLV_ALLELES) 

    # Load and reformat sonar-baldr output table 
    (sonar_data, baldr_bins) = process_baldr_sonar(SONAR_TABLE, cluster_map)

    # Outer-join clustered reads with sonar-baldr output table
    merged_data = sonar_data.merge(clustered_reads, on=['cell_id', 'locus', 'clustered_vcall'], how='outer')
    merged_data.set_index('cell_id', inplace=True)

    # Assign locus
    merged_data['locus'] = ['IGH' if gene.startswith('IGH') 
                            else 'IGK' if gene.startswith('IGK') 
                            else 'IGL' 
                            for gene in merged_data.clustered_vcall]
    
    # Load pseudogenes list
    pseudogenes = {gene.strip() for gene in open(PSEUDOGENE_FILE)}

    ### Note: We never assign productive rearrangements, we only keep the ones in sonar ###
    # Assign rearrangement type
    merged_data['rearrangement_type'] = [
        rearr if pd.notnull(rearr) 
        else 'pseudogene' if gene in pseudogenes 
        else 'passenger'
        for (gene, rearr) in zip(merged_data.clustered_vcall, merged_data.rearrangement_type)
    ]


    # Load the cell metadata
    cols = ['cell_name', 'phenotype', 'epitope', 'specificity', 'mapped_reads']
    meta = pd.read_csv(META_TABLE, sep='\t', index_col='cell_name', usecols=cols)
    meta.index.name = 'cell_id'

    # Compute RPM with mapped reads from metadata table
    merged_data['RPM'] = 1e6 * merged_data.read_count / meta.loc[merged_data.index, "mapped_reads"]
    # Filter out clustered_vcalls with low RPMs (<100)
    merged_data = merged_data[(merged_data.RPM > 100) | (merged_data.rearrangement_type == 'productive')]

    # Aggregate rearrangement types, sum read counts, take average SHM, and count # of v call genes 
    merged_data = (
        merged_data
        .groupby(['cell_id', 'locus', 'rearrangement_type'])
        .agg(dict(
            read_count=sum,
            RPM=sum,
            SHM='mean',
            clustered_vcall=len
        )).rename(columns=dict(clustered_vcall='n_genes'))
    )

    # Add cell info from metadata
    merged_data = merged_data.merge(meta, left_index=True, right_index=True, how='left')

    return merged_data

if __name__ == '__main__':

    # Run pipeline
    merged_data = main()
    
    # Save outputs
    merged_data.to_csv(f'{ROOT_DIR}/appended_meta.csv')









    
