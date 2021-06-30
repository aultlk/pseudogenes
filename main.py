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
SONAR_TABLE = Path(f"{ROOT_DIR}/sonarOutput_withCellBins_20200707.tsv")
PSEUDOGENE_FILE = Path(f"{ROOT_DIR}/pseud.HKL")
META_TABLE = Path(f"{ROOT_DIR}/cellInfo_allMaster.csv")
CLUSTERED_READS = Path(f"{ROOT_DIR}/judah_clustered_gene_counts.tsv") 
IGHV_ALLELES = Path(f"{ROOT_DIR}/IGHV_clusterLookup.txt") 
IGKV_ALLELES = Path(f"{ROOT_DIR}/IGKV_clusterLookup.txt") 
IGLV_ALLELES = Path(f"{ROOT_DIR}/IGLV_clusterLookup.txt")  


def main():

    all_gene_rearrs = {}
    cells_w_gene_duplicates = {}
    unmapped_sonar_gene_cells = {}
    cells_not_in_sonar = []

    # Load all data sets
    pseudogenes = {gene.strip() for gene in open(PSEUDOGENE_FILE)}
    clustered_reads = pd.read_csv(CLUSTERED_READS, sep='\t')
    cluster_map = load_cluster_map(IGHV_ALLELES, IGKV_ALLELES, IGLV_ALLELES) 
    (sonar_data, baldr_bins) = process_baldr_sonar(SONAR_TABLE, cluster_map)

    # Outer-join gene_cluster read counts with sonar data
    merged_data = sonar_data.merge(clustered_reads, left_on=['cell_id', 'clustered_vcall'], 
                                   right_on=['cell_id', 'gene_cluster'], how='outer')
    merged_data.clustered_vcall = merged_data.clustered_vcall.fillna(merged_data.gene_cluster)

    # Assign locus to gene cluster
    merged_data['locus'] = ['IGH' if gene.startswith('IGH') else 'IGK' if gene.startswith('IGK') else 'IGL' 
                            for gene in merged_data.clustered_vcall]
    merged_data = merged_data.set_index('cell_id').drop(columns=['gene_cluster'])
    
    # Adding the meta information
    cols = ['cell_name', 'phenotype', 'epitope', 'specificity', 'mapped_reads']
    meta = pd.read_csv(META_TABLE, sep='\t', index_col='cell_name', usecols=cols)
    meta.index.name = 'cell_id'
    merged_data = merged_data.merge(meta, left_index=True, right_index=True, how='left')

    

    # Drop cells with a max read count < 100; make sure this doesn't include cells with productive assignments 
    # merged_data = merged_data.groupby(level='cell_id').filter(lambda x: x.read_count.max() > 100)

    # Assign large read counts to productive for cells with no sonar/baldr productive bin

    # Add pseudogene information to sonar data. Check w Chaim if I should only be taking the 'productive' classification from sonar
    merged_data['rearrangement_type'] = [
        rearr if pd.notnull(rearr) 
        else 'pseudogene' if gene in pseudogenes 
        else 'passenger'
        for (gene, rearr) in zip(merged_data.clustered_vcall, merged_data.rearrangement_type)]

    # Aggregate rearrangement types, sum read counts and count # of v call genes aggregated
    merged_data = merged_data.groupby(['cell_id', 'locus', 'rearrangement_type']).agg(dict(
        SHM='first',
        read_count='sum',
        clustered_vcall=len
    )).rename(columns=dict(clustered_vcall='n_genes')).reset_index()

    # Merge data with meta cell information and baldr locus bins
    final_data = combine_data(META_TABLE, merged_data, baldr_bins) 
    


    # cleaned_gene_rearrs = clean_data(all_gene_rearrs) 
    # appended_meta = append_meta(META_TABLE, cleaned_gene_rearrs)


    # cells_w_unmapped_sonar_genes = pd.DataFrame.from_dict(unmapped_sonar_gene_cells, orient='index')
    # cells_not_in_sonar = pd.DataFrame(cells_not_in_sonar)

    return final_data

if __name__ == '__main__':

    # Run pipeline
    (final_data) = main()
    
    # Save outputs
    final_data.to_csv(f'{ROOT_DIR}/appended_meta.csv')
    # cells_w_unmapped_sonar_genes.to_csv(f'{ROOT_DIR}/cells_w_unmapped_sonar_genes.csv')
    # cells_not_in_sonar.to_csv(f'{ROOT_DIR}/cells_not_in_sonar.csv')









    
