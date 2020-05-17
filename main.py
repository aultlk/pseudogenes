"""
Pull mapped IG genes and number of reads from `bedtools muticov` output
Merge with master/meta.tsv while also pulling data from filtering2_rearrangements_single-cell.tsv file
"""

from pathlib import Path
import sys 

import pandas as pd

from loader import parse_bedtools_output 
from loader import parse_sonar_output
from loader import parse_pseudogene_list
from loader import parse_dropped_cells
from loader import transform_all_data 
from loader import clean_data
from loader import append_meta

ROOT_DIR = Path(__file__).resolve().parent.parent
SONAR_OUTPUT_FILE = Path(f"{ROOT_DIR}/filtering2_rearrangements_single-cell.tsv")
PSEUDOGENE_FILE = Path(f"{ROOT_DIR}/pseud.HKL")
CELLS_TO_REMOVE_FILE = Path(f"{ROOT_DIR}/cellstoremove.txt")
IGH_IGNORE_FILE = Path(f"{ROOT_DIR}/IGHignore.txt")
IGKL_IGNORE_FILE = Path(f"{ROOT_DIR}/IGKLignore.txt")
META_TABLE = Path(f"{ROOT_DIR}/meta.tsv")


def group_bedtools_output():
    '''
    - iterate paths to bedtools output for each single cell
    - collect single cell id information
    - calls parse_bedtools_output() to extract relevant single cell data 
    -
    '''
    all_sc_info = {}
    sonar_info = parse_sonar_output(SONAR_OUTPUT_FILE)
    pseudogene_info = parse_pseudogene_list(PSEUDOGENE_FILE)
    cells_to_drop = parse_dropped_cells(CELLS_TO_REMOVE_FILE, IGH_IGNORE_FILE, IGKL_IGNORE_FILE)

    bedtools_files = list(Path(f"{ROOT_DIR}/bedtools_mappedIG").glob("**/*.bed"))
    n_cells = len(bedtools_files)
    
    # Troubleshoot: make header for cells w/o sonar gene hits
    with open('../cells_not_in_sonar.csv', 'w') as handle:
        handle.write('cell_id\n')            
    
    with open('../cells_in_sonar_wo_gene-hits.csv', 'w') as handle:
        handle.write('cell_id\tsonar_V_call\tbedtools_mapped_IGs\n')        
        
    for (i, sc_path) in enumerate(bedtools_files):

        # For now, simple progress, real progressbar later
        if i % 10 == 0:
            print("{:,}/{:,}".format(i, n_cells), end='\r')

        cell_id = sc_path.parent.name

        # only collect bedtools output on cells not in cells to filter list
        if cell_id not in cells_to_drop:
            mapped_IGs = parse_bedtools_output(sc_path)
            transformed_info = transform_all_data(cell_id, mapped_IGs, sonar_info, pseudogene_info)

        if transformed_info is not None:
                all_sc_info[cell_id] = transformed_info
    
    return all_sc_info


# allows script to run group_bedtools_output() function, always at the end 
if __name__ == '__main__':
    all_sc_info = group_bedtools_output()
    cleaned_sc_info = clean_data(all_sc_info)
    appended_meta = append_meta(META_TABLE, cleaned_sc_info)

    # Export final merged raw reads output table
    cleaned_sc_info.to_csv(f'{ROOT_DIR}/final_merged_raw_reads.csv')
    
    # Export appended meta table with summary read counts
    appended_meta.to_csv(f'{ROOT_DIR}/appended_meta.csv')
