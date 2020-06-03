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

import pandas as pd

from loader import load_dropped_cells 
from loader import process_raw_reads
from loader import process_baldr_sonar
from loader import agg_gene_duplicates
from loader import load_pseudogenes
from loader import get_gene_rearrs
from loader import clean_data
from loader import append_meta

ROOT_DIR = Path(__file__).resolve().parent.parent
SONAR_TABLE = Path(f"{ROOT_DIR}/meta.cas.tsv")
PSEUDOGENE_LIST = Path(f"{ROOT_DIR}/pseud.HKL")
CELLS_TO_REMOVE_LIST = Path(f"{ROOT_DIR}/cellstoremove.txt")
IGH_IGNORE_LIST = Path(f"{ROOT_DIR}/IGHignore.txt")
IGKL_IGNORE_LIST = Path(f"{ROOT_DIR}/IGKLignore.txt")
META_TABLE = Path(f"{ROOT_DIR}/meta.tsv")


def main():

    all_gene_rearrs = {}
    cells_w_gene_duplicates = {}
    unmapped_sonar_gene_cells = {}

    sonar_db = process_baldr_sonar(SONAR_TABLE)
    pseudogenes = load_pseudogenes(PSEUDOGENE_LIST)
    dropped_cells = load_dropped_cells(CELLS_TO_REMOVE_LIST, IGH_IGNORE_LIST, IGKL_IGNORE_LIST)

    # Progress bar 
    bedtools_beds = list(Path(f"{ROOT_DIR}/bedtools_mappedIG").glob("**/*.bed"))
    n_cells = len(bedtools_beds)
    for sc_path in tqdm(bedtools_beds):        
        cell_id = sc_path.parent.name    

        if cell_id not in dropped_cells:
            sc_sonar = sonar_db[sonar_db.cell_id == cell_id]   

            if sc_sonar.index.duplicated().any():
                cells_w_gene_duplicates[cell_id] = (sc_sonar[sc_sonar.index.duplicated()]
                                                    .index.tolist())
                sc_sonar = agg_gene_duplicates(cell_id, sc_sonar)

            mapped_Vs = process_raw_reads(sc_path)    
            (gene_rearrs, unmapped_sonar_genes) = get_gene_rearrs(cell_id, mapped_Vs, 
                                                                  sc_sonar, pseudogenes)
                
            if gene_rearrs is not None:
                all_gene_rearrs[cell_id] = gene_rearrs
            if unmapped_sonar_genes is not None:
                unmapped_sonar_gene_cells[cell_id] = unmapped_sonar_genes

    cleaned_gene_rearrs = clean_data(all_gene_rearrs)
    appended_meta = append_meta(META_TABLE, cleaned_gene_rearrs)
    cells_w_gene_duplicates = pd.DataFrame.from_dict(cells_w_gene_duplicates, orient='index')
    cells_w_unmapped_sonar_genes = pd.DataFrame.from_dict(unmapped_sonar_gene_cells, orient='index')

    return (appended_meta, cells_w_gene_duplicates, cells_w_unmapped_sonar_genes)


if __name__ == '__main__':
 
    # Run pipeline
    (appended_meta, cells_w_gene_duplicates, cells_w_unmapped_sonar_genes) = main()
    
    # Save outputs
    appended_meta.to_csv(f'{ROOT_DIR}/appended_meta.csv')
    cells_w_gene_duplicates.to_csv(f'{ROOT_DIR}/cells_w_gene_duplicates.csv')
    cells_w_unmapped_sonar_genes.to_csv(f'{ROOT_DIR}/cells_w_unmapped_sonar_genes.csv')










    
