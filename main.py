# Pull mapped IG genes and number of reads from `bedtools muticov` output
# Merge with master/meta.tsv while also pulling data from filtering2_rearrangements_single-cell.tsv file

from pathlib import Path
import sys 

import pandas as pd

from loader import parse_bedtools_output 
from loader import parse_sonar_output
from loader import parse_pseudogene_list
from loader import transform_table_into_cool_info 

ROOT_DIR = Path(__file__).resolve().parent.parent
SONAR_OUTPUT_FILE = Path(f"{ROOT_DIR}/filtering2_rearrangements_single-cell.tsv")
PSEUDOGENE_FILE = Path(f"{ROOT_DIR}/allpseudogenes.list")

def group_bedtools_output():
    '''
    - iterate paths to bedtools output for each single cell
    - collect single cell id information
    - calls parse_bedtools_output() to extract relevant single cell data 
    -
    '''
    all_sc_info = []
    sonar_info = parse_sonar_output(SONAR_OUTPUT_FILE)
    pseudogene_info = parse_pseudogene_list(PSEUDOGENE_FILE)

    for sc_path in Path(f"{ROOT_DIR}/mapped_IG").glob("**/bedtoolsOutput.csv"):
        sc_id = [x.name for x in sc_path.parents][:4]
        print(sc_id, sc_path)

        sc_info = parse_bedtools_output(sc_path) # continue loop to next function
        # sc_info.append(sc_id)
        transformed_info = transform_table_into_cool_info(sc_id, sc_info, sonar_info, pseudogene_info)
        all_sc_info.append(sc_info)

    import ipdb; ipdb.set_trace()


# Next, we will add section that checks if gene is in pseudogene database. 
# Then add ipdb so we can check it if it is still not found. 
# Be careful about parenthesis in pseudogene list (e.g. II, III) 

# allows script to run group_bedtools_output() function, always at the end 
if __name__ == '__main__':
    group_bedtools_output()



