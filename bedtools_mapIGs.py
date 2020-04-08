# This script retreives single cell data, filters dropped cells, and extracts mapped
# immunoglobulin genes and read counts using the bedtools multicov function.   

from pathlib import Path

import os
import pandas as pd

ROOT_DIR = Path(__file__).resolve().parent.parent
SAMPLE_FILE = Path(ROOT_DIR, 'jointSampleSheet.csv')
DROPPED_IGH_FILE = Path(ROOT_DIR, 'IGHignore.txt')
DROPPED_IGKL_FILE = Path(ROOT_DIR, 'IGKLignore.txt')
RNASEQ_DIR = Path('/hpcdata/vrc/vrc1_data/douek_lab/projects/RNASeq/')

# Read sample sheet and create data frame with single cell identifiers
sample_df = pd.read_csv(SAMPLE_FILE,
                        sep=',', 
                        header=None, 
                        usecols=[1, 2, 3, 4, 7], 
                        dtype={'lane': str},
                        names=['project', 'flow_cell', 'lane', 'index', 'cell_id'])


# Runs bedtools multicov function on os command line
# Directs outputs to single cell identifiable directories
for (project, flow_cell, lane, index, cell_id) in sample_df.to_numpy():

    # Create single cell output directory, reformat source_ID
    output_dir = Path(ROOT_DIR, 'bedtools_mappedIG', project, cell_id)
    output_dir.mkdir(parents=True)

    # Locate bam file for single cell 
    bam_file = Path(RNASEQ_DIR, 
                    project,
                    'data', 
                    flow_cell,
                    lane,
                    'STAR', 
                    'recalibrated.bam')

    # Define bed and output file
    bed_file = f"{ROOT_DIR}/genes.sorted.bed"
    output_file = f"{flow_cell}-{lane}-read1_index_{index}.bed"

    # run command in os
    cmd = f"bedtools multicov -bams {bam_file} -bed {bed_file} > {output_dir}/{output_file}"
    os.system(cmd)


 




