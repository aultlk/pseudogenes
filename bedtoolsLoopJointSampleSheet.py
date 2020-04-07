#!/usr/bin/env python

import os
import pandas as pd

## read column extracted joint sample sheet and label columns, default row labels
jointSampleSheetdf = pd.read_csv('/nethome/aultlk/Projects/Pseudogenes/jointSampleSheet.csv', 
                                 sep=',', 
                                 header=None, 
                                 usecols=[1, 2, 3, 4, 7], 
                                 names=['project', 'flow_cell', 'lane', 'index', 'cell_id'])

jointSampleSheetdf['new_col'] = jointSampleSheetdf['project'] + '-' + jointSampleSheetdf['lane']

## for each row of dataframe, iterate across each column and organize directories for each flow cell bedtools output
for (row_names, row_values) in jointSampleSheetdf.iterrows():
    cmd = "mkdir -p /nethome/aultlk/Projects/Pseudogenes/{}/{}/{}/{}/".format(
	  row_values['project'],
	  row_values['flow_cell'],
	  row_values['lane'],
          row_values['index']
	  )
    os.system(cmd)

## for each row of dataframe, iterate across each column and run command on OS system using bedtools multicov function, redirect all outputs to .csv files organ## ized into separate directories per project, flow cell, lane, and index

for (row_names, row_values) in jointSampleSheetdf.iterrows():
    cmd = "bedtools multicov -bams /hpcdata/vrc/vrc1_data/douek_lab/projects/RNASeq/{}/data/{}/{}/{}/STAR/recalibrated.bam -bed /nethome/aultlk/Projects/Pseudogenes/genes.sorted.bed > /nethome/aultlk/Projects/Pseudogenes/{}/{}/{}/{}/bedtoolsOutput.csv".format(
    	  row_values['project'],
	  row_values['flow_cell'],
	  row_values['lane'],
	  row_values['index'],
	  row_values['project'],
          row_values['flow_cell'],
          row_values['lane'],
          row_values['index']
	  )	  
    os.system(cmd)

print('done')

