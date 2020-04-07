import pandas as pd

# load file with duplicate count and v identity single cell info 
data = pd.read_csv('../filtering2_rearrangements_single-cell.tsv', sep='\t')

# split v_call column so that multiple alleles are separated into a list of alleles
data['v_call'] = data['v_call'].str.split(',')

# convert list of alleles into their own row
data = data.explode('v_call')


# subset on 'status' == 'good'
data = data[data['status'] == 'good']


print("Number of duplicated entries with status==good: {}".format(
    data[['v_call', 'source_id', 'cell_id']].duplicated().sum())
)

# remove unknown number tagged onto the end of the source_id after the index
# cannot find this number in jointSampleSheet.csv 
# only identifiers in jointSampleSheet are: 
# cell_id and source_id = flow_cell-lane-index
data['source_id_no_mapped-unmapped'] = data['source_id'].str.replace('IG\-mapped.*', '')

print("Number of duplicated entries with status==good and without unknown identifier : {}".format(
    data[['v_call', 'source_id_no_mapped-unmapped', 'cell_id']].duplicated().sum())
)

