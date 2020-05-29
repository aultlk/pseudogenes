"""
   This script contains code that can be used to bin loci and cells based on detected
   gene rearrangements from mapped raw reads. Can be integrated into loader.py if 
   we decide to compare binning strategies using BALDR and mapped raw reads.

"""



def bin_cells(df):    


    def get_bin(rearr_serie):
        bin_assignment = 'not_detected'
        
        for bin_label in ['pseudogene', 'passenger', 'productive']:
            if (rearr_serie == bin_label).any():
                bin_assignment = bin_label
                break # if condition met, exit the loop
        return pd.Series(bin_assignment, index=rearr_serie.index)

    df_rpm_pos = df[df.RPM > 0]

    # Bin per locus 
    df_rpm_pos['bedtools_locus_bin'] = (
        df_rpm_pos.reset_index()
        .groupby(['cell_id', 'locus'])
        .rearrangement
        .transform(get_bin).tolist())

    # Bin per cell
    df_rpm_pos['bedtools_cell_bin'] = (df_rpm_pos.reset_index()
                              .groupby('cell_id')
                              .rearrangement
                              .transform(get_bin)
                              .tolist())

    return df_rpm_pos
