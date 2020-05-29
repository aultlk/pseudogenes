

# Track cells not in sonar and cells without gene hits
ok_genes = set(mapped_IGs.index).intersection(sonar_info['gene'])

if sonar_info.size == 0:
    with open('../cells_not_in_sonar.csv', 'a') as handle:
        handle.write('{}\n'.format(cell_id))
    return

if not ok_genes:
    with open('../cells_in_sonar_wo_gene-hits.csv', 'a') as handle:
        gene_names_bedtools = ','.join(mapped_IGs['gene'].tolist())
        gene_names_sonar = ','.join(sonar_info['gene'].tolist())
        handle.write(
            '{}\t{}\t{}\n'
            .format(cell_id, genes_sonar, genes_bedtools)
            )
    return

# Troubleshoot: make header for cells w/o sonar gene hits
with open('../cells_not_in_sonar.csv', 'w') as handle:
    handle.write('cell_id\n')

with open('../cells_in_sonar_wo_gene-hits.csv', 'w') as handle:
     handle.write('cell_id\tsonar_V_call\tbedtools_mapped_IGs\n')
