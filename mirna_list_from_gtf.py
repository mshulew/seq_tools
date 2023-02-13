#!/usr/bin/env python3
# coding: utf-8

"""
generates a list of all miRNAs from a gtf file along with gene name of precursor miRNA
"""

import sys
import tqdm
        
if __name__ == "__main__":

    version_number = 1 
    
# parse arguments
    if len(sys.argv) == 2:
        if sys.argv[1] == '-h' or sys.argv[1] == '-help':
            print('*'*100)
            print(os.path.basename(__file__) + ' version ' + str(version_num))
            print('Usage: ' + os.path.basename(__file__) + ' [gtf mature filename] [gtf precursor filename] [output filename]')
            print('gtf source = gtf with mature miRNAs only')
            print('gtf target = gtf with precursor miRNAs only')
            exit()
    elif len(sys.argv) == 4:
        gtf_mature_filename = sys.argv[1]
        gtf_precursor_filename = sys.argv[2]
        output_filename = sys.argv[3]
    else:
        print('missing input and/or output file(s)')
        exit()
        
# import gtf file with precursor miRNAs
    print('importing gtf file with precursor miRNAs...')
    precursor_gtf_data = []
    with open(gtf_precursor_filename, 'r') as precursor_gtf_file:
        for line in precursor_gtf_file:
            if line[0] != '#':
                split_line = line.strip('\n').split('\t')
                precursor_gtf_data.append([split_line[0],split_line[1],split_line[2],split_line[3],split_line[4],split_line[5],split_line[6],split_line[7],split_line[8]])
 
 # import gtf file with mature miRNAs
    print('importing gtf file with mature miRNas...')
    mature_gtf_data = []
    with open(gtf_mature_filename, 'r') as mature_gtf_file:
        for line in mature_gtf_file:
            if line[0] != '#':
                split_line = line.strip('\n').split('\t')
                mature_gtf_data.append([split_line[0],split_line[1],split_line[2],split_line[3],split_line[4],split_line[5],split_line[6],split_line[7],split_line[8]])
                            
# analyze precursor gtf file
    print('analyzing gtf file with precursor miRNAs...')
    output_data = []
    for entry in tqdm.tqdm(precursor_gtf_data):
        gene_name = ''
        gene_product = ''
        if entry[2] == 'exon':
            if 'miRBase' in entry[8]:
                annotations = entry[8].split(';')
                for annotation in annotations:
                    if 'gene "' in annotation:
                        gene_name = annotation.split('gene "')[1].strip('"')
                    if 'product "' in annotation:
                        gene_product = annotation.split('product "')[1].strip('"')
                output_data.append([gene_name,gene_product])
                
# analyze mature gtf file
    print('analyzing gtf file with mature miRNAs...')
    for entry in tqdm.tqdm(mature_gtf_data):
        gene_name = ''
        gene_product = ''
        if entry[2] == 'exon':
            if 'Miltenyi_SKU' in entry[8]:
                annotations = entry[8].split(';')
                for annotation in annotations:
                    if 'Miltenyi_SKU "' in annotation:
                        miltenyi_sku = annotation.split('Miltenyi_SKU "')[1].strip('"')
                    if 'product "' in annotation:
                        gene_product = annotation.split('product "')[1].strip('"')
                for mature_mirna in output_data:
                    if gene_product == mature_mirna[1]:
                        output_data.append([mature_mirna[0],miltenyi_sku])
                        break

# write to file
    print('generating output...')
    output_file = open(output_filename, 'w')
    output_file.write('gene_name\tgene_product\r\n')
    for entry in output_data:
        output_file.write('{}\t{}\r\n'.format(entry[0],entry[1]))
    output_file.close()
            
    
