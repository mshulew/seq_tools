#!/usr/bin/env python3
# coding: utf-8

"""
converts read counts for miRNAs using gene name
consolidates all mature miRNA gene counts and tpm under gene name of precursor
used for comparing read counts generated using miRNA precursor gene to read counts generated using mature miRNAs
"""

import sys
import os
import tqdm
        
if __name__ == "__main__":

    version_number = 1 
    
# parse arguments
    if len(sys.argv) == 2:
        if sys.argv[1] == '-h' or sys.argv[1] == '-help':
            print('*'*100)
            print(os.path.basename(__file__) + ' version ' + str(version_number))
            print('Usage: ' + os.path.basename(__file__) + ' [read counts by gene] [list of miRNA precursor genes and mature miRNA product] [output filename]')
            exit()
    elif len(sys.argv) == 4:
        read_count_filename = sys.argv[1]
        mirna_gene_product_filename = sys.argv[2]
        output_filename = sys.argv[3]
    else:
        print('missing input and/or output file(s)')
        exit()
        
# import read count file 
    print('importing read counts by gene file...')
    read_count_data = []
    with open(read_count_filename, 'r') as read_count_file:
        for line in read_count_file:
            if 'gene_id' not in line:
                split_line = line.strip('\n').split('\t')
                read_count_data.append([split_line[0],split_line[1],split_line[2],split_line[3],split_line[4],split_line[5],split_line[6],split_line[7]])
                
# importing list of miRNA precursor genes and corresponding product
    print('importing list of miRNA precursor genes and corresponding product...')
    mirna_list = {}
    with open(mirna_gene_product_filename, 'r') as mirna_gene_product_file:
        for line in mirna_gene_product_file:
            split_line = line.strip('\n').split('\t')
            mirna_list.setdefault(split_line[1],split_line[0])
            
# processing read counts by gene
    print('processing read counts by gene...')
    processed_data = {}
    for gene in tqdm.tqdm(read_count_data):
        if gene[1] in mirna_list:
            gene_name = mirna_list[gene[1]]
            processed_data.setdefault(gene_name,[gene_name,gene_name,gene[2],gene[3],gene[4],0,0,gene[7]])
            processed_data[gene_name][5] = processed_data[gene_name][5] + int(gene[5]) 
            processed_data[gene_name][6] = processed_data[gene_name][6] + float(gene[6])
        else:
            processed_data.setdefault(gene[1],[gene[0],gene[1],gene[2],gene[3],gene[4],gene[5],gene[6],gene[7]])
            
# convert processed_data to list
    output_data = []
    for gene_name,entry in processed_data.items():
        output_data.append([entry[0],entry[1],entry[2],entry[3],entry[4],entry[5],entry[6],entry[7]])
                        

# write to file
    print('generating output...')
    output_file = open(output_filename, 'w')
    output_file.write('{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\r\n'.format('gene_id','gene_name','transcript_id','avg. transcript length','effective length','read_counts','TPM','RNA type'))
    for entry in output_data:
        output_file.write('{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\r\n'.format(entry[0],entry[1],entry[2],entry[3],entry[4],entry[5],entry[6],entry[7]))
    output_file.close()
            
    
