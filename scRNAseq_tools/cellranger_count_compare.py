#!/usr/bin/env python3
# coding: utf-8

"""
compares count matrix create by cellranger to BAM output file (possorted_genome_bam.bam - found in outs/ directory in cellranger-7.2.0)
attempts to filter alignments for each gene to reproduce cellranger counts
as of 12/9/2023 is almost there, but counts are slightly higher than cellranger for some gene

can reproduce cellranger counts by filtering for xf:i:25 tag and filtering reads in barcodes.tsv file

pre-running: convert Aligned.sortedByCoord.bam to sam (include tags for barcode, umi & genes when running STARsolo)
identifies most of the accepted reads as cellranger
tip: xf:i:25 will filter all accepted reads

usage: python3.11 cellranger_analysis.py [path to cellranger output folder] [path to output files created by this script]

"""

import sys
import os
import time

def createDictionary(filepath, column):
# creates dictionary from tab delimiated file
    dictionary = {}
    entry_num = 0
    with open(filepath, 'r') as file:
        for entry in file:
            entry_num += 1
            dictionary.setdefault(entry_num,entry.splitlines()[0].split('\t')[column])
    return dictionary

def timeFormat(start,current):
    elapsed_time = round(current-start,0)
    return elapsed_time

def extractGeneName(samEntry):
    gene_index = [i for i, s in enumerate(samEntry) if 'GN:Z:' in s][0]
    GeneName = samEntry[gene_index][5:]
    return GeneName
            
def main():
    print('cellranger analysis version 1.0 running...')

    
# start clock
    start_time = time.time()
    
# enter path to STARsolo output firectory
    cellranger_path = sys.argv[1]
    output_path = sys.argv[2]
    matrixfile = cellranger_path + 'outs/filtered_feature_bc_matrix/matrix.mtx'
    barcodefile = cellranger_path + 'outs/filtered_feature_bc_matrix/barcodes.tsv'
    featurefile = cellranger_path + 'outs/filtered_feature_bc_matrix/features.tsv'
    samfile = cellranger_path + 'outs/possorted_genome_bam.sam'
    log_path = output_path + '/log.tsv'

# generate output dir
    if not os.path.isdir(output_path):
        os.makedirs(output_path)

# create dictionaries from barcode and feature files
    bc_dict = createDictionary(barcodefile,0)
    feature_dict = createDictionary(featurefile,1)

# create list of accepted barcodes
    accepted_bc = []
    with open(barcodefile, 'r') as bcfile:
        for line in bcfile:
            accepted_bc.append(line.splitlines()[0])

# convert matrix into [feature_name,barcode,counts]
    matrix = []
    entry_num = 0
    with open(matrixfile, 'r') as mtx:
        for entry in mtx:
            entry_num += 1
            if entry_num > 3:
                entry_list = [int(x) for x in entry.splitlines()[0].split()]
                matrix.append([feature_dict[entry_list[0]],bc_dict[entry_list[1]],entry_list[2]])

# sort matrix by gene and then by barcode
    matrix_by_gene = sorted(matrix, key = lambda x:(x[0],x[1]))

# create list of genes with sublist of barcode and counts
    gene_bc_counts = {}
    for entry in matrix_by_gene:
        gene = entry[0]
        barcode = entry[1]
        count = entry[2]
        gene_bc_counts.setdefault(gene,[])
        gene_bc_counts[gene].append([barcode,count])

# read SAM file
    with open(log_path, 'w') as output_file:
        output_file.write('#start analysis\n')
        output_file.write('#reading sam file\n')
        output_file.write('elapsed time\t{} sec\n'.format(timeFormat(start_time,time.time())))
    
# speed things up by indexing by gene name
    sam = {}
    total_alignments = 0
    with open(samfile, 'r') as samf:
        for alignment in samf:
            total_alignments += 1
            if total_alignments % 1000000 == 0:
                with open(log_path, 'a') as output_file:
                    output_file.write('sam file entries reads:\t{:,}\n'.format(total_alignments))
                    output_file.write('elapsed time\t{} sec\n'.format(timeFormat(start_time,time.time())))
            if 'GN:Z:' in alignment:
                geneName = extractGeneName(alignment.splitlines()[0].split('\t'))
                sam.setdefault(geneName,[])
                sam[geneName].append(alignment.splitlines()[0].split('\t'))

# output to log
    with open(log_path, 'a') as output_file:
        output_file.write('total alignments\t{:,}\n'.format(total_alignments))
        output_file.write('alignments assigned to genes\t{:,}\n'.format(len(sam)))
        output_file.write('cell barcodes\t{:,}\n'.format(len(accepted_bc)))
        output_file.write('total genes assigned\t{:,}\n'.format(len(gene_bc_counts)))
        output_file.write('elapsed time\t{} sec\n'.format(timeFormat(start_time,time.time())))
        output_file.write('---\n')
        output_file.write('{}\t{}\t{}\t{}\t{}\n'.format('#','gene','umis_accepted','matrix_counts','elapsed_time(sec)'))
                          

# iterate through gene_bc_counts and extract reads from SAM file
    genes_processed = 0
    for gene in list(gene_bc_counts.items()):
        genes_processed += 1
        if gene[0] in list(sam.keys()):
            alignments = sam[gene[0]]
            geneName = gene[0]
            accepted_barcodes = {}
            accepted = 0
            rejected = 0
            for alignment in alignments:
                samflag = alignment[1]
                mapq = alignment[4]
                barcode_index = [i for i, s in enumerate(alignment) if 'CR:Z:' in s][0]
                cell_barcode = alignment[barcode_index][5:]
                if any(cell_barcode in x for x in accepted_bc) and "{0:b}".format(int(samflag))[-11:-10] != '1' and int(mapq) > 3:
                    accepted_barcodes.setdefault(cell_barcode,0)
                    accepted_barcodes[cell_barcode] += 1
                    accepted += 1
                elif any(cell_barcode in x for x in accepted_bc) and "{0:b}".format(int(samflag))[-11:-10] != '1' and int(mapq) <= 3:
                    rejected += 1
                elif not any(cell_barcode in x for x in accepted_bc) and "{0:b}".format(int(samflag))[-11:-10] != '1':
                    rejected += 1
                elif any(cell_barcode in x for x in accepted_bc) and "{0:b}".format(int(samflag))[-11:-10] == '1':
                    rejected += 1
                elif not any(cell_barcode in x for x in accepted_bc) and "{0:b}".format(int(samflag))[-11:-10] == '1':
                    rejected += 1

# count total reads from matrix
            total_matrix_counts = 0
            for barcode in gene[1]:
                total_matrix_counts += barcode[1]

# combine accepted barcode list
            combined_bc = [['barcode','accepted','matrix']]
            for a in sorted(list(accepted_barcodes.items())):
                bc_found = False
                for b in gene[1]:
                    if a[0] in b[0]:
                        bc_found = True
                        combined_bc.append([a[0],a[1],b[1]])
                        break
                if not bc_found:
                    combined_bc.append([a[0],a[1],'NA'])
                    
# write to file
# speed things up by writing to file occasionally
            with open(log_path, 'a') as output_file:
                output_file.write('{}\t{}\t{}\t{}\t{}\n'.format(genes_processed,geneName,accepted,total_matrix_counts,timeFormat(start_time,time.time())))
                if abs(total_matrix_counts - accepted) > 2:
                    with open(output_path + '/' + str(genes_processed) + '_' + geneName + '.tsv','w') as gene_output:
                        gene_output.write('barcode\tumis\tmatrix\n')
                        for a in combined_bc:
                            gene_output.write('{}\t{}\t{}\n'.format(a[0],a[1],a[2]))
    

if __name__ == "__main__":
    main()
