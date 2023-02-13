#!/usr/bin/env python3
# coding: utf-8

"""

retrieves nucleotide sequence for a transcript from reference genome
requires transcript_id, gtf file and fasta genome reference

"""

import sys
import os
import subprocess
import tqdm


if __name__ == "__main__":

# parse arguments
    transcript = ''
    gtf_filename = ''
    ref_genome_filename = ''
    
    seq_dict = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N':'N', ' ':''}
    
    if len(sys.argv) < 4:
        print('usage: {} --t [transcript_id] --gtf [gtf file] --genome [reference genome] --o [output_file]'.format(os.path.basename(__file__)))
        print('instead of --t, can use file with names of transcripts with --f [file name]')
        exit()
    else:
        input_list = []
        input_trigger = False
        for a in range(1, len(sys.argv)):
            if input_trigger == True:
                if sys.argv[a][0:2] == '--':
                    input_trigger = False
                else:
                    input_list.append(sys.argv[a])
            if input_trigger == False:
                if sys.argv[a] == '--t':
                    input_trigger = True
                elif sys.argv[a] == '--f':
                    for line in open(sys.argv[a+1], 'r'):
                        input_list.append(line.strip('\n'))
                elif sys.argv[a] == '--gtf':
                    gtf_filename = sys.argv[a+1]
                elif sys.argv[a] == '--genome':
                    ref_genome_filename = sys.argv[a+1] 
                elif sys.argv[a] == '--o':
                    output_filename = sys.argv[a+1] 

# count total lines in ref genome
    result = subprocess.run(['wc','-l',ref_genome_filename], stdout=subprocess.PIPE)
    total_lines = int(str(result.stdout).split()[0].strip("b'"))
    
# import reference genome
    print('importing reference genome...')
    ref_genome = {}
    with open(ref_genome_filename, 'r') as ref_genome_file:
        ref_genome_chromosome = ''
        ref_genome_sequence = ''
        for line in tqdm.tqdm(ref_genome_file,total=total_lines):
            if line[0] == '>':
                if ref_genome_chromosome != '':
                    ref_genome.setdefault(ref_genome_chromosome,ref_genome_sequence)
                    ref_genome_sequence = ''
                    ref_genome_chromosome = line[1:].split()[0]
                elif ref_genome_chromosome == '':
                    ref_genome_chromosome = line[1:].split()[0]
            elif line[0] != '>':
                ref_genome_sequence = ref_genome_sequence + line.strip('\n')
                
        ref_genome.setdefault(ref_genome_chromosome,ref_genome_sequence)

# iterate through gtf file
    print('identifying exons for each transcript...')
    exon_data = [[],[],[],[],[],[]]
    for transcript in tqdm.tqdm(input_list, total=len(input_list)):
        transcript_search = '"{}"'.format(transcript)
        with open(gtf_filename, 'r') as gtf_file:
            transcript_found = False
            for line in gtf_file:
                split_line = line.strip('\n').split('\t')
                if split_line[2] == 'exon':
                    if transcript_found == False:
                        if transcript_search in split_line[8]:
                            exon_data[0].append(transcript)
                            exon_data[1].append(split_line[0])
                            exon_data[2].append(int(split_line[3]))
                            exon_data[3].append(int(split_line[4]))
                            exon_data[4].append('')
                            exon_data[5].append(split_line[6])
                            transcript_found = True
                    elif transcript_found == True:
                        if transcript_search in split_line[8]:
                            exon_data[0].append(transcript)
                            exon_data[1].append(split_line[0])
                            exon_data[2].append(int(split_line[3]))
                            exon_data[3].append(int(split_line[4]))
                            exon_data[4].append('')
                            exon_data[5].append(split_line[6])
                        else:
                            break
                            
# retrieve sequences for each exon
    print('retrieving sequences for each exon...')
    total_lines = len(exon_data[0])
    for a in tqdm.tqdm(range(0,len(exon_data[0]))):
        chromosome=exon_data[1][a]
        exon_start = exon_data[2][a] - 1
        exon_stop = exon_data[3][a]
        sequence = ref_genome[chromosome][exon_start:exon_stop]
        if exon_data[5][a] == '-':
            rev_comp  = ''
            for bp in reversed(sequence):
                rev_comp += seq_dict[bp]
            sequence = rev_comp
        exon_data[4][a] = sequence
                        
# print exon boundaries
    output_file = open(output_filename, 'w')
    for b in range(0,len(exon_data[0])):
        chrom_loc = exon_data[1][b] + ':' + str(exon_data[2][b]) + ':' + str(exon_data[3][b])
        output_file.write('>{}\r\n'.format(exon_data[0][b]))
        output_file.write('>{}\r\n'.format(chrom_loc))
        output_file.write('{}\r\n'.format(exon_data[4][b]))
    output_file.close()
        
                
