#!/usr/bin/env python3
# coding: utf-8

"""

analyzes a fastq illumina sequence and lists all index sequences
searches for perfect matches with index seq

"""

import sys
import os
import tqdm
import regex

if __name__ == "__main__":

# parse arguments
    
    if len(sys.argv) != 4:
        print('usage: {} [input filename] [output filename] [miseq/nextseq/none]'.format(os.path.basename(__file__)))
        exit()

    input_filename = sys.argv[1]
    output_filename = sys.argv[2]
    sequencer = sys.argv[3]
    
    if sequencer == 'miseq':
        i5_sequences = ['TATAGCCT','ATAGAGGC','CCTATCCT','GGCTCTGA','AGGCGAAG','TAATCTTA']
        i5_dictionary = {'TATAGCCT':'i501','ATAGAGGC':'i502','CCTATCCT':'i503','GGCTCTGA':'i504',
                         'AGGCGAAG':'i505','TAATCTTA':'i506','none':'none'}
    elif sequencer == 'nextseq':
        i5_sequences = ['AGGCTATA','GCCTCTAT','AGGATAGG','TCAGAGCC','CTTCGCCT','TAAGATTA']
        i5_dictionary = {'AGGCTATA':'i501','GCCTCTAT':'i502','AGGATAGG':'i503','TCAGAGCC':'i504',
                         'CTTCGCCT':'i505','TAAGATTA':'i506','none':'none'}
        
    i7_sequences = ['ATTACTCG','TCCGGAGA','CGCTCATT','GAGATTCC','ATTCAGAA','GAATTCGT',
                    'CTGAAGCT','TAATGCGC','CGGCTATG','TCCGCGAA','TCTCGCGC','AGCGATAG']
    i7_dictionary = {'ATTACTCG':'i701','TCCGGAGA':'i702','CGCTCATT':'i703','GAGATTCC':'i704','ATTCAGAA':'i705','GAATTCGT':'i706',
                     'CTGAAGCT':'i707','TAATGCGC':'i708','CGGCTATG':'i709','TCCGCGAA':'i710','TCTCGCGC':'i711','AGCGATAG':'i712','none':'none'}

# count total lines in input file
    print('counting total lines in input file...')
    with open(input_filename, 'r') as input_file:
        total_lines = 0
        for line in tqdm.tqdm(input_file):
            total_lines += 1


# iterate through input file
    with open(input_filename, 'r') as input_file:
        sequences = {}
        sequence_count = 0
        line_num = 0
        for line in tqdm.tqdm(input_file, total=total_lines):
            line_num += 1
            if line_num == 1:
                first_line = line.strip('\n')
                obs_index = line.split()[1].split(':')[3].strip('\n')
                obs_i7 = obs_index.split('+')[0]
                obs_i5 = obs_index.split('+')[1]
                
# find matching i7_seq
                if sequencer != 'none':
                    if any(x in obs_i7 for x in i7_sequences) and any(y in obs_i5 for y in i5_sequences):
                        sequences.setdefault(obs_index, 0)
                        sequences[obs_index] += 1
                        sequence_count += 1
                else:
                    sequences.setdefault(obs_index, 0)
                    sequences[obs_index] += 1
                    sequence_count += 1
                
            if line_num == 2:
                second_line = line.strip('\n')
            if line_num == 3:
                third_line = line.strip('\n')
            if line_num == 4:
                fourth_line = line.strip('\n')
                
                line_num = 0

# save counts to to file
        output_file= open(output_filename, 'w')
        output_file.write('input file: {} \r\n'.format(input_filename))
        output_file.write('index\count\r\n')
        for obs_index,count in sequences.items():
            output_file.write('{}\t{}\r\n'.format(obs_index,count))
        output_file.close()

