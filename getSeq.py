#!/usr/bin/env python3
# coding: utf-8

"""
retrieves nucleotide sequence from fasta file
usage: seqFromFasta.py [fasta file] [list of sequences to extract, each as name chromosome strand startpos stoppos] [output file]
"""

import sys
import subprocess
import tqdm
        
if __name__ == "__main__":
    fasta_filename = sys.argv[1]
    list_of_sequences_file = sys.argv[2]
    output_filename = sys.argv[3]
    
    seq_dict = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N':'N'}
    
    print('importing fasta file...')
# count number of lines in fasta file
    result = subprocess.run(['wc','-l',fasta_filename], stdout=subprocess.PIPE)
    total_lines = int(str(result.stdout).split()[0].strip("b'"))

    fasta_data = {}
    entry_name = ''
    entry_seq = ''
    for line in tqdm.tqdm(open(fasta_filename, 'r'),total=total_lines):
        if line[0] == '>':
            if entry_seq == '':
                entry_name = line.strip('\n').split()[0][1:]
            else:
                fasta_data.setdefault(entry_name,entry_seq)
                entry_name = line.strip('\n').split()[0][1:]
                entry_seq = ''
        else:
            entry_seq = entry_seq + line.strip('\n')
            
    print('getting sequences...')
    list_of_sequences = []
    for line in open(list_of_sequences_file, 'r'):
        split_line = line.strip('\n').split('\t')
        list_of_sequences.append([split_line[0],split_line[1],split_line[2],int(split_line[3]),int(split_line[4])])
        
    output_data = []
    for entry in tqdm.tqdm(list_of_sequences,total=len(list_of_sequences)):
        try:
            extracted_sequence = fasta_data[entry[1]][entry[3]:entry[4]]
            if entry[2] == '-':
                input_seq = extracted_sequence.upper()
                extracted_sequence = ''
                for a in reversed(input_seq):
                    extracted_sequence += seq_dict[a]
                       
        except:
            extracted_sequence = 'sequence not extracted'
        
        output_data.append([entry[0],entry[1],entry[2],entry[3],extracted_sequence,entry[4]])
    
    output_data = sorted(output_data)
            
    print('writing to file...')
    output_file = open(output_filename, 'w')
    for entry in output_data:
        output_file.write('{}\t{}\t{}\t{}\t{}\n'.format(entry[0],entry[1],entry[2],entry[3],entry[4]))
    output_file.close()
            
    
    
