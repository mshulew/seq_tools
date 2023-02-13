#!/usr/bin/env python3
# coding: utf-8

"""
retrieves nucleotide sequence from reference genome
usage: seqFromBed.py [reference genome as fasta file] [lbed file] [output file]
"""

import sys
import subprocess
import tqdm
        
if __name__ == "__main__":
    fasta_filename = sys.argv[1]
    bed_filename = sys.argv[2]
    output_filename = sys.argv[3]
    
    seq_dict = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N':'N'}
    
    print('importing reference genome...')
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
    seqs_to_get = []
    for line in open(bed_filename, 'r'):
        line = line.strip('\n').split('\t')
        seqs_to_get.append([line[0],int(line[1]),int(line[2]),line[3],line[5]])
        
    output_data = []
    for entry in tqdm.tqdm(seqs_to_get,total=len(seqs_to_get)):
        try:
            extracted_sequence = fasta_data[entry[0]][entry[1]:entry[2]]
            if entry[4] == '-':
                input_seq = extracted_sequence.upper()
                extracted_sequence = ''
                for a in reversed(input_seq):
                    extracted_sequence += seq_dict[a]
                       
        except:
            extracted_sequence = 'sequence not extracted; max pos: {} start: {} stop: {}'.format(len(fasta_data[entry[0]]),entry[1],entry[2])
        
        output_data.append([entry[0],entry[1],entry[2],entry[3],entry[4],extracted_sequence])
           
    print('writing to file...')
    output_file = open(output_filename, 'w')
    for entry in output_data:
        output_file.write('{}\t{}\t{}\t{}\t{}\t{}\n'.format(entry[0],entry[1],entry[2],entry[3],entry[4],entry[5]))
    output_file.close()
            
    
    
