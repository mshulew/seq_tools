#!/usr/bin/env python3
# coding: utf-8

"""

Aligns short sequences to a reference
Slow
Finds best alignment

python3.6 seq_tools/slowalign.py --F ribodepletion_probes/ref/ribodepletion_probes.fa --R ribodepletion_probes/ref/mm_rRNA_references.fa --outDir slowalign

"""

import os
import sys
import regex

def best_match(query):
    query_name = query[0]
    query_sequence = query[1]
    subject_name = ''
    subject_start = -1
    aligned_seq = ''
    match_type = ''
    
    for ref_seq in ref_sequences:
        for b in range(0,len(query_sequence)):
            if subject_start > -1:
                break 
                
            for ref_seq in ref_sequences:
                try:
                    aligned_seq = regex.search("(" + query_sequence + "){e<=" + str(b) + "}",ref_seq[1]).group()
                    match_type = str(b) + "  mismatch"
                    subject_start = ref_seq[1].find(aligned_seq)
                except:
                    rev_comp_seq = ''
                    for c in reversed(query_sequence):
                        rev_comp_seq += seq_dict[c]
                    try:
                        aligned_seq = regex.search("(" + rev_comp_seq + "){e<=" + str(b) + "}",ref_seq[1]).group()
                        match_type = str(b) + " mismatch, reverse complement"
                        subject_start = ref_seq[1].find(aligned_seq)
                    except:
                        aligned_seq = ''
                    
                if subject_start > -1:
                    subject_name = ref_seq[0]
                    subject_start = subject_start + 1
                    break
                    
    search_results = [query_name,query_sequence,subject_name,match_type,subject_start,aligned_seq]
    return search_results

if __name__ == "__main__":
    
    version_num = 1.0
    
    print('*'*100)
    print('Slow Aligner')
    print('2020 Mark Shulewiz Bio-Rad Laboratories, Inc.')
    print('version ' + str(version_num))
    print('*'*100)
    
    seq_dict = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N':'N'}
    
# parse arguments
    query_fasta_filename = ''
    input_ref_filename = ''
    output_dir = os.getcwd() + '/slowalign_output'

    if len(sys.argv) == 2:
        if sys.argv[1] == '--help' or sys.argv[1] == '--h':
            print('*'*10)
            print('Usage: python3.6 ' + os.path.basename(__file__) + ' --F [fasta file with query sequences] --R [reference FASTA file] --outDir [output directory]')
            print('*'*10)
            exit()
        else:
            print('ERROR with command line: unknown argument {}'.format(sys.argv[1]))
            exit()
    elif len(sys.argv) > 2:
        for a in range(0, len(sys.argv)):
            if sys.argv[a] == '--F':
                query_fasta_filename = sys.argv[a+1]
            elif sys.argv[a] == '--R':
                input_ref_filename = sys.argv[a+1]
            elif sys.argv[a] == '--outDir':
                output_dir = sys.argv[a+1]
    else:
        print('something is missing')
        exit()
       
# check for output directory and create if necessary
    if output_dir != os.getcwd():
        if not os.path.isdir(output_dir):
            os.makedirs(output_dir)
            
    print('Output directory: {}'.format(output_dir))
    print('*'*10)
        
# prepare reference sequences
    print('Importing reference sequences')
    ref_sequences = []
    ref_seq_content = open(input_ref_filename, 'r').read().strip().split('\n')
    for a in range(0, len(ref_seq_content),2):  
        ref_sequences.append([ref_seq_content[a][1:],ref_seq_content[a+1]])
        
# iterate through query sequences
    query_sequences = []
    query_seq_content = open(query_fasta_filename, 'r').read().strip().split('\n')
    for a in range(0, len(query_seq_content),2):
        query_sequences.append([query_seq_content[a][1:],query_seq_content[a+1]])
    for query_seq in query_sequences:
        print(best_match(query_seq))
                
                
