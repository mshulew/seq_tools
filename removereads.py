#!/usr/bin/env python3
# coding: utf-8

"""
Removes reads from a fastq file
"""

import os
import sys
import shutil
import io
import gzip
import random

if __name__ == "__main__":
    
    version_num = 1.0
    
    print('*'*100)
    print('Remove Reads')
    print('2020 Mark Shulewiz Bio-Rad Laboratories, Inc.')
    print('version ' + str(version_num))
    print('*'*100)



    if len(sys.argv) == 2:
        if sys.argv[1] == '--help' or sys.argv[1] == '--h':
            print('*'*10)
            print('Usage: python3.6 ' + os.path.basename(__file__) + ' [input fastq] [list of reads to remove] [fraction to remove] > [output file]')
            print('*'*10)
            exit()
        else:
            print('ERROR with command line: unknown argument {}'.format(sys.argv[1]))
            exit()
    elif len(sys.argv) == 4:
        input_fastq_filename = sys.argv[1]
        reads_to_remove_filename = sys.argv[2]
        fraction_reads_to_keep = sys.argv[3]
    else:
        print('something is missing')
        exit()
    
    input_fastq = []
            
    with io.TextIOWrapper(gzip.open(input_fastq_filename, 'r')) as input_fastq_file:
        input_fastq_content = input_fastq_file.read().strip().split('\n')
    for a in range(0, len(input_fastq_content), 4):
        input_fastq.append([input_fastq_content[a],input_fastq_content[a+1],input_fastq_content[a+2],input_fastq_content[a+3]])
                    
# import list of reads
    reads_to_remove = []
    for line in open(reads_to_remove_filename, 'r'):
       reads_to_remove.append(line.strip('\n'))
    
# shuffle list and remove reads
    number_of_reads = round(int(fraction_reads_to_keep) * len(reads_to_remove))
    reads_to_remove = random.sample(reads_to_remove, number_of_reads)
    
# create output
    for read in input_fastq:
        read_name = read[0][1:].split()[0]
        if read_name not in reads_to_remove:
            print(read[0])
            print(read[1])
            print(read[2])
            print(read[3])
   


    
    
                
        
   
    
    
