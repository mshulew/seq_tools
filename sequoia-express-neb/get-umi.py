#!/usr/bin/env python3
# coding: utf-8

"""

Transfers 11 bp UMI to R2
For use with NEB UMI adapters

"""

import os
import sys
import io
import gzip


if __name__ == "__main__":
    
    version_num = 1.0
    
    print('*'*100)
    print('NEB Add UMI to R2')
    print('2021 Mark Shulewiz Bio-Rad Laboratories, Inc.')
    print('version ' + str(version_num))
    print('*'*100)
    
# parse arguments
    read_filename = ''
    index_filename = ''
    output_filename = ''

    if len(sys.argv) == 2:
        if sys.argv[1] == '--help' or sys.argv[1] == '--h':
            print('*'*10)
            print('Usage: python3.6 ' + os.path.basename(__file__) + ' --read [r3 fastq filename] --umi [r2 fastq filename = umi] --output [output fastq]')
            print('*'*10)
            exit()
        else:
            print('ERROR with command line: unknown argument {}'.format(sys.argv[1]))
            exit()
    elif len(sys.argv) > 2:
        for a in range(0, len(sys.argv)):
            if sys.argv[a] == '--read':
                read_filename = sys.argv[a + 1]
            elif sys.argv[a] == '--umi':
                index_filename = sys.argv[a + 1]
            elif sys.argv[a] == '--output':
                output_filename = sys.argv[a + 1]
    else:
        print('something is missing')
        exit()
        
# process files
    if '.fastq.gz' not in read_filename:
        print('invalid read file - does not end in fastq.gz')
        exit()
    if '.fastq.gz' not in index_filename:
        print('invalid index file - does not end in fastq.gz')
        exit()
    if '.fastq.gz' not in output_filename:
        print('invalid output filename - does not end in fastq.gz')
        exit()
        
    output_fastq = []
    
    with io.TextIOWrapper(gzip.open(read_filename, 'r')) as read_file,io.TextIOWrapper(gzip.open(index_filename, 'r')) as index_file:
        read_file_content = read_file.read().strip().split('\n')
        index_file_content = index_file.read().strip().split('\n')
        for a in range(0, len(read_file_content), 4):
            umi = index_file_content[a + 1][:11]
            umi_q = index_file_content[a + 3][:11]
            read_seq = umi + read_file_content[a + 1]
            read_q = umi_q + read_file_content[a + 3]
            output_fastq.append([read_file_content[a],read_seq,read_file_content[a + 2],read_q])
    
# write to file
    output_file = gzip.open(output_filename, 'wt')
    for entry in output_fastq:
        output_file.write('{}\n'.format(entry[0]))
        output_file.write('{}\n'.format(entry[1]))
        output_file.write('{}\n'.format(entry[2]))
        output_file.write('{}\n'.format(entry[3]))
    output_file.close()
