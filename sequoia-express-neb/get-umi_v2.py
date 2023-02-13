#!/usr/bin/env python3
# coding: utf-8

"""

Transfers 11 bp UMI to R2
For use with NEB UMI adapters
version 2: faster

"""

import os
import sys
import io
import shutil
import subprocess


if __name__ == "__main__":
    
    version_num = 2.0
    
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
        
# check for .gz extension
    if '.fastq.gz' not in read_filename:
        print('invalid read file - does not end in fastq.gz')
        exit()
    if '.fastq.gz' not in index_filename:
        print('invalid index file - does not end in fastq.gz')
        exit()
    if '.fastq.gz' not in output_filename:
        print('invalid output filename - does not end in fastq.gz')
        exit()
        
# create temporary directory
    tmp_dir = os. getcwd() + '/tmp'
    if os.path.isdir(tmp_dir):
        shutil.rmtree(tmp_dir)
    os.makedirs(tmp_dir)
    
# unzip fastq
    print('unzipping file...')
    
    cmd = 'paste <(pigz -d -c {} | paste - - - -) <(pigz -d -c {} | paste - - - -) > {}'.format(read_filename[0],index_filename[1],tmp_dir + '/compressed.tsv')
    process = subprocess.Popen(cmd, shell=True, executable='/bin/bash')
    process.wait()
    
# import reads
    print('adding UMIs to R2...')
    
    output_fastq = []
    for line in open(tmp_dir + '/compressed.tsv', 'r'):
        line = line.strip('\n').split('\t')
        umi = line[5][:11]
        umi_q = line[7][:11]
        read_seq = umi + read_file_content[1]
        read_q = umi_q + read_file_content[3]
        output_fastq.append([line[0],read_seq,line[2],read_q])
    
# write to file
    print('writing to file...')
    
    output_filename = output_filename.split('.gz')[0]
    
    with open(output_filename, 'w') as output_file:
        for entry in output_fastq:
            output_file.write('{}\n'.format(entry[0]))
            output_file.write('{}\n'.format(entry[1]))
            output_file.write('{}\n'.format(entry[2]))
            output_file.write('{}\n'.format(entry[3]))
            
# zip outputfile
    print('zipping output file...')
    cmd = 'pigz {}'.format(output_filename)
    process = subprocess.Popen(cmd, shell=True, executable='/bin/bash')
    process.wait()
    
# delete temp directory
    shutil.rmtree(tmp_dir)
