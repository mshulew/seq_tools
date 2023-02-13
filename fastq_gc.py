#!/usr/bin/env python3
# coding: utf-8

"""
Generates GC content from fastq file
"""

import os
import sys
import subprocess
import io
import gzip
import shutil

def calcPercentGC(sequence):
    gc = 0
    at = 0
    total = 0
    for bp in sequence:
        if bp == 'G' or bp == 'C':
            gc += 1
            total += 1
        elif bp == 'A' or bp == 'T':
            total += 1
    percentgc = (gc/total) * 100
    
    return percentgc

if __name__ == "__main__":
    
    version_num = 1.0
    
    print('*'*100)
    print('FASTQ: GC Content')
    print('2021 Mark Shulewiz Bio-Rad Laboratories, Inc.')
    print('version ' + str(version_num))
    print('*'*100)
    
# parse arguments

    input_filename = ''
    output_filename = ''

    if len(sys.argv) == 2:
        if sys.argv[1] == '--help' or sys.argv[1] == '--h':
            print('*'*10)
            print('Usage: ' + os.path.basename(__file__) + ' --I [input filename] --O [output filename]')
            print('*'*10)
            exit()
        else:
            print('Error 1 command line: unknown argument {}'.format(sys.argv[1]))
            exit()
    elif len(sys.argv) > 2:
        for a in range(0, len(sys.argv)):
            if sys.argv[a] == '--I':
                input_filename = sys.argv[a + 1]
            elif sys.argv[a] == '--O':
                output_filename = sys.argv[a + 1]
    else:
        print('Error 2: something is missing in the commandline')
        exit()
        
    if input_filename == '':
        print('Error 3: no input file supplied')
        exit()
        
    if output_filename == '':
        print('Error 3: no output file supplied')
        exit()
        
    if '.gz' not in input_filename:
        print('Error 4: Input file {} does not end in .gz'.format(input_filename[0]))
        exit()
        
 # create temporary directory
    tmp_dir = 'tmp'
    if os.path.isdir(tmp_dir):
        shutil.rmtree(tmp_dir)
    os.makedirs(tmp_dir)
        
# process fastq
            
# unzip fastq
    print('unzipping file {} ...'.format(input_filename))
    
    cmd = 'pigz -d -c {} | paste - - - - > {}'.format(input_filename,tmp_dir + '/unzipped.fastq')
    process = subprocess.Popen(cmd, shell=True, executable='/bin/bash')
    process.wait()
    
# importing unzipped fastq 
    print('importing unzipped file ...')
    sequence_list = []
    with open(tmp_dir + '/unzipped.fastq', 'r') as input_file:
        for line in input_file:
            line = line.strip('\n').split('\t')
            sequence_list.append(line[1])
            
# count reads by GC content
    print('calculating percent gc per read ...')
    gc_dist = {}
    for sequence in sequence_list:
        if any(bp in ["A","C","G","T"] for bp in sequence):
            rounded_gc = round(calcPercentGC(sequence))
            gc_dist.setdefault(rounded_gc,0)
            gc_dist[rounded_gc] += 1         

# sort read
    print('sorting and writing to file ...')
    gc_dist_as_list = []
    for percent_gc,counts in gc_dist.items():
        gc_dist_as_list.append([percent_gc,counts])
    gc_dist_sorted = sorted(gc_dist_as_list) 
    
# write to file
    with open(output_filename, 'w') as outputfile:
        outputfile.write('{}\t{}\n'.format("percent_gc","reads"))
        for percent_gc in gc_dist_sorted:
            outputfile.write('{}\t{}\n'.format(percent_gc[0],percent_gc[1]))
            
# delete temporary directory
        shutil.rmtree(tmp_dir)
