#!/usr/bin/env python3
# coding: utf-8

"""

determine length of each read in a fastq file

"""

import sys
import io
import gzip


if __name__ == "__main__":

# parse arguments

    if len(sys.argv) == 2:
        if sys.argv[1] == '--help' or sys.argv[1] == '--h':
            print('*'*10)
            print('Usage: python3.6 ' + os.path.basename(__file__) + ' [fastq file] [output file]')
            print('*'*10)
            exit()
        else:
            print('Something is missing')
            exit()
    elif len(sys.argv) == 3:
            fastq_filename = sys.argv[1]
            output_filename = sys.argv[2]
    else:
        print('something is missing')
        exit()
        
    lengthreads = {}
            
    with io.TextIOWrapper(gzip.open(fastq_filename, 'r')) as fastq_file:
        fastq_file_content = fastq_file.read().strip().split('\n')
        for a in range(0, len(fastq_file_content), 4):
            readlength = len(fastq_file_content[a+1])
            lengthreads.setdefault(readlength,0)
            lengthreads[readlength] += 1
            
    sorteddata = []
    for length,reads in lengthreads.items():
        sorteddata.append([length,reads])
    sorteddata = sorted(sorteddata)
    
    output_file = open(output_filename, 'w')
    for entry in sorteddata:
        output_file.write('{}\t{}\n'.format(entry[0],entry[1]))
    output_file.close()

    
    
                
        
   
    
    
