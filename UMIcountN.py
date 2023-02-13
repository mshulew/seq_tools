#!/usr/bin/env python3
# coding: utf-8

"""

determines N (ambiguous base) and length distribution in UMI

"""

import sys
import io
import gzip
from collections import Counter


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
        
    ncontent = {}
    lendist = {}
    umisize = 8
            
    with io.TextIOWrapper(gzip.open(fastq_filename, 'r')) as fastq_file:
        fastq_file_content = fastq_file.read().strip().split('\n')
        for a in range(0, len(fastq_file_content), 4):
            barcode = fastq_file_content[a+1][:umisize]
            numberOfN = Counter(barcode).get("N",0)
            ncontent.setdefault(numberOfN,0)
            ncontent[numberOfN] += 1
            lendist.setdefault(len(barcode),0)
            lendist[len(barcode)] += 1
                
    with open(output_filename, 'w') as output_file:
        output_file.write('{}\t{}\n'.format('Number_of_Ns','UMIs'))
        for a in range(0,umisize+1):
            try:
                output_file.write('{}\t{}\n'.format(a,ncontent[a]))
            except:
                output_file.write('{}\t{}\n'.format(a,0))
        output_file.write('{}\t{}\n'.format('UMI length','UMIs'))
        for a in range(0,umisize+1):
            try:
                output_file.write('{}\t{}\n'.format(a,lendist[a]))
            except:
                output_file.write('{}\t{}\n'.format(a,0))
