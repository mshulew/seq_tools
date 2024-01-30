#!/usr/bin/env python3
# coding: utf-8

"""
Keeps only reads in a pair of FASTQ files that are included in a supplied list
FASTQ files must be uncompressed
"""

import sys
import os

def main():
# parse arguments
    r1_path = 'NONE'
    r2_path = 'NONE'
    keep_path = 'NONE'
    output_path = 'NONE'

    if len(sys.argv) == 2:
        if sys.argv[1] == '--help' or sys.argv[1] == '--h':
            print('*'*10)
            print('Remove reads from FASTQ files')
            print('Usage: ' + os.path.basename(__file__) + ' --r1 [R1 FASTQ] --r2 [R2 FASTQ] --keep [list of reads to keep] --o [output path')
            print('*'*10)
            exit()
        else:
            print('Error 1 command line: unknown argument {}'.format(sys.argv[1]))
            exit()
    elif len(sys.argv) > 2:
        for a in range(len(sys.argv)):
            if sys.argv[a].lower() == '--r1':
                r1_path = sys.argv[a + 1]
            elif sys.argv[a].lower() == '--r2':
                r2_path = sys.argv[a + 1]
            elif sys.argv[a].lower() == '--keep':
                keep_path = sys.argv[a + 1]
            elif sys.argv[a].lower() == '--o':
                output_path = sys.argv[a + 1]
    else:
        print('Error 2: something is missing in the commandline')
        exit()

    if r1_path == 'NONE':
        print('No R1 file was supplied')
        exit()
    if r2_path == 'NONE':
        print('No R2 file was supplied')
        exit()
    if keep_path == 'NONE':
        print('No list of files to keep was supplied')
        exit()
    if output_path == 'NONE':
        print('No output path was supplied')
        exit()
        
    r1_name = r1_path.split('/')[-1].split('.fastq')[0] + '_reads_removed.fastq'
    r2_name = r2_path.split('/')[-1].split('.fastq')[0] + '_reads_removed.fastq'
    r1_output_path = output_path + '/' + r1_name
    r2_output_path = output_path + '/' + r2_name

# import list of reads to keep
    print('importing list of reads to keep...')
    keep = []
    with open(keep_path, 'r') as keep_file:
        for line in keep_file:
            keep.append(line.splitlines()[0])
                
# import fastq files
    print('importing FASTQ files...')
    input_fastqs = {}
    with open(r1_path, 'r') as r1,open(r2_path, 'r') as r2:
        for lane1,lane2 in zip(r1,r2):
            
            r1_name = lane1.strip('\n')
            r1_seq = next(r1).strip('\n')
            r1_sep = next(r1).strip('\n')
            r1_qual = next(r1).strip('\n')
            r2_name = lane2.strip('\n')
            r2_seq = next(r2).strip('\n')
            r2_sep = next(r2).strip('\n')
            r2_qual = next(r2).strip('\n')

            assert r1_name.split()[0] == r2_name.split()[0]
            assert len(r1_seq) == len(r1_qual)
            assert len(r2_seq) == len(r2_qual)
                          
            input_fastqs[r1_name.split()[0][1:]] = [r1_name,r1_seq,r1_sep,r1_qual,r2_name,r2_seq,r2_sep,r2_qual]

# create new fastq
    print('creating new fastq...')
    kept = []
    for i,read in enumerate(keep):
        kept.append(input_fastqs[read])
                    
    print('input reads: {} reads kept {} '.format(len(list(input_fastqs.keys())),len(kept)))

# write to file
    print('writing to file...')
    
    with open(r1_output_path,'w') as r1_output_file:
        for fastq in kept:
            r1_output_file.write('{}\n'.format(fastq[0]))
            r1_output_file.write('{}\n'.format(fastq[1]))
            r1_output_file.write('{}\n'.format(fastq[2]))
            r1_output_file.write('{}\n'.format(fastq[3]))

   
    with open(r2_output_path,'w') as r2_output_file:
        for fastq in kept:
            r2_output_file.write('{}\n'.format(fastq[4]))
            r2_output_file.write('{}\n'.format(fastq[5]))
            r2_output_file.write('{}\n'.format(fastq[6]))
            r2_output_file.write('{}\n'.format(fastq[7]))

if __name__ == "__main__":

    main()
            
            
                                  


    

    

