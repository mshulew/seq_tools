#!/usr/bin/env python3
# coding: utf-8

"""
Extracts UMI from Qiagen miRNeasy libraries and creates a new R2 with UMI alone
"""

import os
import sys
import subprocess
import shutil
import regex

if __name__ == "__main__":
    
    version_num = 1.0
    
    print('*'*100)
    print('Qiagen miRNeasy library processor')
    print('2021 Mark Shulewiz Bio-Rad Laboratories, Inc.')
    print('version ' + str(version_num))
    print('*'*100)
    
# parse arguments
    input_filenames = []
    output_dir = os.getcwd()

    if len(sys.argv) == 2:
        if sys.argv[1] == '--help' or sys.argv[1] == '--h':
            print('*'*10)
            print('Usage: ' + os.path.basename(__file__) + ' --I [input files] --O [output directory]')
            print('*'*10)
            exit()
        else:
            print('Error 1 command line: unknown argument {}'.format(sys.argv[1]))
            exit()
    elif len(sys.argv) > 2:
        input_trigger = False
        for a in range(0, len(sys.argv)):
            if input_trigger == True:
                if sys.argv[a][0:2] == '--':
                    input_trigger = False
                else:
                    input_filenames.append(sys.argv[a])
            if input_trigger == False:
                if a < len(sys.argv) - 1:
                    if sys.argv[a] == '--I':
                        input_trigger = True
                    elif sys.argv[a] == '--O':
                        output_dir = sys.argv[a + 1]
    else:
        print('Error 2: something is missing in the commandline')
        exit()
        
# check for output directory and create if necessary
    if output_dir != os.getcwd():
        if not os.path.isdir(output_dir):
            os.makedirs(output_dir)
            
    print('Output directory: {}'.format(output_dir))
    print('*'*10)
                
# check input file names for .gz extension and 'R1' in name
    for input_filename in input_filenames:
        if '.gz' not in input_filename:
            print('Error 3: Input filename {} does not end in .gz'.format(input_filename[0]))
            exit()
        if 'R1' not in input_filename:
            print('Error 4: Input filename {} does not contain "R1"'.format(input_filename[0]))
            exit()
            
    print('Processing files')
    print('*'*5)
    
    for input_filename in input_filenames:
        
        print('processing {}...'.format(input_filename))
        
 # create temporary directory
        tmp_dir = output_dir + '/tmp'
        if os.path.isdir(tmp_dir):
            shutil.rmtree(tmp_dir)
        os.makedirs(tmp_dir)
            
# unzip fastq
        print('unzipping file...')
    
        base_filename = input_filename.split('.gz')[0]
        if '/' in base_filename:
            base_filename = base_filename.split('/')[len(base_filename.split('/')) - 1]
    
        cmd = 'pigz -d -c {} > {}'.format(input_filename,tmp_dir + '/' + base_filename)
        process = subprocess.Popen(cmd, shell=True, executable='/bin/bash')
        process.wait()
        
# import fastq
        print('importing fastq...')
        
        with open(tmp_dir + '/' + base_filename, 'r') as input_file:
            input_file_content = input_file.read().strip().split('\n')
            
# process reads
        print('processing reads...')
    
# adapter sequence, umi length
        adapter='AACTGTAGGCACCATCAAT' 
        umi_length = 12
        
# R1 and R2 output filenames
        output_r1_filename = output_dir + '/' + base_filename
        output_r2_filename = output_dir + '/' + base_filename.replace('R1','R2')
        
# open outputfiles to write to
        with open(output_r1_filename, 'w') as outputfile_r1, open(output_r2_filename, 'w') as outputfile_r2:
        
# iterate through fastq reads
            for a in range(0, len(input_file_content), 4):
        
# search for adapter (allows up to 2 mismatches)
                adapter_match = regex.search("(?b)(" + adapter + "){e<=2}",input_file_content[a+1],partial=False)
    
                if adapter_match:
            
# get umi sequence
                    umi = input_file_content[a+1][adapter_match.end():(adapter_match.end()+umi_length)]
    
# create R2 read names
                    readname_r1 = input_file_content[a]
                    if readname_r1.split(' ')[1][0] == '1':
                        readname_split = readname_r1.split(' ')
                        readname_split[1] = '2' + readname_split[1][1:]
                        readname_r2 = ' '.join(readname_split)
                        
# write R1 file    
                    outputfile_r1.write('{}\n'.format(readname_r1))
                    outputfile_r1.write('{}\n'.format(input_file_content[a+1][:adapter_match.start()]))
                    outputfile_r1.write('{}\n'.format(input_file_content[a+2]))
                    outputfile_r1.write('{}\n'.format(input_file_content[a+3][:adapter_match.start()]))
                                            
# write R2 file
                    outputfile_r2.write('{}\n'.format(readname_r2))
                    outputfile_r2.write('{}\n'.format(umi))
                    outputfile_r2.write('{}\n'.format(input_file_content[a+2]))
                    outputfile_r2.write('{}\n'.format(input_file_content[a+3][adapter_match.end():(adapter_match.end()+umi_length)]))
            
# zip output files
        print('compressing output files...')
        cmd = 'gzip {};gzip {}'.format(output_r1_filename,output_r2_filename)
        process = subprocess.Popen(cmd, shell=True, executable='/bin/bash')
        process.wait()
                                          
# delete temporary directory
        shutil.rmtree(tmp_dir)
