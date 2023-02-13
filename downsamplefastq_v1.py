#!/usr/bin/env python3
# coding: utf-8

"""

Downsamples fastq file(s)
Accepts single R1 or R1/R2 pairs

"""

import os
import sys
import io
import gzip
import random
import time

def generatename(inputfilename,sizename):
    if '_R1' in inputfilename:
        delimiter = '_R1'
    elif '_R2' in inputfilename:
        delimiter = '_R2'
        
    splitloc = len(inputfilename.split('/')) - 1
    basefilename = inputfilename.split('/')[splitloc].split('.')[0]
    filename = basefilename.split(delimiter)[0] + '_' + sizename + delimiter + basefilename.split(delimiter)[1] +'.fastq.gz'

    return filename

if __name__ == "__main__":
    
    version_num = 1.1
    
    start_time = time.time()
    
    print('*'*100)
    print('FASTQ downsampler')
    print('2020 Mark Shulewiz Bio-Rad Laboratories, Inc.')
    print('version ' + str(version_num))
    print('*'*100)
    
# parse arguments
    cli_input_filenames = []
    input_filenames = []
    output_dir = os.getcwd()
    input_number_of_reads = 0

    if len(sys.argv) == 2:
        if sys.argv[1] == '--help' or sys.argv[1] == '--h':
            print('*'*10)
            print('Usage: python3.6 ' + os.path.basename(__file__) + ' --I [input files] --O [output directory] --size [number of reads]')
            print('*'*10)
            exit()
        else:
            print('ERROR with command line: unknown argument {}'.format(sys.argv[1]))
            exit()
    elif len(sys.argv) > 2:
        input_trigger = False
        for a in range(0, len(sys.argv)):
            if input_trigger == True:
                if sys.argv[a][0:2] == '--':
                    input_trigger = False
                else:
                    cli_input_filenames.append(sys.argv[a])
            if input_trigger == False:
                if a < len(sys.argv) - 1:
                    if sys.argv[a] == '--I':
                        input_trigger = True
                    elif sys.argv[a] == '--O':
                        output_dir = sys.argv[a + 1]
                    elif sys.argv[a].lower() == '--size':
                        input_number_of_reads = sys.argv[a + 1]
    else:
        print('something is missing')
        exit()
        
# process number of reads
    if input_number_of_reads != 0:
        try:
            number_of_reads = int(input_number_of_reads)
        except:
            if input_number_of_reads[len(input_number_of_reads) - 1] == 'M':
                number_of_reads = int(input_number_of_reads[:len(input_number_of_reads) - 1])*1000000
            elif input_number_of_reads[len(input_number_of_reads) - 1] == 'K':
                number_of_reads = int(input_number_of_reads[:len(input_number_of_reads) - 1])*1000
            else:
                print('ERROR: invalid entry for size: {}'.format(input_number_of_reads))
                exit()
    else:
        print('ERROR: invalid entry for size: {}'.format(input_number_of_reads))
        exit()
        
# format output number of reads  
    if str(number_of_reads)[-3:] == '000':
        if str(number_of_reads)[-6:] == '000000':
            readable_reads = str(number_of_reads)[:-6] + 'M'
        else:
            readable_reads = str(number_of_reads)[:-3] + 'K'
    else:
        readable_reads = str(number_of_reads)
       
# check for output directory and create if necessary
    if output_dir != os.getcwd():
        if not os.path.isdir(output_dir):
            os.makedirs(output_dir)
            
    print('Output directory: {}'.format(output_dir))
    print('*'*10)
    
    print('Number of output reads: {}'.format(readable_reads))
    print('*'*10)
        
# process input filenames
    input_filenames = []
    if len(cli_input_filenames) == 0:
        print('missing input files')
        exit()
    elif len(cli_input_filenames) == 1:
        if 'R1' not in cli_input_filenames[0]:
            print('ERROR with input filename: {}'.format(cli_input_filenames[0]))
        else:
            input_filenames.append([cli_input_filenames[0]])
    elif len(cli_input_filenames) > 1:
        for b in range(0, len(cli_input_filenames)):
            if 'R1' in cli_input_filenames[b]:
                if b < len(cli_input_filenames) - 1:
                    if 'R2' in cli_input_filenames[b + 1]:
                        if cli_input_filenames[b].replace('R1','R2') == cli_input_filenames[b + 1]:
                            input_filenames.append([cli_input_filenames[b],cli_input_filenames[b + 1]])
                    else:
                        input_filenames.append([cli_input_filenames[b]])
                else:
                    input_filenames.append([cli_input_filenames[b]])
            elif 'R2' not in cli_input_filenames[b]:
                print('ERROR with input filename: {}'.format(cli_input_filenames[b]))    
                                             
    print('Processing files')
    print('*'*5)
    
    input_fastq = []
    
    for input_filename in input_filenames:
# process single reads
        if len(input_filename) == 1:
            if '.fastq.gz' not in input_filename[0]:
                print('Input file must end in fastq.gz')
                exit()
            print('Processing R1 file: {}'.format(input_filename[0]))
            
            with io.TextIOWrapper(gzip.open(input_filename[0], 'r')) as r1_file:
                r1_file_content = r1_file.read().strip().split('\n')
                for a in range(0, len(r1_file_content), 4):
                    input_fastq.append([r1_file_content[a],r1_file_content[a+1],r1_file_content[a+2],r1_file_content[a+3]])
          
# shuffle list of reads
            if number_of_reads > len(input_fastq):
                print('Sample size greater than number of reads')
                print('Number of reads: {}'.format(len(input_fastq)))
            else:
                if number_of_reads == 0:
                    number_of_reads = len(input_fastq)
                input_number_of_reads = str(len(input_fastq))
                shuffled_fastq = random.sample(input_fastq, number_of_reads)
                
# generate output filename
                r1_output_filename = output_dir + '/' + generatename(input_filename[0],readable_reads)
                print('R1 output filename: {}'.format(r1_output_filename))
                print('*'*10)
    
# write to file
                output_r1_file = gzip.open(r1_output_filename, 'wt')
                for entry in shuffled_fastq:
                    output_r1_file.write('{}\n'.format(entry[0]))
                    output_r1_file.write('{}\n'.format(entry[1]))
                    output_r1_file.write('{}\n'.format(entry[2]))
                    output_r1_file.write('{}\n'.format(entry[3]))
                output_r1_file.close()
            
                
# process read pairs
        elif len(input_filename) == 2:
            if '.fastq.gz' not in input_filename[0] or 'fastq.gz' not in input_filename[1]:
                print('Input file must end in fastq.gz')
                exit()
            print('Processing R1 file: {}'.format(input_filename[0]))
            print('Processing R2 file: {}'.format(input_filename[1]))
            
            with io.TextIOWrapper(gzip.open(input_filename[0], 'r')) as r1_file,io.TextIOWrapper(gzip.open(input_filename[1], 'r')) as r2_file:
                r1_file_content = r1_file.read().strip().split('\n')
                r2_file_content = r2_file.read().strip().split('\n')
                for a in range(0, len(r1_file_content), 4):
                    input_fastq.append([r1_file_content[a],r1_file_content[a+1],r1_file_content[a+2],r1_file_content[a+3],r2_file_content[a],r2_file_content[a+1],r2_file_content[a+2],r2_file_content[a+3]])
                    
 # shuffle list of reads
            if number_of_reads > len(input_fastq):
                print('Sample size greater than number of reads')
                print('Number of reads: {}'.format(len(input_fastq)))
            else:
                if number_of_reads == 0:
                    number_of_reads = len(input_fastq)
                input_number_of_reads = str(len(input_fastq))
           
                shuffled_fastq = random.sample(input_fastq, number_of_reads)
            
# generate output filenames
                r1_output_filename = output_dir + '/' + generatename(input_filename[0],readable_reads)
                r2_output_filename = output_dir + '/' + generatename(input_filename[1],readable_reads)
                print('R1 output filename: {}'.format(r1_output_filename))
                print('R2 output filename: {}'.format(r2_output_filename))
                print('*'*10)
    
# write to file
                output_r1_file = gzip.open(r1_output_filename, 'wt')
                output_r2_file = gzip.open(r2_output_filename, 'wt')
                for entry in shuffled_fastq:
                    output_r1_file.write('{}\n'.format(entry[0]))
                    output_r1_file.write('{}\n'.format(entry[1]))
                    output_r1_file.write('{}\n'.format(entry[2]))
                    output_r1_file.write('{}\n'.format(entry[3]))
                    output_r2_file.write('{}\n'.format(entry[4]))
                    output_r2_file.write('{}\n'.format(entry[5]))
                    output_r2_file.write('{}\n'.format(entry[6]))
                    output_r2_file.write('{}\n'.format(entry[7]))
                output_r1_file.close()
                output_r2_file.close()
                
    elapsed_time_sec = time.time() - start_time
    elapsed_time_minutes = (elapsed_time_sec - (elapsed_time_sec % 60))/60
    print('finished, elapsed time: {} min {} sec'.format(elapsed_time_minutes,int(elapsed_time_sec % 60)))
