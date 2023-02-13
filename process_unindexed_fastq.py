#!/usr/bin/env python3
# coding: utf-8

"""

processes raw illumina sequence data that did not parsed by index
places all sequences with same indices into a unique fastq file
creates R1 & R2 files

"""

import sys
import os
import subprocess
import tqdm

if __name__ == "__main__":

# parse arguments
    
    if len(sys.argv) < 3:
        print('usage: {} [input r1 filename] [input r2 filename] [output dir]'.format(os.path.basename(__file__)))
        exit()

    r1_file = sys.argv[1]
    r2_file = sys.argv[2]
    if len(sys.argv) == 4:
        output_parent_dir = sys.argv[3]
    elif len(sys.argv) == 3:
        output_parent_dir = os.getcwd()
    else:
        print('missing a parameter')
        exit()   

# create output dir
    folder_num = 1
    for foldername in os.listdir(output_parent_dir):
        if 'processed_files_' in str(foldername):
            folder_number = int(foldername.split('processed_files_')[1])
            if folder_number >= folder_num:
                folder_num = folder_number + 1

    output_dir = output_parent_dir + '/processed_files_' + str(folder_num)
    os.makedirs(output_dir)

# count total lines in file
    print('processing {} ...'.format(r1_file))
    result = subprocess.run(['wc','-l',r1_file], stdout=subprocess.PIPE)
    total_lines = int(str(result.stdout).split()[0].strip("b'"))

# iterate through r1 file
    with open(r1_file, 'r') as r1_input_file, open(r2_file, 'r') as r2_input_file:
        sequences = {}
        sequence_count = 0
        line_num = 0
        for line_r1,line_r2 in tqdm.tqdm(zip(r1_input_file,r2_input_file), total=total_lines):
            line_num += 1
            if line_num == 1:
                first_line_r1 = line_r1.strip('\n')
                first_line_r2 = line_r2.strip('\n')
                obs_index = line_r1.split()[1].split(':')[3].strip('\n')               
                
            if line_num == 2:
                second_line_r1 = line_r1.strip('\n')
                second_line_r2 = line_r2.strip('\n')
            if line_num == 3:
                third_line_r1 = line_r1.strip('\n')
                third_line_r2 = line_r2.strip('\n')
            if line_num == 4:
                fourth_line_r1 = line_r1.strip('\n')
                fourth_line_r2 = line_r2.strip('\n')

                sequences.setdefault(obs_index,[])
                sequences[obs_index].append([first_line_r1,second_line_r1,third_line_r1,fourth_line_r1,first_line_r2,second_line_r2,third_line_r2,fourth_line_r2])
                sequence_count += 1
                    
                line_num = 0

                if sequence_count == 10000000:
# save to file every 100,000 entries
                    for seq_index,entries in sequences.items():
                        output_file_r1 = open(output_dir + '/' + seq_index + '_S1_L001_R1_001.fastq', 'a')
                        output_file_r2 = open(output_dir + '/' + seq_index + '_S1_L001_R2_001.fastq', 'a')
                        for entry in entries:
                            first_line_r1 = entry[0]
                            second_line_r1 = entry[1]
                            third_line_r1 = entry[2]
                            fourth_line_r1 = entry[3]
                            first_line_r2 = entry[4]
                            second_line_r2 = entry[5]
                            third_line_r2 = entry[6]
                            fourth_line_r2 = entry[7]
                            output_file_r1.write('{}\r\n'.format(first_line_r1))
                            output_file_r1.write('{}\r\n'.format(second_line_r1))
                            output_file_r1.write('{}\r\n'.format(third_line_r1))
                            output_file_r1.write('{}\r\n'.format(fourth_line_r1))
                            output_file_r2.write('{}\r\n'.format(first_line_r2))
                            output_file_r2.write('{}\r\n'.format(second_line_r2))
                            output_file_r2.write('{}\r\n'.format(third_line_r2))
                            output_file_r2.write('{}\r\n'.format(fourth_line_r2))
                        output_file_r1.close()
                        output_file_r2.close()
                    sequence_count = 0
                    sequences = {}
                    
# save remaining entries to file
        for seq_index,entries in sequences.items():
            output_file_r1 = open(output_dir + '/' + seq_index + '_S1_L001_R1_001.fastq', 'a')
            output_file_r2 = open(output_dir + '/' + seq_index + '_S1_L001_R2_001.fastq', 'a')
            for entry in entries:
                first_line_r1 = entry[0]
                second_line_r1 = entry[1]
                third_line_r1 = entry[2]
                fourth_line_r1 = entry[3]
                first_line_r2 = entry[4]
                second_line_r2 = entry[5]
                third_line_r2 = entry[6]
                fourth_line_r2 = entry[7]
                output_file_r1.write('{}\r\n'.format(first_line_r1))
                output_file_r1.write('{}\r\n'.format(second_line_r1))
                output_file_r1.write('{}\r\n'.format(third_line_r1))
                output_file_r1.write('{}\r\n'.format(fourth_line_r1))
                output_file_r2.write('{}\r\n'.format(first_line_r2))
                output_file_r2.write('{}\r\n'.format(second_line_r2))
                output_file_r2.write('{}\r\n'.format(third_line_r2))
                output_file_r2.write('{}\r\n'.format(fourth_line_r2))
            output_file_r1.close()
            output_file_r2.close()
        
    print('finished')
