#!/usr/bin/env python3
# coding: utf-8

"""

generate reverse complement of fastq file

"""

import sys
import os
import subprocess
import tqdm

if __name__ == "__main__":

# parse arguments
    
    if len(sys.argv) < 3:
        print('usage: {} [input fastq] [output reverse complement fastq]'.format(os.path.basename(__file__)))
        exit()

    input_fastq_filename = sys.argv[1]
    reverse_comp_fastq_filename = sys.argv[2]

# count total lines in file
    print('processing {} ...'.format(input_fastq_filename))
    result = subprocess.run(['wc','-l',input_fastq_filename], stdout=subprocess.PIPE)
    total_lines = int(str(result.stdout).split()[0].strip("b'"))
      
    seq_dict = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N':'N', ' ':''}
    
# iterate through fastq
    fastq_reverse_complement = []
    with open(input_fastq_filename, 'r') as input_fastq_file:
        line_counter = 0
        for line in tqdm.tqdm(input_fastq_file,total=total_lines):
            line_counter += 1
            if line_counter == 1:
                fastq_reverse_complement.append(line.strip('\n'))
            if line_counter == 2:
                reverse_complement_sequence = ''
                for basepair in reversed(line.strip('\n')):
                    reverse_complement_sequence += seq_dict[basepair]
                fastq_reverse_complement.append(reverse_complement_sequence)
            if line_counter == 3:
                fastq_reverse_complement.append(line.strip('\n'))
            if line_counter == 4:
                reverse_quality_score = ''
                for basepair in reversed(line.strip('\n')):
                    reverse_quality_score += basepair
                fastq_reverse_complement.append(reverse_quality_score)
                line_counter = 0
                                                
                
# write reverse complement fastq to file
    output_file = open(reverse_comp_fastq_filename, 'w')
    for line in fastq_reverse_complement:
        output_file.write('{}\n'.format(line))
    output_file.close()
