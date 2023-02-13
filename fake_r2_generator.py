#!/usr/bin/env python3
# coding: utf-8

"""
generates a fake R2 file to fulfill requirement of DBG pipeline for R2

"""
import sys
import os
import tqdm
import random

if __name__ == "__main__":

# parse arguments

    r1_input_filename = sys.argv[1]
    r2_output_filename = sys.argv[2]
    
    umi_random_library = {1:'A', 2:'G', 3:'C', 4:'T'}
    quality_random_library = {0:'!', 1:'"', 2:'#', 3:'$', 4:'%', 5:'&', 6:'(', 7:'*', 8:'+', 9:',', 10:'-', 11:'.', 12:'/',
                              13:'0', 14:'1', 15:'2', 16:'3', 17:'4', 18:'5', 19:'6', 20:'7', 21:'8', 22:'9', 23:':', 24:';',
                              25:'<', 26:'=', 27:'>', 28:'?', 29:'@', 30:'A', 31:'B', 32:'C', 33:'D', 34:'E', 35:'F', 36:'G',
                              37:'H', 38:'I', 39:')'}

# iterate through r1 and create an r2
    print('iterating through r1 ...')
    with open(r1_input_filename, 'r') as r1_input_file:
        total_rows = 0
        for line in r1_input_file:
            total_rows += 1
        
    r2_data = []
    with open(r1_input_filename, 'r') as r1_input_file:
        line_number = 0
        for line in tqdm.tqdm(r1_input_file, total=total_rows):
            line_number += 1
            if line_number == 1 or line_number == 3:
                r2_data.append(line.strip('\n'))
            elif line_number == 2:
                umi = ''
                for a in range(0,8):
                    umi = umi + umi_random_library[random.randint(1,4)]
                r2_data.append(umi)
            if line_number == 4:
                quality = ''
                for a in range(0,8):
                    quality = quality + quality_random_library[random.randint(1,39)]
                r2_data.append(quality)
                
                line_number = 0
                
    print('writing new r2 to file...')
    output_file = open(r2_output_filename, 'w')
    for entry in r2_data:
        output_file.write('{}\n'.format(entry))
    output_file.close()
                    
                

