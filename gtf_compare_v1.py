#!/usr/bin/env python3
# coding: utf-8

"""
compares 2 gtf files and outputs differences
"""

import sys
import tqdm
import multiprocessing

def compare_gtfs(gtf_data):
    dataset = gtf_data[0]
    chromosome = gtf_data[1]
    line_number = gtf_data[2]
    line = gtf_data[3]
    
# compare gtf1 to gtf2 
    if dataset == 'gtf1':
        gtf2_chromosome = gtf2_data[chromosome]
        match_found = False
        for entry in gtf2_chromosome:
            if line == entry[1]:
                match_found = True
                break
                
        if match_found == False:
            comparison_info = []
            comparison_info.append(['1',line_number,'in gtf1 not in gtf2',line])
            return comparison_info
        
# compare gtf2 to gtf1 
    if dataset == 'gtf2':
        gtf1_chromosome = gtf1_data[chromosome]
        match_found = False
        for entry in gtf1_chromosome:
            if line == entry[1]:
                match_found = True
                break
                
        if match_found == False:
            comparison_info = []
            comparison_info.append(['2',line_number,'in gtf2 not in gtf1',line])
            return comparison_info
        
if __name__ == "__main__":

    version_number = 1 
    
# parse arguments
    if len(sys.argv) == 2:
        if sys.argv[1] == '-h' or sys.argv[1] == '-help':
            print('*'*100)
            print(os.path.basename(__file__) + ' version ' + str(version_num))
            print('Usage: ' + os.path.basename(__file__) + ' [gtf filename 1] [gtf filename 2] > [output filename]')
            exit()
    elif len(sys.argv) == 4:
        gtf1_filename = sys.argv[1]
        gtf2_filename = sys.argv[2]
        output_filename = sys.argv[3]
    else:
        print('missing input and/or output file(s)')
        exit()
        
# import first gtf file
    print('indexing gtf 1...')
    gtf1_data = {}
    with open(gtf1_filename, 'r') as gtf1_file:
        line_num = 0
        for line in gtf1_file:
            line_num += 1
            if line[0] != '#':
                chromosome = line.split('\t')[0]
                gtf1_data.setdefault(chromosome, [])
                gtf1_data[chromosome].append([str(line_num),line.strip('\n')])

        
# import second gtf file
    print('indexing gtf 2...')
    gtf2_data = {}
    with open(gtf2_filename, 'r') as gtf2_file:
        line_num = 0
        for line in gtf2_file:
            line_num += 1
            if line[0] != '#':
                chromosome = line.split('\t')[0]
                gtf2_data.setdefault(chromosome, [])
                gtf2_data[chromosome].append([str(line_num),line.strip('\n')])

                                                   
# preparing for multiprocessing
    multiprocessing_input = []
    for chromosome in gtf1_data.keys():
        for line in gtf1_data[chromosome]:
            multiprocessing_input.append(['gtf1',chromosome,line[0],line[1]])
    for chromosome in gtf2_data.keys():
        for line in gtf2_data[chromosome]:
            multiprocessing_input.append(['gtf2',chromosome,line[0],line[1]])
        
# multiprocessing 
    print('comparing gtf1 to gtf2...')
    builtin_outputs = map(compare_gtfs, multiprocessing_input)
    pool_size = multiprocessing.cpu_count()
    pool = multiprocessing.Pool(processes=pool_size)
    
    pool_outputs = list(tqdm.tqdm(pool.imap(compare_gtfs, multiprocessing_input), total=len(multiprocessing_input)))
                
    pool.close()  # no more tasks
    pool.join()  # wrap up current tasks 
    
# prepare output data
    output_data = []
    for top_level_entry in pool_outputs:
        if top_level_entry is not None:
            for entry in top_level_entry:
                output_data.append([entry[0],entry[1],entry[2],entry[3]])
            
    output_data = sorted(output_data, key=lambda x:(x[0], x[1]))

# write to file
    output_file = open(output_filename, 'w')
    for entry in output_data:
        output_file.write('{}\t{}\t{}\t{}\r\n'.format(entry[0],entry[1],entry[2],entry[3]))
    output_file.close()
