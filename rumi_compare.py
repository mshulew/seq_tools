#!/usr/bin/env python3
# coding: utf-8

"""
Compares pre- and post-rumi deduplication BAM files
Identifies removed alignments and check that they are valid PCR duplicates
"""

import sys
import os
import shutil
import subprocess
import textdistance
import multiprocessing
import tqdm

def process_alignments(input_entry):
    entry_number = input_entry[0]
    entry_data = input_entry[1]
    alignment_data = []
# entry_data = [chromosome, startpost [alignment(s)]]    
# alignments = name umi retained/removed
# output should be: chromosome start name umi retained/removed best_match_umi hamming distance 

# iterate through each position in entry_data
    for position in entry_data:
        position_chromosome = position[0]
        position_start = position[1]
        if len(position[2]) == 1:
            alignment = position[2][0]
            if alignment[2] == 'retained':
                alignment_data.append([position_chromosome,position_start,alignment[0],alignment[1],alignment[2],'N/A','N/A'])
            else:
                alignment_data.append([position_chromosome,position_start,alignment[0],alignment[1],alignment[2],'N/A','single alignment at site'])
        else:
            for a in range(0, len(position[2])):
                alignment = position[2][a]
                if alignment[2] == 'retained':
                    alignment_data.append([position_chromosome,position_start,alignment[0],alignment[1],alignment[2],'N/A','N/A'])
                else:
                    query_umi = alignment[1]
                    hamming_results = ()
# if alignment is removed, compared against retained alignments
                    for b in range(0,len(position[2])):
                        if a != b:
                            if position[2][b][2] == 'retained':
                                test_umi = position[2][b][1]
                                hamming_distance = textdistance.hamming(query_umi,test_umi)
# if hamming distance is 0 keep result and move along
                                if hamming_distance == 0:
                                    hamming_results = (test_umi,hamming_distance)
                                    break
                                else:
                                    if len(hamming_results) == 0:
                                        hamming_results = (test_umi,hamming_distance)
                                    else:
                                        if hamming_distance < hamming_results[1]:
                                            hamming_results = (test_umi,hamming_distance)
                    if hamming_results == ():
                        hamming_results = ('N/A','all reads at site removed')
                    alignment_data.append([position_chromosome,position_start,alignment[0],query_umi,alignment[2],hamming_results[0],hamming_results[1]])

# write chromosome data to file
    with open(outputdir + '/' + str(entry_number) + '.tsv','w') as outputfile:
        for entry in alignment_data:
            outputfile.write('{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(entry[0],entry[1],entry[2],entry[3],entry[4],entry[5],entry[6]))
    
    return alignment_data


if __name__ == "__main__":
    
    version_num = 0.1

    print('*'*100)
    print('rumi deduplication checker')
    print('2022 Mark Shulewiz Bio-Rad Laboratories, Inc.')
    print('version ' + str(version_num))
    print('*'*100)

    inputdir = ''
    outputdir = 'rumiDedupCheckerOutput'

    if len(sys.argv) == 2:
        if sys.argv[1] == '--h' or sys.argv[1] == '--help':
            print('Usage: ' + os.path.basename(__file__) + ' --i [input directory] --o [output directory]')
            print('*' * 100)
            exit()
    elif len(sys.argv) > 2:
        for a in range(0, len(sys.argv)):
            if sys.argv[a] == '--i':
                inputdir = sys.argv[a+1]
            elif sys.argv[a] == '--o':
                outputdir = sys.argv[a+1]
    else:
        print('missing parameters/options')
        exit()

# create output directory
    if not os.path.isdir(outputdir):
        os.makedirs(outputdir)

# create temp directory
    tmpdir = outputdir + '/tmp'
    if not os.path.isdir(tmpdir):
        os.makedirs(tmpdir)

# convert BAM files to SAM and copy to temp directory
    print('converting BAM files to SAM...')
    if os.path.isfile(inputdir + '/umiTagging/Aligned.sortedByCoord.tagged.bam'):
        cmd = 'samtools view {} > {}'.format(inputdir + '/umiTagging/Aligned.sortedByCoord.tagged.bam',tmpdir + '/pre.sam')
        process = subprocess.Popen(cmd, shell=True, executable='/bin/bash')
        process.wait()
    else:
        print('{} not found'.format(inputdir + '/umiTagging/Aligned.sortedByCoord.tagged.bam'))
        exit()
    if os.path.isfile(inputdir + '/dedup/rumi_dedup.sort.bam'):
        cmd = 'samtools view {} > {}'.format(inputdir + '/dedup/rumi_dedup.sort.bam',tmpdir + '/post.sam')
        process = subprocess.Popen(cmd, shell=True, executable='/bin/bash')
        process.wait()
    else:
        print('{} not found'.format(inputdir + '/dedup/rumi_dedup.sort.bam'))
        exit()

# import pre-deduplicated BAM file
    print('importing pre-deduplicated BAM file...')
    prededup_data = []
    for line in open(tmpdir + '/pre.sam', 'r'):
        line = line.strip('\n').split('\t')
        name = line[0]
        chrom = line[2]
        start = int(line[3])
        umi = ''
        for element in line:
            if 'XU:Z:' in element:
                umi = element.split('XU:Z:')[1]
        prededup_data.append([name,chrom,start,umi])

# import post-deduplication BAM file
    print('importing post-deduplication BAM file...')
    postdedup_data = []
    for line in open(tmpdir + '/post.sam', 'r'):
        line = line.strip('\n').split('\t')
        name = line[0]
        chrom = line[2]
        start = int(line[3])
        umi = ''
        for element in line:
            if 'XU:Z:' in element:
                umi = element.split('XU:Z:')[1]
        postdedup_data.append([name,chrom,start,umi])

# combine pre and post deduplication data
    print('combining pre and post deduplication data...')
    prededup_data = sorted(prededup_data,key=lambda x: (x[1], x[0]))
    postdedup_data = sorted(postdedup_data,key=lambda x: (x[1], x[0]))
    combined_data = []
    a = 0
    b = 0
    while a < len(prededup_data):
        if b == len(postdedup_data):
            combined_data.append([prededup_data[a][0],prededup_data[a][1],prededup_data[a][2],prededup_data[a][3],'removed'])
            a += 1
        else:
            if prededup_data[a][0] == postdedup_data[b][0]:
                combined_data.append([prededup_data[a][0],prededup_data[a][1],prededup_data[a][2],prededup_data[a][3],'retained'])
                a += 1
                b += 1
            else:
                combined_data.append([prededup_data[a][0],prededup_data[a][1],prededup_data[a][2],prededup_data[a][3],'removed'])
                a += 1

    combined_data = sorted(combined_data,key=lambda x: (x[1], x[2]))
    
# print summary of alignments inputed
    print('{} pre-deduplication alignments'.format(a))
    print('{} post-deduplication alignments'.format(b))

# group alignments by chromosome - max 10,000 alignments by start position
    print('preparing alignments for multiprocessing...')
    multi_input = []
    single_entry = []
    position_entry = []
    entry_num = 0
    alignments_in_entry = 0
    current_chromosome = ''
    current_position = ''
    for a in range(0,len(combined_data)):
        alignment = combined_data[a]
        if alignments_in_entry > 10000:
            if current_chromosome == alignment[1] and current_position == alignment[2]:
                position_entry[2].append([alignment[0],alignment[3],alignment[4]])
                alignments_in_entry += 1
            else:
                entry_num += 1
                single_entry.append(position_entry)
                multi_input.append([entry_num,single_entry])
                single_entry = []
                position_entry = [alignment[1],alignment[2],[[alignment[0],alignment[3],alignment[4]]]]
                current_chromosome = alignment[1]
                current_position = alignment[2]
                alignments_in_entry = 1
        elif a == len(combined_data) - 1:
            entry_num += 1
            if current_chromosome == alignment[1] and current_position == alignment[2]: 
                position_entry[2].append([alignment[0],alignment[3],alignment[4]])
                single_entry.append(position_entry)
                multi_input.append([entry_num,single_entry]) 
            else: 
                single_entry.append(position_entry)
                position_entry = [alignment[1],alignment[2],[[alignment[0],alignment[3],alignment[4]]]]
                single_entry.append(position_entry)
                entry_num += 1
                multi_input.append([entry_num,single_entry])
        else:
            if current_chromosome == alignment[1] and current_position == alignment[2]:
                position_entry[2].append([alignment[0],alignment[3],alignment[4]])
                alignments_in_entry += 1
            else:
                if len(position_entry) > 0:
                    single_entry.append(position_entry)
                    position_entry = [alignment[1],alignment[2],[[alignment[0],alignment[3],alignment[4]]]]
                    alignments_in_entry += 1
                    current_chromosome = alignment[1]
                    current_position = alignment[2]
                else:
                    position_entry = [alignment[1],alignment[2],[[alignment[0],alignment[3],alignment[4]]]]
                    alignments_in_entry += 1
                    current_chromosome = alignment[1]
                    current_position = alignment[2]

# write multiprocessing input to file
    input_out = []
    for entry in multi_input:
        entry_num = str(entry[0])
        for position in entry[1]:
            chromosome = position[0]
            startpos = str(position[1])
            for alignment in position[2]:
                input_out.append([entry_num,chromosome,startpos,alignment[0],alignment[1],alignment[2]])

    with open(outputdir + '/multiprocessing_input.tsv', 'w') as outputfile:
        for entry in input_out:
            outputfile.write('{}\n'.format('\t'.join(entry)))

# multiprocessing
    print('multiprocessing...')
    print('{} entries to process'.format(len(multi_input)))

    builtin_outputs = map(process_alignments,multi_input)
    pool_size = multiprocessing.cpu_count()
    pool = multiprocessing.Pool(processes=pool_size)

    pool_outputs = list(tqdm.tqdm(pool.imap(process_alignments,multi_input), total=len(multi_input)))

    pool.close()
    pool.join()

    multi_output = []
    for pool_output in pool_outputs:
        multi_output.append(pool_output)

# rumi deduplication stats
    print('generating stats...')
    removal_stats = {}
    output_data = []

    for entry in multi_output:
        for alignment in entry:
            if alignment[4] == 'removed':
                removal_stats.setdefault(alignment[6],0)
                removal_stats[alignment[6]] += 1
            else:
                removal_stats.setdefault(alignment[4],0)
                removal_stats[alignment[4]] += 1
            output_data.append([alignment[0],str(alignment[1]),alignment[2],alignment[3],alignment[4],alignment[5],str(alignment[6])])
# write output
    with open(outputdir + '/alignment_data.tsv', 'w') as outputfile:
        for entry in output_data:
            outputfile.write('{}\n'.format('\t'.join(entry)))

    with open(outputdir + '/rumi_deduplication_stats.tsv', 'w') as outputfile:
        for hamming_distance,occurrences in removal_stats.items():
            outputfile.write('{}\t{}\n'.format(hamming_distance,occurrences))

# remove tmp dir
    shutil.rmtree(tmpdir)
